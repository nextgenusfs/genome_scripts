#!/usr/bin/env python

import sys, os, subprocess, argparse, shutil, warnings, csv, multiprocessing
from Bio import SeqIO
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SearchIO

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='findHMM.py',
    description='''Find HMM model in GenBank flatfile genome.''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='Genome (GBK) format')
parser.add_argument('-m','--hmm', required=True, help='HMM file')
parser.add_argument('-o','--out', required=True, help='Protein FASTA output')
parser.add_argument('--evalue', default='1e-50', help='HMM evalue')
parser.add_argument('--cpus', default=2, type=int, help='Number of CPUs')
parser.add_argument('--debug', action='store_true', help='Debug')
parser.add_argument('--maxIntron', default='3000', help='Maximum intron length')
args=parser.parse_args()
FNULL = open(os.devnull, 'w')

def getSize(filename):
    st = os.stat(filename)
    return st.st_size

def group_by_heading( some_source ):
    buffer= []
    for line in some_source:
        if line.startswith( ">" ):
            if buffer: yield buffer
            buffer= [ line ]
        else:
            buffer.append( line )
    yield buffer

def gb2output(input, output1, output2, output3):
    with open(output1, 'w') as proteins:
        with open(output2, 'w') as transcripts:
            with open(output3, 'w') as scaffolds:
                with open(input, 'rU') as gbk:
                    SeqRecords = SeqIO.parse(gbk, 'genbank')
                    for record in SeqRecords:
                        scaffolds.write(">%s\n%s\n" % (record.id, record.seq))
                        for f in record.features:
                            if f.type == "CDS":
                                proteins.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], f.qualifiers['translation'][0]))
                            if f.type == "mRNA":
                                feature_seq = f.extract(record.seq)
                                transcripts.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], feature_seq))

def tblastnFilter(input, query, cpus, output):
    global HitList
    HitList = []
    #start by formatting blast db/dustmasker filtered format
    subprocess.call(['dustmasker', '-in', input, '-infmt', 'fasta', '-parse_seqids', '-outfmt', 'maskinfo_asn1_bin', '-out', os.path.join(output,'genome_dust.asnb')], stdout = FNULL, stderr = FNULL)
    subprocess.call(['makeblastdb', '-in', input, '-dbtype', 'nucl', '-parse_seqids', '-mask_data', os.path.join(output, 'genome_dust.asnb'), '-out', os.path.join(output, 'genome')], stdout = FNULL, stderr = FNULL)
    #okay, now run tblastn using uniprot proteins
    subprocess.call(['tblastn', '-num_threads', str(args.cpus), '-db', os.path.join(output, 'genome'), '-query', query, '-max_target_seqs', '1', '-db_soft_mask', '11', '-threshold', '999', '-max_intron_length', args.maxIntron, '-evalue', '1e-10', '-outfmt', '6', '-out', os.path.join(output,'filter.tblastn.tab')], stdout = FNULL, stderr = FNULL)
    #now parse through results, generating a list for exonerate function
    with open(os.path.join(output, 'filter.tblastn.tab')) as input:
        reader = csv.reader(input, delimiter='\t')
        for cols in reader:
            hit = cols[0] + '::' + cols[1]
            if hit not in HitList:
                HitList.append(hit)

def runExonerate(input):
    global Missing
    s = input.split('::')
    if s[0].startswith('sp|'):
        name = s[0].split("|")[1] + '_' + s[1]
    else:
        name = s[0].split()[0] + '_' + s[1]
    query = os.path.join(tmpdir, name+'.fa')
    with open(query, 'w') as output:
        rec = record_dict[s[0]]
        output.write(">%s\n%s\n" % (rec.id, rec.seq))
    scaffold = s[1] + '.fa'
    scaffold = os.path.join(tmpdir, scaffold)
    exonerate_out = 'exonerate_' + name + '.out'
    exonerate_out = os.path.join(tmpdir, exonerate_out)
    ryo = ">%qi|pident=%pi|%ti:%tcb-%tce|Exonerate-Partial\n%tcs\n"
    with open(exonerate_out, 'w') as output4:
        subprocess.call(['exonerate', '--model', 'p2g', '--showvulgar', 'no', '--showalignment', 'no', '--showquerygff', 'no', '--showtargetgff', 'no', '--maxintron', args.maxIntron, '--percent', '25', '--ryo', ryo , query, scaffold], stdout = output4)
    os.remove(query)


#set up tmpdir
tmpdir = 'tmp_'+str(os.getpid())
os.makedirs(tmpdir)

#Split GBK into parts
base = args.input.rsplit('.', 1)[0]
if '/' in base:
    base = base.split('/') [-1]
Proteins = os.path.join(tmpdir, base+'.proteins.fa')
Transcripts = os.path.join(tmpdir, base+'.transcripts.fa')
Genome = os.path.join(tmpdir, base+'.genome.fa')
gb2output(args.input, Proteins, Transcripts, Genome)

#load proteins into dictionary
protein_dict = SeqIO.to_dict(SeqIO.parse(Proteins, 'fasta'))

#check number of HMMer models
HMMstat = os.path.join(tmpdir, 'hmmstat.txt')
with open(HMMstat, 'w') as output:
    subprocess.call(['hmmstat', args.hmm], stdout = output, stderr = FNULL)
HMMmodels = []
with open(HMMstat, 'rU') as input:
    for line in input:
        if line.startswith('\n'):
            continue
        if not line.startswith('#'):
            cols = line.split(' ')
            cols = filter(None, cols)
            if not cols[1] in HMMmodels:
                HMMmodels.append(cols[1])
print "Looking for %i protein HMM model(s)" % len(HMMmodels)
#do hmmer search of proteins
print "Scanning proteome using HMMsearch"
HMM = os.path.join(tmpdir, base+'.hmmsearch.txt')
subprocess.call(['hmmsearch', '-o', HMM, '--cpu', str(args.cpus), '-E', args.evalue, args.hmm, Proteins], stdout = FNULL, stderr = FNULL)
Results = {}
with open(HMM, 'rU') as results:
    for qresult in SearchIO.parse(results, "hmmer3-text"):
        query_length = qresult.seq_len #length of HMM model
        hits = qresult.hits
        num_hits = len(hits)
        if num_hits > 0:
            query = hits[0].id
            hit = hits[0].query_id
            score = hits[0].bitscore
            evalue = hits[0].evalue
            num_hsps = len(hits[0].hsps)
            aln_length = 0
            for x in range(0,num_hsps):
                aln_length += hits[0].hsps[x].aln_span
            if hit not in Results:
                Results[hit] = (query, score, evalue, aln_length)
 
with open(args.out, 'w') as output:
    for k,v in Results.items():
        description = k+"|"+v[0]+"|evalue="+str(v[2])+"|HMMer3-Complete"
        rec = protein_dict[v[0]]
        rec.id = description
        rec.description = ''
        rec.name = ''
        SeqIO.write(rec, output, 'fasta')
        
notfound = []
for i in HMMmodels:
    if not i in Results:
        notfound.append(i)

if len(notfound) > 0: #have to do some more work here for these to be sure they really don't exist
    #get consensus from hmm model
    print "%i missing models [%s]" % (len(notfound), ', '.join(notfound))
    Consensus = os.path.join(tmpdir, 'missing.consensi.tmp')
    Consensi = os.path.join(tmpdir, 'missing.consensi.fa')
    with open(Consensus, 'w') as output1:
        subprocess.call(['hmmemit', '-c', args.hmm], stdout = output1, stderr = FNULL)
    with open(Consensi, 'w') as output2:
        with open(Consensus, 'rU') as input:
            for rec in SeqIO.parse(input, 'fasta'):
                rec.id = rec.id.replace('-consensus', '')
                rec.name = ''
                rec.description = ''
                if rec.id in notfound:
                    SeqIO.write(rec, output2, 'fasta')

    #now run tblastn against genome with those notfound
    Blast = os.path.join(tmpdir, 'tblastn.blast.tab')
    print "Try to recover models using tBlastn/Exonerate"
    tblastnFilter(Genome, Consensi, args.cpus, tmpdir)
    print "found %i preliminary tBlastn alignments" % (len(HitList))

    #split genome fasta into individual scaffolds
    with open(Genome, 'rU') as input:
        for record in SeqIO.parse(input, "fasta"):
            SeqIO.write(record, os.path.join(tmpdir, record.id + ".fa"), "fasta")

    #Now run exonerate on hits
    print "Polishing alignments with Exonerate"
    record_dict = SeqIO.to_dict(SeqIO.parse(Consensi, 'fasta'))
    p = multiprocessing.Pool(args.cpus)
    rs = p.map_async(runExonerate, HitList)
    p.close()
    while (True):
        if (rs.ready()): break

    #now collect all exonerate results into one
    print "Saving all hits to %s" % args.out
    Exonerate = os.path.join(tmpdir, 'exonerate.output.txt')
    skip = ['Command line', '%', ' ','tmp_', '\n', '--', 'Hostname']
    with open(Exonerate, 'w') as output5:
        for root, dirs, files in os.walk(tmpdir):
            for file in files:
                if file.endswith('.out'):
                    filename = os.path.join(root, file)
                    with open(filename, 'rU') as readfile:
                        for line in group_by_heading(readfile):
                            for i in line:
                                if not any(i.startswith(x) for x in skip):
                                    output5.write(i)
    with open(args.out, 'ab') as output6:
        with open(Exonerate, 'rU') as input:
            for rec in SeqIO.parse(input, 'fasta'):
                output6.write('>%s\n%s\n' % (rec.id, rec.seq.translate()))
            

if not args.debug:
    shutil.rmtree(tmpdir)
