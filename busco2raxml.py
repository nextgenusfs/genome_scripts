#!/usr/bin/env python

import sys, argparse, os, subprocess, shutil
from Bio import SeqIO

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

parser=argparse.ArgumentParser(prog='busco2raxml.py',
    description='''Script runs BUSCO2 on several GBK files, then parses results, finds all in common, and draws tree using RAxML''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--input', nargs='+', required=True, help='Input GBK files')
parser.add_argument('-b','--busco', required=True, help='Path to BUSCO2')
parser.add_argument('-d','--busco_db', required=True, help='BUSCO2 DB to use')
parser.add_argument('-o','--out', required=True, help='Output name')
parser.add_argument('-n','--num', type=int, help='Number of BUSCO models to use for tree')
parser.add_argument('-c','--cpus', type=int, default=1, help='Number of CPUs')
args=parser.parse_args()

def gb2name(input):
    with open(input, 'rU') as gbk:
        for record in SeqIO.parse(gbk, 'genbank'):
            for f in record.features:
                if f.type == "source":
                    organism = f.qualifiers.get("organism", ["???"])[0]
                    isolate = f.qualifiers.get("isolate", ["???"])[0]
                    if isolate == "???":
                        isolate = f.qualifiers.get("strain", ["???"])[0]
    if organism == '???': #default to file name
        organism = os.path.basename(input).split('.',-1)[0]
    else:
        organism = organism.replace(' ', '_').lower()
    if isolate == '???':
        return organism
    else:  
        return organism+'_'+isolate

def gb2prots(input, tmpdir):
    outname = gb2name(input)+'.prots.fa'
    with open(os.path.join(tmpdir, outname), 'w') as proteins:
        with open(input, 'rU') as gbk:
            for record in SeqIO.parse(gbk, 'genbank'):
                for f in record.features:
                    if f.type == "CDS":
                        proteins.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], f.qualifiers['translation'][0]))

def parseBUSCO(input):
    #now parse output
    results = {}
    with open(input, 'rU') as busco:
        for line in busco:
            if line.startswith('#'):
                continue
            col = line.split('\t')
            if col[1] == 'Complete':
                ProtID = col[2]
                BuscoHit = col[0]
                if not BuscoHit in results:
                    results[BuscoHit] = ProtID
                else: #don't want any multiple hits, not sure that can ever happen but prevent here
                    del results[BuscoHit]
    return results

'''
first parse the input, if given proteins, i.e. ends with .fa, .fasta, .fna - then just copy to tempfolder
if passed GBK files, then convert to protein fasta files and put in tempfolder
then grab all protein fasta files in tempfolder and run BUSCO protein search
'''
#first parse the arguments
tmpdir = 'busco2raxml_'+str(os.getpid())
os.makedirs(tmpdir)

for i in args.input:
    if i.endswith('.gbk') or i.endswith('.gbf') or i.endswith('.gbff'):
        gb2prots(i, tmpdir)
    else:
        shutil.copyfile(i, os.path.join(tmpdir, i))

#now only files in tmpdir are protein fasta files, grab them all
file_list = []
for file in os.listdir(tmpdir):
    file_list.append(file)

#now loop through each and run BUSCO, collect complete results in list of dictionaries
AllResults = []
for x in file_list:
    name = os.path.basename(x).split('.',-1)[0]
    cmd = [args.busco, '-i', x, '-m', 'proteins', '-l', os.path.abspath(args.busco_db), '-o', name, '-c', str(args.cpus), '-f']
    subprocess.call(cmd, cwd=tmpdir)
    BuscoResults = parseBUSCO(os.path.join(tmpdir, 'run_'+name, 'full_table_'+name+'.tsv'))
    AllResults.append(BuscoResults)

#now get a list of all
AllBuscos = []
for x in AllResults:
    for k,v in x.items():
        if not k in AllBuscos:
            AllBuscos.append(k)

#loop through all buscos and determine if present in every dictionary
BadBuscos = []
for x in AllBuscos:
    for y in AllResults:
        if not x in y:
            BadBuscos.append(x)

#now find those that are found in all results
BadBuscos = set(BadBuscos)
Keepers = [x for x in AllBuscos if x not in BadBuscos]

#now loop through data and pull out the proteins for each
with open(args.out, 'w') as output:
    for i in range(0,len(args.input)):
        SpeciesName = os.path.basename(file_list[i]).split('.',-1)[0]
        Proteins = SeqIO.index(os.path.join(tmpdir, file_list[i]), 'fasta')
        Seq = []
        for y in Keepers:
            ID = AllResults[i].get(y)
            Seq.append(str(Proteins[ID].seq))
        output.write('>%s\n%s\n' % (SpeciesName, ''.join(Seq)))

#finalize
print "Found %i BUSCOs conserved in all genomes" % len(Keepers)
print("%s\n" %  ', '.join(Keepers))
print "Concatenated protein sequences for all %i genomes located in: %s" % (len(args.input), args.out)

