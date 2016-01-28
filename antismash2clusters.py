#!/usr/bin/env python

import sys, os, shutil, subprocess, argparse
from Bio import SeqIO
from natsort import natsorted

#can make this a bit more complicated but simplify output
#would require bedtools, bwa, samtools

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='antismash2clusters.py',
    description='''Script that can map paired end jumping library sequences to genome, parse antiSMASH SM cluster results, and determine with jumping library hits contain secondary metabolite clusters.''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='Genome in FASTA format')
parser.add_argument('-g','--gff', required=True, help='GFF feature file for genome')
parser.add_argument('-a','--antismash', required=True, help='antiSMASH results in GenBank format (.gbk)')
parser.add_argument('-o','--out', required=True, help='Output files base name')
parser.add_argument('--tmpdir', default='antismash2clusters', help='Folder name to collect intermediate files')
parser.add_argument('-f','--fwd_reads', help='Forward reads from jumping library (BAC forward sequences)')
parser.add_argument('-r','--rev_reads', help='Reverse reads from jumping library (BAC reverse sequences)')
parser.add_argument('--cpus', default='1', help='Number of CPUs to use with BWA)')
parser.add_argument('--lib_range', default='5000-200000', help='Range of insert sizes in jumping library')
parser.add_argument('--cluster_padding', default=15000, type=int, help='Number of bp to pad each side of predicted cluster for summary output')
args=parser.parse_args()

def which(name):
    try:
        with open(os.devnull) as devnull:
            subprocess.Popen([name, '--version'], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True

def CheckDependencies(input):
    missing = []
    for p in input:
        if which(p) == False:
            missing.append(p)
    if missing != []:
        error = ", ".join(missing)
        log.error("Missing Dependencies: %s.  Please install missing dependencies and re-run script" % (error))
        sys.exit(1)

def countfasta(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if line.startswith (">"):
                count += 1
    return count

def countfastq(input):
    lines = sum(1 for line in open(input))
    count = int(lines) / 4
    return count

def ParseAntiSmash(input, tmpdir, output):
    print("Now parsing antismash genbank result, finding clusters and backbone enzymes")
    global BackBone, SMCOGs, bbSubType, bbDomains, Offset, smProducts, InterPro, PFAM
    BackBone = {}; SMCOGs = {}; bbSubType = {}; bbDomains = {}; Offset = {}; smProducts = {}; InterPro = {}; PFAM = {}
    backboneCount = 0; clusterCount = 0; cogCount = 0
    #parse antismash genbank to get clusters in bed format and slice the record for each cluster prediction
    with open(output, 'w') as antibed:
        with open(input, 'rU') as input:
            SeqRecords = SeqIO.parse(input, 'genbank')
            for record in SeqRecords:
                for f in record.features:
                    if f.type == "source":
                        record_start = f.location.start
                        record_end = f.location.end
                    if f.type == "cluster":
                        clusterCount += 1
                        chr = record.id
                        start = f.location.start
                        end = f.location.end
                        clusternum = f.qualifiers.get("note")[0].replace("Cluster number: ", "")
                        antibed.write("%s\t%s\t%s\tCluster_%s\t0\t+\n" % (chr, start, end, clusternum))
                        sub_start = start - args.cluster_padding
                        sub_end = end + args.cluster_padding
                        if sub_start < 1:
                            sub_start = 1
                        if sub_end > record_end:
                            sub_end = record_end
                        sub_record = record[sub_start:sub_end]
                        Offset['Cluster_'+clusternum] = sub_start
                        sub_record_name = os.path.join(tmpdir, 'Cluster_'+clusternum+'.gbk')
                        with open(sub_record_name, 'w') as clusterout:
                            SeqIO.write(sub_record, clusterout, 'genbank')
                    Domains = []
                    if f.type == "CDS":
                        ID = f.qualifiers.get('locus_tag')[0]
                        #if not ID.endswith('-T1'):
                            #ID = ID + '-T1'                       
                        if f.qualifiers.get('sec_met'):            
                            for k, v in f.qualifiers.items():
                                if k == 'sec_met':
                                    for i in v:
                                        if i.startswith('Type:'):
                                            type = i.replace('Type: ', '')
                                            backboneCount += 1
                                            BackBone[ID] = type
                                        if i.startswith('NRPS/PKS subtype:'):
                                            subtype = i.replace('NRPS/PKS subtype: ', '')
                                            bbSubType[ID] = subtype
                                        if i.startswith('NRPS/PKS Domain:'):
                                            doms = i.replace('NRPS/PKS Domain: ', '')
                                            doms = doms.split('. ')[0]
                                            Domains.append(doms)
                                bbDomains[ID] = Domains
                        for k,v in f.qualifiers.items():
                            if k == 'note':
                                for i in v:
                                    if i.startswith('smCOG:'):
                                        COG = i.replace('smCOG: ', '')
                                        COG = COG.split(' (')[0]
                                        SMCOGs[ID] = COG
                                        cogCount += 1
                                    elif not i.startswith('smCOG tree'):
                                        notes = i
                                        smProducts[ID] = notes
                            if k == 'db_xref':
                                for i in v:
                                    if i.startswith('InterPro:'):
                                        interpro = i.replace('InterPro:', '')
                                        if not ID in InterPro:
                                            InterPro[ID] = [interpro]
                                        else:
                                            InterPro[ID].append(interpro)
                                    if i.startswith('PFAM:'):
                                        pfam = i.replace('PFAM:', '')
                                        if not ID in PFAM:
                                            PFAM[ID] = [pfam]
                                        else:
                                            PFAM[ID].append(pfam)
                            
            print("Found %i clusters, %i biosynthetic enyzmes, and %i smCOGs predicted by antiSMASH" % (clusterCount, backboneCount, cogCount))

def GetClusterGenes():
    global dictClusters
    #pull out genes in clusters from GFF3, load into dictionary
    print("Finding genes in each cluster")
    with open(GenesInClusters, 'w') as output:
        subprocess.call(['bedtools', 'intersect','-wo', '-a', AntiSmashBed, '-b', GFF], stdout = output)
    dictClusters = {}
    with open(GenesInClusters, 'rU') as input:
        for line in input:
            cols = line.split('\t')
            if cols[8] != 'gene':
                continue
            gene = cols[14].replace('ID=', '')
            ID = cols[3]
            if ID not in dictClusters:
                dictClusters[ID] = [gene]
            else:
                dictClusters[ID].append(gene)

#run some checks of dependencies first
programs = ['bwa', 'samtools', 'bedtools']
CheckDependencies(programs)

#constraints on the BAC sizes, sometimes some map outside of a realistic range
if not '-' in args.lib_range:
    print("Error, you need to supply a range for BAC insert size, i.e. 5000-200000")
    os._exit(1)
else:
    range = args.lib_range.split('-')
    smallest = int(range[0])
    largest = int(range[1])

if args.fwd_reads:
    if not args.rev_reads:
        print("You supplied forward reads, but not reverse reads")
        os._exit(1)
    #check if fasta or fastq
    with open(args.fwd_reads, 'rU') as input:
        for line in input:
            if line.startswith('>'):
                format = 'fasta'
                break
            elif line.startswith('@'):
                format = 'fastq'
                break
            else:
                print("Reads are not properly formatted, need to be either FASTA or FASTQ format")
                os._exit(1)
    #now make sure properly paired
    if format == 'fasta':
        forwardCount = countfasta(args.fwd_reads)
        reverseCount = countfasta(args.rev_reads)
    elif format == 'fastq':
        forwardCount = countfastq(args.fwd_reads)
        reverseCount = countfastq(args.rev_reads)
    if forwardCount != reverseCount:
        print("Reads do not appear to be properly paired, exiting")
        os._exit(1)

#define files/input/output
genome = args.input
GFF = args.gff
f_reads = args.fwd_reads
r_reads = args.rev_reads
antismash = args.antismash
outputDir = args.tmpdir
cpus = args.cpus
FinalOut = args.out + '.bac.overlap.txt'
ClustersOut = args.out + '.secmet.clusters.txt'
sam_out = os.path.join(outputDir, 'bwa.sam')
bam_out = os.path.join(outputDir, 'bwa.bam')
sortbam = os.path.join(outputDir, 'bwa.sort')
realsortedbam = sortbam + '.bam'
BedPE = os.path.join(outputDir, 'bedpe.bed')
bacBED = os.path.join(outputDir, 'reads.bed')
AntiSmashBed = os.path.join(outputDir, 'antismash.bed')
ClusterBed = os.path.join(outputDir, 'full.clusters.bed')
GenesInBacs = os.path.join(outputDir, 'genes.in.bacs.txt')
GenesInClusters = os.path.join(outputDir, 'genes.in.clusters.txt')
FNULL = open(os.devnull, 'w')

#first would be to create output folder, create bwa index and then map reads
if not os.path.isdir(outputDir):
    os.makedirs(outputDir)

if args.fwd_reads:
    if not os.path.isfile(os.path.join(outputDir, 'bwa.sort.bam')):
        #move genome file into output folder to make index
        tmp_genome = os.path.join(outputDir, 'genome.fasta')
        shutil.copyfile(genome, os.path.join(outputDir, 'genome.fasta'))

        #create bwa index
        print("Creating BWA index...")
        subprocess.call(['bwa', 'index', tmp_genome], stderr = FNULL)

        #now map reads to indexed genome
        print("Mapping reads to indexed genome...")
        with open(sam_out, 'w') as output:
            subprocess.call(['bwa', 'mem', '-t', cpus, tmp_genome, f_reads, r_reads], stdout = output, stderr = FNULL)
        with open(bam_out, 'w') as output:
            subprocess.call(['samtools', 'view', '-bS', sam_out], stdout = output)
        subprocess.call(['samtools', 'sort', '-n', bam_out, sortbam])
    else:
        print("BWA mapping has already been run, moving on...")

    #now convert sorted bam to bedpe
    print("Converting BAM output to PE bed file")
    with open(BedPE, 'w') as output:
        subprocess.call(['bedtools', 'bamtobed', '-bedpe', '-i', realsortedbam], stdout = output, stderr = FNULL)

    #convert bedpe to single interval bam
    print("Converting PE Bed file to BAC interval bed")
    with open(bacBED, 'w') as output:
        with open(BedPE, 'rU') as input:
            for line in input:
                if line.startswith('.'):
                    continue
                cols = line.split('\t')
                if cols[0] == cols[3] and cols[8] != cols[9]:
                    chr = cols[0]
                    start = cols[1]
                    end = cols[5]
                    name = cols[6]
                    length = int(end) - int(start)
                    output.write("%s\t%s\t%s\t%s\t0\t+\n" % (chr, start, end, name))

    #parse antismash and get genes in cluster
    ParseAntiSmash(antismash, outputDir, AntiSmashBed)
    GetClusterGenes()

    #do an intersection of the bacs with the gff3 file to eventually get the genes overlapping in each bac - load into dictionary
    print("Finding genes located in each BAC")
    with open(GenesInBacs, 'w') as output:
        subprocess.call(['bedtools', 'intersect','-wo', '-a', bacBED, '-b', GFF], stdout = output)
    #now convert to dictionary
    dictBAC = {}
    with open(GenesInBacs, 'rU') as input:
        for line in input:
            cols = line.split('\t')
            if cols[8] != 'gene':
                continue
            gene = cols[14].replace('ID=', '')
            ID = cols[3]
            if ID not in dictBAC:
                dictBAC[ID] = [gene]
            else:
                dictBAC[ID].append(gene)

    #now intersect these bed files finding BACs that contain even partial clusters
    print("Determining overlaps between BACs and SM clusters")
    with open(ClusterBed, 'w') as output:
        subprocess.call(['bedtools', 'intersect', '-wo', '-a', AntiSmashBed, '-b', bacBED], stdout = output)

    #finally parse all information to create a tab delimited output of results
    with open(FinalOut, 'w') as output:
        output.write("#Cluster\tCluster Location\tCluster length\tBackbone enzymes(s)\tNumber of cluster genes\tBAC name\tBAC location\tBAC length\tNumber of genes in BAC\tCluster Coverage\tCluster genes\tBAC genes\n")
        with open(ClusterBed, 'rU') as input:
            for line in input:
                cols = line.split('\t')
                cluster = cols[3]
                cluster_genes = dictClusters[cluster]#this should be list of genes in cluster
                backbone_enzyme = []
                for k, v in BackBone.items(): #loop through backbone genes and grab those from cluster
                    if k in cluster_genes:
                        f_enz = k + ' (' + v + ')'
                        backbone_enzyme.append(f_enz)
                enzymes = ", ".join([str(x) for x in backbone_enzyme])
                scaffold = cols[0]
                cluster_start = cols[1]
                cluster_end = cols[2]
                cluster_length = int(cluster_end) - int(cluster_start)
                BAC = cols[9]
                BAC_start = cols[7]
                BAC_end = cols[8]
                BAC_length = int(BAC_end) - int(BAC_start)
                if BAC_length < smallest or BAC_length > largest:
                    continue
                coverage = int(cols[12]) / float(cluster_length) * 100.0
                gene_list = dictBAC.get(BAC)
                if gene_list:
                    genes = ", ".join([str(x) for x in gene_list])
                else:
                    genes = 'none'
                Cluster = ", ".join([str(x) for x in cluster_genes])
                #subtract SM cluster list from BAC list to get which genes are missing
                gene_list = set(gene_list)
                miss_genes = [x for x in cluster_genes if x not in gene_list]
                if miss_genes:
                    missing = ", ".join([str(x) for x in miss_genes])
                else:
                    missing = 'none'
                output.write("%s\t%s:%s-%s\t%s\t%s\t%i\t%s\t%s:%s-%s\t%s\t%i\t%s%%\t%s\t%s\t%s\n" % (cluster, scaffold, cluster_start, cluster_end, cluster_length, enzymes, len(cluster_genes), BAC, scaffold, BAC_start, BAC_end, BAC_length, len(gene_list), coverage, Cluster, genes, missing))

else: #do not map reads, but just reformat antiSMASH results in same way as above

    ParseAntiSmash(antismash, outputDir, AntiSmashBed)
    GetClusterGenes()

#so if you got here, that means you have 1) already created your BAC file and/or 2) have clusters split into gbk files
#now lets pull out information for each cluster into a tab delimited file
'''#dictionaries available now 
BackBone --> ID = type
SMCOGs --> ID = SMCOG
bbSubType --> ID = NRPS/PKS subtype
dictClusters --> ClusterID = [genes]
bbDomains --> ID = [domains]
'''
#print SMCOGs
#print BackBone
#print bbSubType
#print bbDomains
#print dictClusters
#print Offset
#print smProducts

for file in os.listdir(outputDir):
    if file.endswith('.gbk'):
        base = file.replace('.gbk', '')
        base = base.capitalize()
        outputName = os.path.join(outputDir, base+'.secmet.cluster.txt')
        file = os.path.join(outputDir, file)
        with open(outputName, 'w') as output:
            output.write("#%s\n" % base)
            output.write("#GeneID\tChromosome:start-stop\tStrand\tClusterPred\tBackbone Enzyme\tDomains\tProduct\tsmCOGs\tInterPro\tPFAM\tNote\tProtein Seq\tDNA Seq\n")
            with open(file, 'rU') as input:
                SeqRecords = SeqIO.parse(input, 'genbank')
                for record in SeqRecords:
                    for f in record.features:
                        if f.type == "CDS":
                            name = f.qualifiers["locus_tag"][0]
                            prot_seq = f.qualifiers['translation'][0]
                            cluster_genes = dictClusters.get(base) #this should be list of genes in each cluster
                            start = f.location.nofuzzy_start
                            actualStart = int(start) + int(Offset.get(base)) + 1 #account for python numbering shift?
                            end = f.location.nofuzzy_end
                            actualEnd = int(end) + int(Offset.get(base))
                            strand = f.location.strand
                            if strand == 1:
                                strand = '+'
                                DNA_seq = record.seq[start:end]
                            elif strand == -1:
                                strand = '-'
                                DNA_seq = record.seq[start:end].reverse_complement()
                            chr = record.id
                            product = f.qualifiers["product"][0]
                            if name in cluster_genes:
                                location = 'cluster_pred'
                            else:
                                location = 'flanking'
                            if name in bbSubType:
                                enzyme = bbSubType.get(name)
                            else:
                                if name in BackBone:
                                    enzyme = BackBone.get(name)
                                else:
                                    enzyme = '.'
                            if name in SMCOGs:
                                cog = SMCOGs.get(name)
                                #cog = enzyme.split(":")[-1]
                            else:
                                cog = '.'
                            if name in bbDomains:
                                domains = ";".join(bbDomains.get(name))
                            else:
                                domains = '.'
                            if name in smProducts:
                                note = smProducts.get(name)
                            else:
                                note = '.'
                            if name in InterPro:
                                IP = ";".join(InterPro.get(name))
                            else:
                                IP = '.'
                            if name in PFAM:
                                PF = ":".join(PFAM.get(name))
                            else:
                                PF = '.'
                            output.write("%s\t%s:%i-%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (name, chr, actualStart, actualEnd, strand, location, enzyme, domains, product, cog, IP, PF, note, prot_seq, DNA_seq))
                                              
#now put together into a single file
finallist = []
for file in os.listdir(outputDir):
    if file.endswith('secmet.cluster.txt'):
        file = os.path.join(outputDir, file)
        finallist.append(file)
with open(ClustersOut, 'w') as output:
    for file in natsorted(finallist):
        with open(file, 'rU') as input:
            output.write(input.read())
            output.write('\n\n')

os._exit(1)