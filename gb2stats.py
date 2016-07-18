#!/usr/bin/env python

import sys
from Bio import SeqIO

def genomeStats(input):
    from Bio.SeqUtils import GC
    lengths = []
    GeeCee = []
    Genes = 0
    tRNA = 0
    Prots = 0
    with open(input, 'rU') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'genbank')
        for record in SeqRecords:
            lengths.append(len(record.seq))
            GeeCee.append(str(record.seq))
            for f in record.features:
                if f.type == "source":
                    organism = f.qualifiers.get("organism", ["???"])[0]
                    isolate = f.qualifiers.get("isolate", ["???"])[0]
                if f.type == "CDS":
                    Prots += 1
                if f.type == "gene":
                    Genes += 1
                if f.type == "tRNA":
                    tRNA += 1
    GenomeSize = sum(lengths)
    LargestContig = max(lengths)
    ContigNum = len(lengths)
    AvgContig = int(round(GenomeSize / ContigNum))
    pctGC = round(GC("".join(GeeCee)), 2)
    
    #now get N50
    lengths.sort()
    nlist = []
    for x in lengths:
        nlist += [x]*x
    if len(nlist) % 2 == 0:
        medianpos = int(len(nlist) / 2)
        N50 = int((nlist[medianpos] + nlist[medianpos-1]) / 2)
    else:
        medianpos = int(len(nlist) / 2)
        N50 = int(nlist[medianpos])
        
    #get L50, number of scaffods
    LfiftyLen = GenomeSize / 2
    cumulativelength = 0
    Lfifty = 0
    for i in sorted(lengths, reverse=True):
        if cumulativelength < LfiftyLen:
            cumulativelength += i
            Lfifty += 1
    
    #return values in a list
    return [organism, isolate, "{0:,}".format(GenomeSize)+' bp', "{0:,}".format(LargestContig)+' bp', "{0:,}".format(AvgContig)+' bp', "{0:,}".format(ContigNum), "{0:,}".format(N50)+' bp', "{0:,}".format(Lfifty),"{:.2f}".format(pctGC)+'%', "{0:,}".format(Genes), "{0:,}".format(Prots), "{0:,}".format(tRNA)]

stats = genomeStats(sys.argv[1])

print "-----------------------------------------"
print("%s %s" % (stats[0], stats[1]))
print "-----------------------------------------"
print ("Genome size:\t%s" % (stats[2]))
print ("Contigs:\t%s" % (stats[5]))
print ("Largest Contig:\t%s" % (stats[3]))
print ("Average Contig:\t%s" % (stats[4]))
print ("Contig N50:\t%s" % (stats[6]))
print ("Contig L50:\t%s" % (stats[7]))
print ("GC Content:\t%s" % (stats[8]))
print ("Genes:\t\t%s" % (stats[9]))
print ("Protein Coding:\t%s" % (stats[10]))
print ("tRNA:\t\t%s" % (stats[11]))

