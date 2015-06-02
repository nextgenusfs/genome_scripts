#!/usr/bin/env python
#This script generates some basic genome stats from a genbank flat file format (.gbk)
#Run like: gb_stats.py file.gbk

#The N50 calculation is from Cara Magnabosco https://caramagnabosco.wordpress.com/2014/02/20/calculate-the-n50-of-assembled-contigs/


import sys
import numpy
from Bio import SeqIO

gbk_filename = sys.argv[1]

gene_length = 0
avg_gene_length = 0

num_records = list(SeqIO.parse(gbk_filename, "genbank"))
contigs = len(num_records)

percent = "%"
ending = "Mb"
ending2 = "bp"
cds_count = 0
gene_count = 0
trna_count = 0
total = 0
C_count = 0
G_count = 0
seqlength = []
for record in SeqIO.parse(gbk_filename, "genbank"):
    total += len(record.seq)
    bp = len(record.seq)
    seqlength.append(bp)
    mb_total = float(total) / 1000000
    C_count += record.seq.count("C")
    G_count += record.seq.count("G")
    GC_total = C_count + G_count
    GC_count = GC_total / float(total) * 100  
    for f in record.features:
        if f.type == "source":
            organism = f.qualifiers.get("organism", ["???"])[0]
            isolate = f.qualifiers.get("strain", ["???"])[0]
        if f.type == "CDS":
            cds_count = cds_count + 1
        if f.type == "gene":
        	gene_count = gene_count + 1
        if f.type == "tRNA":
            trna_count = trna_count + 1
        if f.type == "gene":
	        gene_length = gene_length + len(f) 
	        
avg_gene_length = gene_length / gene_count 
         
seqlength = sorted(seqlength)
unique = []
for entry in seqlength:
    if not entry in unique:
        unique.append(entry)

n50 = []
for entry in unique:
    multiplier = seqlength.count(entry) * entry
    for i in range(multiplier):
        n50.append(entry)
index = len(n50)/2
avg = []

if index % 2==0:
    first = n50[index -1]
    second = n50[index]
    avg.append(first)
    avg.append(second)
    n50 = numpy.mean(avg)
    n50 = n50
else:
    n50 = n50[index -1]
   
print "-----------------------------------------"
print("%s %s" % (organism, isolate))
print "-----------------------------------------"
print ("Scaffolds:\t{:,}".format(contigs))
print ("Total bp:\t{:.2f} {}".format (mb_total, ending))
print ("Largest Contig:\t{:,} {}".format(seqlength[-1], ending2))
print ("Scaffold N50:\t{:,g} {}".format(n50, ending2))
print ("GC Content:\t%.2f %s" % (GC_count,percent))
print ("Genes:\t\t{:,}".format(gene_count))
print ("Protein Coding:\t{:,}".format(cds_count))
print ("tRNA:\t\t{:,}".format(trna_count))
print ("Avg gene len:\t{:,g} {}".format (avg_gene_length, ending2))


