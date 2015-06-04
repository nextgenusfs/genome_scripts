#!/usr/bin/env python

#This script converts genbank format into tab delimited SMURF format
#written by Jon Palmer palmer.jona at gmail dot com

import sys
from Bio import SeqIO
import argparse

parser=argparse.ArgumentParser(
    description='''Script that generates SMURF input from Genbank Flat File ''',
    epilog="""Written by Jon Palmer (2015)  palmer.jona@gmail.com""")
parser.add_argument('--gb', help='Genbank file')
parser.add_argument('--prot', default="smurf.proteins.fasta", help='Fasta output')
parser.add_argument('--dna', default="smurf.scaffolds.fasta", help='Fasta output')
parser.add_argument('--smurf', default="smurf.annotations.txt", help='Smurf output')
args=parser.parse_args()

if len(sys.argv) < 2:
    parser.print_usage()
    sys.exit(1)

gbk_file = args.gb
smurf_file = args.smurf
fasta_file = args.prot
contig_file = args.dna

smurf_file_handle = open(smurf_file, "w")
fasta_file_handle = open(fasta_file, "w")
contig_file_handle = open(contig_file, "w")

for record in SeqIO.parse(gbk_file, "genbank"):
    scaffold_name = record.name[6:]
    contig_file_handle.write(">%s\n%s\n" % (
        scaffold_name.lstrip("0"), 
        record.seq))
    for f in record.features:
        name = record.name[6:]
        if f.type == "CDS":
            fasta_file_handle.write(">%s\n%s\n" % (
                   f.qualifiers['locus_tag'][0],
                   f.qualifiers['translation'][0]))
            locus_tag = f.qualifiers.get("locus_tag", ["No ID"])[0]
            product_name = f.qualifiers.get("product", ["No Description"])[0]
            mystart = f.location.start
            myend = f.location.end
            strand = f.location.strand
            if strand == 1:
                smurf_file_handle.write("%s\t%s\t%s\t%s\t%s\n" % (locus_tag, name.lstrip("0"), int(mystart), int(myend), product_name))
            else:
                smurf_file_handle.write("%s\t%s\t%s\t%s\t%s\n" % (locus_tag, name.lstrip("0"), int(myend), int(mystart), product_name))

contig_file_handle.close()
fasta_file_handle.close()
smurf_file_handle.close()
