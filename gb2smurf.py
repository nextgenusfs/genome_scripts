#!/usr/bin/env python

#This script converts genbank format into tab delimited SMURF format
#written by Jon Palmer palmer.jona at gmail dot com

import sys
import re
from Bio import SeqIO
import argparse

parser=argparse.ArgumentParser(
    description='''Script that generates SMURF input from Genbank Flat File ''',
    epilog="""Written by Jon Palmer (2015)  palmer.jona@gmail.com""")
parser.add_argument('gbk', help='Genbank file (Required)')
parser.add_argument('-p','--prot', default="smurf.proteins.fasta", help='Fasta output')
parser.add_argument('-g','--dna', default="smurf.scaffolds.fasta", help='Fasta output')
parser.add_argument('-s','--smurf', default="smurf.annotations.txt", help='Smurf output')
parser.add_argument('--ncbi', action='store_true', help='NCBI chromosome names')
parser.add_argument('--jgi', action='store_true', help='JGI naming')
args=parser.parse_args()

if len(sys.argv) < 2:
    parser.print_usage()
    sys.exit(1)

gbk_file = args.gbk
smurf_file = args.smurf
fasta_file = args.prot
contig_file = args.dna

smurf_file_handle = open(smurf_file, "w")
fasta_file_handle = open(fasta_file, "w")
contig_file_handle = open(contig_file, "w")

for record in SeqIO.parse(gbk_file, "genbank"):
    if args.ncbi:
        scaffold_name = record.name[6:]
    else:
        scaffold_name = re.sub('[^0-9]','', record.name)
    contig_file_handle.write(">%s\n%s\n" % (
        scaffold_name.lstrip("0"), 
        record.seq))
    for f in record.features:
        if args.ncbi:
            name = record.name[6:]
        else:
            name = re.sub('[^0-9]','', record.name)
        if f.type == "CDS":
            if args.jgi:
                fasta_file_handle.write(">%s\n%s\n" % (
                   f.qualifiers['gene'][0],
                   f.qualifiers['translation'][0]))
                locus_tag = f.qualifiers.get("gene", ["No ID"])[0]
                product_name = f.qualifiers.get("product", ["No Description"])[0]
                mystart = f.location.start
                myend = f.location.end
                strand = f.location.strand
                if strand == 1:
                    smurf_file_handle.write("%s\t%s\t%s\t%s\t%s\n" % (locus_tag, name.lstrip("0"), int(mystart), int(myend), product_name))
                else:
                    smurf_file_handle.write("%s\t%s\t%s\t%s\t%s\n" % (locus_tag, name.lstrip("0"), int(myend), int(mystart), product_name))
            else:
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
