#!/usr/bin/env python

import sys
from Bio import SeqIO
import argparse
 
parser = argparse.ArgumentParser()
parser.add_argument('gbk', help='Genbank file')
parser.add_argument('fasta', help='fasta file')
parser.add_argument('--maker', action="store_true", help='Maker -T1 ending to names')
parser.add_argument('--jgi', action="store_true", help='JGI use name instead of locus tag')
args = parser.parse_args()

gbk_filename = args.gbk
faa_filename = args.fasta

input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")
ending = "-T1"
for seq_record in SeqIO.parse(input_handle, "genbank") :
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            assert len(seq_feature.qualifiers['translation'])==1
            if args.maker:
                output_handle.write(">%s%s\n%s\n" % (
                    seq_feature.qualifiers['locus_tag'][0],
                    ending,
                    seq_feature.qualifiers['translation'][0]))
            else:
                if args.jgi:
                    output_handle.write(">%s\n%s\n" % (
                        seq_feature.qualifiers['gene'][0],
                        seq_feature.qualifiers['translation'][0]))
                else:
                    output_handle.write(">%s\n%s\n" % (
                        seq_feature.qualifiers['locus_tag'][0],
                        seq_feature.qualifiers['translation'][0]))

output_handle.close()
input_handle.close()
print "Proteins extracted from " + gbk_filename + " saved to " + faa_filename
