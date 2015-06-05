#!/usr/bin/env python

import sys
from Bio import SeqIO
import argparse
 
parser = argparse.ArgumentParser()
parser.add_argument('gbk', help='Genbank file')
parser.add_argument('fasta', help='fasta file')
args = parser.parse_args()

gbk_filename = args.gbk
faa_filename = args.fasta

input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")
ending = "-T1"
for seq_record in SeqIO.parse(input_handle, "genbank") :
    print "Dealing with GenBank record %s" % seq_record.id
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            assert len(seq_feature.qualifiers['translation'])==1
            output_handle.write(">%s%s\n%s\n" % (
                   seq_feature.qualifiers['locus_tag'][0],
                   ending,
                   seq_feature.qualifiers['translation'][0]))

output_handle.close()
input_handle.close()
print "Done"
