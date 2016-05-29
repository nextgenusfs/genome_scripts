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
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            try:
                ID = seq_feature.qualifiers['locus_tag'][0]
            except KeyError:
                try:
                    ID = seq_feature.qualifiers['protein_id'][0]
                except KeyError:
                    pass
            try:
                Seq = seq_feature.qualifiers['translation'][0]
            except KeyError:
                pass
            
            if ID:
                if Seq:
                    output_handle.write(">%s\n%s\n" % (ID, Seq))
                else:
                    print "%s has no sequence" % ID

output_handle.close()
input_handle.close()
print "Proteins extracted from " + gbk_filename + " saved to " + faa_filename
