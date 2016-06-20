#!/usr/bin/env python

import sys
from Bio import SeqIO
import argparse
 
parser = argparse.ArgumentParser()
parser.add_argument('--id', default='locus_tag', choices=['locus_tag', 'protein_id'], help='Genbank file')
parser.add_argument('-i','--input', required=True, help='Genbank file')
parser.add_argument('-o', '--output', required=True, help='fasta file')
args = parser.parse_args()

input_handle  = open(args.input, "r")
output_handle = open(args.output, "w")
for seq_record in SeqIO.parse(input_handle, "genbank") :
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            try:
                locusID = seq_feature.qualifiers['locus_tag'][0]
            except KeyError:
                pass            
            try:
                protID = seq_feature.qualifiers['protein_id'][0]
            except KeyError:
                pass
            try:
                Seq = seq_feature.qualifiers['translation'][0]
            except KeyError:
                pass
            
            if Seq:
                if args.id == 'locus_tag':
                    if locusID:
                        ID = locusID
                    else:
                        if protID:
                            ID = protID
                        else:
                            continue                   
                elif args.id == 'protein_id':
                    if protID:
                        ID = protID
                    else:
                        if locusID:
                            ID = locusID
                        else:
                            continue                        
                output_handle.write(">%s\n%s\n" % (ID, Seq))
            else:   
                print "%s has no sequence" % ID

output_handle.close()
input_handle.close()
print "Proteins extracted from " + args.input + " saved to " + args.output
