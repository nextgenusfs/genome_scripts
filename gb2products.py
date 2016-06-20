#!/usr/bin/env python

import sys
from Bio import SeqIO
import argparse
 
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', required=True, help='Genbank file')
args = parser.parse_args()

with open(args.input, 'rU') as gbk:
    for record in SeqIO.parse(gbk, 'genbank'):
        for seq_feature in record.features:
            if seq_feature.type == "CDS":
                ID = seq_feature.qualifiers['locus_tag'][0]
                protID = seq_feature.qualifiers['protein_id'][0]
                product = seq_feature.qualifiers['product'][0]
            
                sys.stdout.write('%s\t%s\t%s\n' % (ID, protID, product))