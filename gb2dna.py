#!/usr/bin/env python

import sys
from Bio import SeqIO
import argparse
 
parser = argparse.ArgumentParser()
parser.add_argument('gbk', help='Genbank file')
args = parser.parse_args()

gbk_filename = args.gbk

input_handle  = open(gbk_filename, "rU")

for seq_record in SeqIO.parse(input_handle, "genbank"):
    sys.stdout.write(">%s\n%s\n" % (seq_record.id, seq_record.seq))

input_handle.close()

