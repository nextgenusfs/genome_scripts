#!/usr/bin/env python

#This script concatenates input fasta file a

import sys, os
from Bio import SeqIO

if len(sys.argv) < 3:
    print "Usage: cat_fasta.py input.fasta header_name"
    os._exit(1)

with open(sys.argv[1], 'rU') as input:
    sys.stdout.write('>%s\n' % sys.argv[2])
    for record in SeqIO.parse(input, 'fasta'):
        sys.stdout.write('%s' % str(record.seq))
    sys.stdout.write('\n')


