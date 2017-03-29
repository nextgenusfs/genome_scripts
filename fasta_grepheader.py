#!/usr/bin/env python

#This script removes fasta files from list in file

import sys
from Bio import SeqIO

if len(sys.argv) < 3:
    print "Usage: fasta_grepheader.py input.fasta header"
    sys.exit(1)

with open(sys.argv[1], 'rU') as input:
	for rec in SeqIO.parse(input, 'fasta'):
		if sys.argv[2] in rec.id:
			SeqIO.write(rec, sys.stdout, 'fasta')



