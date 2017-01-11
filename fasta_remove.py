#!/usr/bin/env python

#This script removes fasta files from list in file

import sys
from Bio import SeqIO

if len(sys.argv) < 3:
    print "Usage: fasta_remove.py input.fasta list2remove.txt"
    sys.exit(1)

#get list of names from file
remove = []
with open(sys.argv[2], 'rU') as filein:
	for line in filein:
		line = line.replace('\n', '').rstrip()
		if not line in remove:
			remove.append(line)

with open(sys.argv[1], 'rU') as input:
	for rec in SeqIO.parse(input, 'fasta'):
		if not rec.id in remove:
			SeqIO.write(rec, sys.stdout, 'fasta')



