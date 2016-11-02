#!/usr/bin/env python

#This script sorts fasta file by header alphabetically

import sys, fileinput, os
from Bio import SeqIO
from natsort import natsorted

#hack to accept stdin, write to temp file, then parse to dictionary, delete original
with open('sort_fasta_tmp.fa', 'w') as output:
    for line in fileinput.input():
        output.write(line)

with open('sort_fasta_tmp.fa', 'rU') as input:
    seqs_dict = SeqIO.to_dict(SeqIO.parse(input, 'fasta'))

sorted_list = natsorted(seqs_dict)

for i in sorted_list:
    SeqIO.write(seqs_dict[i], sys.stdout, 'fasta')

os.remove('sort_fasta_tmp.fa')


