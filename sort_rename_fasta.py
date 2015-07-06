#!/usr/bin/env python

#This script sorts fasta file and then renames headers
#written by Jon Palmer palmer.jona at gmail dot com

import sys
import os
from Bio import SeqIO
import argparse

parser=argparse.ArgumentParser(
    description='''Script that sorts fasta be length and renames headers ''',
    epilog="""Jon Palmer (2015)  palmer.jona@gmail.com""")
parser.add_argument('fasta', help='Fasta file (Required)')
parser.add_argument('-o','--out', default="sort.renamed.fasta", help='Fasta output')
parser.add_argument('-b','--base', default="Supercontig", help='Base name for fasta headers')
parser.add_argument('-l','--length', default="1", help='min length to keep seq')
args=parser.parse_args()

if len(sys.argv) < 2:
    parser.print_usage()
    sys.exit(1)

fasta_in = args.fasta
base_name = args.base
length = int(args.length)

#get rid of sequences shorter than min length
seq_pass = []
for record in SeqIO.parse(fasta_in, "fasta"):
    if len(record.seq) > length :
        seq_pass.append(record)
print "Found %i sequences" % len(seq_pass)
tempout = open("fasta.tmp", "w")
SeqIO.write(seq_pass, tempout, "fasta")


#sort records and write temp file
records = list(SeqIO.parse("fasta.tmp", "fasta"))
records.sort(cmp=lambda x,y: cmp(len(y),len(x)))
SeqIO.write(records, "fasta.tmp2", "fasta")

#load temp fasta back in and change names, prob unnecessary but couldn't figure out how to do with SeqIO
fasta_tmp = open("fasta.tmp2")
fasta_out = open(args.out, 'wb')
count = 1
for line in fasta_tmp.readlines():
	if line.startswith('>'):
		contig_id = '>' + base_name + '_' + str(count) + '\n'
		fasta_out.write(contig_id)
		count += 1
	else:
		fasta_out.write(line)
tempout.close()
fasta_out.close()
fasta_tmp.close()
os.remove("fasta.tmp")
os.remove("fasta.tmp2")
