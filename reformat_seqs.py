#!/usr/bin/env python

from Bio import SeqIO
import sys

counter = 1
master = 1
with open(sys.argv[3], 'w') as out2:
	with open(sys.argv[2], 'w') as output:
		with open(sys.argv[1], 'rU') as input:
			for rec in SeqIO.parse(input, "fasta"):
				if 'scaffold' in rec.id:
					rec.id.replace('scaffold', 'unplaced')
					rec.name = ''
					rec.description = ''
				SeqIO.write(rec, output, 'fasta')
				out2.write('%s\t%i\n' % (rec.id, master))
				master += 1