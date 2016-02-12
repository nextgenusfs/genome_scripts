#!/usr/bin/env python

from Bio import SeqIO
import sys

if len(sys.argv) < 2:
    print("Usage: gb_remove_records.py genbank.gbk removal_list.txt > output.gbk")
    sys.exit(1)

#load in list to set for removal
removal = []
with open(sys.argv[2], 'rU') as remove:
    for line in remove:
        line = line.replace('\n', '')
        removal.append(line)

removal = set(removal)
sys.stderr.write("Removing %i contigs from genbank file...\n" % len(removal))

Count = 0
with open('filtered_out_proteins.fa', 'w') as output:
    with open(sys.argv[1], 'rU') as input:
        SeqRecords = SeqIO.parse(input, 'genbank')
        for record in SeqRecords:
            if record.id in removal:
                for f in record.features:
                    if f.type == 'CDS':
                        Count += 1
                        output.write(">%s|%s\n%s\n" % (f.qualifiers['locus_tag'][0], record.id, f.qualifiers['translation'][0]))
            else:
                SeqIO.write(record, sys.stdout, 'genbank')
sys.stderr.write("%i protein models filtered out, written to filtered_out_proteins.fa\n" % Count)