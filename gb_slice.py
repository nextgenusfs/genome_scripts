#!/usr/bin/env python

import sys
from Bio import SeqIO

if len(sys.argv) < 2:
    print "Usage: %s File.gbk Rec_ID start end Output.gbk <reverse_comp>" % sys.argv[0]
    sys.exit(1)


input = open(sys.argv[1], 'rb')
output = open(sys.argv[5], 'wb')
find_rec = sys.argv[2]
SeqRecords = SeqIO.parse(input, 'genbank')
start = int(sys.argv[3])
end = int(sys.argv[4])
for record in SeqRecords:
    if record.id == find_rec:
        print record.id
        sub_record = record[start:end]
        try:
            if sys.argv[6] == 'reverse_comp':
                print "reverse-complementing record"
                sub_record = sub_record.reverse_complement(id=record.id)
        except IndexError:
            sub_record = sub_record
        SeqIO.write(sub_record, output, 'genbank')
input.close()
output.close()