#!/usr/bin/env python

import sys, argparse
from Bio import SeqIO

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='gbslice.py', usage="%(prog)s [options] -i genbank -o sliced_genbank",
    description='''Various ways to slice a GenBank file''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='GenBank genome file')
parser.add_argument('-o','--out', required=True, help='Name of output GenBank file')
parser.add_argument('-c','--contig', help='Name of Contig or RecID to grab')
parser.add_argument('-l','--location', nargs='+', help='Coordinates of slice; start end')
parser.add_argument('-r','--reverse_comp', action='store_true', help='reverse complement file')
args=parser.parse_args()


SeqRecords = SeqIO.parse(args.input, 'genbank')
with open(args.out, 'w') as output:
	for record in SeqRecords:
		if args.contig:
			if not args.contig == record.id:
				continue
			print "Found %s contig" % args.contig
			if args.location:
				if len(args.location) == 2:
					start = int(args.location[0])
					end = int(args.location[1])
				else: #throw error
					print "You did not specify two locations for slicing, i.e. -l 100 3000"
					sys.exit(1)
				if end > len(record.seq):
					print "%i end location is out of bounds: %s is %i bp" % (end, record.id, len(record.seq))
					end = len(record.seq)
				sub_record = record[start:end]
			else:
				sub_record = record
			if args.reverse_comp:
				sub_record = sub_record.reverse_complement(id=record.id)
			SeqIO.write(sub_record, output, 'genbank')
		else:
			if args.reverse_comp:
				sub_record = record.reverse_complement(id=record.id)
				SeqIO.write(sub_record, output, 'genbank')
			else:
				print "nothing to do, exiting"
				sys.exit(1)
				