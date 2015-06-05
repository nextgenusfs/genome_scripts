#!/usr/bin/env python

#This script converts genbank format into GFF Format
#requires BCBio gff parser, install: pip install bcbio-gff
#script keeps stupid genbank ridiculously long chromosome/scaffold names

import sys
from Bio import SeqIO
from BCBio import GFF
import argparse

parser=argparse.ArgumentParser(
    description='''Script that converts genbank to GFF3 format ''',
    epilog="""Written by Jon Palmer (2015)  palmer.jona@gmail.com""")
parser.add_argument('gbk', help='Genbank file')
parser.add_argument('-g','--gff', default="gb2gff.gff3", help='GFF3 output')
args=parser.parse_args()

if len(sys.argv) < 2:
    parser.print_usage()
    sys.exit(1)
 
in_file = args.gbk
out_file = args.gff
in_handle = open(in_file)
out_handle = open(out_file, "w")
 
GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)
 
in_handle.close()
out_handle.close()
