#!/usr/bin/env python

#This script converts genbank format into tab delimited SMURF format

import sys
from Bio import SeqIO

gbk_filename = sys.argv[1]

for record in SeqIO.parse(gbk_filename, "genbank"):
    for f in record.features:
        name = record.name
        if f.type == "CDS":
            locus_tag = f.qualifiers.get("locus_tag", ["???"])[0]
            product_name = f.qualifiers.get("product", ["???"])[0]
            mystart = f.location.start
            myend = f.location.end
            strand = f.location.strand
            if strand == 1:
                print("%s\t%s\t%s\t%s\t%s" % (locus_tag, name, mystart, myend, product_name))
            else:
                print("%s\t%s\t%s\t%s\t%s" % (locus_tag, name, myend, mystart, product_name))
