#!/usr/bin/env python

#convert snippy VCF to 100 bp interval bed file

import os, sys

table = []
with open(sys.argv[1], 'rU') as input:
    for line in input:
        line = line.replace('\n', '')
        if line.startswith('#'):
            continue
        cols = line.split('\t')
        scaffold = cols[0]
        start = int(cols[1]) - 50
        end = int(cols[1]) + 50
        sys.stdout.write('%s\t%i\t%i\n' % (scaffold, start, end))

