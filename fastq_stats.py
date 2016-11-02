#!/usr/bin/env python

import sys
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('fastq', help='fastq file')
args = parser.parse_args()

lengths = []
count = 0
with open(args.fastq, 'rU') as fastqin:
    for title, sequence, qual in FastqGeneralIterator(fastqin):
        lengths.append(len(sequence))
        count += 1

#convert to numpy array
a = np.array(lengths)
minLen = np.amin(a)
maxLen = np.amax(a)
avgLen = np.mean(a)
totalLen = np.sum(a)

#output stats
print "Reads:    "+'{0:,}'.format(count)
print "AvgLen:   "+'{0:,}'.format(int(avgLen))+ ' bp'
print "Shortest: "+'{0:,}'.format(minLen)+ ' bp'
print "Longest:  "+'{0:,}'.format(maxLen)+ ' bp'
print "Total:    "+'{0:,}'.format(totalLen)+ ' bp'