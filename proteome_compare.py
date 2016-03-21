#!/usr/bin/env python

import sys, subprocess
from natsort import natsorted

#script to run global alignment using usearch8 and then parse the output

if len(sys.argv) < 2:
    print "Usage: proteome_compare.py query.faa reference.faa singletons.output.txt okay.hits.txt"
    sys.exit(1)

def compare(input):
    global singletons
    global good
    singletons = []
    good = []
    cantsithereseatstaken = []
    with open(input, 'rU') as input:
        identical = 0; great = 0; okay = 0; partial = 0; unique = 0; others = 0; bad = 0;
        for line in input:
            col = line.split('\t')
            if col[1] == '*':
                singletons.append([col[0]])
                unique += 1
                continue
            if float(col[2]) == 100.0 and int(col[3]) == int(col[4]) and int(col[5]) == int(col[6]):
                identical += 1
                if not col[1] in cantsithereseatstaken:
                    good.append([col[0],col[1]])
                    cantsithereseatstaken.append(col[1])
                continue
            elif float(col[2]) >= 95.0 and int(col[5]) >= 95 and int(col[6]) >= 95:
                great += 1
                if not col[1] in cantsithereseatstaken:
                    good.append([col[0],col[1]])
                    cantsithereseatstaken.append(col[1])
                continue
            elif float(col[2]) >= 80.0 and int(col[5]) >= 80 and int(col[6]) >= 80:
                okay += 1
                if not col[1] in cantsithereseatstaken:
                    good.append([col[0],col[1]])
                    cantsithereseatstaken.append(col[1])
                continue
            elif float(col[2]) >= 60.0 and int(col[5]) >= 60 and int(col[6]) >= 60:
                others +=1
                if not col[1] in cantsithereseatstaken:
                    good.append([col[0],col[1]])
                    cantsithereseatstaken.append(col[1])
                continue
            else:
                bad += 1
                singletons.append([col[0]])
    return (identical, great, okay, others, unique, bad)

subprocess.call(['usearch8', '-usearch_global', sys.argv[1], '-db', sys.argv[2], '-id', '0.6', '-userout', 'prot_compare1.txt', '-userfields', 'query+target+id+ql+tl+qcov+tcov', '-output_no_hits', '-top_hit_only', '-notrunclabels'])

subprocess.call(['usearch8', '-db', sys.argv[1], '-usearch_global', sys.argv[2], '-id', '0.6', '-userout', 'prot_compare2.txt', '-userfields', 'query+target+id+ql+tl+qcov+tcov', '-output_no_hits', '-top_hit_only', '-notrunclabels'])

#parse output
qvsr = compare('prot_compare1.txt')

with open(sys.argv[3], 'w') as output:
    for hit in singletons:
        output.write("%s\n" % (hit[0]))

with open(sys.argv[4], 'w') as output:
    for hit in good:
        output.write("%s\t%s\n" % (hit[0], hit[1]))

rvsq = compare('prot_compare2.txt')

print"--------------------------------"
print"Query vs Ref:"
print"--------------------------------"
print"Identical (100%% id,cov): %i" % (qvsr[0])
print"Great (>95%% id,cov):     %i" % (qvsr[1])
print"Good (>80%% id,cov):      %i" % (qvsr[2])
print"Others (>60%% id,cov):     %i" % (qvsr[3])
print"Bad (>60%%, cov < 60%%):    %i" % (qvsr[5])
print"Unique (to the query):    %i" % (qvsr[4])
print"Missing (unique to ref):  %i" % (rvsq[4])
print"--------------------------------"
