#!/usr/bin/env python

import sys, csv, collections

def flatten(l):
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el

#load in iprscan results, eggnot results, blast results, proteinortho and add annotation to proteinortho data

if len(sys.argv) < 4:
    print "Usage: %s proteinortho eggnog blast iprscan"
    sys.exit(1)

protortho = sys.argv[1]
eggnog = sys.argv[2]
blast = sys.argv[3]
iprscan = sys.argv[4]


EggNog = {}
Blast = {}
IprScan = {}


with open(eggnog, 'rU') as eggnog:
    reader = csv.reader(eggnog, delimiter='\t')
    for line in reader:
        EggNog[line[0]] = line[1]
#print EggNog

with open(blast, 'rU') as blast:
    reader = csv.reader(blast, delimiter='\t')
    for line in reader:
        Blast[line[0]] = line[2]
#print Blast

with open(iprscan, 'rU') as iprscan:
    reader = csv.reader(iprscan, delimiter='\t')
    for line in reader:
        if line[0] in IprScan:
            IprScan[line[0]].append(line[1])
        else:
            IprScan[line[0]] = [line[1]]
#print IprScan

with open(protortho, 'rU') as protortho:
    reader = csv.reader(protortho, delimiter='\t')
    for line in reader:
        if line[0].startswith("#"):
            line.append('nr_Blast')
            line.append('EggNog4.1')
            line.append('IprScan')
            print "\t".join(line)+'\n'
            continue
        #now deal with blast hits, then EggNog, then IprScan
        bhit = []
        ehit = []
        ihit = []
        for i in line[3:]:
            if i != "*":
                hit1 = Blast.get(i)
                hit2 = EggNog.get(i)
                hit3 = IprScan.get(i)
                if not hit1 in bhit:
                    bhit.append(hit1)
                if not hit2 in ehit:
                    ehit.append(hit2)
                if not hit3 in ihit:
                    ihit.append(hit3)
        ihit = flatten(ihit)
        bhit = list(filter(lambda x: x!= None, bhit))
        ehit = list(filter(lambda x: x!= None, ehit))
        ihit = filter(lambda x: x!= None, ihit)
        fhit1 = '; '.join(bhit)
        fhit2 = '; '.join(ehit)
        fhit3 = '; '.join(ihit)
        line.append(fhit1)
        line.append(fhit2)
        line.append(fhit3)
        print "\t".join(line)+'\n'


