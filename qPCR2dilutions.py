#!/usr/bin/env python

import sys, itertools, math
from natsort import natsorted

if len(sys.argv) < 2:
    print("Usage: qPCR2dilutions.py StepOnePlusout.txt length qPCRdilution OutputConcentration(pM)\nExample: qPCR2dilutions.py run1.txt 400 6000 1000")
    sys.exit(1)

InFile = sys.argv[1]
Length = int(sys.argv[2])
Dilution = int(sys.argv[3])
outputconc = int(sys.argv[4])

#1 pM = 6.022e5 how about standardizing to 1000 pmol or 6.022e8 copies/ul; c1v1=c2v2
conc = 6.022e5 * outputconc

factor = 183 / float(Length)

results = {'#Sample_name' : ('Copies/ul', 'volume', 'water', 'total', 'note')}
with open(InFile, 'rU') as input:
    for line in itertools.islice(input, 9, None):
        cols = line.split('\t')
        if cols[1]:
            well = cols[0]
            name = cols[1].replace(' ', '_')
            copies = cols[9]
            orig_copies = float(copies) * Dilution
            adjusted = float(orig_copies) * factor

            volume = 5
            if adjusted < conc:
                dil = volume
                water = 0
                note = '[amplicons] is < than %i pM' % outputconc
            else:
                numerator = conc * volume
                dil = numerator / adjusted
                if dil < 2:
                    adj = 2 / dil
                    adj = math.ceil(adj)
                    volume = volume * adj
                else:
                    adj = 1
                dil = round(dil * adj, 1)
                water = volume - dil

                note = 'none'
            results[name] = ("{:.3e}".format(adjusted), dil, water, int(volume), note)

print("#Calcuations based on %i pM for each sample, then combine equally for OT2 reaction, i.e. 5 ul from each" % outputconc) 
print("#Then dilute to 100 pM, then mix 6.25 ul of 100 pM amplicons with 18.75 ul water for final 26 pM at 25 ul")
for k,v in natsorted(results.items()):
    print "%s\t%s\t%s\t%s\t%s\t%s" % (k,v[0],v[1],v[2],v[3],v[4])