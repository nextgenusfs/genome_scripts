#!/usr/bin/env python

import sys, os, argparse, subprocess, shutil
from Bio import SeqIO
from natsort import natsorted

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='cluster2utax.py',
    description='''Take output from bold2utax and cluster to 97 pct, then use LCA to rename taxonomy.''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='FASTA file from bold2utax.py')
parser.add_argument('-o','--out', required=True, help='Output basename')
parser.add_argument('--min_tax', default='k', choices=['k','p','c','o','f','g','s'], help='Minimum levels of taxonomy to keep seq')
parser.add_argument('--remove_dups', action='store_true', help='remove duplicate species')
parser.add_argument('--debug', action='store_true', help='Keep intermediate files')
args=parser.parse_args()

def countfasta(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if line.startswith (">"):
                count += 1
    return count

TxLv = {'k': 0, 'p': 1, 'c': 2, 'o': 3, 'f': 4, 'g': 5, 's': 6}
Lv = TxLv.get(args.min_tax)
#set up tmpdir
tmpdir = 'tmp_'+str(os.getpid())
os.makedirs(tmpdir)
print "----------------------------------"
#run dereplicatoin using vsearch
print "Dereplication"
print "----------------------------------"
Derep = os.path.join(tmpdir, 'derep.fa')
subprocess.call(['vsearch', '--derep_fulllength', args.input, '--sizeout', '--output', Derep, '--notrunclabels'])
print "----------------------------------"

#now cluster at 97%
print "Clustering"
print "----------------------------------"
OTUs = os.path.join(tmpdir, 'otus.fa')
UC = os.path.join(tmpdir, 'otus.uc')
subprocess.call(['vsearch', '--cluster_size', Derep, '--id', '0.97', '--centroids', OTUs, '--uc', UC, '--notrunclabels', '--xsize'])
print "----------------------------------"
#load uclust into dictionary for LCA taxonomy adjustments
print 'Determining Last Common Ancestor (LCA) for clustered seqs'
LCA = {}
with open(UC, 'rU') as input:
    for line in input:
        line = line.replace('\n', '')
        cols = line.split('\t')
        if cols[3] != '*':
            centroid = cols[9].rsplit(';', 2)[0]
            hit = cols[8].rsplit(';', 2)[0]
            centroidTax = centroid.split('tax=')[-1]
            centroidID = centroid.split('tax=')[0].rstrip()
            hitTax = hit.split('tax=')[-1]
            centroidLCA = centroidTax.split(',')
            hitLCA = hitTax.split(',')
            if not centroidLCA == hitLCA:
                if len(centroidLCA) == len(hitLCA):
                    #tax levels are the same, so now get LCA
                    for x in range(0,len(centroidLCA)):
                        if not centroidLCA[x] == hitLCA[x]:
                            Tax = centroidLCA[:x]
                else:
                    #if not the same length, then determine if they are identical to same level of taxonomy
                    if len(centroidLCA) > len(hitLCA):
                        for x in range(0,len(hitLCA)):
                            if not centroidLCA[x] == hitLCA[x]:
                                Tax = centroidLCA[:x]
                            else:
                                Tax = centroidLCA
                    else:
                        for x in range(0,len(centroidLCA)):
                            if not hitLCA[x] == centroidLCA[x]:
                                Tax = hitLCA[:x]
                            else:
                                Tax = hitLCA
                #generate name change dictionary called LCA
                LCA[centroid] = '%stax=%s' % (centroidID, ','.join(Tax))
print "----------------------------------"
#ok, now load in the centroid fasta file, rename if needed, and filter based on duplicates
#they have been clustered, so sequences with the same taxonomy down to species level need to be resolved
#also for utax training, likely not worth keeping seqs only classified to class
print "Renaming centroids and filtering data"
TaxSeen = []
UtaxOut = args.out+'.utax.fa'
UsearchOut = args.out+'.usearch.fa'
total = 0
okayTax = 0
duplicateSp = 0
with open(UtaxOut, 'w') as utax:
    with open(UsearchOut, 'w') as usearch:
        with open(OTUs, 'rU') as input:
            for rec in SeqIO.parse(input, 'fasta'):
                if not 'tax=' in rec.description:
                    continue
                if 'Pending' in rec.description:
                    continue
                #entire ID is in rec.description
                if rec.description in LCA:
                    ID = LCA.get(rec.description)
                else:
                    ID = rec.description
                #for usearch out, we want all files, so output here
                rec.id = ID
                rec.name = ''
                rec.description = ''
                SeqIO.write(rec, usearch, 'fasta')
                total += 1

                tax = ID.split('tax=')[-1]
                taxLevels = tax.count(',') #0 = kingdom, 1 = phylum, 2 = class, 3 = order, 4 = family, 5 = genus, 6 = species
                if taxLevels >= Lv:
                    okayTax += 1
                    if taxLevels > 5:
                        if not tax in TaxSeen:
                            TaxSeen.append(tax)
                        else:
                            if args.remove_dups:
                                duplicateSp += 1
                                continue  #skip if duplicated, just taking first one, hopefully this is okay...built as option.
                    SeqIO.write(rec, utax, 'fasta')
print "----------------------------------"
print '%i input Seqs' % countfasta(args.input)
print '%i Seqs after 97%% clustering' % total
print '%i Seqs saved to %s' % (total, UsearchOut)
print '%i Seqs removed low taxonomy' % (total - okayTax)
print '%i Seqs removed duplicate species' % (duplicateSp)
print '%i Seqs saved to %s' % (countfasta(UtaxOut), UtaxOut)
print "----------------------------------"
if not args.debug:
    shutil.rmtree(tmpdir)
