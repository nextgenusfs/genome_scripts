#!/usr/bin/env python

#script to process RNAseq output from kallisto/sleuth and funannotate annotation table

import os, sys, argparse
from natsort import natsorted

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='sleuth2table.py',
    description='''Script that converts sleuth tsv file, funnotate annotation, and sleuth TSP file.''',
    epilog="""Written by Jon Palmer (2018) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='Sleuth output')
parser.add_argument('-a','--annotation', help='Funannotate annotation.txt (tsv file)')
parser.add_argument('-t','--tsp', required=True, help='Kallisto/Sleuth tsp expression values for each sample')
args=parser.parse_args()

#functions from https://stackoverflow.com/questions/15389768/standard-deviation-of-a-list
def mean(data):
    """Return the sample arithmetic mean of data."""
    n = len(data)
    if n < 1:
        raise ValueError('mean requires at least one data point')
    return sum(data)/float(n) # in Python 2 use sum(data)/float(n)

def _ss(data):
    """Return sum of square deviations of sequence data."""
    c = mean(data)
    ss = sum((x-c)**2 for x in data)
    return ss

def stddev(data, ddof=0):
    """Calculates the population standard deviation
    by default; specify ddof=1 to compute the sample
    standard deviation."""
    n = len(data)
    if n < 2:
        raise ValueError('variance requires at least two data points')
    ss = _ss(data)
    pvar = ss/(n-ddof)
    return pvar**0.5

#load the output of kallisto_table to append to sleuth tsv
'''
geneID : {condition: [(sample1, est_counts, tmp, eff_len, length), (sample2, est_counts, tmp, eff_len, length)]
		  condition: [(sample1, est_counts, tmp, eff_len, length), (sample2, est_counts, tmp, eff_len, length)]}
'''
ExpValues = {}
with open(args.tsp, 'rU') as infile:
	for line in infile:
		line = line.rstrip()
		if ',' in line:
			cols = line.split(',')
			cols = [x.replace('"', '') for x in cols]
			if cols[1].startswith('target_id'):
				continue
			num,id,sample,est_counts,tmp,eff_len,length,condition = cols
			if not id in ExpValues:
				ExpValues[id] = {condition: [(sample, est_counts, tmp, eff_len, length)]}
			else:
				if condition in ExpValues[id]:
					ExpValues[id][condition].append((sample, est_counts, tmp, eff_len, length))
				else:
					ExpValues[id][condition] = [(sample, est_counts, tmp, eff_len, length)]

#loop through and calcuate average, stdev of each condition
ExpFinal = {}
for k,v in natsorted(ExpValues.items()):
	for cond,sample in natsorted(v.items()):
		AllC = []
		AllT = []
		name = []
		effLen = sample[0][3]
		length = sample[0][4]
		for s in sample:
			AllC.append(float(s[1]))
			AllT.append(float(s[2]))
			name.append(s[0])
		AvgCount = (mean(AllC), stddev(AllC))
		AvgTmp = (mean(AllT), stddev(AllT))
		if not k in ExpFinal:
			ExpFinal[k] = [{'condition': cond, 'name': name, 'counts': AllC, 'tmp': AllT, 'avgcounts': AvgCount, 'avgtmp': AvgTmp, 'length': length, 'efflen': effLen}]
		else:
			ExpFinal[k].append({'condition': cond, 'name': name, 'counts': AllC, 'tmp': AllT, 'avgcounts': AvgCount, 'avgtmp': AvgTmp, 'length': length, 'efflen': effLen})

for key,v in natsorted(ExpFinal.items()):
	outlist = []
	newlist = sorted(v, key=lambda k: k['condition'])
	for x in newlist:
		outlist.append()
#load in funannotate annotations to append to sleuth data
'''
#dumbass should unify this in funannotate....
#header from compare looks like
GeneID	scaffold:start-end	strand	length	description	Ortho Group	EggNog	BUSCO	Secreted	TransMembrane	Protease family	CAZyme family	Transcription factor	InterPro Domains	PFAM Domains	GO terms	SecMet Cluster	SMCOG

#header from annotate looks like
GeneID	Feature	Contig	Start	Stop	Strand	Name	Product	BUSCO	PFAM	InterPro	EggNog	COG	GO Terms	Secreted	Membrane	Protease	CAZyme	Notes	Translation
'''
Annotations = {}
if args.annotation:
	with open(args.annotation, 'rU') as annot:
		for line in annot:
			ID,feature,start,strand,length,product,ortho,eggnog,busco,secreted,membrane,protease,cazyme,tf,ipr,pfam,go,secmet,smcog,cog,notes,protein = (None,)*22
			line = line.rstrip()
			if 'Feature' in line:
				ID,feature,contig,start,stop,strand,name,product,busco,pfam,ipr,eggnog,cog,go,secreted,membrane,protease,cazyme,notes,protein = line.split('\t')
			else:
				ID,start,strand,length,product,ortho,eggnog,busco,secreted,membrane,protease,cazyme,tf,ipr,pfam,go,secmet,smcog = line.split('\t')
			if ID == 'GeneID': #header
				continue
			if not ID in Annotations:
				Annotations[ID] = {'feature': feature, 'contig': contig, 'product': product, 'busco': busco,
								   'name': name, 'pfam': pfam, 'iprscan': ipr, 'eggnog': eggnog, 'cog': cog,
								   'secreted': secreted, 'membrane': membrane, 'protease': protease, 'cazyme': cazyme, 
								   'note': notes, 'secmet': secmet, 'smcog': smcog, 'tf':tf, 'ortholog': ortho}
			
#okay, now lets parse sleuth table, append tpm values, and append annotation if exists
		
			