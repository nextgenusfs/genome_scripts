#!/usr/bin/env python

import sys, argparse
import pandas as pd
from Bio import SeqIO

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='gb2tab.py', usage="%(prog)s [options] -i genbank -o tsv_gene_output",
    description='''Genbank to tsv file for gene models''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='GenBank genome file')
parser.add_argument('-o','--out', required=True, help='Name of output basename file')
parser.add_argument('-p','--proteinortho', help='proteinOrtho POFF file to annotate orthologs')
args=parser.parse_args()

def remove_values_from_list(the_list, val):
   return [value for value in the_list if value != val]

def protortho2dict(poff, basename):
	#function takes a proteinortho poff file and a locus_tab basename to generate a dictionary of orthologs
	#poff files look likethis
	#load into pandas dataframe
	poffdf = pd.read_csv(poff, delimiter='\t', header=0)
	numspecies = len(poffdf.columns) - 3
	slice = poffdf[(poffdf["# Species"] == numspecies) & (poffdf["Genes"] == numspecies)].head(1)
	column = ''
	for col in slice:
		for x in slice[col]:
			if basename in str(x):
				column = col
	df = poffdf.drop(["# Species", "Genes", "Alg.-Conn."], axis=1)
	df.set_index(column, inplace=True)
	poffdict = {}
	for index, row in df.iterrows():
		values = []
		for x in range(0, len(row)):
			if row[x] != '*':
				hits = row[x].split(',')
				for y in hits:
					values.append(y)
		if ',' in index:
			keys = index.split(',')
			for i in range(0, len(keys)):
				others = keys[:i] + keys[i+1 :]
				values = values + others
				poffdict[keys[i]] = values
		else:
			poffdict[index] = values
	return poffdict

def getGBbaseName(input):
	basename = ''
	with open(input, 'rU') as gbk:
		for record in SeqIO.parse(gbk, 'genbank'):
			for f in record.features:
				if f.type == 'CDS':
					ID = f.qualifiers['locus_tag'][0]
					basename = ID.split('_')[0]
					break
	return basename
	
def gb2tab(input, Proteins, Transcripts, tabbedout, lookup):
    #this will not output any UTRs for gene models, don't think this is a problem right now....
    with open(tabbedout, 'w') as tabout:
        tabout.write("GeneID\tType\tScaffold\tStart\tEnd\tStrand\tProduct\tOrthologs\n")
        with open(Transcripts, 'w') as tranout:
			with open(Proteins, 'w') as proteins:
				with open(input, 'rU') as gbk:
					for record in SeqIO.parse(gbk, 'genbank'):
						for f in record.features:
							if f.type == "mRNA":
								feature_seq = f.extract(record.seq)
								tranout.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], feature_seq))
							if f.type == 'CDS':
								proteins.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], f.qualifiers['translation'][0]))
								chr = record.id
								ID = f.qualifiers['locus_tag'][0]
								try:
									product = f.qualifiers['product'][0]
								except KeyError:
									product = "hypothetical protein"
								start = f.location.nofuzzy_start + 1
								end = f.location.nofuzzy_end
								strand = f.location.strand
								if strand == 1:
									strand = '+'
								elif strand == -1:
									strand = '-'
								if lookup:
									if ID in lookup:
										orthologs = ','.join(lookup.get(ID))
									else:
										orthologs = 'NA'
								else:
									orthologs = 'NA'
								tabout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (ID, 'CDS', chr, str(start), str(end), strand, product, orthologs))
							if f.type == 'tRNA':
								ID = f.qualifiers['locus_tag'][0]
								start = f.location.nofuzzy_start
								end = f.location.nofuzzy_end
								strand = f.location.strand
								if strand == 1:
									strand = '+'
								elif strand == -1:
									strand = '-'
								try:
									product = f.qualifiers['product'][0]
								except KeyError:
									product = "tRNA-XXX"
								chr = record.id
								tabout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\tNA\n" % (ID, 'tRNA', chr, str(start), str(end), strand, product))

#run function
TransOut = args.out+'.transcripts.fasta'
ProtOut = args.out+'.proteins.fasta'
TabOut = args.out+'.tsv'
if args.proteinortho:
	genename = getGBbaseName(args.input)
	Orthologs = protortho2dict(args.proteinortho, genename)
	gb2tab(args.input, ProtOut, TransOut, TabOut, Orthologs)
else:
	gb2tab(args.input, ProtOut, TransOut, TabOut, False)

				