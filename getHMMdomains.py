#!/usr/bin/env python

#script takes input protein fasta and HMM model to pull out those domains
#written by Jon Palmer palmer.jona at gmail dot com
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SearchIO
import sys
import os
from Bio import SeqIO
import argparse
from os.path import expanduser
home = expanduser("~")

#define HMM db here
hmm_db = home + '/projects/DB/HMM'

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)
parser=argparse.ArgumentParser(prog='getHMMdomains.py', usage="%(prog)s [options] proteins.fasta > domains.fasta\n%(prog)s -h for help menu",
    description='''Script to identify HMM domains from proteins and pull out only sequence corresponding to the domain in multi-fasta format.  HMMER3 is a dependency''',
    formatter_class=MyFormatter,
    epilog="""Written by Jon Palmer (2015) palmer.jona@gmail.com""")
parser.add_argument('fasta', help='Fasta input (aa seq)')
parser.add_argument('-m','--hmm', default='KS_domain.hmm', help='HMM model (Can also specificy full path)')
parser.add_argument('-e','--evalue', default="1e-10", help='HMM Evalue cutoff')
parser.add_argument('-l','--length', default="50", help='min length of domain to keep')
parser.add_argument('-c','--cpus', default="1", help='HMM model')
args=parser.parse_args()
    
fasta = args.fasta
if not "/" in args.hmm:
    hmm = hmm_db + '/' + args.hmm
else:
    hmm = args.hmm
evalue = args.evalue
length = args.length
cpus = args.cpus

#run hmmerscan search
os.system("%s %s %s %s %s %s %s" % ('hmmscan --cpu', cpus, '-E', evalue, hmm, fasta, '> hmmscan.temp.txt'))

#format results and print fasta seq based on coordinates
hmmer_results = open('hmmscan.temp.txt', 'r')
inputSeqFile = open(fasta, "rU")
SeqRecords = SeqIO.to_dict(SeqIO.parse(inputSeqFile, "fasta"))
for qresult in SearchIO.parse(hmmer_results, "hmmer3-text"):
    hits = qresult.hits
    if len(hits) > 0:
        beste = hits[0].hsps[0].evalue
        query = hits[0].query_id
        hit = hits[0].id.replace('.hmm', '')
        start = hits[0].hsps[0].env_start
        end = hits[0].hsps[0].env_end
        hit_length = end - start
        if int(hit_length) >= int(length) and float(beste) <= float(evalue):
            subseq = SeqRecords[query][start:end].seq
            print("%s%s|%s|%s-%s\n%s" % ('>', query, hit, start, end, subseq))

#close and clean
hmmer_results.close()
inputSeqFile.close()
os.remove('hmmscan.temp.txt')