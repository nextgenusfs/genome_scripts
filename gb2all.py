#!/usr/bin/env python

#script takes a GBK input and outputs parts

import os,sys, argparse
from Bio import SeqIO

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='gb2all.py',
    description='''Script that converts GenBank flatfile to it's "parts".''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='GenBank flatfile (annotated)')
parser.add_argument('-o','--out', required=True, help='Basename output')
args=parser.parse_args()


def gb2allout(input, GFF, Proteins, Transcripts, DNA):
    #this will not output any UTRs for gene models, don't think this is a problem right now....
    with open(GFF, 'w') as gff:
        gff.write("##gff-version 3\n")
        with open(Proteins, 'w') as proteins:
            with open(Transcripts, 'w') as transcripts:
                with open(DNA, 'w') as scaffolds:
                    with open(input, 'rU') as gbk:
                        for record in SeqIO.parse(gbk, 'genbank'):
                            scaffolds.write(">%s\n%s\n" % (record.id, record.seq))
                            for f in record.features:
                                if f.type == "mRNA":
                                    feature_seq = f.extract(record.seq)
                                    transcripts.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], feature_seq))
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
                                    num_exons = len(f.location.parts)
                                    current_phase = 0
                                    gff.write("%s\tGenBank\tgene\t%s\t%s\t.\t%s\t.\tID=%s\n" % (chr, start, end, strand, ID))
                                    gff.write("%s\tGenBank\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s-T1;Parent=%s;product=%s\n" % (chr, start, end, strand, ID, ID, product))


#setup outputs
genome = args.out+'.genome.fa'
transcript = args.out+'.transcripts.fa'
proteins = args.out+'.proteins.fa'
gff = args.out+'.gff3'

gb2allout(args.input, gff, proteins, transcript, genome)

