#!/usr/bin/env python

#script takes a GBK input and outputs proteins, transcripts, DNA, and minimalistic GFF file

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
                                	if not 'pseudo' in f.qualifiers:
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
										if num_exons < 2: #only a single exon
											ex_start = str(f.location.nofuzzy_start + 1)
											ex_end = str(f.location.nofuzzy_end)
											gff.write("%s\tGenBank\texon\t%s\t%s\t.\t%s\t.\tID=%s-T1.exon1;Parent=%s-T1\n" % (chr, ex_start, ex_end, strand, ID, ID))
											gff.write("%s\tGenBank\tCDS\t%s\t%s\t.\t%s\t0\tID=%s-T1.cds;Parent=%s-T1\n" % (chr, ex_start, ex_end, strand, ID, ID))
										else: #more than 1 exon, so parts are actually in correct orientation, so loop through
											for i in range(0,num_exons):
												ex_start = str(f.location.parts[i].nofuzzy_start + 1)
												ex_end = str(f.location.parts[i].nofuzzy_end)
												ex_num = i + 1
												gff.write("%s\tGenBank\texon\t%s\t%s\t.\t%s\t.\tID=%s-T1.exon%i;Parent=%s-T1\n" % (chr, ex_start, ex_end, strand, ID, ex_num, ID))
												gff.write("%s\tGenBank\tCDS\t%s\t%s\t.\t%s\t%i\tID=%s-T1.cds;Parent=%s-T1\n" % (chr, ex_start, ex_end, strand, current_phase, ID, ID))
												current_phase = (current_phase - (int(ex_end) - int(ex_start) + 1)) % 3
												if current_phase == 3:
													current_phase = 0


#setup outputs
genome = args.out+'.genome.fa'
transcript = args.out+'.transcripts.fa'
proteins = args.out+'.proteins.fa'
gff = args.out+'.gff3'

gb2allout(args.input, gff, proteins, transcript, genome)

