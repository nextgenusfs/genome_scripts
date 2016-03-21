#!/usr/bin/env python

from Bio import SeqIO
import sys

sys.stdout.write("##gff-version 3\n")
with open(sys.argv[1], 'rU') as gbk:
    for record in SeqIO.parse(gbk, 'genbank'):
        for f in record.features:
            if f.type == 'CDS':
                chr = record.id
                ID = f.qualifiers['locus_tag'][0]
                product = f.qualifiers['product'][0]
                start = f.location.nofuzzy_start + 1
                end = f.location.nofuzzy_end + 1
                strand = f.location.strand
                if strand == 1:
                    strand = '+'
                elif strand == -1:
                    strand = '-'
                num_exons = len(f.sub_features)
                current_phase = 0
                sys.stdout.write("%s\tGenBank\tgene\t%s\t%s\t.\t%s\t.\tID=%s\n" % (chr, start, end, strand, ID))
                sys.stdout.write("%s\tGenBank\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s-T1;Parent=%s;product=%s\n" % (chr, start, end, strand, ID, ID, product))
                if num_exons < 1: #only a single exon
                    ex_start = str(f.location.nofuzzy_start + 1)
                    ex_end = str(f.location.nofuzzy_end + 1)
                    sys.stdout.write("%s\tGenBank\texon\t%s\t%s\t.\t%s\t.\tID=%s-T1.exon1;Parent=%s-T1\n" % (chr, ex_start, ex_end, strand, ID, ID))
                    sys.stdout.write("%s\tGenBank\tCDS\t%s\t%s\t.\t%s\t%i\tID=%s-T1.cds;Parent=%s-T1\n" % (chr, ex_start, ex_end, strand, current_phase, ID, ID))
                else: #more than 1 exon, so parts sub_features
                    if f.location.strand == 1:
                        for i in range(0,num_exons):
                            ex_start = str(f.sub_features[i].location.nofuzzy_start + 1)
                            ex_end = str(f.sub_features[i].location.nofuzzy_end + 1)
                            ex_num = i + 1
                            sys.stdout.write("%s\tGenBank\texon\t%s\t%s\t.\t%s\t.\tID=%s-T1.exon%i;Parent=%s-T1\n" % (chr, ex_start, ex_end, strand, ID, ex_num, ID))
                            sys.stdout.write("%s\tGenBank\tCDS\t%s\t%s\t.\t%s\t%i\tID=%s-T1.cds;Parent=%s-T1\n" % (chr, ex_start, ex_end, strand, current_phase, ID, ID))
                            current_phase = (current_phase - (int(ex_end) - int(ex_start) + 1)) % 3
                            if current_phase == 3:
                                current_phase = 0
                    else:
                        for i in reversed(range(0,num_exons)):
                            ex_start = str(f.sub_features[i].location.nofuzzy_start + 1)
                            ex_end = str(f.sub_features[i].location.nofuzzy_end + 1)
                            ex_num = num_exons - i
                            sys.stdout.write("%s\tGenBank\texon\t%s\t%s\t.\t%s\t.\tID=%s-T1.exon%i;Parent=%s-T1\n" % (chr, ex_start, ex_end, strand, ID, ex_num, ID))
                            sys.stdout.write("%s\tGenBank\tCDS\t%s\t%s\t.\t%s\t%i\tID=%s-T1.cds;Parent=%s-T1\n" % (chr, ex_start, ex_end, strand, current_phase, ID, ID))
                            current_phase = (current_phase - (int(ex_end) - int(ex_start) + 1)) % 3
                            if current_phase == 3:
                                current_phase = 0

            if f.type == 'tRNA':
                ID = f.qualifiers['locus_tag'][0]
                start = f.location.nofuzzy_start
                end = f.location.nofuzzy_end
                strand = f.location.strand
                if strand == 1:
                    strand = '+'
                elif strand == -1:
                    strand = '-'
                product = f.qualifiers['product'][0]
                chr = record.id
                sys.stdout.write("%s\tGenBank\tgene\t%s\t%s\t.\t%s\t.\tID=%s\n" % (chr, start, end, strand, ID))
                sys.stdout.write("%s\tGenBank\ttRNA\t%s\t%s\t.\t%s\t.\tID=%s-T1;Parent=%s;product=%s\n" % (chr, start, end, strand, ID, ID, product))
                sys.stdout.write("%s\tGenBank\texon\t%s\t%s\t.\t%s\t.\tID=%s-T1.exon1;Parent=%s-T1\n" % (chr, start, end, strand, ID, ID))