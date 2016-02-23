#!/usr/bin/env python
import sys, argparse, os
from Bio import SeqIO

parser=argparse.ArgumentParser(
    description='''Script that takes Genbank file input and outputs ProteinOrtho input ''',
    epilog="""Jon Palmer (2015)  palmer.jona@gmail.com""")
parser.add_argument('-i','--input', required=True, help='Genbank file (Required)')
parser.add_argument('-n','--name', required=True, help='Base name for output files')
parser.add_argument('-o','--out', default='protortho', help='Output folder')
args=parser.parse_args()

def gb2proteinortho(input, folder, name):
    history = []
    gffOut = os.path.join(folder, name+'.gff')
    FastaOut = os.path.join(folder, name+'.faa')
    with open(gffOut, 'w') as gff:
        with open(FastaOut, 'w') as fasta:
            with open(input, 'rU') as input:
                SeqRecords = SeqIO.parse(input, 'genbank')
                for record in SeqRecords:
                    for f in record.features:
                        if f.type == 'CDS':
                            protID = f.qualifiers['protein_id'][0]
                            locusID = f.qualifiers['locus_tag'][0]
                            start = f.location.nofuzzy_start
                            end = f.location.nofuzzy_end
                            strand = f.location.strand
                            if strand == 1:
                                strand = '+'
                            elif strand == -1:
                                strand = '-'
                            translation = f.qualifiers['translation'][0]
                            product = f.qualifiers['product'][0]
                            chr = record.id
                            if '.' in chr:
                                chr = chr.split('.')[0]
                            if not protID in history:
                                history.append(protID)
                                gff.write("%s\tNCBI\tCDS\t%s\t%s\t.\t%s\t.\tID=%s;Alias=%s;Product=%s;\n" % (chr, start, end, strand, locusID, protID, product))
                                fasta.write(">%s\n%s\n" % (locusID, translation))
                                
if not os.path.isdir(args.out):
    os.makedirs(args.out)
gb2proteinortho(args.input, args.out, args.name)
print "Script complete.  %s.gff and %s.faa saved in %s" % (args.name, args.name, args.out)
