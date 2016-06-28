#!/usr/bin/env python

import sys, os, subprocess, argparse

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='snippy2stats.py', usage="%(prog)s [options] genome1 genome2",
    description='''Run BCFTOOLS stats on snippy output to generate some summary stats.''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', nargs='+', required=True, help='List of snippy folders')
args=parser.parse_args()

def bcftools_stats(folder):
    file = os.path.join(folder, 'snps.vcf.gz')
    tmpdata = 'snippy2stats.tmp'
    total = ''
    snps = ''
    mnps = ''
    indels = ''
    others = ''
    with open(tmpdata, 'w') as output:
        subprocess.call(['bcftools', 'stats', file], stdout = output)
    with open(tmpdata, 'rU') as input:
        for line in input:
            line = line.replace('\n', '')
            if line.startswith('SN'):
                if 'number of records:' in line:
                    total = line.split('\t')[-1]
                if 'number of SNPs:' in line:
                    snps = line.split('\t')[-1]
                if 'number of MNPs:' in line:
                    mnps = line.split('\t')[-1]
                if 'number of indels:' in line:
                    indels = line.split('\t')[-1]
                if 'number of others:' in line:
                    others = line.split('\t')[-1]
    os.remove(tmpdata)
    return (folder, total, snps, mnps, indels, others)

#write header
sys.stdout.write('Strain\tTotal\tSNPs\tMNPs\tIndels\tOthers\n')
num_input = len(args.input)
for i in range(0, num_input):
    stats = bcftools_stats(args.input[i])
    sys.stdout.write('\t'.join(stats)+'\n')
