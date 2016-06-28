#!/usr/bin/env python

import sys, os, subprocess, argparse
from natsort import natsorted

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='snippy2indelavg.py', usage="%(prog)s [options] -i genome1 genome2",
    description='''Average the distribution of snippy indels .''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', nargs='+', required=True, help='List of funannotate genome folders')
args=parser.parse_args()

def bcftools_avg(folder):
    global data
    file = os.path.join(folder, 'snps.vcf.gz')
    tmpout = 'snippy2indel.tmp'
    with open(tmpout, 'w') as output:
        subprocess.call(['bcftools', 'stats', file], stdout = output)
    with open(tmpout, 'rU') as input:
        for line in input:
            line = line.replace('\n', '')
            if line.startswith('IDD'):
                cols = line.split('\t')
                location = int(cols[2])
                if not location in data:
                    data[location] = int(cols[3])
                else:
                    num = data.get(location) + int(cols[3])
                    data[location] = num
data = {}
num_input = len(args.input)
for i in range(0, num_input):
    bcftools_avg(args.input[i])

#divide by number of input
for k,v in data.items():
    data[k] = v / num_input
    
for k,v in natsorted(data.iteritems()):
    sys.stdout.write('%s\t%s\n' % (k,v))