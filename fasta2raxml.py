#!/usr/bin/env python

import sys, pylab, argparse, os, subprocess, shutil
from Bio import Phylo
from Bio.Phylo.Consensus import get_support

#script takes multi-fasta file, runs MAFFT/trimAl/RAxML and then draw tree in pdf format
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

parser=argparse.ArgumentParser(prog='fasta2raxml.py',
    description='''Script runs Mafft -> trimAl -> RAxML. takes multi-fasta as input''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--input', required=True, help='Input multi-fasta file')
parser.add_argument('--outgroup', help='Fasta header name of outgroup for RAxML')
parser.add_argument('--bootstrap', default=100, type=int, help='Num bootstrap replicates to run')
parser.add_argument('--debug', action='store_true', help='Keep all intermediate files')
parser.add_argument('--cpus', default=1, type=int, help='Number of cpus')
parser.add_argument('-m', '--method', default='prot', choices=['prot', 'nucl'], help='Protein or Nucleotide aligment')
parser.add_argument('-o','--out', required=True, help='Output phylogeny tree name (PDF)')
args=parser.parse_args()

FNULL = open(os.devnull, 'w')

#make tmpdir
tmp = 'fasta2raxml_tmp'
if os.path.isdir(tmp):
    shutil.rmtree(tmp)
os.makedirs(tmp)

print "Aligning sequences using MAFFT"
mafft_out = os.path.join(tmp, 'alignment.fa')
with open(mafft_out, 'w') as output:
    subprocess.call(['mafft', '--thread', str(args.cpus), args.input], stdout = output, stderr = FNULL)

print "Trimming alignment using trimAl"
trimal_out = os.path.join(tmp, 'trimal.phylip')
trimal_out = os.path.abspath(trimal_out)
subprocess.call(['trimal', '-in', mafft_out, '-out', trimal_out, '-automated1', '-phylip'], stderr = FNULL, stdout = FNULL)

print "Running RAxML"
if args.method == 'prot':
    method = 'PROTGAMMAAUTO'
elif args.method == 'nucl':
    method = 'GTRGAMMA'
if args.cpus == 1:
    if not args.outgroup:
        subprocess.call(['raxmlHPC-PTHREADS', '-f', 'a', '-m', method, '-p', '12345', '-x', '12345', '-#', str(args.bootstrap), '-s', trimal_out, '-n', 'nwk'], cwd = tmp)
    else:
        subprocess.call(['raxmlHPC-PTHREADS', '-f', 'a', '-m', method, '-o', args.outgroup, '-p', '12345', '-x', '12345', '-#', str(args.bootstrap), '-s', trimal_out, '-n', 'nwk'], cwd = tmp)
else:
    if not args.outgroup:
        subprocess.call(['raxmlHPC-PTHREADS', '-T', str(args.cpus), '-f', 'a', '-m', method, '-p', '12345', '-x', '12345', '-#', str(args.bootstrap), '-s', trimal_out, '-n', 'nwk'], cwd = tmp)
    else:
        subprocess.call(['raxmlHPC-PTHREADS', '-T', str(args.cpus), '-f', 'a', '-m', method, '-o', args.outgroup, '-p', '12345', '-x', '12345', '-#', str(args.bootstrap), '-s', trimal_out, '-n', 'nwk'], cwd = tmp)

print "Drawing tree inferred from RAxML"
#parse with biopython and draw
trees = list(Phylo.parse(os.path.join(tmp, 'RAxML_bootstrap.nwk'), 'newick'))
best = Phylo.read(os.path.join(tmp, 'RAxML_bestTree.nwk'), 'newick')
support_tree = get_support(best, trees)
Phylo.write(support_tree, args.out.split('.')[0]+'.nwk', 'newick')
Phylo.draw(support_tree, do_show=False)
pylab.axis('off')
pylab.savefig(args.out, format='pdf', bbox_inches='tight', dpi=1000) 

if not args.debug:
    shutil.rmtree(tmp)
