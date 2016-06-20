#!/usr/bin/env python

#script to process snippy output to phylogenetic tree

import os, sys, argparse, subprocess, shutil, random, pylab
from Bio import SeqIO
from Bio import Phylo
from Bio.Phylo.Consensus import get_support
from natsort import natsorted

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='snippy2tree.py',
    description='''Script that converts snippy-core output to phylogeny via several different methods.''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', default='core', help='snippy prefix name (basename)')
parser.add_argument('-m','--method', required=True, choices=['all', 'binary', 'snps', 'window'], default='snps', help='method to use for reconstructing phylogeny')
parser.add_argument('--bootstrap', default=100, type=int, help='Number of bootstraps to run with RAxML')
parser.add_argument('--outgroup', default='Reference', help='Outgroup for RAxML')
parser.add_argument('--window', default=25, type=int, help='flanking region for each SNP in window method')
parser.add_argument('--cpus', default=2, type=int, help='number of cpus')
args=parser.parse_args()

FNULL = open(os.devnull, 'w')

def getheaders(input):
    #get the sample names
    with open(input+'.tab', 'rU') as info:
        header_line = next(info)
        cols = header_line.split('\t')
        headers = cols[2:-4]
    return headers

#snps method - take alignment directly from snippy and run RAxML
def runsnps(input, outgroup):
    print "Running concatenated SNP method."
    file_in = input+'.aln'
    tmpdir = 'snps_tmp'
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir)
    if not outgroup:
        subprocess.call(['raxmlHPC-PTHREADS', '-T', str(args.cpus), '-f', 'a', '-m', 'GTRGAMMA', '-p', '12345', '-x', '12345', '-#', str(args.bootstrap), '-s', os.path.abspath(file_in), '-n', 'snps.nwk'], cwd = tmpdir, stdout = FNULL, stderr = FNULL)
    else:
        subprocess.call(['raxmlHPC-PTHREADS', '-T', str(args.cpus), '-f', 'a', '-m', 'GTRGAMMA', '-p', '12345', '-x', '12345', '-#', str(args.bootstrap), '-s', os.path.abspath(file_in), '-o', outgroup, '-n', 'snps.nwk'], cwd = tmpdir, stdout = FNULL, stderr = FNULL)

     #parse with biopython and draw
    trees = list(Phylo.parse(os.path.join(tmpdir, 'RAxML_bootstrap.snps.nwk'), 'newick'))
    best = Phylo.read(os.path.join(tmpdir,'RAxML_bestTree.snps.nwk'), 'newick')
    support_tree = get_support(best, trees)
    Phylo.write(support_tree, args.input+'.snps.phylogeny.nwk', 'newick')
    Phylo.draw(support_tree, do_show=False)
    pylab.axis('off')
    pylab.savefig(args.input+'.snps.phylogeny.pdf', format='pdf', bbox_inches='tight', dpi=1000)
    shutil.rmtree(tmpdir)

#take snps, pull out window around snp, concatenate and then run RAxML
def runwindow(input, outgroup):
    print "Running sliding-window method, using %i bp flanking each SNP." % args.window
    tmpdir = 'window_tmp'
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir)
    snpdb = {}
    count = 0
    with open(input+'.tab', 'rU') as data:
        next(data) #skip header
        for line in data:
            count += 1
            cols = line.split('\t')
            scaffold = cols[0]
            pos = cols[1]
            name = 'snp'+str(count)
            if not name in snpdb:
                snpdb[name] = [scaffold, int(pos) - args.window-1, int(pos) + args.window]
    sequencedb = {}
    for num, id in enumerate(samples):
        if id == 'Reference':
            fastafile = os.path.join(samples[num+1], 'reference', 'ref.fa')
        else:
            fastafile = os.path.join(id, 'snps.consensus.subs.fa')
        sequence = []
        with open(fastafile, 'rU') as fasta:
            sequence_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
            for k,v in natsorted(snpdb.iteritems()):
                rec = sequence_dict.get(v[0])
                sequence.append(str(rec.seq[v[1]:v[2]]))
        if id not in sequencedb:
            sequencedb[id] = ''.join(sequence)
    with open(os.path.join(tmpdir, 'windows.fa'), 'w') as output:
        for k,v in sequencedb.iteritems():
            output.write('>%s\n%s\n' % (k,v))
    #run alignment on output
    with open(os.path.join(tmpdir, 'windows.mafft.fa'), 'w') as output:
        subprocess.call(['mafft', '--thread', str(args.cpus), '--auto', os.path.join(tmpdir, 'windows.fa')], stdout = output, stderr = FNULL)
    #now trim alignment
    subprocess.call(['trimal', '-in', os.path.join(tmpdir, 'windows.mafft.fa'), '-out', os.path.join(tmpdir, 'windows.trimal.fa'), '-automated1'])

    if not outgroup:
        subprocess.call(['raxmlHPC-PTHREADS', '-T', str(args.cpus), '-f', 'a', '-m', 'GTRGAMMA', '-p', '12345', '-x', '12345', '-#', str(args.bootstrap), '-s', 'windows.trimal.fa', '-n', 'windows.nwk'], cwd = tmpdir, stdout = FNULL, stderr = FNULL)
    else:
        subprocess.call(['raxmlHPC-PTHREADS', '-T', str(args.cpus), '-f', 'a', '-m', 'GTRGAMMA', '-p', '12345', '-x', '12345', '-#', str(args.bootstrap), '-s', 'windows.trimal.fa', '-o', outgroup, '-n', 'windows.nwk'], cwd = tmpdir, stdout = FNULL, stderr = FNULL)

     #parse with biopython and draw
    trees = list(Phylo.parse(os.path.join(tmpdir, 'RAxML_bootstrap.windows.nwk'), 'newick'))
    best = Phylo.read(os.path.join(tmpdir,'RAxML_bestTree.windows.nwk'), 'newick')
    support_tree = get_support(best, trees)
    Phylo.write(support_tree, args.input+'.windows.phylogeny.nwk', 'newick')
    Phylo.draw(support_tree, do_show=False)
    pylab.axis('off')
    pylab.savefig(args.input+'.windows.phylogeny.pdf', format='pdf', bbox_inches='tight', dpi=1000)
    shutil.rmtree(tmpdir)

#parse the snippy vcf file, pull out the "binary" alleles, concatenate, and run RAxML
def runbinary(input, outgroup):
    print "Runnning binary (presence/absence) method."
    tmpdir = 'binary_tmp'
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir)
    datadb = []
    with open(input+'.binary.fa', 'w') as output:
        with open(input+'.vcf', 'rU') as vcf:
            for line in vcf:
                line = line.replace('\n', '')
                if line.startswith('##'):
                    continue
                data = line.split('\t')[9:]
                if '3' in data:
                    data1 = ['0' if x=='3' else x for x in data]
                    data1 = ['0' if x=='2' else x for x in data1]
                    datadb.append(data1)
                    data2 = ['0' if x=='1' else x for x in data]
                    data2 = ['1' if x=='2' else x for x in data2]
                    data2 = ['0' if x=='3' else x for x in data2]
                    data3 = ['0' if x=='1' else x for x in data]
                    data3 = ['0' if x=='2' else x for x in data3]
                    data3 = ['1' if x=='3' else x for x in data3]
                    datadb.append(data3)
                elif '2' in data:
                    data1 = ['0' if x=='2' else x for x in data]
                    datadb.append(data1)
                    data2 = ['0' if x=='1' else x for x in data]
                    data2 = ['1' if x=='2' else x for x in data2]
                    datadb.append(data2)
                else:
                    datadb.append(data)
        binarydata = [list(x) for x in zip(*datadb)]
        for i in binarydata:
            output.write(">%s\n%s\n" % (i[0], ''.join(i[1:])))

    if not outgroup:
        subprocess.call(['raxmlHPC-PTHREADS', '-T', str(args.cpus), '-f', 'a', '-m', 'BINGAMMA', '-p', '12345', '-x', '12345', '-#', str(args.bootstrap), '-s', os.path.abspath(input+'.binary.fa'), '-n', 'binary.nwk'], cwd = tmpdir, stdout = FNULL, stderr = FNULL)
    else:
        subprocess.call(['raxmlHPC-PTHREADS', '-T', str(args.cpus), '-f', 'a', '-m', 'BINGAMMA', '-p', '12345', '-x', '12345', '-#', str(args.bootstrap), '-s', os.path.abspath(input+'.binary.fa'), '-o', outgroup, '-n', 'binary.nwk'], cwd = tmpdir, stdout = FNULL, stderr = FNULL)

     #parse with biopython and draw
    trees = list(Phylo.parse(os.path.join(tmpdir, 'RAxML_bootstrap.binary.nwk'), 'newick'))
    best = Phylo.read(os.path.join(tmpdir,'RAxML_bestTree.binary.nwk'), 'newick')
    support_tree = get_support(best, trees)
    Phylo.write(support_tree, args.input+'.binary.phylogeny.nex', 'nexus')
    Phylo.draw(support_tree, do_show=False)
    pylab.axis('off')
    pylab.savefig(args.input+'.binary.phylogeny.pdf', format='pdf', bbox_inches='tight', dpi=1000)
    #shutil.rmtree(tmpdir)

#first check the input and make sure files are present
if not os.path.isfile(args.input+'.tab'):
    print "Error: %s.tab not found, check --input parameter" % args.input
    os._exit(1)
if not os.path.isfile(args.input+'.vcf'):
    print "Error: %s.vcf not found, check --input parameter" % args.input
    os._exit(1)
if not os.path.isfile(args.input+'.aln'):
    print "Error: %s.aln not found, check --input parameter" % args.input
    os._exit(1)

#determine which samples were run with snippy-core
global samples
samples = getheaders(args.input)
if not args.outgroup in samples:
    print "%s not found in snippy-core ouput: [%s]" % (args.outgroup, ", ".join(samples))
    print "Running RAxML with no outgroup."
    outgroup = False
else:
    print "Running RAxML with %s as outgroup." % args.outgroup
    outgroup = args.outgroup

#run specified method
if args.method == 'snps':
    runsnps(args.input, outgroup)
elif args.method == 'window':
    runwindow(args.input, outgroup)
elif args.method == 'binary':
    runbinary(args.input, outgroup)
elif args.method == 'all':
    runsnps(args.input, outgroup)
    runbinary(args.input, outgroup)
    runwindow(args.input, outgroup)

print "Script has completed successfully!"