#!/usr/bin/env python

#This script takes folder of GBK files, reads taxonomy, and moves them into folders based on phylum
#written by Jon Palmer palmer.jona at gmail dot com

import sys, argparse, Bio, glob, os
from Bio import SeqIO
class bcolors:
    OKGREEN = '\033[92m'
    BLUE = '\033[36m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

parser=argparse.ArgumentParser(
    description='''Script that reads folder of GBK files and splits them by taxonomy''',
    epilog="""Written by Jon Palmer (2015)  palmer.jona@gmail.com""")
parser.add_argument('-d', '--directory', default="current", help='folder')
args=parser.parse_args()

if args.directory == "current":
    dir = os.getcwd()
else:
    dir = os.chdir(args.directory)
os.chdir(dir)
for file in glob.glob("*.gbk"):
    print bcolors.WARNING + "Working on: " + bcolors.ENDC + file
    file_path = os.path.abspath(file)
    dir_path = os.path.abspath(dir)
    record = next(SeqIO.parse(file, "genbank"))
    tax = record.annotations["taxonomy"]
    try:
        if 'Ascomycota' in tax:
            level = "Ascomycota"
        if 'Basidiomycota' in tax:
            level = "Basidiomycota"
        if 'Microsporidia' in tax:
            level = "Microsporidia"
        if 'Zygomycota' in tax:
            level = "Zygomycota"
        if 'Mucormycotina' in tax:
            level = "Zygomycota"
        if 'Chytridiomycota' in tax:
            level = "Chytridiomycota"
        if 'Rozella' in tax:
            level = "Chytridiomycota"
        if 'Glomeromycota' in tax:
            level = "Glomeromycota"
    except IndexError:
        level = "no_taxonomy"  
    if not os.path.exists(level):
        os.makedirs(level)
    new_path = dir_path + "/" + level + "/" + file
    os.rename(file_path, new_path)