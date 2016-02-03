#!/usr/bin/env python
import sys

'''
Written by Jon Palmer, nextgen.usfs@gmail.com
last updated 2/3/2016
This script takes the output of ITSx to parse the locations of the SSU, ITS1, 5.8S, ITS2, and LSU regions.

Run ITSx as follows:
ITSx -i multi-fasta_ITS.fa -o output -t F --allow_single_domain

This will produce an output file called, 'output.positions.txt', you can then run the following:

python ITSx2tbl.py output.positions.txt > genbank.tbl
'''

with open(sys.argv[1], 'rU') as input:
    for line in input:
        if line.startswith('\n'):
            continue
        cols = line.split('\t')
        ID = cols[0]
        length = cols[1].replace(' bp.', '')
        SSU = cols[2].replace('SSU: ', '')
        if SSU.startswith('No'):
            SSU = ''
        else:
            SSUstart = SSU.split('-')[0]
            SSUend = SSU.split('-')[1]
        ITS1 = cols[3].replace('ITS1: ', '')
        if ITS1.startswith('No'):
            ITS1 = ''
        else:
            ITS1start = ITS1.split('-')[0]
            ITS1end = ITS1.split('-')[1]
        S = cols[4].replace('5.8S: ', '')
        if S.startswith('No'):
            S = ''
        else:
            Sstart = S.split('-')[0]
            Send = S.split('-')[1]
        ITS2 = cols[5].replace('ITS2: ', '')
        if ITS2.startswith('No'):
            ITS2 = ''
        else:
            ITS2start = ITS2.split('-')[0]
            ITS2end = ITS2.split('-')[1]
        LSU = cols[6].replace('LSU: ', '')
        if LSU.startswith('No'):
            LSU = ''
        else:
            LSUstart = LSU.split('-')[0]
            LSUend = LSU.split('-')[1]
        #okay, so all 4 regions are parsed, if empty they do not exist
        sys.stdout.write(">Feature %s\n" % cols[0])
        #sys.stdout.write("1\t%s\trRNA\n" % length)
        if SSU and ITS1:
            sys.stdout.write("<%s\t%s\trRNA\n" % (SSUstart, SSUend))
            sys.stdout.write("\t\t\tproduct\tpartial small subunit (18S) rRNA\n")
            sys.stdout.write("%s\t%s\trRNA\n" % (ITS1start, ITS2end))
            sys.stdout.write("\t\t\tproduct\tinternal transcribed spacer 1 (ITS1) of rRNA\n")
        elif ITS1 and int(ITS1start) != 1:
            SSUend = int(ITS1start) - 1
            sys.stdout.write("<1\t%i\trRNA\n" % (SSUend))
            sys.stdout.write("\t\t\tproduct\tpartial small subunit (18S) rRNA\n")
            sys.stdout.write("%s\t%s\trRNA\n" % (ITS1start, ITS2end))
            sys.stdout.write("\t\t\tproduct\tinternal transcribed spacer 1 (ITS1) of rRNA\n")
        elif ITS1:
            sys.stdout.write("%s\t%s\trRNA\n" % (ITS1start, ITS2end))
            sys.stdout.write("\t\t\tproduct\tinternal transcribed spacer 1 (ITS1) of rRNA\n")           
        if S and ITS2:
            sys.stdout.write("%s\t%s\trRNA\n" % (Sstart, Send))
            sys.stdout.write("\t\t\tproduct\t5.8S rRNA\n")        
            sys.stdout.write("%s\t%s\trRNA\n" % (ITS2start, ITS2end))
            sys.stdout.write("\t\t\tproduct\tinternal transcribed spacer 2 (ITS2) of rRNA\n")
        elif ITS2 and int(ITS2start) != 1:
            Send = int(ITS2start) - 1
            sys.stdout.write("<1\t%s\trRNA\n" % (Send))
            sys.stdout.write("\t\t\tproduct\tpartial 5.8S rRNA\n")        
            sys.stdout.write("%s\t%s\trRNA\n" % (ITS2start, ITS2end))
            sys.stdout.write("\t\t\tproduct\tinternal transcribed spacer 2 (ITS2) of rRNA\n")
        elif ITS2:
            sys.stdout.write("%s\t%s\trRNA\n" % (ITS2start, ITS2end))
            sys.stdout.write("\t\t\tproduct\tinternal transcribed spacer 2 (ITS2) of rRNA\n")        
        if LSU:
            sys.stdout.write("%s\t>%s\trRNA\n" % (LSUstart, LSUend))
            sys.stdout.write("\t\t\tproduct\tpartial large subunit (28S) rRNA\n")
        elif ITS2end < length:
            LSUstart = int(ITS2end) + 1
            sys.stdout.write("%s\t>%s\trRNA\n" % (LSUstart, length))
            sys.stdout.write("\t\t\tproduct\tpartial large subunit (28S) rRNA\n")   