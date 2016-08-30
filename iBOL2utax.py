#!/usr/bin/env python

import sys, re

if len(sys.argv) < 2:
    print("Purpose: Script will parse Bold DB txt output to fasta file, keeping only records deposited in GenBank\nUsage: bold2utax.py bold_data.txt > formatted.fa")
    sys.exit(1)

with open(sys.argv[1], 'rU') as input:
    for line in input:
        line = line.replace('\n', '')
        if line.startswith('processid'):
            header = line.split('\t')
            pid = header.index('phylum_reg')
            cid = header.index('class_reg')
            oid = header.index('order_reg')
            fid = header.index('family_reg')
            gid = header.index('genus_reg')
            sid = header.index('species_reg')
            seqid = header.index('nucraw')
            boldid = header.index('seqdataid')
            gbid = header.index('accession')
            continue
        col = line.split('\t')
        K = 'k:Animalia'
        P = 'p:'+col[pid].strip()
        C = 'c:'+col[cid].strip()
        O = 'o:'+col[oid].strip()
        F = 'f:'+col[fid].strip()
        G = 'g:'+col[gid].strip()
        S = 's:'+col[sid].strip().replace('.', '')
        if ' sp ' in S: #remove those that have sp. in them
            S = ''
        ID = col[boldid].strip()
        GB = col[gbid].strip()
        Seq = col[seqid].replace('-', '')
        Seq = re.sub('N*$', '', Seq)
        Seq = re.sub('^N*', '', Seq)
        tax = []
        for i in [K,P,C,O,F,G,S]:
            if not i.endswith(':'):
                tax.append(i)
        tax_fmt = ','.join(tax)
        if tax_fmt.endswith(','):
            tax_fmt = tax_fmt.rsplit(',',1)[0]
        sys.stdout.write('>%s;tax=%s\n%s\n' % (ID, tax_fmt, Seq))
