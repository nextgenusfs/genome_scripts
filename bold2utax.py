#!/usr/bin/env python

import sys, re

with open(sys.argv[1], 'rU') as input:
    for line in input:
        line = line.replace('\n', '')
        if line.startswith('processid'):
            continue
        col = line.split('\t')
        K = 'k:Animalia'
        P = 'p:'+col[8].strip()
        C = 'c:'+col[10].strip()
        O = 'o:'+col[12].strip()
        F = 'f:'+col[14].strip()
        G = 'g:'+col[18].strip()
        S = 's:'+col[20].strip()
        ID = col[39].strip()
        GB = col[41].strip()
        Seq = col[42].replace('-', '')
        Seq = re.sub('N*$', '', Seq)
        Seq = re.sub('^N*', '', Seq)
        tax = []
        for i in [K,P,C,O,F,G,S]:
            if not i.endswith(':'):
                tax.append(i)
        tax_fmt = ','.join(tax)
        if not GB == '':
            sys.stdout.write('>%s|%s;tax=%s\n%s\n' % (ID, GB, tax_fmt, Seq))
        else:
            sys.stdout.write('>%s;tax=%s\n%s\n' % (ID, tax_fmt, Seq))
