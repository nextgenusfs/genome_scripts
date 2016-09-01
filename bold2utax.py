#!/usr/bin/env python

import sys, re, argparse
from natsort import natsorted

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='bold2utax.py',
    description='''Parse BOLD DB TSV data dump into FASTA with UTAX compatible labels.''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='Bold data dump TSV format')
parser.add_argument('-o','--out', required=True, help='UTAX formated FASTA output')
parser.add_argument('--primer', default='GGTCAACAAATCATAAAGATATTGG', help='Forward Primer Sequence')
parser.add_argument('--trim_length', default='250', help='Distance to Trim sequences to')
parser.add_argument('--primer_mismatch', default='4', help='Mismatches allowed in primer')
parser.add_argument('--require_genbank', action='store_true', help='Require output to have GenBank Accessions')
args=parser.parse_args()

LetterToSet = {}
LetterToSet['A'] = "A"
LetterToSet['C'] = "C"
LetterToSet['G'] = "G"
LetterToSet['T'] = "T"
LetterToSet['M'] = "AC"
LetterToSet['R'] = "AG"
LetterToSet['W'] = "AT"
LetterToSet['S'] = "CG"
LetterToSet['Y'] = "CT"
LetterToSet['K'] = "GT"
LetterToSet['V'] = "ACG"
LetterToSet['H'] = "ACT"
LetterToSet['D'] = "AGT"
LetterToSet['B'] = "CGT"
LetterToSet['X'] = "GATC"
LetterToSet['N'] = "GATC"

def MatchLetter(a, b):
	global LetterToSet
	try:
		sa = LetterToSet[a.upper()]
	except:
		return False
	try:
		sb = LetterToSet[b.upper()]
	except:
		return False
	for ca in sa:
		if ca in sb:
			return True
	return False

def MatchPrefix(Seq, Primer):
	L = len(Seq)
	PrimerLength = len(Primer)
	n = PrimerLength
	if L < n:
		n = L
	Diffs = 0
	for i in range(0, n):
		if not MatchLetter(Seq[i], Primer[i]):
			Diffs += 1
	return Diffs

def BestMatch(Seq, Primer):
	L = len(Seq)
	PrimerLength = len(Primer)
	BestDiffs = PrimerLength
	BestPos = -1
	for Pos in range(0, L-PrimerLength+1):
		d = MatchPrefix(Seq[Pos:], Primer)
		if d < BestDiffs:
			BestDiffs = d
			BestPos = Pos
	return BestPos, BestDiffs

def TrimPrimer(Sequence, primer, trimlen, mismatch):
    #find primer location
    Diffs = BestMatch(Sequence, primer)
    if Diffs[1] > int(mismatch):
        #assume that primer was trimmed off
        Seq = Sequence[:int(trimlen)]
    else:
        Seq = Sequence[Diffs[0]+len(primer):int(trimlen)]
    #do a length check, if shorter than 250 bp, discard
    if len(Seq) >= int(trimlen):
        return Seq

Total = -1
NoPrimer = 0
PrimerFound = 0
TooShort = 0
UniqueBIN = 0
DupBIN = 0

bins_seen = {}
with open(args.input, 'rU') as input:
    for line in input:
        Total += 1
        line = line.replace('\n', '')
        if line.startswith('processid'):
            header = line.split('\t')
            pid = header.index('phylum_name')
            cid = header.index('class_name')
            oid = header.index('order_name')
            fid = header.index('family_name')
            gid = header.index('genus_name')
            sid = header.index('species_name')
            seqid = header.index('nucleotides')
            boldid = header.index('sequenceID')
            gbid = header.index('genbank_accession')
            bin = header.index('bin_uri')
            idby = header.index('identification_provided_by')
            continue
        col = line.split('\t')
        #some idiots have collector names in these places in the DB, this DB is kind of a mess, doing the best I can....
        bd = col[idby].strip()
        badnames = bd.split(' ')
        badfiltered = ['1','2','3','4','5','6','7','8','9','0']
        for y in badnames:
            if len(y) > 2:
                if y != 'Art':
                    if y != 'Eric':
                        badfiltered.append(y)
        K = 'k:Animalia'
        P = 'p:'+col[pid].strip()
        C = 'c:'+col[cid].strip()
        O = 'o:'+col[oid].strip()
        F = 'f:'+col[fid].strip()
        G = 'g:'+col[gid].strip()
        S = 's:'+col[sid].strip().replace('.', '')
        if ' sp ' in S: #remove those that have sp. in them
            S = ''
        if S.endswith(' sp'):
            S = ''
        if badfiltered:
            if any(bad in G for bad in badfiltered):
                G = ''
            if any(bad in S for bad in badfiltered):
                S = ''
        ID = col[boldid].strip()
        GB = col[gbid].strip()
        if args.require_genbank:
            if GB: #if there is a GB accession
                if 'Pending' in GB:
                    continue
                else:
                    pass
            else:
                continue
        Seq = col[seqid].replace('-', '')
        Seq = re.sub('N*$', '', Seq)
        Seq = re.sub('^N*', '', Seq)
        tax = []
        for i in [K,P,C,O,F,G,S]:
            if not i.endswith(':'):
                tax.append(i)
        tax_fmt = ','.join(tax)
        if tax_fmt.endswith(',') or tax_fmt.endswith(', '):
            tax_fmt = tax_fmt.rsplit(',',1)[0]
        BIN = col[bin].strip()
        if BIN == '':
            continue
        if not BIN in bins_seen:
            TrimSeq = TrimPrimer(Seq, args.primer, args.trim_length, args.primer_mismatch)
            if TrimSeq:
                bins_seen[BIN] = (tax_fmt, GB, TrimSeq)
                UniqueBIN += 1
            else:
                TooShort +=1
        else:
            DupBIN += 1
            oldtax = bins_seen.get(BIN)[0]
            oldlevels = oldtax.count(',')
            if tax_fmt.count(',') > oldlevels:
                TrimSeq = TrimPrimer(Seq, args.primer, args.trim_length, args.primer_mismatch)
                if TrimSeq:
                    bins_seen[BIN] = (tax_fmt, GB, TrimSeq)

with open(args.out, 'w') as output:
    for k,v in natsorted(bins_seen.items()):
        if v[1] != '':
            output.write('>%s_%s;tax=%s\n%s\n' % (k, v[1], v[0], v[2]))
        else:
            output.write('>%s;tax=%s\n%s\n' % (k, v[0], v[2]))

print "%i total records processed" % Total
print "%i records duplicate BIN" % DupBIN
print "%i seq too short (%s bp)" % (TooShort, args.trim_length)
print "%i records written to %s" % (UniqueBIN, args.out)

