#!/usr/bin/env python

#script to filter down a genbank tbl annotation file that match a list

import sys, argparse, re

parser=argparse.ArgumentParser(
    description='''Script that does filters an NCBI feature file (.tbl) based on scaffolds/contigs in a list''',
    epilog="""Jon Palmer (2016) nextgenusfs@gmail.com""")
parser.add_argument('-i', '--input', required=True, help='Genbank .tbl annotation file(Required)')
parser.add_argument('-o','--out', required=True, help='Output .tbl file (filtered)')
parser.add_argument('-l','--list', required=True, help='List of scaffolds (csv)')
parser.add_argument('-r','--rename', action='store_true', help='Rename scaffolds with 2nd column in list')
parser.add_argument('-p','--print', dest='show',  action='store_true', help='Print to STDOUT')
parser.add_argument('-a','--annotation', action='store_true', help='only scaffolds with annotation')
args=parser.parse_args()

def group_by_heading( some_source ):
    buffer= []
    for line in some_source:
        if line.startswith( ">" ):
            if buffer: yield buffer
            buffer= [ line ]
        else:
            buffer.append( line )
    yield buffer

#parse list of scaffolds and put in dictionary if two columns
filter = []
rename = {}
with open(args.list, 'rU') as inputlist:
    for line in inputlist:
        line = line.replace('\n', '')
        cols = line.split(',')
        filter.append(cols[0])
        if not cols[0] in rename:
            if len(cols) > 1:
                rename[cols[0]] = cols[1]
            else:
                print("Error: you did not pass a 2 column CSV file")
                sys.exit(1)

#turn list into regex
find = re.compile(r'\b(?:%s)\b' % '|'.join(filter))

#now parse genbank tbl file
with open(args.out, 'w') as output:
    with open(args.input, 'rU') as tbl:
        for line in group_by_heading(tbl):
            if find.search(line[0]):
                if args.rename:
                    match_name = find.search(line[0]).group(0)
                    line[0] = line[0].replace(match_name, rename.get(match_name))
                if not args.show:
                    if args.annotation:
                        if len(line) > 3:
                            output.write(''.join(line))
                    else:
                        output.write(''.join(line))
                else:
                    if args.annotation:
                        if len(line) > 3:
                            sys.stdout.write(''.join(line))
                    else:
                        sys.stdout.write(''.join(line))

