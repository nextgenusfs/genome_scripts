#!/bin/bash

# script to format FungiDB GFF and Fasta files to GBK and SMURF input
# Jon Palmer: palmer.jona@gmail.com

#run test if no files given, print usage
if [ -z "$1" ]; then
	echo "Run this script as follows:
	$0 <GFF file> <Fasta file> <Genome Name>"
	exit
else
    #check to see if all programs are installed and in path
    command -v gsed >/dev/null 2>&1 || { echo "I require gsed but it's not installed.  Aborting." >&2; exit 1; }
    command -v gag.py >/dev/null 2>&1 || { echo "I require gag.py but it's not installed.  Aborting." >&2; exit 1; }
    command -v tbl2asn >/dev/null 2>&1 || { echo "I require tbl2asn but it's not installed.  Aborting." >&2; exit 1; }
    command -v gb2proteins.py >/dev/null 2>&1 || { echo "I require gb2proteins.py but it's not installed.  Aborting." >&2; exit 1; }
    command -v gsort >/dev/null 2>&1 || { echo "I require gsort but it's not installed.  Aborting." >&2; exit 1; }
    command -v gawk >/dev/null 2>&1 || { echo "I require gawk but it's not installed.  Aborting." >&2; exit 1; }
    command -v gb2smurf.py >/dev/null 2>&1 || { echo "I require gb2smurf.py but it's not installed.  Aborting." >&2; exit 1; }
    command -v python >/dev/null 2>&1 || { echo "I require python but it's not installed.  Aborting." >&2; exit 1; }
    echo "------------------------"
    echo "Formatting input files."
    echo "------------------------"
    #get current directory for writing results to
    dir=$PWD
    #make folder
    mkdir -p $3
    #get name of fasta headers to grep gff file with as they have fasta sequence file at bottom
    header=$(gsed 's/ |.*//g' $2 | head -n1 | gsed 's/>//g')
    if [[ $header == *"_SC"* ]]; then
        echo "Fasta header contains _SC:  ($header)"
        echo "Running extra steps to format correctly"
        gsed 's/ |.*//g' $2 | gsed 's/>.*_SC/>SC/g' > $3/genome.fasta
        short_head=${header:0:2}
        #first need to reformat GFF - alter some problematic regions
        gsed 's/^.*_SC/SC/g' $1 | grep -v $'\tsupercontig\t' | grep "^SC.*\tFungi" | gsed 's/;/;\t/g' | gawk -F"\t" 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9$13;}' | gsed 's/web_id.*$//g' > $3/genome.gff
    else
        gsed 's/ |.*//g' $2 > $3/genome.fasta
        short_head=${header:0:3}
        #first need to reformat GFF - alter some problematic regions
        grep "^$short_head" $1 | grep -v $'\tsupercontig\t' | gsed 's/;/;\t/g' | gawk -F"\t" 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9$13;}' | gsed 's/web_id.*$//g' > $3/genome.gff
    fi    
    #now get descriptions for annotation file
    grep "^$short_head" $1 | grep -v $'\tsupercontig\t' | gsed 's/;/;\t/g' | grep $'\tgene\t' | gawk -F"\t" '{ print $9"\t""product""\t"$11;}' | gsed 's/;//g' | gsed 's/ID=/rna_/g' | gsed 's/\tproduct/-1\tproduct/g' | gsed 's/description=//g' | gsed 's/+/ /g'  > $3/annotation.tmp
    #remove strange formatting and remove empty lines as gag chokes on that
    python -c 'import sys, urllib; print urllib.unquote(sys.stdin.read())' < $3/annotation.tmp | gsed '/^$/d' > $3/annotations.txt
    #try running command line gag
    cd $dir/$3
    echo "Running GAG."
    echo "------------------------"
    gag.py fasta=genome.fasta gff=genome.gff anno=annotations.txt
    cd $dir/$3/gag_output
    cp genome.fasta genome.fsa
    echo "------------------------"
    echo "Running tbl2asn."
    echo "------------------------"
    tbl2asn -p . -t $HOME/test.sbt -M n -Z discrep -a r10u -l paired-ends -j "[organism=$3]" -V b -c fx
    cp genome.gbf $dir/$3.gbk
    cd $dir
    echo "------------------------"
    echo "Running gb2smurf.py"
    echo "------------------------"
    gb2smurf.py $3.gbk -p $3.proteins.fasta -g $3.scaffolds.fasta -s $3.smurf.txt
    echo "
    SMURF output is complete, script finished!
    
    Results:
        $3.gbk
        $3.proteins.fasta
        $3.scaffolds.fasta
        $3.smurf.txt
        
    Now you can submit $3.proteins.fasta and $3.smurf.txt to http://jcvi.org/smurf/upload.php
    
    And you can submit $3.gbk to AntiSmash (http://antismash.secondarymetabolites.org)
    
    You likely want to delete the folder $dir/$3, unless you need the intermediate files for something else.
    
    "
     echo "------------------------"
fi