#!/bin/bash

# script to format JGI GFF and Fasta files to GBK
# Jon Palmer: nextgenusfs@gmail.com

#run test if no files given, print usage
if [ -z "$1" ]; then
	echo "Run this script as follows:
	$0 <GFF file> <Fasta file> <KOG tab> <Genome Name>
GFF file downloaded from JGI, Genome FASTA file, and KOG tab file"
	exit
else
    #check to see if all programs are installed and in path
    command -v gsed >/dev/null 2>&1 || { echo "I require gsed but it's not installed.  Aborting." >&2; exit 1; }
    command -v gag.py >/dev/null 2>&1 || { echo "I require gag.py but it's not installed.  Aborting." >&2; exit 1; }
    command -v tbl2asn >/dev/null 2>&1 || { echo "I require tbl2asn but it's not installed.  Aborting." >&2; exit 1; }
    command -v blastp >/dev/null 2>&1 || { echo "I require blastp but it's not installed.  Aborting." >&2; exit 1; }
    command -v gsort >/dev/null 2>&1 || { echo "I require gsort but it's not installed.  Aborting." >&2; exit 1; }
    command -v gawk >/dev/null 2>&1 || { echo "I require gawk but it's not installed.  Aborting." >&2; exit 1; }

    #get directory of script to find genbank template
    script_dir=$(command -v format_jgi_kog.sh | gsed 's,/format_jgi_kog.sh,,g')

    #get current directory for writing results to
    dir=$PWD
    #make folder
    mkdir -p $4
    #get name of fasta headers from human input, get first 5 headers and print out.
    header_count=$(grep -c "^>" $2)
    header=$(grep "^>" $2 | head -n5 | gsed 's/>//g')
    printf "\nThere are $header_count fasta files: first 5 headers look like this:\n\n$header\n\n"
    printf "Do you need to edit headers? (you need to simplify complicated headers) [y/n]: "
    read header_question
    if [ $header_question == "n" ]; then
        cp $2 $4/genome.fasta
        gsed 's/jgi.p|//g' $1 | gsed 's/|/_/g' > $4/genome.gff
    else
        printf "Enter leading portion of header to remove (copy and paste works well): "
        read header_trim
        #now reformat fasta headers and save
        gsed "s/>$header_trim/>/g" $2 > $4/genome.fasta
        #print back new headers
        new_header=$(grep "^>" $4/genome.fasta | head -n5 | gsed 's/>//g')
        echo "New fasta headers look like this:

$new_header"
        #now reformat GFF to match fasta headers
        gsed "s/^$header_trim//g" $1 | gsed 's/jgi.p|//g' | gsed 's/|/_/g' > $4/genome.gff
    fi
    #now get annotations from KOG file, first get example gene name from gff
    name=$(grep $'\tgene\t' $4/genome.gff | head -n1 | gsed 's/;/\t/g' | gsed 's/Name=//g' | cut -f10)
    name_stem=${name/_*/}
    #now need to make a 2 column "association" file
    grep $'\tmRNA\t' $1 | cut -f9 | gsed 's/ID=//g' | gsed 's/;/\t/g' | gsed 's/proteinId=//g' | gawk -F"\t" -v n="$name_stem" '{ print n"_"$4"\t"$1;}' > $4/annotations.map
    grep -v "#transcriptId" $3 | gawk -F"\t" -v n="$name_stem" '{ print n"_"$2"\t""product""\t"$4;}' > $4/annotations.tmp
    gawk -F"\t" 'NR==FNR{a[$1]=$2;next}$1 in a{print a[$1]"\t"$2"\t"$3;}' $4/annotations.map $4/annotations.tmp > $4/annotations.txt
    #try running command line gag
    cd $dir/$4
    echo "------------------------"
    echo "Running GAG."
    echo "------------------------"
    gag.py -f genome.fasta -g genome.gff -a annotations.txt
    cd $dir/$4/gag_output
    cp genome.fasta genome.fsa
    echo "------------------------"
    echo "Running tbl2asn."
    echo "------------------------"
    tbl2asn -p . -t $script_dir/test.sbt -M n -Z discrep -a r10u -l paired-ends -j "[organism=$4]" -V b -c fx
    cp genome.gbf $dir/$4.gbk
    cd $dir
    
fi