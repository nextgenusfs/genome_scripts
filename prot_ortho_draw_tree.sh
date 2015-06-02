#!/bin/bash
#script for taking proteinorthov5.11 results, pulling orthologs, randomizing to desired number, concatenating all sequences together, and then drawing phylogentic tree

simulated_readarray ()
{
somearray="$1"
local file="$2"

[ -f "$file" ] || return 1

local i=0
while read line; do
eval $somearray[$i]='"$line"'
local i=$((i+1))
done <$file

return 0
}

if [ -z "$1" ]; then
	echo "Run this script as follows:
	$0 test.poff number_prots outnamestem cpus"
	exit
else

	#Script for running and formatting proteinortho results into a concatenated alignment for phylogeny
	#run like sh protortho_draw_tree.sh test.poff 100 outnamestem cpus
	num_proteins=$2

	#make directory to store fasta files
	mkdir -p fasta_storage
	mkdir -p prot_tmp

	grep '^#' $1 | cut -d$'\t' -f 4- | gsed 's/.faa/ /g' | gsed 's/\t//g' > prot_tmp/column_names.txt
	simulated_readarray species prot_tmp/column_names.txt
	num=$(tr -cd ' ' <prot_tmp/column_names.txt | wc -c | gawk '{print $1;}')
	echo "I've detected $num species in this file and they are:
	${species[@]}
	"
	#ask a single question about outgrouping
	echo "Enter species name for outgroup in RaXML, or type none for no outgroup"
	read outgroup
	#how many bootstraps
	echo "How many bootstrap replicates do you want to run?"
	read boot
	#process the prot ortho file to save space
	gawk 'NR == 1 {print}' $1 > subsample.poff
	grep $"^$num\t$num\t1" $1 >> subsample.poff

	#then run following in directory to pull fasta files
	perl /usr/local/proteinortho_v5.11/tools/grab_proteins.pl subsample.poff

	#move .fasta files into folder
	mv *.fasta fasta_storage/

	echo "Ok, now finding desired number of random orthologous groups of proteins."

	#count number of orthologs
	orthologs=$(grep $"^$num\t$num\t1" subsample.poff | wc -l | gawk '{print $1;}')

	#determine number of random orthologs to use.
	if [ "$orthologs" -lt "$num_proteins" ]; then
		echo "$num_proteins is out of range, setting to largest value of $orthologs"
		num_proteins=$orthologs
	else
		echo "Test passed: $num_proteins is in range"
	fi

	#select random from the first column of fasta protein names
	grep $"^$num\t$num\t1" subsample.poff | cut -f4 | gsort -R | head -n $num_proteins > prot_tmp/$num_proteins.random.txt
	#now search through directory above to pull out file names
	grep -Rlf prot_tmp/$num_proteins.random.txt fasta_storage |  gsed 's/fasta_storage\///g' | gsed 's/.fasta//g' > prot_tmp/$num_proteins.random.list.txt

	#now run a loop to rename
	while read name; do gsed "s/^>/>group_"$name"_/g" fasta_storage/$name.fasta; done <prot_tmp/$num_proteins.random.list.txt > prot_tmp/$num_proteins.random.renamed.fa

	#now need to split apart each one based on the header name from the .proteinortho file.
	for species in ${species[@]}
	do
		fasta_tool --grep_header $species prot_tmp/$num_proteins.random.renamed.fa | fasta_tool --seq_only | tr -d '\n' | gsed "s/^/>$species\n/g" | fasta_tool --print --wrap 60 | gsed 's/*//g' > prot_tmp/$species.$num_proteins
	done

	#now put into one fasta file
	cat prot_tmp/*.$num_proteins > $3.fasta

	#run mafft on this file to generate an alignment
	mafft $3.fasta > $3.align.fasta

	#Run TrimAl to remove gaps from alignment
	trimal -in $3.align.fasta -out $3.trimal.fasta -automated1

	#Convert to PHYLIP
	seqmagick convert --input-format fasta --output-format phylip $3.trimal.fasta $3.trimal.phylip

	#Run fasttree
	fasttree < $3.trimal.phylip > $3.trimal.fasttree.nwk
	#run raxml
	if [ $outgroup == none ]; then
		raxmlHPC-PTHREADS -T $4 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -\# $boot -s $3.trimal.phylip -n $3.raxml.nwk
	else
		raxmlHPC-PTHREADS -T $4 -f a -m PROTGAMMAAUTO -o $outgroup -p 12345 -x 12345 -\# $boot -s $3.trimal.phylip -n $3.raxml.nwk
	fi

	echo "You are finished, yay!"
	#clean up your mess
	rm subsample.poff
	rm -R fasta_storage/
fi





