#!/bin/bash

#Script for HMMer searching the EggNog4.1 database and annotating the protein families	

if [ -z "$1" ]; then
	echo "Run this script as follows:
	$0 yourproteins.fasta hmmer_base_name Evalue num_cpus
	Example:  $0 Anid.proteins.fasta Anid.eggnog 1e-10 6"
	exit
else
	#first run Hmmer, saving sig domains less than 1e-10
	hmmscan --cpu $4 --domtblout $2.out -E $3 $HMMDB/fuNOG.hmm $1
	#parse results with awk, just get values we care about, length 50% - 150% of target, and then get best match
	grep -v "^#" $2.out | gawk '{if($3/$6 >= 0.5 && $3/$6 <= 1.5) print $1"\t"$4"\t"$7"\t"$3"\t"$6"\t"$3/$6;}' | gsort -k2,2 -k3,3g | gsort -u -k2,2 > $2.filtered
	#now grab descriptions
	gawk 'BEGIN {FS=OFS="\t"} NR==FNR {v[$1]=$3;next} {print $2,$1": "v[$1],$3,$4,$5,$6}' $HMMDB/fuNOG.mapping.txt $2.filtered > $2.results.txt
	#clean up
	rm $2.out
	rm $2.filtered
	echo "Finished.  Ran HMMscan against fuNOG database with Evalue of $3, and then filtered results."
fi
