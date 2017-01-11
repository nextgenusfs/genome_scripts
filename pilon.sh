#!/bin/bash

#iterative pilon

#genome in $1, forward reads in $2, reverse reads in $3, cpus in $4, final output in $5
mkdir -p pilon
cd pilon
cp $1 .
#number of loops is hardcoded to 10 runs of pilon
for i in {1..10}; do
	if [[ $i == '1' ]]; then
		GENOME=$1
	else
		x=$(($i - 1))
		GENOME=pilon$x.fasta
	fi
	if [[ ! -f pilon$1.fasta ]]; then
		echo "Working on $GENOME"
		bwa index $GENOME
		bwa mem -t $4 $GENOME $2 $3 | samtools view -@ $4 -bS - | samtools sort -@ $4 -o pilon$i.bam -
		samtools index pilon$i.bam
		#now run pilon
		pilon --genome $GENOME --frags pilon$i.bam --output pilon$i --threads $4 --changes
	else
		echo "$GENOME has already been run, moving onto next iteration"
	fi
done

cp pilon10.fasta ../$5
cp pilon10.changes ../$5.changes.txt
