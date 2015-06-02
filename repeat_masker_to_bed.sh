#!/bin/bash

#Run like RM_to_bed.sh options repeatmasker.out fileout.bed

if [ $1 == "all" ];then
	#reformat RepeatMasker .out file, write to temp.bed
	tail -n +4 $2 | gsed 's/^\s*//' | gsed -E 's/\s+/\t/g' | grep -v 'Low_complexity' | cut -f5-7 | gawk 'BEGIN { OFS="\t" } { print $1, $2-1, $3}' | gsort -k1,1 -k2,2n | mergeBed > $3
	counts=$(wc -l $3 | gawk '{print $1;}')
	echo "Script is complete, bed file is composed of $counts regions."
fi

if [ $1 == "nosimple" ];then
	#reformat RepeatMasker .out file, write to temp.bed
	tail -n +4 $2 | gsed 's/^\s*//' | gsed -E 's/\s+/\t/g' | grep -v 'Low_complexity' | grep -v 'Simple_repeat' | cut -f5-7 | gawk 'BEGIN { OFS="\t" } { print $1, $2-1, $3}' | gsort -k1,1 -k2,2n | mergeBed > $3
	counts=$(wc -l $3 | gawk '{print $1;}')
	echo "Script is complete, bed file is composed of $counts regions."
fi

if [ $1 == "onlysimple" ];then
	#reformat RepeatMasker .out file, write to temp.bed
	tail -n +4 $2 | gsed 's/^\s*//' | gsed -E 's/\s+/\t/g' | grep -v 'Low_complexity' | grep 'Simple_repeat' | cut -f5-7 | gawk 'BEGIN { OFS="\t" } { print $1, $2-1, $3}' | gsort -k1,1 -k2,2n | mergeBed > $3
	counts=$(wc -l $3 | gawk '{print $1;}')
	echo "Script is complete, bed file is composed of $counts regions."
fi
