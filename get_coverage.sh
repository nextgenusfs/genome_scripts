#!/bin/bash

if [ -z "$2" ]; then
    echo "Usage: get_coverage.sh list_bam_files.txt location.bed "
    exit
fi

#now loop through each file, make bam index and run bedtools coverage

echo "Sample,Coverage,MatingType"
while read line; do if [[ ! -f "$line.bai" ]];then samtools index $line;fi; echo -n "$line,"; mat=$(bedtools multicov -bams $line -bed $2 | cut -f5); echo -n $mat","; if [[ $mat > 10 ]];then echo -n "MAT1-1"; else echo -n "MAT1-2";fi; echo ""; done < $1
