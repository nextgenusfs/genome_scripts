#!/bin/bash
#script to draw karyotype map after rearrangement with CEA10 and Af293 anchors

if [ -z "$1" ]; then
    echo "Usage: redraw_karyotype.sh isolate_name"
    exit
fi

isolate=$1
cd $isolate

#just need to update the seqids to draw all chromosomes, layout is same as last time you drew
cut -f1 Af293.bed | gsort | uniq | paste -d',' -s - > seqids
cut -f1 ReMap.bed | gsort | uniq | paste -d',' -s - >> seqids
cut -f1 CEA10.bed | gsort | uniq | paste -d',' -s - >> seqids

#then run script
python -m jcvi.graphics.karyotype seqids layout

#rename
mv karyotype.pdf $isolate.all_contigs.pdf