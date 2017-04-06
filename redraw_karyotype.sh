#!/bin/bash
#script to draw karyotype map after rearrangement with CEA10 and Af293 anchors

if [ -z "$1" ]; then
    echo "Usage: redraw_karyotype.sh isolate_name"
    exit
fi

isolate=$1
cd $isolate

#script to draw the SM clusters


#just need to update the seqids to draw all chromosomes, layout is same as last time you drew
cut -f1 Af293.bed | gsort | uniq | paste -d',' -s - > seqids
cut -f1 ReMap.bed | gsort | uniq | paste -d',' -s - >> seqids
cut -f1 CEA10.bed | gsort | uniq | paste -d',' -s - >> seqids

#then run script
python -m jcvi.graphics.karyotype seqids layout

#rename
mv karyotype.pdf $isolate.all_contigs.pdf

#Microsynteny
#python -m jcvi.compara.synteny mcscan grape.bed grape.peach.lifted.anchors --iter=1 -o grape.peach.i1.blocks
python -m jcvi.compara.synteny mcscan Af293.bed Af293.ReMap.lifted.anchors --iter=1 -o Af293.ReMap.i1.blocks
python -m jcvi.compara.synteny mcscan ReMap.bed ReMap.CEA10.lifted.anchors --iter=1 -o ReMap.CEA10.i1.blocks

#need to select region of interest
gawk '/AFUA_3G02400/ {p=1}; p; /AFUA_3G02770/ {p=0}' Af293.ReMap.i1.blocks > Af293.blocks
#gawk '/AFUB_045860/ {p=1}; p; /AFUB_045570/ {p=0}' ReMap.CEA10.i1.blocks
gawk '/AFUC_03224/ {p=1}; p; /AFUC_03260/ {p=0}' ReMap.CEA10.i1.blocks > CEA10.blocks
python -m jcvi.graphics.synteny blocks grape_peach.bed blocks.layout