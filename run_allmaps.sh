#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: run_allmaps.sh isolate_name GeneName"
    exit
fi

isolate=$1
name=$2

mkdir -p $isolate
cd $isolate

#symlink the necessary files into this working folder
ln -s /Volumes/LinuxHD/Afumi24genomes/assemblies/$isolate\_annotate/predict_results/aspergillus_fumigatus.gff3 $isolate.gff
ln -s /Volumes/LinuxHD/Afumi24genomes/assemblies/$isolate\_annotate/predict_results/aspergillus_fumigatus.scaffolds.fa $isolate.scaffolds.fa
ln -s /Volumes/LinuxHD/Afumi24genomes/assemblies/$isolate\_annotate/predict_results/aspergillus_fumigatus.transcripts.fa $isolate.transcripts.fa
ln -s /usr/local/test_jcvi/Af293.bed .
ln -s /usr/local/test_jcvi/Af293.cds .

#get BED file from GFF
python -m jcvi.formats.gff bed --type=mRNA --key=Parent -o $isolate.bed $isolate.gff

#get FASTA file from transcripts
python -m jcvi.formats.fasta format --sep=" " --nodesc $isolate.transcripts.fa $isolate.cds

#generate overlaps
python -m jcvi.compara.catalog ortholog Af293 $isolate --cscore=.99
python -m jcvi.compara.synteny screen --minspan=30 --simple Af293.$isolate.anchors Af293.$isolate.anchors.new

#create seqids files
cut -f1 Af293.bed | gsort | uniq | paste -d',' -s - > seqids
cut -f1 $isolate.bed | gsort | uniq | paste -d',' -s - >> seqids

#create layout file
gsed "s/ReMap/$isolate/g" /usr/local/test_jcvi/layout_chr3 > layout

#draw first graphic
python -m jcvi.graphics.karyotype seqids layout
mv karyotype.pdf $isolate.af293.original.pdf

#convert anchors to synteny
python -m jcvi.assembly.syntenypath bed Af293.$isolate.anchors --switch --scale=100000 -o $isolate.anchors.bed

#one could use multiple bed files, but need in correct format for allmaps
python -m jcvi.assembly.allmaps mergebed $isolate.anchors.bed -o $isolate.allmaps.bed

#now run allmaps
python -m jcvi.assembly.allmaps path $isolate.allmaps.bed $isolate.scaffolds.fa

#liftover annotations
liftOver -gff $isolate.gff $isolate.allmaps.chain $isolate.updated.gff $isolate.unplaced.gff

#filter GFF file for renaming gene models in new assembly
grep -v $'_codon\t' $isolate.updated.gff | gsed 's/scaffold_/unplaced_/g' | python -m jcvi.formats.gff sort - > $isolate.clean.gff
reformat_seqs.py $isolate.allmaps.fasta $isolate.final.fasta $isolate.order.txt
/usr/local/opt/funannotate/libexec/util/maker_map_ids.pl --prefix $name\_ --justify 5 --suffix -T --iterate 1 --sort_order $isolate.order.txt $isolate.clean.gff > mapping.ids
/usr/local/opt/funannotate/libexec/util/map_gff_ids.pl mapping.ids $isolate.clean.gff

#now run GAG development version on resulting data
/usr/local/GAG/gag.py -f $isolate.final.fasta -g $isolate.clean.gff -o gag --fix_start_stop
cp gag/genome.fasta gag/genome.fsa
tbl2asn -p gag -t /usr/local/opt/funannotate/libexec/lib/test.sbt -M n -Z discrep -j "[organism=Aspergillus fumigatus] [isolate=$isolate]" -V b -c fx
cp gag/genome.gbf $isolate.gbk

#create tabbed output
gb2tab.py -i $isolate.gbk -o $isolate

#get new anchors for mapping to Chr3 region
python -m jcvi.formats.gff bed --type=mRNA --key=Parent -o ReMap.bed $isolate.clean.gff
python -m jcvi.formats.fasta format --sep=" " --nodesc $isolate.transcripts.fasta ReMap.cds
python -m jcvi.compara.catalog ortholog Af293 ReMap --cscore=.99
python -m jcvi.compara.synteny screen --minspan=30 --simple Af293.ReMap.anchors Af293.ReMap.anchors.new
#also map to CEA10
cp /usr/local/test_jcvi/CEA10.bed .
cp /usr/local/test_jcvi/CEA10.cds .
python -m jcvi.compara.catalog ortholog ReMap CEA10 --cscore=.99
python -m jcvi.compara.synteny screen --minspan=30 --simple ReMap.CEA10.anchors ReMap.CEA10.anchors.new
cp /usr/local/test_jcvi/layout_3way layout
cp /usr/local/test_jcvi/seqids_3way seqids
python -m jcvi.graphics.karyotype seqids layout
mv karyotype.pdf $isolate.chr3.pdf
