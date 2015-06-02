#!/bin/bash

#This script filters repetetive elements generated from RepeatModeler, TransposonPSI, and RepARK

#############set up TransposonPSI run for Pseudogymnoascus species#############################
##start
genome=$1;
#make a temp directory to store all intermediate files
mkdir -p tmp;

#run blastp to get protein matches to atypical transposons.
echo "Running blastp against atypical transposons"
echo "---------------------------------------------"
blastp -num_threads 5 -query ../$genome.proteins.fasta -db /usr/local/TransposonPSI_08222010/transposon_ORF_lib/transposonDB -outfmt "6 qseqid sseqid evalue pident slen qlen" -seg yes -evalue 1e-5 -max_target_seqs 1 -out tmp/$genome.blastp.transposonDB
echo "Now filtering blast results"
echo "---------------------------------------------"
#Now get coordinates from file and print DNA sequences
cut -f1 tmp/$genome.blastp.transposonDB > tmp/$genome.transposonDB.txt
cat tmp/$genome.transposonDB.txt | sort | uniq | gsed 's/-T1//g' > tmp/$genome.PSI.list.txt
grep -f tmp/$genome.PSI.list.txt ../$genome.gff | grep $'\tgene\t' | gsed 's/locus_tag=/locus_tag=\t/g' | gawk -F"\t" '{ if($7 == "+") print $10"\t"$1":"$4"-"$5; else print $10"\t"$1":"$5"-"$4;}' > tmp/$genome\_prot_coords.txt
#format genome scaffold file so doesn't cause problems
mv ../$genome.scaffolds.fasta ../$genome.scaffolds.fasta.bak
fasta_tool --print --wrap 60 ../$genome.scaffolds.fasta.bak > ../$genome.scaffolds.fasta
#Run perl script to print out fasta file of TransposonPSI DNA sequence
get_fasta_from_gff.pl -i tmp/$genome\_prot_coords.txt -g ../$genome.scaffolds.fasta > tmp/$genome.prot.transposonPSI.fa

############################ get transposonPSI nucleotide DNA sequences in fasta format ################################
gsed 's/ID=/ID=\t/g' $genome.scaffolds.fasta.TPSI.allHits.chains.bestPerLocus.gff3 | gsed 's/; Target=/; Target=\t/g' | gsed 's/;/\t;/g' | gawk -F"\t" '{ if ($7 == "+") print $10"|"$12"\t"$1":"$4"-"$5; else print $10"|"$12"\t"$1":"$5"-"$4; }' > tmp/$genome\_nuc_coords.txt 
get_fasta_from_gff.pl -i tmp/$genome\_nuc_coords.txt -g ../$genome.scaffolds.fasta > tmp/$genome.nuc.transposonPSI.fa
############################  END ####################

#rename RepeatModeler folder to something more reasonable and cuz we are going to make another one
find . -depth -type d -name 'RM_*' -exec mv {} RepeatModeler \;

#get RpMod data into temp folder
cp RepeatModeler/consensi.fa.classified tmp/RpMod.fa

#get RepArk data into temp folder
cp $genome\_repark.fa tmp/$genome\_repark.fa

echo "Comparing/Combining repeats from the different software"
echo "---------------------------------------------"
#now want to start with repeatmodeler and determine any overlap with transposonPSI
makeblastdb -in tmp/RpMod.fa -dbtype nucl -title $genome\_rpmod -out tmp/$genome\_rpmod
blastn -query tmp/$genome.nuc.transposonPSI.fa -db tmp/$genome\_rpmod -outfmt "6 qseqid sseqid evalue pident slen qlen" -dust yes -evalue 1e-150 -max_target_seqs 1 -out tmp/$genome\_PSI.blast
#filter blast results to 100% ID
gawk '{ if ($4 >= 99.00) print $1; }' tmp/$genome\_PSI.blast > tmp/grep_PSI_list.txt
cat tmp/RpMod.fa tmp/$genome.nuc.transposonPSI.fa > tmp/$genome.all.PSI.fa
fasta_tool --remove tmp/grep_PSI_list.txt tmp/$genome.all.PSI.fa > tmp/$genome.all.filtered.PSI.fa
#now make new blastdb with filtered PSI results
makeblastdb -in tmp/$genome.all.filtered.PSI.fa -input_type fasta -dbtype nucl -title $genome\_all_PSI -out tmp/$genome\_all
#now blast RpMod against PSI database
blastn -query tmp/$genome\_repark.fa -db tmp/$genome\_all -outfmt "6 qseqid sseqid evalue pident slen qlen" -dust yes -evalue 1e-150 -max_target_seqs 1 -out tmp/$genome\_RepARK.blast
#filter results to a lower threshold because this was done de novo and could have a few errors
awk '{ if ($4 >= 97.00) print $1; }' tmp/$genome\_RepARK.blast > tmp/grep_RepARK_list.txt
cat tmp/$genome.all.filtered.PSI.fa tmp/$genome\_repark.fa > tmp/$genome.all.PSI.RpMod.RepARK.fa
fasta_tool --wrap 80 --remove tmp/grep_RepARK_list.txt tmp/$genome.all.PSI.RpMod.RepARK.fa > tmp/$genome.all.repeats.fa
total_repeats=$(grep -c "^>" tmp/$genome.all.repeats.fa)
echo "---------------------------------------------"
echo "Done. We've run RepeatModeler, TransposonPSI, and RepARK.  Then we filtered, cleaned, and combined these results for $genome.  This has resulted in identification of $total_repeats repetitive elements and now I will run RepeatClassifier to classify these sequences."
echo "---------------------------------------------"
#Run RepeatClassifier on combined elements
perl /usr/local/Cellar/repeatmodeler/1.0.8/RepeatClassifier -engine ncbi -consensi tmp/$genome.all.repeats.fa
mv tmp/$genome.all.repeats.fa.classified $genome\_repeats_classified.fa
echo "---------------------------------------------"
#Now lets take these repeats and run RepeatMasker to determine the composition in the genome
RepeatMasker -pa 4 -dir . -lib $genome\_repeats_classified.fa ../$genome.scaffolds.fasta
echo "Repeat Pipeline is finished."
echo "---------------------------------------------"
