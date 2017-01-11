#!/bin/bash

GENOME=$1
CONTAM=$2

#drop these scaffolds, now down to 14, now need to run NCBI contaminant screen
wget -c --tries=0 --read-timeout=20 ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/contam_in_euks.fa.gz
gunzip contam_in_euks.fa.gz
wget -c --tries=0 --read-timeout=20 ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/contam_in_prok.fa
wget -c --tries=0 --read-timeout=20 ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/mito.nt.gz
gunzip mito.nt.gz
makeblastdb -in contam_in_euks.fa -dbtype nucl -title NCBI_contam -out CONTAM_EUKS
makeblastdb -in contam_in_prok.fa -dbtype nucl -title NCBI_contam -out CONTAM_PROKS
makeblastdb -in mito.nt -dbtype nucl -out MITO


#screen euks
echo "#EUKARYOTIC Contamination Screen" >> $CONTAM
echo "#--------------------------------" >> $CONTAM
blastn -query $GENOME -db CONTAM_EUKS -num_threads 6 -dust yes -soft_masking true -perc_identity 90.0 -lcase_masking -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | gawk '($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)' > euk.screen.txt
cat euk.screen.txt | grep -v "^#" >> $CONTAM
echo "#--------------------------------" >> $CONTAM

#screen proks
echo "#PROKARYTOIC Contamination Screen" >> $CONTAM
echo "#--------------------------------" >> $CONTAM
blastn -query $GENOME -db CONTAM_PROKS -num_threads 6 -dust yes -soft_masking true -perc_identity 90.0 -lcase_masking -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | gawk '($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)' > prok.screen.txt
cat prok.screen.txt | grep -v "^#" >> $CONTAM
echo "#--------------------------------" >> $CONTAM

#now check for mitochondrial contamination
echo "#Mitochondrial Contamination Screen" >> $CONTAM
echo "#--------------------------------" >> $CONTAM
blastn -query $GENOME -num_threads 6 -db MITO -dust yes -soft_masking true -perc_identity 98.6 -outfmt 7 | gawk '$4>=120' >mito.screen.txt
cat mito.screen.txt | grep -v "^#" >> $CONTAM
echo "#--------------------------------" >> $CONTAM

