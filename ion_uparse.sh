#!/bin/bash

#########################################################
#    Script for OTU clustering/taxonomy using UPARSE	#
#########################################################

#You will need the following software installed
#USEARCH8 - should be named usearch8, i.e. with softlink
#Robert's python scripts - in the scripts folder
#RDP Classfier 2.10.1 - also included in the folder
#GNUplot - install with sudo apt-get or brew install gnuplot
#maker script: fasta_tool
#GNU sed - install with brew install gsed gawk
#GNU awk
#GNU grep


#global variables, set for your installation
python_scripts=$HOME/scripts;
rdp_directory=$HOME/scripts;
DB_directory=$HOME/projects/DB;

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
echo "---------------------------------------------------------------------"
echo "Welcome to script for ITS filtering/clustering for Ion Torrent Data." 
echo "It looks like you are in this directory: $PWD"
echo "---------------------------------------------------------------------"
echo "Is this the directory that you will be working in? [y/n]:  \c"
read question1
if [ $question1 == y ]; then
	d=$PWD
else
	echo "Enter directory (leave off last backslash): " 
	read d
fi
echo "---------------------------------------------------------------------"
echo "Ok, we've set working directory to: $d"
mkdir -p $d/tmp
echo "---------------------------------------------------------------------"
echo "Enter a base name for this project (i.e. FACE):  \c"
read project
echo "---------------------------------------------------------------------"
echo "Enter the base name(s) of your FASTQ file(s) (leave off .fastq): \c"
read -a fastq_array
echo "${fastq_array[@]}" > $d/tmp/file_names
echo "---------------------------------------------------------------------"
echo "The default behavior of this script is to find all 45 PGM barcodes from your
dataset.  This is often what you want to do, however, you should only specifiy the
barcodes that you've actually used, so in the next question you can either answer 'all'
which will specifiy all 45 barcodes, 'add' which will prompt you to type in which 
barcodes that you have used, or 'subtract' which will prompt you to type in which 
of the 45 barcodes you want to remove, i.e. this is easier than typing 40 of the 
barcodes that you used.  You simply need to list the barcodes by number."
echo "---------------------------------------------------------------------"
echo "Will you be processing multiple PGM chips together? [y/n] \c"
read question4
if [ $question4 == y ]; then
	echo "---------------------------------------------------------------------"
	echo "Ok, now we will find barcodes, primers, and then trim reads to a desired length"
	echo "Did you use fITS7 forward primer? (AGTGARTCATCGAATCTTTG) [y/n]:  \c"
	read question3
	if [ $question3 == y ]; then
		fwd_primer=AGTGARTCATCGAATCTTTG
	else
		echo "Enter forward primer sequence:  "
		read fwd_primer
	fi
	echo "---------------------------------------------------------------------"
	echo "Did you use ITS4 reverse primer? (TCCTCCGCTTATTGATATGC) [y/n]:  \c"
	read question3
	if [ $question3 == y ]; then
		rev_primer=TCCTCCGCTTATTGATATGC
	else
		echo "Enter reverse primer sequence:  "
		read rev_primer
	fi
	echo "---------------------------------------------------------------------"
	echo "What length would you like to trim the reads to (recommended is 175 bp):  \c"
	read trim_length
	echo "---------------------------------------------------------------------"
	echo "What expected error would you like to quality trim data (recommended < 0.5):  \c"
	read quality_score
	simulated_readarray files $d/tmp/file_names
	for files in ${files[@]} 
		do
			file_bytes=$(stat -f '%z' $d/$files.fastq)
			file_mega_bytes=$((file_bytes / 10**6))
			num_reads=$(grep -c "^@" $files.fastq)
			echo "---------------------------------------------------------------------"
			echo "Now processing: $files.fastq"
			echo "file size: $file_mega_bytes MB"
			echo "Number of reads: $num_reads"
			echo "Enter a name for this chip:  \c"
			read chip_name
			echo "---------------------------------------------------------------------"
			echo "Which barcodes did you use: [all, add, subtract]  \c"
			read question2
			if [ $question2 == all ]; then
				cp $DB_directory/pgm_barcodes.fa $d/tmp/pgm_barcodes.fa
			fi
			if [ $question2 == add ]; then
				echo "Ok, list the barcodes that you used (i.e. 1 4 5 6 7):  "
				read -a barcode_list
				echo "${barcode_list[@]}" > $d/tmp/barcode_list
				gsed -i 's/ /\n/g' $d/tmp/barcode_list
				gsed -i 's/^/BC_/g' $d/tmp/barcode_list
				fasta_tool --select $d/tmp/barcode_list $DB_directory/pgm_barcodes.fa > $d/tmp/pgm_barcodes.fa
			fi
			if [ $question2 == subtract ]; then
				echo "Ok, list the barcodes that you want to remove  (i.e. 1 4 5 6 7):  "
				read -a barcode_list
				echo "${barcode_list[@]}" > $d/tmp/barcode_list
				gsed -i 's/ /\n/g' $d/tmp/barcode_list
				gsed -i 's/^/BC_/g' $d/tmp/barcode_list
				fasta_tool --remove $d/tmp/barcode_list $DB_directory/pgm_barcodes.fa > $d/tmp/pgm_barcodes.fa
			fi
			echo "---------------------------------------------------------------------"
			echo "Now running demultiplexing script, this may take some time."
			python $python_scripts/fastq_strip_barcode_relabel7.py $files.fastq $fwd_primer $rev_primer $d/tmp/pgm_barcodes.fa Reads_ $chip_name 50 $trim_length > $d/$chip_name\_$trim_length\_trim.fastq 
			
			usearch8 -fastq_filter $d/$chip_name\_$trim_length\_trim.fastq -fastq_qmax 45 -fastq_maxee $quality_score -fastaout $d/tmp/$chip_name\_$trim_length\_filter.fa.tmp
		done
	find tmp -name '*.tmp' -exec cat {} \; > $d/$project\_$trim_length\_ee$quality_score.fa
else
	simulated_readarray files $d/tmp/file_names
	echo "---------------------------------------------------------------------"
	echo "Which barcodes did you use: [all, add, subtract]  \c"
	read question2
	if [ $question2 == all ]; then
		cp $DB_directory/pgm_barcodes.fa $d/tmp/pgm_barcodes.fa
	fi
	if [ $question2 == add ]; then
		echo "Ok, list the barcodes that you used (i.e. 1 4 5 6 7):  "
		read -a barcode_list
		echo "${barcode_list[@]}" > $d/tmp/barcode_list
		gsed -i 's/ /\n/g' $d/tmp/barcode_list
		gsed -i 's/^/BC_/g' $d/tmp/barcode_list
		fasta_tool --select $d/tmp/barcode_list $DB_directory/pgm_barcodes.fa > $d/tmp/pgm_barcodes.fa
	fi
	if [ $question2 == subtract ]; then
		echo "Ok, list the barcodes that you want to remove  (i.e. 1 4 5 6 7):  "
		read -a barcode_list
		echo "${barcode_list[@]}" > $d/tmp/barcode_list
		gsed -i 's/ /\n/g' $d/tmp/barcode_list
		gsed -i 's/^/BC_/g' $d/tmp/barcode_list
		fasta_tool --remove $d/tmp/barcode_list $DB_directory/pgm_barcodes.fa > $d/tmp/pgm_barcodes.fa
	fi
	echo "---------------------------------------------------------------------"
	echo "Ok, now we will find barcodes, primers, and then trim reads to a desired length"
	echo "Did you use fITS7 forward primer? (AGTGARTCATCGAATCTTTG) [y/n]:  \c"
	read question3
	if [ $question3 == y ]; then
		fwd_primer=AGTGARTCATCGAATCTTTG
	else
		echo "Enter forward primer sequence:  "
		read fwd_primer
	fi
	echo "---------------------------------------------------------------------"
	echo "Did you use ITS4 reverse primer? (TCCTCCGCTTATTGATATGC) [y/n]:  \c"
	read question3
	if [ $question3 == y ]; then
		rev_primer=TCCTCCGCTTATTGATATGC
	else
		echo "Enter reverse primer sequence:  "
		read rev_primer
	fi
	echo "---------------------------------------------------------------------"
	echo "What length would you like to trim the reads to (recommended is 175 bp):  \c"
	read trim_length
	echo "---------------------------------------------------------------------"
	echo "What expected error would you like to quality trim data (recommended < 0.5):  \c"
	read quality_score
	echo "---------------------------------------------------------------------"
	echo "Now running demultiplexing script, this may take some time."
	python $python_scripts/fastq_strip_barcode_relabel6.py $files.fastq $fwd_primer $rev_primer $d/tmp/pgm_barcodes.fa Reads_ 50 $trim_length > $d/$files\_$trim_length\_trim.fastq
	usearch8 -fastq_filter $d/$files\_$trim_length\_trim.fastq -fastq_qmax 45 -fastq_maxee $quality_score -fastaout $d/$project\_$trim_length\_ee$quality_score.fa
fi
echo "---------------------------------------------------------------------"
echo "Now running dereplication using USEARCH8"
usearch8 -derep_fulllength $d/$project\_$trim_length\_ee$quality_score.fa -sizeout -fastaout $d/$project\_$trim_length\_ee$quality_score\_derep.fa
echo "---------------------------------------------------------------------"
echo "Now removing singletons using USEARCH8"
usearch8 -sortbysize $d/$project\_$trim_length\_ee$quality_score\_derep.fa -minsize 2 -fastaout $d/$project\_$trim_length\_ee$quality_score\_sort.fa
echo "---------------------------------------------------------------------"
echo "Now clustering at 97% using USEARCH8"
usearch8 -cluster_otus $d/$project\_$trim_length\_ee$quality_score\_sort.fa -sizein -sizeout -relabel OTU_ -otus $d/$project\_$trim_length\_ee$quality_score\_otus1.fa
echo "---------------------------------------------------------------------"
echo "Now removing sequencing padding"
gsed 's/N.*$//g' $d/$project\_$trim_length\_ee$quality_score\_otus1.fa > $d/$project\_$trim_length\_ee$quality_score\_otus1_clean.fa
echo "---------------------------------------------------------------------"
echo "Now filtering chimeric sequences using UCHIME-ITS2"
usearch8 -uchime_ref $d/$project\_$trim_length\_ee$quality_score\_otus1_clean.fa -db $DB_directory/uchime_ITS2.fasta -strand plus -nonchimeras $d/$project\_$trim_length\_ee$quality_score\_otus.fa
#Run Blast search here to remove 'no blast hit' OTUs
echo "---------------------------------------------------------------------"
echo "Now mapping reads back to OTUs using usearch_global"
usearch8 -usearch_global $d/$project\_$trim_length\_ee$quality_score.fa -db $d/$project\_$trim_length\_ee$quality_score\_otus.fa -id 0.97 -strand plus -uc $d/$project\_map.uc
echo "---------------------------------------------------------------------"
echo "Now creating OTU table"
python $python_scripts/uc2otutab.py $d/$project\_map.uc > $d/$project\_otu_table.txt
echo "---------------------------------------------------------------------"
echo "Converting to BIOM format"
biom convert --table-type="OTU table" -i $d/$project\_otu_table.txt -o $d/$project\_otu_table.biom --to-hdf5
echo "---------------------------------------------------------------------"
echo "Assigning taxonomy with UTAX"
usearch8 -utax $d/$project\_$trim_length\_ee$quality_score\_otus.fa -db $DB_directory/utax/unite_its2.fa -taxconfs $DB_directory/utax/its2.tc -tt $DB_directory/utax/unite.tt -utaxout $d/$project\_utax.txt
echo "---------------------------------------------------------------------"
echo "Assigning taxonomy with RDP UNITE training set"
java -Xmx6000m -jar $rdp_directory/rdp_classifier_2.10.1/dist/classifier.jar classify -g fungalits_unite -f fixrank -o $d/$project\_RDP_unite.txt $d/$project\_$trim_length\_ee$quality_score\_otus.fa
#reformat classfication file
gsed 's/\tdomain\t/\t/g' $d/$project\_RDP_unite.txt | gsed 's/\tphylum\t/\t/g' | gsed 's/\tclass\t/\t/g' | gsed 's/\torder\t/\t/g' | gsed 's/\tfamily\t/\t/g' | gsed 's/\tgenus\t/\t/g' | gsed 's/\tspecies\t/\t/g' | gsed 's/\t\tFungi/\tFungi/g' | gawk -F"\t" '{ print $1"\t"$2"("$3");"$4"("$5");"$6"("$7");"$8"("$9");"$10"("$11");"$12"("$13");"$14"("$15")" ; }' > $d/$project\_RDP_unite_taxonomy.txt
echo "---------------------------------------------------------------------"
echo "Assigning taxonomy with RDP Warcup training set"
java -Xmx6000m -jar $rdp_directory/rdp_classifier_2.10.1/dist/classifier.jar classify -g fungalits_warcup -f fixrank -o $d/$project\_RDP_warcup.txt $d/$project\_$trim_length\_ee$quality_score\_otus.fa
#reformat classfication file
gsed 's/\tdomain\t/\t/g' $d/$project\_RDP_warcup.txt | gsed 's/\tphylum\t/\t/g' | gsed 's/\tclass\t/\t/g' | gsed 's/\torder\t/\t/g' | gsed 's/\tfamily\t/\t/g' | gsed 's/\tgenus\t/\t/g' | gsed 's/\tspecies\t/\t/g' | gsed 's/\t\tFungi/\tFungi/g' | gawk -F"\t" '{ print $1"\t"$2"("$3");"$4"("$5");"$6"("$7");"$8"("$9");"$10"("$11");"$12"("$13");"$14"("$15")" ; }' > $d/$project\_RDP_warcup_taxonomy.txt
echo "---------------------------------------------------------------------"
echo "Assigning taxonomy with iterative Blast search"
echo "---------------------------------------------------------------------"
otus=$project\_$trim_length\_ee$quality_score\_otus.fa
num_threads=12
total_otus=$(ggrep -c "^>" $d/$otus)
date
echo "Ok, it looks like you have $total_otus sequences to Blast.
Now first running MegaBlast on $d/$otus 
against UNITE database"
mkdir $d/blast_tmp
ggrep ">" $d/$otus | gsed 's/>//g' > $d/blast_tmp/all_otus_ids.txt
blastn -query $d/$otus -db $DB_directory/BlastDB/UNITE -max_target_seqs 1 -task megablast -num_threads $num_threads -out $d/blast_tmp/first_pass.tab -outfmt "6 qseqid qlen sseqid slen length pident evalue"
echo "---------------------------------------------------------------------"
date	
echo "First search against UNITE db is finished, now filtering results..."
gawk -F"\t" '{ if($5/$2>=0.99 && $6>=97.00) print $0; }' $d/blast_tmp/first_pass.tab | cut -f1 > $d/blast_tmp/first_pass.ids
gawk -F"\t" '{ if($5/$4>=0.99 && $6>=97.00) print $0; }' $d/blast_tmp/first_pass.tab | cut -f1 >> $d/blast_tmp/first_pass.ids
#add semi-colon back onto ids, as it is getting stripped by blast
gsed -i 's/$/;/g' $d/blast_tmp/first_pass.ids
fasta_tool --remove $d/blast_tmp/first_pass.ids $d/$otus > $d/blast_tmp/otus_first_pass.fa
otus_first_pass=$(ggrep -c "^>" $d/blast_tmp/otus_first_pass.fa)
found_first=$(wc -l $d/blast_tmp/first_pass.ids | gawk '{print $1;}')
echo "---------------------------------------------------------------------"
date
echo "I found $found_first hit(s) greater than 97% from UNITE db
Now searching remaining $otus_first_pass sequences against UNITE+INSDC database"
blastn -query $d/blast_tmp/otus_first_pass.fa -db $DB_directory/BlastDB/UNITE_INSDC -max_target_seqs 1 -num_threads $num_threads -task megablast -out $d/blast_tmp/second_pass.tab -outfmt "6 qseqid qlen sseqid slen length pident evalue"
echo "---------------------------------------------------------------------"
date	
echo "Second search against UNITE+INSDC db is finished, now filtering results..."
gawk -F"\t" '{ if($5/$2>=0.99 && $6>=97.00) print $0; }' $d/blast_tmp/second_pass.tab | cut -f1 > $d/blast_tmp/second_pass.ids
gawk -F"\t" '{ if($5/$4>=0.99 && $6>=97.00) print $0; }' $d/blast_tmp/second_pass.tab | cut -f1 >> $d/blast_tmp/second_pass.ids
#add semi-colon back onto ids, as it is getting stripped by blast
gsed -i 's/$/;/g' $d/blast_tmp/second_pass.ids
fasta_tool --remove $d/blast_tmp/second_pass.ids $d/blast_tmp/otus_first_pass.fa > $d/blast_tmp/otus_second_pass.fa
otus_second_pass=$(ggrep -c "^>" $d/blast_tmp/otus_second_pass.fa)
found_second=$(wc -l $d/blast_tmp/second_pass.ids | gawk '{print $1;}')
echo "---------------------------------------------------------------------"
date
echo "I found $found_second hits greater than 97% from UNITE+INSDC
Now searching remaining $otus_second_pass sequences against NCBI nt database"
blastn -query $d/blast_tmp/otus_second_pass.fa -remote -db nt -max_target_seqs 1 -task megablast -out $d/blast_tmp/third_pass.tab -outfmt "6 qseqid qlen sseqid slen length pident evalue stitle"
echo "---------------------------------------------------------------------"
date	
echo "Third search against NCBI nt db is finished, now cleaning and combining results"
#if total alignment length is not at least 95%, run again?
gawk -F"\t" '{ if($5/$2>=0.97 && $6>=97.00) print $0; }' $d/blast_tmp/third_pass.tab | cut -f1 > $d/blast_tmp/third_pass.ids
#add semi-colon back onto ids, as it is getting stripped by blast
gsed -i 's/$/;/g' $d/blast_tmp/third_pass.ids
fasta_tool --remove $d/blast_tmp/third_pass.ids $d/blast_tmp/otus_second_pass.fa > $d/blast_tmp/otus_third_pass.fa
otus_third_pass=$(ggrep -c "^>" $d/blast_tmp/otus_third_pass.fa)
found_third=$(wc -l $d/blast_tmp/third_pass.ids | gawk '{print $1;}')
echo "I found $found_third hits greater than 97% from NCBI nt.  There are still
$otus_third_pass that have not been identified at the 97% threshold.  I will now go back
through the data and pull out those we can identify to lower levels of certainty."
echo "---------------------------------------------------------------------"
date
#total re-write, this is unesecarily complicated
#take blast hits and assign db name first
gawk -F"\t" '{ print $0"\t""UNITE"; }' $d/blast_tmp/first_pass.tab > $d/blast_tmp/first_pass.blast
gawk -F"\t" '{ print $0"\t""UNITE-INSDC"; }' $d/blast_tmp/second_pass.tab > $d/blast_tmp/second_pass.blast
gawk -F"\t" '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""NCBI nt""\t"$8; }' $d/blast_tmp/third_pass.tab > $d/blast_tmp/third_pass.blast.tax
#get taxonomy for UNITE and UNITE-INSDC
gawk -F"\t" 'NR==FNR {h[$3] = $2; next} {print $0"\t"h[$3]}' $DB_directory/iterative_tax.txt $d/blast_tmp/first_pass.blast > $d/blast_tmp/first_pass.blast.tax
gawk -F"\t" 'NR==FNR {h[$3] = $2; next} {print $0"\t"h[$3]}' $DB_directory/iterative_tax.txt $d/blast_tmp/second_pass.blast > $d/blast_tmp/second_pass.blast.tax
#filter, name, sort in one command to get highest Evalue and then rename via identity settings and 98% of length
#split into two files, one greater than, one less than 97% over length of sequence
cat $d/blast_tmp/first_pass.blast.tax $d/blast_tmp/second_pass.blast.tax $d/blast_tmp/third_pass.blast.tax | gsort -V -k1,1 -k7,7g | gsort -u -V -k1,1 | gawk -F"\t" '{ if($5/$2>=0.97) print $0; }' > $d/blast_tmp/blast_pass_97.txt
cat $d/blast_tmp/first_pass.blast.tax $d/blast_tmp/second_pass.blast.tax $d/blast_tmp/third_pass.blast.tax | gsort -V -k1,1 -k7,7g | gsort -u -V -k1,1 | gawk -F"\t" '{ if($5/$2<0.97) print $0; }' > $d/blast_tmp/blast_fail_97.txt
#now rename failing ones
gawk -F"\t" '{ if($6>=97.00) print $0"\t""*species*"; if($6>=90.00 && $6<97.00) print $0"\t""*genus*"; if($6>=85.00 && $6<90.00) print $0"\t""*family*"; if($6>=80.00 && $6<85.00) print $0"\t""*order*"; if($6>=75.00 && $6<80.00) print $0"\t""*class*"; if($6<75.00) print $0"\t""*not significant*"; }' $d/blast_tmp/blast_fail_97.txt > $d/blast_tmp/fail_97.formatted
#rename passing
gawk -F"\t" '{ if($6>=97.00) print $0"\t""species"; if($6>=90.00 && $6<97.00) print $0"\t""genus"; if($6>=85.00 && $6<90.00) print $0"\t""family"; if($6>=80.00 && $6<85.00) print $0"\t""order"; if($6>=75.00 && $6<80.00) print $0"\t""class"; if($6<75.00) print $0"\t""not significant"; }' $d/blast_tmp/blast_pass_97.txt > $d/blast_tmp/pass_97.formatted
#combine and sort again
cat $d/blast_tmp/pass_97.formatted $d/blast_tmp/fail_97.formatted | gsort -t $'\t' -V -k1,1 -k10,10r | gsort -t $'\t' -u -V -k1,1 > $d/blast_tmp/all_results.txt
#print header to file
echo "QueryID	QueryLength	SubjectID	SubjectLength	AlignmentLength	PercentID	Evalue	Database	Taxonomy	IdLevel" > $d/blast_tmp/header.txt
#check to see if we didn't get any blast hits from some OTUs
cut -f1 $d/blast_tmp/all_results.txt | uniq > $d/blast_tmp/final_hits.txt
ggrep -Fwvf $d/blast_tmp/final_hits.txt $d/blast_tmp/all_otus_ids.txt > $d/blast_tmp/no_blast_hit.list
fasta_tool --select $d/blast_tmp/no_blast_hit.list $d/$otus > $d/otus_no_blast_hit.fa
#make taxonomy file to append to OTU table
gawk -F"\t" '{ print $1"\t"$9" ""("$10"-"$6"%)"; }' $d/blast_tmp/all_results.txt | gsed 's/18S.*(/(/g' | gsed 's/5\.8S.*(/(/g' | gsed 's/internal.*(/(/g' | gsed 's/partial.*(/(/g' | gsed 's/small.*(/(/g' > $d/blast_tmp/Blast_tax.txt
gawk -F"\t" '{ print $1"\t""No blast hit"; }' $d/blast_tmp/no_blast_hit.list > $d/blast_tmp/no_hits.txt
cat $d/blast_tmp/Blast_tax.txt $d/blast_tmp/no_hits.txt | gsort -V -k1,1 > $d/$project\_Blast_taxonomy_table.txt
species=$(ggrep -c $'\tspecies' $d/blast_tmp/all_results.txt)
genus=$(ggrep -c $'\tgenus' $d/blast_tmp/all_results.txt)
family=$(ggrep -c $'\tfamily' $d/blast_tmp/all_results.txt)
order=$(ggrep -c $'\torder' $d/blast_tmp/all_results.txt)
class=$(ggrep -c $'\tclass' $d/blast_tmp/all_results.txt)
no_sig=$(ggrep -c $'\t\*not significant\*' $d/blast_tmp/all_results.txt)
species1=$(ggrep -c $'\t\*species\*' $d/blast_tmp/all_results.txt)
genus1=$(ggrep -c $'\t\*genus\*' $d/blast_tmp/all_results.txt)
family1=$(ggrep -c $'\t\*family\*' $d/blast_tmp/all_results.txt)
order1=$(ggrep -c $'\t\*order\*' $d/blast_tmp/all_results.txt)
class1=$(ggrep -c $'\t\*class\*' $d/blast_tmp/all_results.txt)
no_sig1=$(ggrep -c $'\t\*not significant\*' $d/blast_tmp/all_results.txt)
no_hits=$(wc -l $d/blast_tmp/no_blast_hit.list | gawk '{print $1;}')
species2="$((species + species1))"
genus2="$((genus + genus1))"
family2="$((family + family1))"
order2="$((order + order1))"
class2="$((class + class1))"
no_sig2="$((no_sig + no_sig1))"
echo "---------------------------------------------------------------------"
echo "IterativeOTU Blast is finished, here are some stats on your run:
	
	Total sequences blasted:	$total_otus
	
					Align > 97%	Align < 97%	Total
	Hits to species level (97%):	   $species		   $species1	 	 $species2
	Hits to genus level (90%):	   $genus		   $genus1	 	 $genus2
	Hits to family level (85%):	   $family		   $family1	 	 $family2
	Hits to order level (80%):	   $order		   $order1	  	 $order2
	Hits to class level (75%):	   $class		   $class1	  	 $class2
	Low quality hits (<75%):	   $no_sig		   $no_sig1	  	 $no_sig2
	No blast hits at all:		   $no_hits

	Sequences without blast hits are located in: 
	$d/otus_no_blast_hit.fa

	Blast taxonomy file for appending to OTU table is located in:
	$d/Blast_taxonomy_table.txt

	Final blast results are located in: 
	$d/IterativeOTU_blast.txt" 2>&1 | tee $d/IterativeOTU_run_stats.txt;
cat $d/blast_tmp/header.txt $d/blast_tmp/all_results.txt > $d/IterativeOTU_blast.txt
echo "---------------------------------------------------------------------"
date
echo "---------------------------------------------------------------------"










