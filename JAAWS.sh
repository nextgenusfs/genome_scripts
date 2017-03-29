#!/bin/bash

#Jon's Awesome Annotation Wrapper Script (JAAWS)

#point here is to run trimmomatic to trim PE Illumina reads, then Spades, followed by some automated cleanup steps
DIR=$1
R1=$2
R2=$3
QR1=qtrimR1pe.fq.gz
QR2=qtrimR2pe.fq.gz
QCR1=clean.1.fq.gz
QCR2=clean.2.fq.gz
SG1=spades.scaffolds.fasta
CPU=$4
PHYLUM=Ascomycota
OUTPUT=$5

#move into directory, make sure reads are passed as full paths....
mkdir -p $1
cd $1

if [ ! -f clean.1.fq.gz ]; then
	#run trimmomatic, trim a few bp off 5' end and remove adatpers, don't quality trim
	trimmomatic PE -threads $CPU $R1 $R2 $QR1 qtrimR1se.fq.gz $QR2 qtrimR2se.fq.gz ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 HEADCROP:8 MINLEN:36

	#map to phix genome, remove any paired reads that mapß
	wget -c --tries=0 --read-timeout=20 ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/Enterobacteria_phage_phiX174_sensu_lato_uid14015/NC_001422.fna
	bowtie2-build NC_001422.fna phiX
	#map to phiX and recover unmapped paired reads
	bowtie2 -x phiX -p $CPU -1 $QR1 -2 $QR2 --un-conc clean.fq > /dev/null 2>&1
	#re-compress output PE files
	bgzip -@ $CPU clean.1.fq
	bgzip -@ $CPU clean.2.fq
fi

if [ ! -f spades.scaffolds.clean2.fasta ]; then
	#check for spades folder, run on continue if found
	if [ ! -d spades ]; then
		#Now run spades
		spades.py -k 21,33,55,77,99,127 -t $CPU -m 64 --careful -o spades --pe1-1 $QCR1 --pe1-2 $QCR2
	else
		spades.py -k 21,33,55,77,99,127 --continue -t $CPU -m 64 --careful -o spades --pe1-1 $QCR1 --pe1-2 $QCR2
	fi

	#grab scaffolds from output folder
	cp spades/scaffolds.fasta $SG1

	#mapped reads to assembly using BWA
	bwa index $SG1
	bwa mem -t $CPU $SG1 $QCR1 $QCR2 | samtools view -@ $CPU -bS - > spades.bam

	#resulting assembly was blasted against NCBI nt database
	blastn -task megablast -query $SG1 -db $BLASTDB/nt -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -culling_limit 5 -num_threads $CPU -evalue 1e-25 -out spades.nt.megablast.out

	#run Blobtools to look for contamination
	blobtools create -i $SG1 -y spades -t spades.nt.megablast.out -b spades.bam
	blobtools view -i blobDB.json
	grep -v '^#' blobDB.table.txt > blob.summary.txt
	blobtools blobplot -i blobDB.json

	#calculate average coverage of largest 50 scaffolds
	AVGCOV=$(egrep "$PHYLUM|no-hit" blob.summary.txt | head -n 50 | gawk '{ sum += $7; n++ } END { if (n > 0) print sum / n; }')
	MINCOV=$(bc <<< $AVGCOV*0.05)
	MINCOV=${MINCOV%.*}
	grep -v '^#' blob.summary.txt | gawk -v cov=$MINCOV '{ if($7 < cov) print $1}' > blob.low-coverage.txt

	#parse the blob summary output and filter for taxonomy and coverage, get rid of low coverage scaffolds, coverage less than 5% of average
	grep -v '^#' blob.summary.txt | egrep -v "$PHYLUM|no-hit" | gawk '{print$1;}' > blob.contamination.txt

	#now run the NCBI contamination screen
	contamination_screen.sh $SG1 ncbi_contamination.txt
	grep -v '^#' ncbi_contamination.txt | cut -f1 > ncbi_remove-list.txt

	#combine all lits of scaffolds to remove
	cat blob.contamination.txt blob.low-coverage.txt ncbi_remove-list.txt | sort | uniq > all.remove.txt

	#remove contaminant scaffolds
	fasta_remove.py $SG1 all.remove.txt > spades.scafffolds.clean1.fasta

	#now run funannotate clean, dropping scaffolds less than 1 kb
	funannotate clean -i spades.scafffolds.clean1.fasta -o spades.scaffolds.clean2.fasta -m 1000
fi

if [ ! -f $OUTPUT.changes.txt ]; then
	#Clean up assembly using Pilon, 5 iterations should be sufficient
	mkdir -p pilon
	cd pilon
	cp ../spades.scaffolds.clean2.fasta .
	#number of loops is hardcoded to 5 runs of pilon
	for i in {1..5}; do
		if [[ $i == '1' ]]; then
			GENOME=spades.scaffolds.clean2.fasta
		else
			x=$(($i - 1))
			GENOME=pilon$x.fasta
		fi
		if [[ ! -f pilon$i.fasta ]]; then
			echo "Working on $GENOME"
			bwa index $GENOME
			bwa mem -t $CPU $GENOME ../$QCR1 ../$QCR2 | samtools view -@ $CPU -bS - | samtools sort -@ $CPU -o pilon$i.bam -
			samtools index pilon$i.bam
			#now run pilon
			pilon --genome $GENOME --frags pilon$i.bam --output pilon$i --threads $CPU --changes
		else
			echo "$GENOME has already been run, moving onto next iteration"
		fi
	done

	cp pilon5.fasta pilon_final.fasta
	cp pilon5.changes ../$OUTPUT.changes.txt
	cd ..
	#sort by size and rename fasta headers
	funannotate sort -i pilon/pilon_final.fasta -o $OUTPUT
fi

#do a little cleanup as these files take up a lot of space, you don't need all these pilon BAM files, also remove the intermediate cleanup up fastqs
rm MITO*
rm CONTAM_PROKS*
rm CONTAM_EUKS*
rm -r pilon
rm $QR1
rm $QR2

