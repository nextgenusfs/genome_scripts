#!/bin/bash

#Jon's Awesome Annotation Wrapper Script (JAAWS)

if [ -z "$5" ]; then
    echo "Usage: JAAWS_ion.sh Output_Dir IonTorrent_reads NUM_CPUS PHYLUM OUTPUT REF_Genome(optional)"
    exit
fi

#point here is to run Spades, followed by some automated cleanup steps
DIR=$1
R1=$PWD/$2
SG1=spades.contigs.fasta
CPU=$3
PHYLUM=$4
OUTPUT=$5
REF=$PWD/$6

#move into directory, make sure reads are passed as full paths....
mkdir -p $1
cd $1

#since adapters and Quality trimming is done by default on Ion Server, move directly to spades

#check for spades folder, run on continue if found
if [ ! -d spades ]; then
	#Now run spades
	spades.py --iontorrent -t $CPU -m 64 --careful -o spades -s $R1
else
	spades.py --iontorrent --continue -t $CPU -m 64 --careful -o spades -s $R1
fi

#grab contigs from output folder
if [ -n "$6" ]; then
	cp spades/contigs.fasta $SG1
else
	cp spades/scaffolds.fasta $SG1
fi

if [ ! -f spades.bam ]; then
	#mapped reads to assembly using BWA
	bwa index $SG1
	bwa mem -t $CPU $SG1 $R1 | samtools view -@ $CPU -bS - > spades.bam
fi

if [ ! -f spades.nt.megablast.out ]; then
	#resulting assembly was blasted against NCBI nt database
	blastn -task megablast -query $SG1 -db $BLASTDB/nt -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -culling_limit 5 -num_threads $CPU -evalue 1e-25 -out spades.nt.megablast.out
fi

if [ ! -f blob.summary.txt ]; then
	#run Blobtools to look for contamination
	blobtools create -i $SG1 -y spades -t spades.nt.megablast.out -b spades.bam
	blobtools view -i blobDB.json
	grep -v '^#' blobDB.table.txt > blob.summary.txt
	blobtools blobplot -i blobDB.json
fi

if [ ! -f spades.clean3.fasta ]; then
	#calculate average coverage of largest 50 scaffolds
	AVGCOV=$(egrep "$PHYLUM|no-hit" blob.summary.txt | head -n 50 | gawk '{ sum += $7; n++ } END { if (n > 0) print sum / n; }')
	MINCOV=$(bc <<< $AVGCOV*0.05)
	MINCOV=${MINCOV%.*}
	MAXCOV=$(bc <<< $AVGCOV*10)
	grep -v '^#' blob.summary.txt | gawk -v cov=$MINCOV '{ if($7 < cov) print $1}' > blob.low-coverage.txt

	#get high coverage, > 10X average
	grep -v '^#' blob.summary.txt | gawk -v cov=$MAXCOV '{ if($7 > cov) print $1}' > blob.high-coverage.txt

	#parse the blob summary output and filter for taxonomy and coverage, get rid of low coverage scaffolds, coverage less than 5% of average
	grep -v '^#' blob.summary.txt | egrep -v "$PHYLUM|no-hit" | gawk '{print$1;}' > blob.contamination.txt

	#now run the NCBI contamination screen
	contamination_screen.sh $SG1 ncbi_contamination.txt
	grep -v '^#' ncbi_contamination.txt | cut -f1 > ncbi_remove-list.txt

	#combine all lits of scaffolds to remove
	cat blob.contamination.txt blob.low-coverage.txt blob.high-coverage.txt ncbi_remove-list.txt | sort | uniq > all.remove.txt

	#remove contaminant scaffolds
	fasta_remove.py $SG1 all.remove.txt > spades.clean1.fasta

	#now run funannotate clean, dropping scaffolds less than 1 kb and if they are contained in another scaffold
	funannotate clean -i spades.clean1.fasta -o spades.clean2.fasta -m 1000

	#rename for cleaner downstream processing
	funannotate sort -i spades.clean2.fasta -b contig -o spades.clean3.fasta
fi

#run ragout if a reference genome is given
if [ -n "$6" ]; then
	echo -e ".references = ref
.target = $DIR

ref.fasta = $REF
$DIR.fasta = ../spades.clean3.fasta

.blocks = small
.naming_ref = ref" > ragout.spades.rcp
	#now run ragout
	ragout.py -o ragout --overwrite --repeats --threads $CPU ragout.spades.rcp
	#combine unplaced and scaffolded outputs
	cat ragout/$DIR\_scaffolds.fasta ragout/$DIR\_unplaced.fasta > ragout.combined.fasta
	funannotate sort -i ragout.combined.fasta -o $DIR.cleaned.fasta
else
	cp spades.clean3.fasta $DIR.cleaned.fasta
fi	

if [ ! -f $OUTPUT.changes.txt ]; then
	#Clean up assembly using Pilon, 10 iterations should be sufficient
	mkdir -p pilon
	cd pilon
	cp ../$DIR.cleaned.fasta .
	#number of loops is hardcoded to 5 runs of pilon
	for i in {1..10}; do
		if [[ $i == '1' ]]; then
			GENOME=$DIR.cleaned.fasta
		else
			x=$(($i - 1))
			GENOME=pilon$x.fasta
		fi
		if [[ ! -f pilon$i.fasta ]]; then
			echo "Working on $GENOME"
			bwa index $GENOME
			bwa mem -t $CPU $GENOME $R1 | samtools view -@ $CPU -bS - | samtools sort -@ $CPU -o pilon$i.bam -
			samtools index pilon$i.bam
			#now run pilon
			pilon --genome $GENOME --unpaired pilon$i.bam --output pilon$i --threads $CPU --changes
		else
			echo "$GENOME has already been run, moving onto next iteration"
		fi
	done

	cp pilon10.fasta pilon_final.fasta
	cp pilon10.changes ../$OUTPUT.changes.txt
	cd ..
	#sort by size and rename fasta headers
	funannotate sort -i pilon/pilon_final.fasta -o $OUTPUT
fi

#do a little cleanup as these files take up a lot of space, you don't need all these pilon BAM files, also remove the intermediate cleanup up fastqs
rm MITO*
rm CONTAM_PROKS*
rm CONTAM_EUKS*
rm -r pilon
