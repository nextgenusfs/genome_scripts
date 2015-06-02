#!/bin/bash

#################################################
#                                               #
#	The genome annotation script                #
#                                               #
#################################################

#source bashrc
source ~/.bashrc

simplifyFasta ()
{
echo "############################################################################################"
echo "Now Re-naming Fasta Headers to avoid problems"
simplifyFastaHeaders.pl $d/$genome.fasta Supercontig_ $d/$genome.simple.fa $d/header.genome.map
}

TrinityPGMorIllumina ()
{
#Launch VirtualBox now so it is ready by the time you get there
vboxmanage startvm UbuntuServer --type headless
max_mem=50
echo "############################################################################################"
echo -n "Are your RNA-seq reads from Illumina or Ion Torrent? [illumina or ion]  "
read illumina_or_ion
if [ $illumina_or_ion == illumina ]; then
	echo -n "Are your RNA-seq reads single or paired-end [single/pe]?  "
	read illumina_question
	if [ $illumina_question == pe ]; then
		echo -n "Enter the file name of your Illumina Paired-End RNA-seq data, leave off othe _R1.fastq & R2.fastq: "
    	read RNA_reads
    	mv $d/$RNA_reads\_R1.fastq $d/$genome.RNA_R1.fastq
    	mv $d/$RNA_reads\_R2.fastq $d/$genome.RNA_R2.fastq
    	echo -n "Enter number of bp to trim from 5' end: "
    	read five_trim
    	echo -n "Enter the maximum intron length: "
    	read MAX_INTRON_LENGTH
    	echo -n "Enter the distance between PE reads (default is 200 bp):  "
    	read mate_distance
    	echo -n "How many CPUs do you want to use for Trinity? (max=10, recommended=6)  "
    	read trinity_cpus
    	#trim reads
    	trim_galore -o $d -q 20 --phred33 --paired --trim1 --clip_R1 $five_trim --clip_R2 $five_trim $d/$genome.RNA_R1.fastq $d/$genome.RNA_R2.fastq
    	#rename validated reads
    	mv $d/$genome.RNA_R1_val_1.fq $d/$genome.RNA_val_R1.fq
    	mv $d/$genome.RNA_R2_val_2.fq $d/$genome.RNA_val_R2.fq
    	#gzip the originals to save space, do this in background to keep script moving, should only take 1 CPU each and limited RAM
    	nohup gzip -1 $d/$genome.RNA_R1.fastq &
    	nohup gzip -1 $d/$genome.RNA_R2.fastq &
    	#estimate the memory required for Jellyfish/Trinity, i.e. 1 GB per 1 million reads
    	num_reads=$(grep -c "^@" $d/$genome.RNA_val_R1.fq)
    	num_reads2="$((num_reads * 2))"  #multiply by two because of paired end reads
    	memory_req="$(($num_reads2 / 1000000 + 1))"
    	if [ "$memory_req" -gt "$max_mem" ]; then
        	memory_req=50
    	else
        	memory_req=$memory_req
    	fi
    	#run trinity de novo
    	TRINITY=/home/ngs/programs/trinityrnaseq_r20140717
   		directory=/mnt/shared/$d
   		echo "############################################################################################"
    	echo -e "Ok, setting Jellyfish memory to $memory_req GB and running Trinity de novo. "
		#Now send ssh Trinity command to VirtualBox
		ssh -p 2222 ngs@localhost "$TRINITY/Trinity --seqType fq --JM $memory_req\G --output /mnt/shared/$d/trinity --jaccard_clip --left $d/$genome.RNA_val_R1.fq --right $d/$genome.RNA_val_R2.fq --CPU $trinity_cpus"
    	echo "############################################################################################"
   		echo "Trinity de novo has finished successfully, running stats and starting Trinity Genome Guided."	
   		#First run tophat alignment of RNA to genome
    	cd $d
    	mkdir -p trinity_gg
    	#symlinks won't work because path is different from host to guest OS, so move files (clean it up later)
    	mv $genome.RNA_val_R1.fq trinity_gg/reads_R1.fq
    	mv $genome.RNA_val_R2.fq trinity_gg/reads_R2.fq
    	cp $genome.simple.fa trinity_gg/genome.fa
		cd trinity_gg   
    	bowtie2-build genome.fa genome
    	tophat -i 10 -I $MAX_INTRON_LENGTH -p $cpus genome reads_R1.fq reads_R2.fq
    	samtools sort tophat_out/accepted_hits.bam coordSorted
    	#Now launch Trinity Genome Guided via ssh
    	ssh -p 2222 ngs@localhost "$TRINITY/Trinity --genome $directory/genome.fa --genome_guided_use_bam $directory/trinity_gg/coordSorted.bam --genome_guided_max_intron $MAX_INTRON_LENGTH --genome_guided_sort_buffer 4G --genome_guided_CPU 2 --seqType fq --JM 2G --left $directory/trinity_gg/reads_1.fq --right $directory/trinity_gg/reads_R2.fq --CPU $trinity_cpus --jaccard_clip --output $directory/trinity_gg" 
    	#Shutdown Virtual box since we will need the RAM for PASA
		vboxmanage controlvm UbuntuServer poweroff
    	#now we want to concatenate results and use as input to the PASA pipeline to find high quality mRNAs
    	cd /Volumes/MacHD3/shared
    	cat $d/trinity/$species.Trinity.fasta $d/trinity_gg/Trinity-GG.fasta > $d/$species.all.trinity.fasta
    	mkdir $d/pasa
    	cd $d/pasa
		#Generate a list of transcript accession numbers from Trinitiy de novo assembly for comparison purposes - might use this later?
		$PASAHOME/misc_utilities/accession_extractor.pl < $species.Trinity.fasta > tdn.accs
		#Now Run PASA RNA-seq using the combined transcriptome. First we will make a directory and populate it with the necessary files.
		cp $PASAHOME/pasa_conf/pasa.alignAssembly.Template.txt alignAssembly.config
		gsed -i "s,<__MYSQLDB__>,$species," alignAssembly.config
		gsed -i "s,<__MIN_PERCENT_ALIGNED__>,75," alignAssembly.config
		gsed -i "s,<__MIN_AVG_PER_ID__>,95," alignAssembly.config
		$PASAHOME/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g ../$genome.simple.fa -t ../$species.all.trinity.fasta --ALIGNERS blat,gmap --CPU $cpus
		#mail -s "PASA Pipline has now finished." $email
   fi
fi
if [ $illumina_or_ion == ion ]; then
    echo -n "Enter the file name of your Ion Torrent RNA-seq data: "
    read RNA_reads
    mv $d/$RNA_reads $d/$genome.RNA.fastq
    echo -n "Enter number of bp to trim from 5' end: "
    read five_trim
    echo -n "Enter the maximum intron length: "
    read MAX_INTRON_LENGTH
    echo -n "How many CPUs do you want to use for Trinity? (max=10, recommended=6)  "
    read trinity_cpus
    #trim reads
    trim_galore -o $d -q 20 --phred33 --adapter ATCACCGACTGCCCATAGAGAGGCTGAGAC --clip_R1 $five_trim $d/$genome.RNA.fastq
    #estimate the memory required for Jellyfish/Trinity, i.e. 1 GB per 1 million reads
    num_reads=$(grep -c "^@" $d/$genome.RNA_trimmed.fq)
    memory_req="$(($num_reads / 1000000 + 1))"
    if [ "$memory_req" -gt "$max_mem" ]; then
        memory_req=50
    else
        memory_req=$memory_req
    fi
    #run trinity de novo
    TRINITY=/home/ngs/programs/trinityrnaseq_r20140717
    directory=/mnt/shared/$d
    echo "############################################################################################"
    echo -e "Ok, setting Jellyfish memory to $memory_req GB and running Trinity de novo. "
	#Now send ssh Trinity command to VirtualBox
	ssh -p 2222 ngs@localhost "$TRINITY/Trinity --seqType fq --JM $memory_req\G --SS_lib_type F --output $directory/$species --single $directory/$genome.RNA_trimmed.fq --CPU $trinity_cpus"
	#for some reason the symlink or move didn't work, so mv the file just in case from mac os
	mv $d/$species/Trinity.fasta.tmp $d/$species.Trinity.fasta
	ssh -p 2222 ngs@localhost "$TRINITY/util/TrinityStats.pl $directory/$species.Trinity.fasta > $directory/$species.Trinity.stats"
    echo "############################################################################################"
    echo "Trinity de novo has finished successfully, running stats and starting Trinity Genome Guided."
    #First run tophat alignment of RNA to genome
    cd $d
    mkdir -p trinity_gg
    #symlinks won't work because path is different from host to guest OS, so copy files (clean it up later)
    cp $genome.RNA_trimmed.fq trinity_gg/reads.fq
    cp $genome.simple.fa trinity_gg/genome.fa
	cd trinity_gg   
    bowtie2-build genome.fa genome
    tophat -i 10 -I $MAX_INTRON_LENGTH -p $cpus genome reads.fq
    samtools sort tophat_out/accepted_hits.bam coordSorted
    #Now launch Trinity Genome Guided via ssh
    ssh -p 2222 ngs@localhost "$TRINITY/Trinity --genome $directory/genome.fa --genome_guided_use_bam $directory/trinity_gg/coordSorted.bam --genome_guided_max_intron $MAX_INTRON_LENGTH --genome_guided_sort_buffer 4G --genome_guided_CPU 2 --seqType fq --JM 2G --single $directory/trinity_gg/reads.fq --SS_lib_type F --CPU $trinity_cpus --output $directory/trinity_gg" 
fi
echo "############################################################################################"
echo "Trinity Genome-Guided is finished.  Now combining results and running PASA to find likely ORF"
#concatenate the two Trinity results
cd /Volumes/MacHD3/shared
cat $d/$species.Trinity.fasta $d/trinity_gg/Trinity-GG.fasta > $d/$species.all.trinity.fasta
#Shutdown VirtualBox
vboxmanage controlvm UbuntuServer poweroff
#now we want to concatenate results and use as input to the PASA pipeline to find high quality mRNAs
mkdir $d/pasa
cd $d/pasa
#Generate a list of transcript accession numbers from Trinitiy de novo assembly for comparison purposes - might use this later?
$PASAHOME/misc_utilities/accession_extractor.pl < ../$species.Trinity.fasta > tdn.accs
#Now Run PASA RNA-seq using the combined transcriptome. First we will make a directory and populate it with the necessary files.
cp $PASAHOME/pasa_conf/pasa.alignAssembly.Template.txt alignAssembly.config
gsed -i "s,<__MYSQLDB__>,$species," alignAssembly.config
gsed -i "s,<__MIN_PERCENT_ALIGNED__>,75," alignAssembly.config
gsed -i "s,<__MIN_AVG_PER_ID__>,95," alignAssembly.config
$PASAHOME/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -r -R -g ../$genome.simple.fa -t ../$species.all.trinity.fasta --ALIGNERS blat,gmap --CPU $cpus --transcribed_is_aligned_orient --stringent_alignment_overlap 30.0 -I $MAX_INTRON_LENGTH
#now run Transdecoder script, not sure where it is looking for Pfam data in?  I guess it is just looking for ORFs which is apparently fine.
$PASAHOME/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta $species.assemblies.fasta --pasa_transcripts_gff3 $species.pasa_assemblies.gff3
trinity_genes=$(grep -c "^>" ../$species.all.trinity.fasta)
pasa_genes=$(grep -c "^>" $species.assemblies.fasta)
trans_genes=$(grep -c $'\tgene\t' $species.assemblies.fasta.transdecoder.gff3)
echo "Trinity, Trinity genome-guided, and PASA have finished.  Trinity resulted in $trinity_genes transcript assemblies and those were passed to PASA, which assembled them into $pasa_genes gene models.  We then ran TransDecoder to determine likely coding regions to be output, $species genome had $trinity_genes that passed and of those there are $trans_genes best models to use for training." | mail -s "Trinity-PASA done for $species" $email
simplifyFastaHeaders.pl $species.assemblies.fasta.transdecoder.cds pasa_ $species.pasa.transdecoder.fa header.pasa.map
#might need to do some cleaning up here, wait until finished to see the results.

}


RepeatMod ()
{
mkdir $d/RpMod
cd $d/RpMod
BuildDatabase -name $species $d/$genome.simple.fa
RepeatModeler -database $species
cd $d
RPMOD="$(find RpMod -type d -name 'RM_*')"
RepeatMasker -lib $d/$RPMOD/consensi.fa.classified $d/$genome.simple.fa
mail -s "RepeatMasking is finished for $species" $email < $d/$genome.simple.fa.tbl
gsed -i 's,#,_,g' $d/$RPMOD/consensi.fa.classified
perl /usr/local/maker/exe/RepeatMasker/util/rmOutToGFF3.pl $genome.simple.fa.out > $genome.repeatmasker.gff3
}

Cegma ()
{
mkdir $d/cegma
cd $d/cegma
cegma -T $cpus -g $d/$genome.simple.fa
mail -s "CEGMA has finished for $species" $email < $d/cegma/output.completeness_report
cegma2zff output.cegma.gff $d/$genome.simple.fa
fathom $d/cegma/genome.ann $d/cegma/genome.dna -categorize 1000
fathom -export 1000 -plus $d/cegma/uni.ann $d/cegma/uni.dna
forge $d/cegma/export.ann $d/cegma/export.dna
hmm-assembler.pl $genome . > $d/$genome.cegmasnap.hmm
cd $d
}

GeneMark ()
{
mkdir $d/genemark
cd $d/genemark
gmes_petap.pl --ES --fungus --cores $cpus --sequence $d/$genome.simple.fa
cp $d/genemark/output/gmhmm.mod $d/$genome\_gm.mod
cd $d
}

GeneMarkTopHat ()
{
mkdir $d/genemark-ET
#now convert to proper format
bet_to_gff.pl  --bet $d/Trinity_GG/tophat_out/junctions.bed   --gff $d/genemark-ET/introns.gff  --label TopHat2

#now run GeneMark Training
cd $d/genemark-ET
gmes_petap.pl --sequence $d/$genome.simple.fa --ET $d/genemark-ET/introns.gff  --et_score 4 --fungus --ES --cores $cpus

#get the training file
cp $d/genemark-ET/output/gmhmm.mod $d/$genome\_gm.mod
cd $d
}


Rfam ()
{
rfam_scan.pl -blastdb /home/ngs/BlastDB/Rfam.fasta -o $d/$genome.Rfam.gff3 /home/ngs/BlastDB/Rfam.cm $d/$genome.simple.fa
perl -i -pne 's/id=(.+?);/id=$1;Name=$1;/;s/rfam-id=(.+)/rfam-id=$1;Note=$1/;' $d/$genome.Rfam.gff3
}

Maker ()
{
cd $d
maker -CTL
#ok, lets look for the trinity transdecoder output file, if find it, then use it as EST evidence, otherwise ask
if [ -f "$d/$species.trinity.transdecoder.fa" ]; then
    echo "############################################################################################"
    echo "I've found Trinity derived transcripts, using those as EST evidence for Maker."
    sed -i "s,est= #set of ESTs,est=$d/$species.trinity.transdecoder.fa #set of ESTs," $d/maker_opts.ctl
else
    echo -n "Do you have EST data for this genome? [y/n] "
    read EST_data_question
    if [ $EST_data_question == y ]; then
        echo -n "Enter full path of EST data to use with Maker:  "
        read EST_data
        sed -i "s,est= #set of ESTs,est=$EST_data #set of ESTs," $d/maker_opts.ctl
    else
        echo -n "Ok, then enter EST data from a closely related species (full path):  "
        read alt_EST_data
        sed -i "s,altest= #EST/cDNA,altest=$alt_EST_data #EST/cDNA," $d/maker_opts.ctl
    fi
fi
#Now edit the maker opts file
sed -i "s,genome= #genome sequence,genome=$genome.simple.fa #genome sequence," $d/maker_opts.ctl
sed -i "s,protein=  #protein sequence,protein=$BLASTDB/swissprot_fungi.fa #protein sequence," $d/maker_opts.ctl
sed -i "s,snaphmm=,snaphmm=$genome.cegmasnap.hmm," $d/maker_opts.ctl
sed -i "s,gmhmm= #GeneMark,gmhmm=$genome\_gm.mod #GeneMark," $d/maker_opts.ctl
sed -i "s,rmlib= #provide,rmlib=$d/$RPMOD/consensi.fa.classified #provide," $d/maker_opts.ctl
sed -i 's,est2genome=0 #infer,est2genome=1 #infer,' $d/maker_opts.ctl
sed -i 's,protein2genome=0,protein2genome=1,' $d/maker_opts.ctl
sed -i "s,cpus=1,cpus=$cpus," $d/maker_opts.ctl
sed -i 's,split_hit=10000,split_hit=$MAX_INTRON_LENGTH,' $d/maker_opts.ctl

#Ok, now that options file has been edited, run first iteration of maker
maker -base maker1
echo "############################################################################################"
echo "The first iteration of maker  has finished, now retreiving results"
cd $d/maker1.maker.ouput
gff3_merge -o $d/$genome.maker1.gff3 -d maker1_master_datastore_index.log
cd $d

#Now re-train SNAP with output
mkdir $d/maker1snaptrain
cd $d/maker1snaptrain
maker2zff $d/$genome.maker1.gff3
fathom genome.ann genome.dna -categorize 1000
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl $genome . > $d/$genome.maker1.snap.hmm

#Convert to GFF3 to train augustus
zff2gff3.pl genome.ann | perl -plne 's/\t(\S+)$/\t\.\t$1/' > $d/$genome.augustus.gff3
cd $d

#Send an email to summarize results of first pass of Maker
grep_gff3="$(grep -c $'\tgene\t' $genome.maker1.gff3)"
maker1_qual="$(grep 'Name' $genome.augustus.gff3 | cut -f9 | sort | uniq | grep -c 'Name')"
echo "The first round of Maker has finished using GeneMark, Cegma-trained SNAP, protein2genome, cdna2genome, and BLAST searches using Trinity-Transdecoder RNA-seq or EST evidence from a closely related species in addition to using proteins from the swissprot fungal database.  This has produced $grep_gff3 gene models produced by Maker.  There are $maker1_qual high quality gene models that will be used to train Augustus and re-train SNAP." | mail -s "MAKER1 has finished for $species" $email

#Now Train Augustus
autoAug.pl --genome=$genome.simple.fa --species=$species --trainingset=$genome.augustus.gff3 --cdna=$species.trinity.transdecoder.fa singleCPU --useexisting --noutr -v

echo "Augustus training done. Now running second iteration of MAKER to finish the annotation pipeline for $species." | mail -s "Augustus training for $species finished" $email


}

echo "############################################################################################"
echo -e "Welcome to the Fungal Genome Annotation Script.  There are a number of options you can
choose for different analysis.  Please make a selection below.

    1) Run test Trinity RNA-seq assembly
"
echo -n "Enter selection from above:  "
read menu_selection
if [ $menu_selection == 1 ]; then
    echo "############################################################################################"
    echo "Ok, now I will need some additional information from you in order to get started."
    echo "You are currently in the $PWD directory."
    echo -n "Enter the directory that you will be working in (enter the full path - /Users/ngs/projects/myproject ):  "
    read d
    echo "############################################################################################"
    echo -n "Enter the species short name that you are annotating, i.e. 'Pdest':  "
    read species
    echo "############################################################################################"
    echo -n "Enter the genome variable, i.e. the base name of your fasta file, i.e. 'Psp3VT5':  "
    read genome
    echo "############################################################################################"
    echo -n "Enter your email address to receive updates and results:  "
    read email
    echo "############################################################################################"
    echo -n "How many CPUs (cores) do you want to use for annotation?  "
    read cpus
    simplifyFasta
    TrinityPGMorIllumina
fi

exit $?


