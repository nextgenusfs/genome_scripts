#!/bin/bash
#script to run minimal GeneMark and Augustus gene prediction for A.fumigatus

#run script error message
if [ -z "$1" ]; then
	echo "Run this script as follows:
	$0 <Fasta file>"
	exit
else
    #get current directory and save for future
    dir=$PWD
    #get script directory
    script_dir=$(command -v Afumi_pred.sh | gsed 's,/Afumi_pred.sh,,g')
    #take input and make new name for folder
    genome=${1%%.*}
    #get last 3 characters of string to use for Gene naming so we have unique names later on, write a little loop for CEA10
    if [ $genome == "CEA10" ]; then
        gene_name=$genome
    else
        gene_name=${genome: -3}
        gene_name="AF$gene_name"
    fi

    mkdir -p $genome
    cd $genome

    #set up maker control files
    maker -CTL

    #edit opts file
    gsed -i "s,genome= #genome sequence,genome=$dir/$1 #genome sequence," maker_opts.ctl
    gsed -i "s,est= #set of ESTs,est=$dir/Af293_cDNA.fasta #set of ESTs ," maker_opts.ctl
    gsed -i "s,gmhmm= #GeneMark,gmhmm=$dir/Afumi.mod #GeneMark," maker_opts.ctl
    gsed -i "s,augustus_species= #Augustus,augustus_species=aspergillus_fumigatus #Augustus," maker_opts.ctl
    gsed -i "s,protein=  #protein sequence,protein=$BLASTDB/uniprot_sprot.trinotate.pep #protein sequence," maker_opts.ctl
    gsed -i 's,model_org=all #select,model_org=fungi #select,' maker_opts.ctl
    gsed -i "s,cpus=1,cpus=20," maker_opts.ctl
    gsed -i "s,min_protein=0 #require,min_protein=50 #require," maker_opts.ctl
    gsed -i 's,split_hit=10000,split_hit=1500,' maker_opts.ctl
    #gsed -i 's,clean_up=0 #removes,clean_up=1 #removes,' maker_opts.ctl
    gsed -i 's,keep_preds=0 #Concordance,keep_preds=1 #Concordance,' maker_opts.ctl

    #run maker
    maker -base $genome

    #get maker generated GFF3
    cd $dir/$genome/$genome.maker.output
    gff3_merge -n -d $genome\_master_datastore_index.log
    maker_map_ids --prefix $gene_name\_ --justify 5 $genome.all.gff > $genome.all.id.map
    map_gff_ids $genome.all.id.map $genome.all.gff

    #make GBK file for genome via GAG
    mkdir -p $dir/$genome/GAG
    cp $dir/$1 $dir/$genome/GAG/genome.fasta
    cp $dir/$genome/$genome.maker.output/$genome.all.gff $dir/$genome/GAG/genome.gff
    cd $dir/$genome/GAG
    gag.py fasta=genome.fasta gff=genome.gff
    cd $dir/$genome/GAG/gag_output
    mv genome.fasta genome.fsa
    tbl2asn -p . -t $script_dir/test.sbt -M n -Z discrep -a r10u -l paired-ends -j "[organism=$genome]" -V b -c fx
    cp genome.gbf $dir/$genome/$genome.gbk
    cd $dir/$genome
    #run some stats with your script and see how we did
    gb_stats.py $genome.gbk > $genome.stats.txt
    #output SMURF data
    gb2smurf.py $genome.gbk -p $genome.proteins.fasta -g $genome.scaffolds.fasta -s $genome.smurf.txt
    #add ACCESSION numbers to work with MGB
    perl -i -pe "s/^ACCESSION   /$& .$gene_name\_Contig_.++$n/ge" $genome.gbk
    #done
    cd $dir
    echo "Script is finished"

fi

