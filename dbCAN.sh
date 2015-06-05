#!/bin/bash

#Script to run dbCAN HMM search locally
#usage - sh dbCAN.sh proteins.fasta root_name cpus 

if [ -z "$1" ]; then
	echo "Run this script as follows:
	$SCRIPTS/dbCAN.sh yourproteins.fasta root_name num_cpus"
	exit
else
	#first run HMMer with HMM models
	hmm=/usr/local/dbCAN/dbCAN-fam-HMMs.txt

	hmmscan --cpu $3 $hmm $1 > $2.hmmout

	#hmmparser.sh cut and pasted below
	cat $2.hmmout| perl -e 'while(<>){if(/^\/\//){$x=join("",@a);($q)=($x=~/^Query:\s+(\S+)/m);while($x=~/^>> (\S+.*?\n\n)/msg){$a=$&;@c=split(/\n/,$a);$c[0]=~s/>> //;for($i=3;$i<=$#c;$i++){@d=split(/\s+/,$c[$i]);print $q."\t".$c[0]."\t$d[6]\t$d[7]\t$d[8]\t$d[10]\t$d[11]\n" if $d[6]<1;}}@a=();}else{push(@a,$_);}}' \
	| gsort -k 1,1 -k 6,6n -k 7,7n | uniq \
        | perl -e 'while(<>){chomp;@a=split;next if $a[-1]==$a[-2];push(@{$b{$a[0]}},$_);}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-1]-$b[-2];$len2=$c[-1]-$c[-2];$len3=$b[-1]-$c[-2];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[2]<$c[2]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_."\n";}}' \
        | uniq | perl -e 'open(IN,"/usr/local/dbCAN/all.hmm.ps.len.ps");while(<IN>){chomp;@a=split;$b{$a[0]}=$a[1];}while(<>){chomp;@a=split;$r=($a[4]-$a[3])/$b{$a[1]};print $_."\t".$r."\n";}' \
	| perl -e 'while(<>){@a=split(/\t/,$_);if(($a[-2]-$a[-3])>80){print $_ if $a[2]<1e-5;}else{print $_ if $a[2]<1e-3;}}' | gawk '$NF>0.3' | gsed 's/.hmm//g' > $2.dbCAN.results
	#clean up
	rm $2.hmmout
	exit
fi
