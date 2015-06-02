#!/bin/bash

#Script for getting dBCAN.sh script results into comparative matrix table


#grab list of models and keep running tally
for file in *.results; do cut -f2 $file >> all_hmm.list; done
#now clean up this list for uniqs
sort all_hmm.list | uniq > hmm_uniq.list
#clean up files with spaces
for file in *.results; do gsed -i 's/ //g' $file; done
#let's also represent these data as 1 domain per protein, so we get enzyme numbers
for file in *.results; do gsort -k1,1 -k2,2 -k3,3nr $file | gsort -u -k1,1 -k2,2 > $file.enzymes; done
#now loop through each results (which is domains in genome file
for file in *.results; do echo "Now Processing: $file"; while read line; do echo "$line\t\c"; cut -f2 $file | grep -Fxc "$line"; done < hmm_uniq.list > $file.temp; done
#now get the names into file for list
echo "CAZYme\t\c" > dbCAN_domains_table.txt
printf '%s\t' *.temp | gsed 's/.fasta.dbCAN.results.temp//g' >> dbCAN_domains_table.txt
echo "" >> dbCAN_domains_table.txt
#finally count up the numbers of each domain 
gawk 'NF > 1{ a[$1] = a[$1] "\t" $2} END {for( i in a ) print i, a[i]}' *.temp | gsort >> dbCAN_domains_table.txt

#clean up your mess
rm *.temp

#now loop through each results (which is domains in genome file
for file in *.enzymes; do echo "Now Processing: $file"; while read line; do echo "$line\t\c"; cut -f2 $file | grep -Fxc "$line"; done < hmm_uniq.list > $file.temp2; done
#now get the names into file for list
echo "CAZYme\t\c" > dbCAN_enzymes_table.txt
printf '%s\t' *.temp2 | gsed 's/.fasta.dbCAN.results.enzymes.temp2//g' >> dbCAN_enzymes_table.txt
echo "" >> dbCAN_enzymes_table.txt
#finally count up the numbers of each domain 
gawk 'NF > 1{ a[$1] = a[$1] "\t" $2} END {for( i in a ) print i, a[i]}' *.temp2 | gsort >> dbCAN_enzymes_table.txt

#clean up your mess
rm *.temp2
