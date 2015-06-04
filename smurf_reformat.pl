#!/usr/bin/env perl

=head1 DESCRIPTION

	This script takes the Secondary Metabolite Cluster prediction text file from SMURF
	(Secondary-Metabolite-Clusters.txt) to print amino acid and DNA sequence.

=head1 USAGE

	smurf_reformat.pl -i Secondary-Metabolite-Clusters.txt -g genome.fasta -p proteins.fasta

=head2 Options
 -i|--in
 -g|--genome
 -p|--proteins
 -v|--version
	

=head1 Author

	Jon Palmer
	palmerjona at gmail dot com

=cut

use strict;
use warnings;
use Getopt::Long;
use autodie;
use Bio::SeqIO;
use Bio::DB::Fasta;

#set command line arguments
my ($infi, $scaffolds, $proteins) = @ARGV;
my $version="smurf_reformat.pl\tv0.1";
my $helpAsked;
GetOptions(
	'i|in:s'=>\$infi,
	'g|genome:s'=>\$scaffolds,
	'p|proteins:s'=>\$proteins, #protein fasta file
	'h|help:s'=>\$helpAsked,
	'v|version'=>sub{print $version."\n"; exit;},
);

if(defined($helpAsked)) {
	prtUsage();
	exit;
	}
sub prtUsage {
print "\n";
	print "Options:\n";
	print "-------------------------------------------------------------\n";
	print "-i (--in):  SMURF Cluster file \n";
	print "-g (--genome):  gDNA sequences fasta format\n";
	print "-p (--proteins):  protein sequences fasta format \n";
	print "-h (--help): Print this message\n";
	print "-------------------------------------------------------------\n";
	print "example:  smurf_reformat.pl -i Secondary-Metabolite-Clusters.txt -g genome.fasta -p proteins.fasta > output.txt\n\n"
	
}

if(!defined($infi)) {
	print "You did not specify a SMURF fasta file: type $0 --help for usage instructions\n";
	exit;
}

#Reformat genome input to make sure it is wrapped as if single line can be problematic
my $seqIN = Bio::SeqIO->new('-format' => 'fasta', '-file' => "$scaffolds");
my $seqOUT = Bio::SeqIO->new('-format' => 'fasta', '-file' => ">smurf.scaffolds.tmp");
while (my $inseq = $seqIN->next_seq) {
    $seqOUT->write_seq($inseq);
}

#reformat protein input to make sure it is properly wrapped
my $seqIN2 = Bio::SeqIO->new('-format' => 'fasta', '-file' => "$proteins");
my $seqOUT2 = Bio::SeqIO->new('-format' => 'fasta', '-file' => ">smurf.proteins.tmp");
while (my $inseq2 = $seqIN2->next_seq) {
    $seqOUT2->write_seq($inseq2);
}

open my $fh, "<", $infi;
print "#Cluster\tGene_ID\tChromosome-Scaffold\t5' end\t3' end\tDescription\tProtSeq\tDNASeq\n";
my %names;
my $count;
while (<$fh>) {
    next if  /^Cluster/;
    next if /^Backbone/;
    next if /^\n/;
    my @row = split(/\t/,$_);
    if (not exists $names{$row[0]}) {
        $names{$row[0]} = "Cluster_" . ++$count;
    }
    my $db = Bio::DB::Fasta->new('smurf.scaffolds.tmp');
    my $subseq = $db->seq($row[3], $row[5], $row[6]);
    my $db2 = Bio::DB::Fasta->new('smurf.proteins.tmp');
    my $protseq = $db2->seq($row[1]);
    print "$names{$row[0]}\t$row[1]\t$row[3]\t$row[5]\t$row[6]\t$row[9]\t$protseq\t$subseq\n";
    
    }

close $fh;
#clean up temp files
unlink 'smurf.scaffolds.tmp';
unlink 'smurf.scaffolds.tmp.index';
unlink 'smurf.proteins.tmp';
unlink 'smurf.proteins.tmp.index';
exit;
