#!/usr/bin/env perl

=head1 DESCRIPTION

	This script takes the Secondary Metabolite Cluster prediction text file from SMURF
	(Secondary-Metabolite-Clusters.txt) to print amino acid and DNA sequence.

=head1 USAGE

	smurf_reformat.pl -i Secondary-Metabolite-Clusters.txt -g genome.fasta -p proteins.fasta

=head2 Options
 -c|--cluster
 -b|--backbone
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
my ($infi, $backbone, $scaffolds, $proteins) = @ARGV;
my $version="smurf_reformat.pl\tv0.1";
my $out="smurf";
my $helpAsked;
GetOptions(
	'c|cluster:s'=>\$infi,
	'b|backbone:s'=>\$backbone,
	'g|genome:s'=>\$scaffolds,
	'p|proteins:s'=>\$proteins, #protein fasta file
	'o|out:s'=>\$out, #out name stem
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
	print "-c (--cluster):  SMURF Cluster file \n";
	print "-b (--backbone):  SMURF Backbone file\n";
	print "-g (--genome):  gDNA sequences fasta format\n";
	print "-p (--proteins):  protein sequences fasta format \n";
	print "-o (--out): output name stem (out_clusters.txt) \n";
	print "-h (--help): Print this message\n";
	print "-------------------------------------------------------------\n";
	print "example:  smurf_reformat.pl -i Secondary-Metabolite-Clusters.txt -g genome.fasta -p proteins.fasta -o outfile  \n\n"

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

my $output .= "$out\_clusters.txt";
open(OUTPUT, '>'.$output) or die;
open my $fh, "<", $infi;
print OUTPUT "#Cluster\tGene_ID\tChromosome-Scaffold\t5' end\t3' end\tDescription\tProtSeq\tDNASeq\n";
my %names;
my $count;
while (<$fh>) {
    next if  /^Cluster/;
    next if /^Backbone/;
    next if /^\n/;
    next if /^\r$/;
    chomp;
    my @row = split(/\t/,$_);
    if (not exists $names{$row[0]}) {
        $names{$row[0]} = "Cluster_" . ++$count;
    }
    my $db = Bio::DB::Fasta->new('smurf.scaffolds.tmp');
    my $subseq = $db->seq($row[3], $row[5], $row[6]);
    chomp($subseq);
    my $db2 = Bio::DB::Fasta->new('smurf.proteins.tmp');
    my $protseq = $db2->seq($row[1]);
    chomp($protseq);
    print OUTPUT "$names{$row[0]}\t$row[1]\t$row[3]\t$row[5]\t$row[6]\t$row[9]\t$protseq\t$subseq\n";


    }

close $fh;

my $output2 .= "$out\_backbone.txt";
open(OUTPUT2, '>'.$output2) or die;
open my $fh2, "<", $backbone;
print OUTPUT2 "#Gene_ID\tPrediction\tChromosome-Scaffold\t5' end\t3' end\tDescription\tProtSeq\tDNASeq\n";
while (<$fh2>) {
    next if /^Backbone/;
    next if /^\n/;
    next if /^\r$/;
    chomp;
    my @field = split(/\t/,$_);
    my $db3 = Bio::DB::Fasta->new('smurf.scaffolds.tmp');
    my $subseq2 = $db3->seq($field[2], $field[4], $field[5]);
    chomp($subseq2);
    my $db4 = Bio::DB::Fasta->new('smurf.proteins.tmp');
    my $protseq2 = $db4->seq($field[0]);
    chomp($protseq2);
    print OUTPUT2 "$field[0]\t$field[6]\t$field[2]\t$field[4]\t$field[5]\t$field[1]\t$protseq2\t$subseq2\n";


    }
close $fh2;


#clean up temp files
unlink 'smurf.scaffolds.tmp';
unlink 'smurf.scaffolds.tmp.index';
unlink 'smurf.proteins.tmp';
unlink 'smurf.proteins.tmp.index';
exit;
