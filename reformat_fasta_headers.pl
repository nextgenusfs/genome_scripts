#!/usr/bin/env perl -w
=usage 

reformat_fasta_headers.pl -f fasta_file -a annotation file (2 columns tab delimited: find col 1 and replace with col2)

=cut

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

#set command line arguments
my ($fasta, $annot) = @ARGV;
my $version="reformat_fasta_headers.pl\tv0.0.1";
GetOptions(
	'f|fasta:s'=>\$fasta,
	'a|annot:s'=>\$annot,
	'v|version'=>sub{print $version."\n"; exit;},
);


open my $fh, '<', $annot or die $!;
my %annot = map { /(\S+)\s+(.+)/; $1 => $2 } <$fh>;
close $fh;

my $in = Bio::SeqIO->new( -file => $fasta, -format => 'Fasta' );

while ( my $seq = $in->next_seq() ) {
    my $seqID = $annot{ $seq->id } // $seq->id;
    print ">$seqID\n" . $seq->seq . "\n";
}
