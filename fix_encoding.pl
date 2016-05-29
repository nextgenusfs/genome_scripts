#!/usr/bin/env perl
use strict;
use warnings;
use Encode qw(encode decode);
use Encode::Detect::Detector;

my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $counter = 0;

open (OUTFILE, ">", $outfile);
open (INFILE, $infile);
while (<INFILE>) {
    $counter++;
    my $encoding = Encode::Detect::Detector::detect($_);
    if (defined($encoding) and ($encoding eq "windows-1252")) {
        my $tmp_line = decode($encoding, $_); $_ = encode("utf-8", $tmp_line);
        # $_ =~ s/ë/e/g;
        print "Encoding found on line " . $counter . ": " . $encoding . ", will try to encode to utf-8.\n"; 
    }
    print OUTFILE $_;
}

close INFILE;
close OUTFILE;

exit;

