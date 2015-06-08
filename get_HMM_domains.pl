#!/usr/bin/env perl

#Wrapper script to run Hmmer3 to find coordinates of domain, reformat, and write just the domain for each to new file

use strict;
use warnings;
use Getopt::Long;

#set command line arguments
my ($proteins, $hmm, $length) = @ARGV;
my $version="get_HMM_domains.pl\tv0.1.1";
my $osname = $^O;
my $helpAsked;
my $hmm_db="/Users/jon/projects/DB/HMM"; #need to set this to your local folder containing HMM models
my $cpus = 1;
my $evalue = 1e-10;
my $outfile = "default_HMM_output.fasta";
GetOptions(
    'p|proteins:s'=>\$proteins, #protein fasta file
    'm|hmm:s'=>\$hmm, #name of HMM model to search
    'e|evalue:s'=>\$evalue, #directory you are working in
    'l|length:s'=>\$length, #min length required to keep domain
    'c|cpus:i'=>\$cpus, #cpus to use for hmmer search, default is one
    'o|out:s'=>\$outfile, # final multi-fasta file
    'h|help:s'=>\$helpAsked,
    'v|version'=>sub{print $version."\n"; exit;},
);

if(defined($helpAsked)) {
    prtUsage();
    exit;
}
if(!defined($proteins)) {
    print "You did not specify a protein fasta file, type get_HMM_domains.pl -h for help\n";
    exit;
}
print "-------------------------------------------------------------\n";
#check to see if required external programs are in path
my $check_hmm = `hmmscan -h | grep "# HMMER" | sed "s/# //" `;
my $check_bed = `bedtools --version`;

if (defined($check_hmm) && $check_hmm ne "") {
  print "\nFound HMMer: $check_hmm";
}
else {
  print "HMMer installation not found in PATH, exiting\n";
  exit;
}
if (defined($check_bed) && $check_bed ne "") {
  print "Found Bedtools: $check_bed\n";
}
else {
  print "BedTools installation not found in PATH, exiting\n";
  exit;
}
#Sub routine for running hmmscan and creating a bed file
sub runHMMer{
print "Running hmmerscan:\n   hmmscan --domtblout hmmscan.temp.out --cpu $cpus -E $evalue $hmm_db/$hmm $proteins\n";
my $hmm_results = `hmmscan --domtblout hmmscan.temp.out --cpu $cpus -E $evalue $hmm_db/$hmm $proteins`;
}
sub formatHMM{
my $output = 'hmmscan_domains.bed';
open(OUTPUT, '>'.$output) or die;
open(my $fh, "<:encoding(UTF-8)", 'hmmscan.temp.out') or die;
open(OUTPUT, '>>'.$output) or die;
while (<$fh>) {
    next if /^#/;
    my @row = split(/\s+/,$_);
    my $len = $row[20]-$row[19];
    if ($len > $length) {
        print OUTPUT "$row[3]\t$row[19]\t$row[20]\n";
    }
}
close $fh;
close $output;
}

sub getSeqs{
print "\nGetting Fasta coordinates:\n   bedtools getfasta -fullHeader -fi $proteins -bed hmmscan_domains.bed -fo $outfile\n";
my $fasta_seqs = `bedtools getfasta -fullHeader -fi $proteins -bed hmmscan_domains.bed -fo $outfile`
}

#sub routine for help
sub prtHelp {
    print "\n";
    print "All Options are Required:\n";
    print "-------------------------------------------------------------\n";
    print "--proteins flag:  fasta protein file \n";
    print "--hmm flag:  HMM model to use (KS or A)\n";
    print "--out flag:  Output file \n";
    print "--length flag:  Minimum length of domain\n";
    print "--evalue flag:  Evalue cutoff to use \n";
    print "--cpu flag:  Number of threads to use\n";
    print "-------------------------------------------------------------\n";
    print "example:  get_HMM_domains.pl --proteins proteins.fasta --hmm KS.hmm --length 300 --out domains_out.fasta 1\n";
}
#print usage subroutine
sub prtUsage {
    print "OS version: $osname\n";
    print "Usage: perl $0 -h\n";
    prtHelp();
}

#Run sub routines
runHMMer();
formatHMM();
getSeqs();


print "\n";
print "Script has successfully run, your $hmm domain proteins are in $outfile\n";
print "-------------------------------------------------------------\n";
unlink 'hmmscan_domains.bed';
unlink 'hmmscan.temp.out';
