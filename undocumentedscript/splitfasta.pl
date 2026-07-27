#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $flag_needhelp  = 0;
my $flag_equalbase = 0;
my $flag_force     = 0;

GetOptions('help'      => \$flag_needhelp,
	   'equalbase' => \$flag_equalbase,
	   'force'     => \$flag_force
	   );

unless(@ARGV == 2 && !$flag_needhelp){
    print STDERR "Usage: splitfasta.pl [options] <fasta file> <split size>\n";
    print STDERR "    FASTA file is split into FASTA files containing at most\n";
    print STDERR "    <split size> sequences on default. When --equalbase option\n";
    print STDERR "    is specified, splitfasta.pl first counts the number of\n";
    print STDERR "    bases in the given FASTA file first. Then it splits \n";
    print STDERR "    given FASTA sequences into FASTA files containing\n";
    print STDERR "    approximately (total bases) / (split size) bases.\n";
    print STDERR "    --force overwrites existing output files; without it\n";
    print STDERR "    splitfasta.pl refuses to clobber '<fasta file>.<n>'.\n";
    exit 1;
}

my $fastafiletobesplit = shift;
my $splitsize = shift;

unless($splitsize =~ /^\d+$/ && $splitsize >= 1){
    die "split size must be a positive integer (got '$splitsize')\n";
}

print STDERR "Target FASTA file = '$fastafiletobesplit'\n";
print STDERR "split size = $splitsize ";
if($flag_equalbase) {
    print STDERR "files";
} else {
    print STDERR "sequences";
}
print STDERR "\n\n";

my $number_of_bases_to_split = 0;

if($flag_equalbase) {
    print STDERR "Equal #base mode\n";
    print STDERR "  Counting total #bases\n";
    open(my $cfh, '<', $fastafiletobesplit) or die "cannot open '$fastafiletobesplit': $!";
    my $total_number_of_bases = 0;
    my $total_number_of_reads = 0;
    while(<$cfh>){
	if(/^>/){ # definition line
	    $total_number_of_reads ++;
	} else {
	    # baselength() is used here too, so that this pre-count and the
	    # count made while splitting agree on CRLF input.
	    $total_number_of_bases += baselength($_);
	}
    }
    print STDERR "  #sequences = $total_number_of_reads\n";
    print STDERR "  #bases     = $total_number_of_bases\n";
    $number_of_bases_to_split = int($total_number_of_bases / $splitsize);
    print STDERR "  split around ${number_of_bases_to_split}bp\n";
    close $cfh;
} else {
    print STDERR "Equal #sequence mode\n";
}

open(my $fh, '<', $fastafiletobesplit) or die "cannot open '$fastafiletobesplit': $!";
print STDERR "Splitting...\n";

my $defline;
my $seqs    = '';
my $fnumber = 1;
my $numseq  = 0;
my $numbase = 0;
my $ofh;                      # undef until the first record is actually written
my $current_output_filename;

while(<$fh>){
    if(/^>/){ # definition line
	flushifany();
	$defline = $_;
    } else {
	$seqs .= $_;
    }
}
flushifany();
closeoutputf();
close $fh;

sub flushifany
{
    return unless(defined $defline);
    # The counters are updated AFTER the record has been written, so that the
    # numbers reported by closeoutputf() describe what really went into the
    # file that is being closed.
    my $len = baselength($seqs);
    if($flag_equalbase) {         # base split mode
	if($numbase > $number_of_bases_to_split && $fnumber < $splitsize) {
	    closeoutputf();
	}
    } else {                      # sequence split mode
	if($numseq >= $splitsize) {
	    closeoutputf();
	}
    }
    openoutputf() unless(defined $ofh);
    print $ofh $defline or die "cannot write to '$current_output_filename': $!";
    print $ofh $seqs    or die "cannot write to '$current_output_filename': $!";
    $numseq++;
    $numbase += $len;
    $defline = undef;
    $seqs = '';
}

sub baselength
{
    my $seq = shift;
    return 0 unless(defined $seq);
    $seq =~ s/[\r\n]//g;
    return length($seq);
}

sub openoutputf
{
    my $fname = "$fastafiletobesplit.$fnumber";
    if(-e $fname && !$flag_force) {
	die "'$fname' already exists. Remove the previous output (or give --force) before running again.\n";
    }
    print STDERR "Output to '$fname'\n";
    open($ofh, '>', $fname) or die "cannot open '$fname' for output: $!";
    $current_output_filename = $fname;
}

sub closeoutputf
{
    return unless(defined $ofh);
    print STDERR "  $numseq sequences, $numbase bases\n";
    close($ofh) or die "error while writing '$current_output_filename': $!";
    $ofh = undef;
    $numseq  = 0;
    $numbase = 0;
    $fnumber++;
}
