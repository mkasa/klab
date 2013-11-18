#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $flag_needhelp  = 0;
my $flag_equalbase = 0;

GetOptions('help'      => \$flag_needhelp,
	   'equalbase' => \$flag_equalbase
	   );

unless(@ARGV == 2){
    print STDERR "Usage: splitfasta.pl [options] <fasta file> <split size>\n";
    print STDERR "    FASTA file is split into FASTA files containing at most\n";
    print STDERR "    <split size> sequences on default. When --equalbase option\n";
    print STDERR "    is specified, splitfasta.pl first counts the number of\n";
    print STDERR "    bases in the given FASTA file first. Then it splits \n";
    print STDERR "    given FASTA sequences into FASTA files containing\n";
    print STDERR "    approximately (total bases) / (split size) bases.\n";
    exit 1;
}

my $fastafiletobesplit = shift;
my $splitsize = shift;

if($splitsize <= 1){
    die 'splitsize must be greater than one';
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
    open FH, "< $fastafiletobesplit" or die "cannot open '$fastafiletobesplit'";
    my $total_number_of_bases = 0;
    my $total_number_of_reads = 0;
    while(<FH>){
	if(/^>/){ # definition line
	    $total_number_of_reads ++;
	} else {
	    chomp;
	    $total_number_of_bases += length($_);
	}
    }
    print STDERR "  #sequences = $total_number_of_reads\n";
    print STDERR "  #bases     = $total_number_of_bases\n";
    $number_of_bases_to_split = int($total_number_of_bases / $splitsize);
    print STDERR "  split around ${number_of_bases_to_split}bp\n";
    close FH;
} else {
    print STDERR "Equal #sequence mode\n";
}

open FH, "< $fastafiletobesplit" or die "cannot open '$fastafiletobesplit'";
print STDERR "Splitting...\n";

my $defline;
my $seqs;
my $fnumber = 1;
my $numseq  = 0;
my $numbase = 0;
openoutputf();

while(<FH>){
    if(/^>/){ # definition line
	flushifany();
	$defline = $_;
    } else {
	$seqs .= $_;
    }
}
flushifany();
closeoutputf();
close FH;

sub flushifany()
{
    if(defined $defline){
	if($flag_equalbase) {         # base split mode
	    if($numbase > $number_of_bases_to_split) {
		if($fnumber < $splitsize) {
		    closeoutputf();
		    openoutputf();
		}
	    }
	}
	$numseq++;
	$numbase += baselength($seqs);
	unless($flag_equalbase) {   # sequence split mode
	    if($numseq > $splitsize) {
		closeoutputf();
		$numseq++;
		openoutputf();
	    }
	}
	print OFH "$defline";
	print OFH "$seqs";
	$defline = undef;
	$seqs = '';
    }
}

sub baselength(){
    my $seq = shift;
    $seq =~ s/[\r\n]//g;
    return length($seq);
}

sub openoutputf(){
    print STDERR "Output to '$fastafiletobesplit.$fnumber'\n";
    open OFH, "> $fastafiletobesplit.$fnumber" or die "cannot open '$fastafiletobesplit.$fnumber' for output";
}

sub closeoutputf(){
    print STDERR "  $numseq sequences, $numbase bases\n";
    close OFH;
    $numseq  = 0;
    $numbase = 0;
    $fnumber++;
}
