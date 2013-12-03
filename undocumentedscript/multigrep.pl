#!/usr/bin/env perl

use strict;

use Getopt::Long;
use Pod::Usage;

my $flag_man     = 0;
my $flag_help    = 0;
my $flag_invert  = 0;
my $debug        = 0;
my @param_keywords;
my @param_keywordfilename;

GetOptions( 'help|?'  => \$flag_help,
	    'man'     => \$flag_man,
	    'v'       => \$flag_invert,
	    'debug'   => \$debug,
	    'key=s'   => \@param_keywords,
	    'file=s'  => \@param_keywordfilename,
	    ) or pod2usage(2);
pod2usage(1) if $flag_help;
pod2usage(-verbose => 2) if $flag_man; 

for my $keywordfilename (@param_keywordfilename) {
    my $numberOfKeywordsInTheFile = 0;
    open FH, "< $keywordfilename" or die "Cannot open '$keywordfilename'";
    while(<FH>) {
	chomp;
	push(@param_keywords, $_);
	$numberOfKeywordsInTheFile++;
    }
    close FH;
    if($numberOfKeywordsInTheFile == 0) {
	print STDERR "WARNING: '$keywordfilename' has no keywords\n";
    }
}

if(@param_keywords <= 0) {
    print STDERR "ERROR: No keywords are specified.\n";
    exit 1;
}

while(<>){
    my $hasFoundKeyword = 0;
    for my $keyword (@param_keywords) {
	if(/$keyword/) {
	    $hasFoundKeyword = 1; last;
	}
    }
    print if($hasFoundKeyword ^ $flag_invert);
}

=pod

=head1 NAME

multigrep.pl - grep with multikeywords

=head1 SYNOPSIS

multigrep.pl [options] [file ...]

Options:
   -help            brief help message
   -man             full documentation
   -k=keyword       keyword
   -f=keywordfile   keywordfile

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-v>

Outputs only lines that do not contain keywords. It is similar to grep's '-v' option.

=back

=head1 DESCRIPTION

B<multigrep.pl> will grep by multi keywords.
Keywords can be specified at the command line options or by the file.

=head1 EXAMPLES

This example outputs lines containing 'abc' or 'def' in hogehoge.txt
    multigrep.pl -k=abc -k=def hogehoge.txt

This exmaple outputs lines containing keywords specified in keywords.txt. The search target is hogehoge.txt
    multigrep.pl -f=keywords.txt hogehoge.txt
Each line of keywords describes one keyword.


=cut

