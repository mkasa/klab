#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

my $flag_man     = 0;
my $flag_help    = 0;
my $flag_invert  = 0;
my $flag_regex   = 0;
my $debug        = 0;
my @param_keywords;
my @param_keywordfilename;

GetOptions( 'help|?'  => \$flag_help,
	    'man'     => \$flag_man,
	    'v'       => \$flag_invert,
	    'regex'   => \$flag_regex,
	    'debug'   => \$debug,
	    'key=s'   => \@param_keywords,
	    'file=s'  => \@param_keywordfilename,
	    ) or pod2usage(2);
pod2usage(1) if $flag_help;
pod2usage(-verbose => 2) if $flag_man;

for my $keywordfilename (@param_keywordfilename) {
    my $numberOfKeywordsInTheFile = 0;
    open(my $kfh, '<', $keywordfilename) or die "Cannot open '$keywordfilename': $!";
    while(<$kfh>) {
	chomp;
	s/\r\z//;   # chomp alone leaves the CR of a CRLF keyword file
	next if($_ eq '');
	push(@param_keywords, $_);
	$numberOfKeywordsInTheFile++;
    }
    close $kfh or print STDERR "WARNING: error while closing '$keywordfilename': $!\n";
    if($numberOfKeywordsInTheFile == 0) {
	print STDERR "WARNING: '$keywordfilename' has no keywords\n";
    }
}

# Compile all the patterns up front: an empty pattern would mean "reuse the
# last successful match" (i.e. match everything), and an invalid one used to
# blow up in the middle of the input after part of the output was written.
my @patterns;
my %seen_keyword;
for my $keyword (@param_keywords) {
    next unless(defined $keyword);
    $keyword =~ s/[\r\n]+\z//;
    if($keyword eq '') {
	print STDERR "WARNING: an empty keyword is ignored\n";
	next;
    }
    next if($seen_keyword{$keyword}++);
    my $re;
    if($flag_regex) {
	$re = eval { qr/$keyword/ };
	unless(defined $re) {
	    my $errmsg = defined $@ ? $@ : '';
	    $errmsg =~ s/\s+\z//;
	    print STDERR "ERROR: invalid regular expression '$keyword'\n";
	    print STDERR "  $errmsg\n" if($errmsg ne '');
	    exit 1;
	}
    } else {
	$re = qr/\Q$keyword\E/;   # keywords are literal unless --regex is given
    }
    print STDERR "DEBUG: keyword '$keyword' => $re\n" if($debug);
    push(@patterns, $re);
}

if(@patterns <= 0) {
    print STDERR "ERROR: No keywords are specified.\n";
    exit 1;
}

while(<>){
    my $hasFoundKeyword = 0;
    for my $pattern (@patterns) {
	if(/$pattern/) {
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

=item B<-regex>

Treat each keyword as a Perl regular expression. By default keywords are
matched literally (as with grep -F), so characters such as '*' or '(' have no
special meaning. All patterns are validated before any input is read.

=back

=head1 DESCRIPTION

B<multigrep.pl> will grep by multi keywords.
Keywords can be specified at the command line options or by the file.
Empty keywords (e.g. a blank or trailing line in a keyword file) are ignored,
and duplicated keywords are used only once.

=head1 EXAMPLES

This example outputs lines containing 'abc' or 'def' in hogehoge.txt
    multigrep.pl -k=abc -k=def hogehoge.txt

This exmaple outputs lines containing keywords specified in keywords.txt. The search target is hogehoge.txt
    multigrep.pl -f=keywords.txt hogehoge.txt
Each line of keywords describes one keyword.


=cut

