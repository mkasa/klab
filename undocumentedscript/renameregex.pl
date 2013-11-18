#!/usr/bin/env perl

use strict;

use Getopt::Long;
use Pod::Usage;

my $flag_man     = 0;
my $flag_help    = 0;
my $flag_do      = 0;
my $debug        = 0;

GetOptions( 'help|?'  => \$flag_help,
	    'man'     => \$flag_man,
	    'debug'   => \$debug,
	    'do'      => \$flag_do,
	    ) or pod2usage(2);
pod2usage(1) if $flag_help;
pod2usage(-verbose => 2) if $flag_man; 

my $expr = shift;
my $to   = shift;

$expr = "^$expr\$";

if($expr eq '' || $to eq '') {
    print STDERR " usage: renameregex.pl [options...] <from expr> <toexpr>\n\n";
    print STDERR "   Quote by \"'\" in order not to expand regular expression by your shell.\n";
    print STDERR "   For safety, you have to add -do option to actually execute renaming.\n";
    print STDERR "   renameregex.pl just prints renaming commands without -do option.\n\n";
    print STDERR " Example 1. (Rename all *.xml => *.html) \n";
    print STDERR "   % renameregex.pl -do '(.*)\\.xml' '\$1.html'\n\n";
    exit 1;
}

my @files = ((glob "*"), (glob ".*"));
my $count = 0;
for my $file (@files) {
    if($file =~ /$expr/) {
	my $newfname = $file;
	eval "\$newfname =~ s/\$expr/$to/;";
	if($@) {
	    die "Invalid regular expression '$expr' or '$to'";
	}
	print "mv $file\t$newfname\n";
	if($flag_do) {
	    if(rename($file, $newfname)) {
		$count++;
	    } else {
		print STDERR "Failed to mv $file\t$newfname\n";
	    }
	}
    }
}
unless($flag_do) {
    print STDERR "\nIf you are satisfied with the renaming commands shown above,\njust add -do option to actually execute renaming\n\n";
} else {
    print "$count file(s) renamed\n";
}

=pod

=head1 NAME

renameregex.pl - Rename by regular expression

=head1 SYNOPSIS

renameregex.pl [options..] <Source regexp> <Dest regexp>

Options:
   -do              actually execute renaming commands
   -help            brief help message
   -man             full documentation

=head1 OPTIONS

=over 8

=item B<-do>

Without -do option, renameregex.pl will just print out the series of 'mv' commands.
This is for design because we often give a wrong regular expression that may destroy many files.
-do option must be added only when you are satisfied with the series of commands shown.

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will rename files by the rule specified by regular expression.
There are several styles in regular expression, but you must give to this program
regular expression that comforms to Perl's style.

=head1 EXAMPLES

	Example 1. (Rename all *.xml => *.html)
	    % renameregex.pl -do '(.*)\\.xml' '$1.html'

	Example 2. (Rename all abc1.* => abc2.*)
	    % renameregex.pl -do 'abc1\.(.*)' 'abc2\.$1'

=cut
