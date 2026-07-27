#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

my $flag_man     = 0;
my $flag_help    = 0;
my $flag_do      = 0;
my $flag_force   = 0;
my $debug        = 0;

GetOptions( 'help|?'  => \$flag_help,
	    'man'     => \$flag_man,
	    'debug'   => \$debug,
	    'do'      => \$flag_do,
	    'force'   => \$flag_force,
	    ) or pod2usage(2);
pod2usage(1) if $flag_help;
pod2usage(-verbose => 2) if $flag_man;

my $expr = shift;
my $to   = shift;

# NOTE: the guard must come BEFORE the anchoring, otherwise "^$expr\$" is
#       never the empty string and the 'no source pattern given' check is dead.
if(!defined $expr || $expr eq '' || !defined $to || $to eq '') {
    print STDERR " usage: renameregex.pl [options...] <from expr> <toexpr>\n\n";
    print STDERR "   Quote by \"'\" in order not to expand regular expression by your shell.\n";
    print STDERR "   For safety, you have to add -do option to actually execute renaming.\n";
    print STDERR "   renameregex.pl just prints renaming commands without -do option.\n\n";
    print STDERR " Example 1. (Rename all *.xml => *.html) \n";
    print STDERR "   % renameregex.pl -do '(.*)\\.xml' '\$1.html'\n\n";
    exit 1;
}

my $anchored_expr = "^$expr\$";
# Compile the source pattern up front so that a bad pattern is reported
# with a friendly message instead of blowing up at the first match.
my $compiled_expr = eval { qr/$anchored_expr/ };
unless(defined $compiled_expr) {
    my $errmsg = $@;
    $errmsg = '' unless(defined $errmsg);
    $errmsg =~ s/\s+\z//;
    print STDERR "Invalid source regular expression '$expr'\n";
    print STDERR "  $errmsg\n" if($errmsg ne '');
    exit 1;
}

# Expand $1..$9 / ${1}..${nn} in the replacement string.
#
# SECURITY: the replacement used to be interpolated into a string that was
#           then eval'ed as Perl source, so a replacement such as
#           '$1.html/; system("rm -rf /"); #' executed arbitrary code - and it
#           did so even without -do, i.e. in the mode documented as safe.
#           Nothing here is eval'ed and no shell is involved.
sub expand_replacement($$)
{
    my ($template, $groupsref) = @_;
    my $out = $template;
    $out =~ s{ \\(.) | \$\{([1-9][0-9]*)\} | \$([1-9]) }
             { defined $1 ? $1 : group_value($groupsref, defined $2 ? $2 : $3) }gex;
    return $out;
}

sub group_value($$)
{
    my ($groupsref, $n) = @_;
    return '' if($n < 1 || $n > @$groupsref);
    my $v = $groupsref->[$n - 1];
    return defined $v ? $v : '';
}

# Quote a string so that the printed 'mv' command can actually be run.
sub shell_quote($)
{
    my $s = shift;
    return "''" unless(defined $s && $s ne '');
    return $s if($s =~ m|^[A-Za-z0-9_./-]+$|);
    $s =~ s/'/'\\''/g;
    return "'" . $s . "'";
}

my @files = ((glob "*"), (glob ".*"));
my @plan;
for my $file (@files) {
    next if($file eq '.' || $file eq '..');   # glob ".*" returns these too
    next unless($file =~ $compiled_expr);
    my @groups;
    for(my $n = 1; $n < scalar(@-); $n++) {
	push(@groups, defined $-[$n] ? substr($file, $-[$n], $+[$n] - $-[$n]) : undef);
    }
    my $newfname = expand_replacement($to, \@groups);
    print STDERR "DEBUG: '$file' -> '$newfname'\n" if($debug);
    push(@plan, [$file, $newfname]);
}

# Refuse to destroy data: check the whole batch before touching anything.
my @problems;
my %destination_count;
for my $p (@plan) {
    $destination_count{$p->[1]}++;
}
for my $p (@plan) {
    my ($file, $newfname) = @$p;
    if($newfname eq '') {
	push(@problems, "'$file' would be renamed to an empty name");
	next;
    }
    if($destination_count{$newfname} > 1) {
	push(@problems, "several files would be renamed to '$newfname'");
	next;
    }
    if($newfname ne $file && -e $newfname) {
	push(@problems, "'$newfname' already exists (would be overwritten by '$file')");
    }
}
{
    my %seen;
    @problems = grep { !$seen{$_}++ } @problems;
}

for my $p (@plan) {
    my ($file, $newfname) = @$p;
    print "mv ", shell_quote($file), " ", shell_quote($newfname), "\n";
}

if(@problems && !$flag_force) {
    print STDERR "\nERROR: refusing to rename anything because:\n";
    print STDERR "  - $_\n" for(@problems);
    print STDERR "Nothing has been renamed. Add --force if you really mean to overwrite.\n\n";
    exit 1;
}

my $count = 0;
if($flag_do) {
    for my $p (@plan) {
	my ($file, $newfname) = @$p;
	next if($newfname eq $file);
	if(rename($file, $newfname)) {
	    $count++;
	} else {
	    print STDERR "Failed to mv ", shell_quote($file), " ", shell_quote($newfname), " : $!\n";
	}
    }
    print "$count file(s) renamed\n";
} else {
    print STDERR "\nIf you are satisfied with the renaming commands shown above,\njust add -do option to actually execute renaming\n\n";
}

=pod

=head1 NAME

renameregex.pl - Rename by regular expression

=head1 SYNOPSIS

renameregex.pl [options..] <Source regexp> <Dest regexp>

Options:
   -do              actually execute renaming commands
   -force           allow overwriting existing files / colliding destinations
   -help            brief help message
   -man             full documentation

=head1 OPTIONS

=over 8

=item B<-do>

Without -do option, renameregex.pl will just print out the series of 'mv' commands.
This is for design because we often give a wrong regular expression that may destroy many files.
-do option must be added only when you are satisfied with the series of commands shown.
Without -do nothing on disk is touched at all.

=item B<-force>

By default renameregex.pl refuses to run (even with -do) when a destination
file already exists, or when two source files would be renamed to the same
name, because either would silently destroy data. -force disables that check.

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
