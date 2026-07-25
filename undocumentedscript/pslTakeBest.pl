#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

my $flag_man        = 0;
my $flag_help       = 0;
my $debug           = 0;

my $param_length    = 50;
my $param_neartop   = 0.03;
my $param_coverage  = 0;
my $param_identity  = 0;

GetOptions( 'help|?'           => \$flag_help,
	    'man'              => \$flag_man,
	    'debug'            => \$debug,
	    'length=i'         => \$param_length,
	    'neartop=f'        => \$param_neartop,
	    'identity=f'       => \$param_identity,
	    'coverage=f'       => \$param_coverage
    ) or pod2usage(2);
pod2usage(1) if $flag_help;
pod2usage(-verbose => 2) if $flag_man;
pod2usage(2) if(@ARGV != 1);

my $pslfilename = shift;

unless(open FH, "<", $pslfilename) {
    print STDERR "Could not open '$pslfilename'\n";
    exit 1;
}

my $is_first_line = 1;
my $prev_qname    = undef;
my %seen_qname    = ();
my $warned_unsorted = 0;
my @queries = ();
while(<FH>) {
    chomp;
    next if(/^\s*$/);
    my ($match, $mismatch, $repmatch, $num_n, $qgap_count, $qgap_bases, $tgap_count, $tgap_bases, $strand, $qname, $qsize, $qstart, $qend, $tname, $tsize, $tstart, $tend, $block_count, $block_sizes, $qstarts, $tstarts)
	= split(/\t/);
    if($is_first_line && defined $match && $match =~ /^psLayout/) {
	print; print "\n";
	for(my $i = 0; $i < 4; $i++) {
	    my $l = <FH>;
	    unless(defined $l) {
		print STDERR "WARNING: '$pslfilename' ends in the middle of the psl header\n";
		last;
	    }
	    print $l;
	}
	$is_first_line = 0;
	next;
    }
    $is_first_line = 0;
    unless(defined $qname && defined $qsize && defined $qstart && defined $qend) {
	print STDERR "WARNING: skipping malformed psl line $.\n";
	next;
    }
    if(!defined $prev_qname || $prev_qname ne $qname) {
	flush_queries(@queries) if(defined $prev_qname);
	# The best-hit filter works on one contiguous run of a query name at a
	# time, so the input must already be grouped by query name. A psl
	# sorted by target (or concatenated from several runs) would otherwise
	# be filtered independently per run and emit duplicated hits.
	if($seen_qname{$qname} && !$warned_unsorted) {
	    print STDERR "WARNING: '$pslfilename' is not grouped by query name (e.g. '$qname' appears again after another query).\n";
	    print STDERR "WARNING: each contiguous run of a query is filtered independently, so duplicated hits may be emitted.\n";
	    print STDERR "WARNING: sort the psl by query name (e.g. 'sort -k10,10') first.\n";
	    $warned_unsorted = 1;
	}
	$seen_qname{$qname} = 1;
	$prev_qname = $qname;
	@queries    = ();
    }
    # qStart/qEnd are 0-based half-open, so the length is qEnd - qStart.
    my $alignment_length = $qend - $qstart;
    # Real identity: matched bases over the bases that are actually aligned
    # (gaps and Ns excluded). The old 'match / qsize' was a coverage measure,
    # which -coverage already tests below.
    my $aligned_bases = $match + $mismatch + $repmatch;
    push(@queries, {line       => $_,
		    length     => $alignment_length,
		    qlength    => $qsize,
		    match      => $match,
		    matchratio => $aligned_bases > 0 ? ($match + $repmatch) / $aligned_bases : 0,
		    mismatch   => $mismatch}
	);
}
flush_queries(@queries) if(defined $prev_qname);
close FH or print STDERR "WARNING: error while closing '$pslfilename': $!\n";

sub flush_queries
{
    my @filtered;
    for(@_) {
	next if($_->{length}     < $param_length);
	next if($_->{length}     < $_->{qlength} * $param_coverage);
	next if($_->{matchratio} < $param_identity);
	push(@filtered, $_);
    }
    my @sorted = sort {$b->{match} <=> $a->{match}} @filtered;
    return unless(@sorted);
    my $best_match      = $sorted[0]->{match};
    # With no matched base at all there is nothing worth reporting; without
    # this guard the threshold would be 0 and every zero-match record passes.
    return if($best_match <= 0);
    my $threshold_match = $best_match * ( 1 - $param_neartop );
    for(@sorted) {
	last if($_->{match} < $threshold_match);
	print $_->{line}, "\n";
    }
}

=pod

=head1 NAME

pslTakeBest.pl - take the best(s) alignments in psl file

=head1 SYNOPSIS

pslTakeBest.pl [options...] <psl file>

Options:
   -help             brief help message
   -man              full documentation

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-length>

Specify the minimum alignment length (qEnd - qStart; psl coordinates are
0-based half-open). The defaut is 50.

=item B<-neartop>

Specify the threshold value such that alignments with the number of matched bases worse by this value than the best alignment are kept in the output.
The default value is 0.03.

=item B<-identity>

Specify the minimum identity required to be output.
Identity is (matches + repMatches) / (matches + repMatches + misMatches),
i.e. the fraction of the actually aligned bases that agree; gaps and Ns are
not counted. The default value is 0.

=item B<-coverage>

Specify the minimum coverage required to be output.
Alignments that are not covered to this threshold are discarded.
The default value is 0.

=back

=head1 DESCRIPTION

B<pslTakeBest.pl> will first process the input psl file.

=cut

