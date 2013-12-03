#!/usr/bin/env perl

use strict;

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
my @queries = ();
while(<FH>) {
    chomp;
    my ($match, $mismatch, $repmatch, $num_n, $qgap_count, $qgap_bases, $tgap_count, $tgap_bases, $strand, $qname, $qsize, $qstart, $qend, $tname, $tsize, $tstart, $tend, $block_count, $block_sizes, $qstarts, $tstarts)
	= split(/\t/);
    if($is_first_line && $match =~ /^psLayout/) {
	print; print "\n";
	for(my $i = 0; $i < 4; $i++) {
	    my $l = <FH>; print $l;
	}
    } else {
	$is_first_line = 0;
	if($prev_qname ne $qname) {
	    flush_queries(@queries);
	    $prev_qname = $qname;
	    @queries    = ();
	}
	push(@queries, {line       => $_,
			length     => $qend - $qstart + 1,
			qlength    => $qsize,
			match      => $match,
			matchratio => $qsize == 0 ? 0 : $match / $qsize,
			mismatch   => $mismatch}
	    );
    }
}
flush_queries(@queries) if(defined $prev_qname);
close FH;

sub flush_queries(@)
{
    my @filtered;
    for(@_) {
	next if($_->{length}     < $param_length);
	next if($_->{length}     < $_->{qlength} * $param_coverage);
	next if($_->{matchratio} < $param_identity);
	push(@filtered, $_);
    }
    my @sorted = sort {$b->{match} <=> $a->{match}} @filtered;
    my $best_match      = scalar(@sorted) > 0 ? $sorted[0]->{match} : 0;
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

Specify the minimum alignment length. The defaut is 50.

=item B<-neartop>

Specify the threshold value such that alignments with the number of matched bases worse by this value than the best alignment are kept in the output.
The default value is 0.03.

=item B<-identity>

Specify the minimum identity required to be output.
The default value is 0.

=item B<-coverage>

Specify the minimum coverage required to be output.
Alignments that are not covered to this threshold are discarded.
The default value is 0.

=back

=head1 DESCRIPTION

B<pslTakeBest.pl> will first process the input psl file.

=cut

