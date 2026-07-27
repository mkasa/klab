#!/usr/bin/env perl
# NOTE: '#!/usr/bin/env perl -w' passes "perl -w" to env as a SINGLE argument
#       on Linux and fails with "env: 'perl -w': No such file or directory";
#       'use warnings' below does the same job portably.

use strict;
use warnings;

use FindBin;
use lib $FindBin::Bin;   # Primer3.pm sits next to this script and '.' is no
                         # longer in @INC since perl 5.26.

use Getopt::Long;
use Pod::Usage;
use Primer3;

my $flag_man      = 0;
my $flag_help     = 0;
my $flag_noheader = 0;
my $debug         = 0;
my $param_target        = '';
my $param_product_range = '70-110';
my $param_maxn          = 0;
my $param_mintm         = undef;
my $param_maxtm         = undef;
my $param_opttm         = undef;
my $param_mingc         = undef;
my $param_maxgc         = undef;
my $param_optgc         = undef;
my $param_minsize       = undef;
my $param_maxsize       = undef;
my $param_optsize       = undef;

GetOptions( 'help|?'         => \$flag_help,
		    'man'            => \$flag_man,
		    'noheader'       => \$flag_noheader,
		    'debug'          => \$debug,
		    'target=s'       => \$param_target,
		    'productrange=s' => \$param_product_range,
		    'maxn=i'         => \$param_maxn,
		    'mintm=f'        => \$param_mintm,
		    'maxtm=f'        => \$param_maxtm,
		    'opttm=f'        => \$param_opttm,
		    'mingc=f'        => \$param_mingc,
		    'maxgc=f'        => \$param_maxgc,
		    'optgc=f'        => \$param_optgc,
		    'minsize=i'      => \$param_minsize,
		    'maxsize=i'      => \$param_maxsize,
		    'optsize=i'      => \$param_optsize
		    ) or pod2usage(2);

my $sequence = shift;
# "GTAGTCAGTAGACNATGACNACTGACGATGCAGACNACACACACACACACAGCACACAGGTATTAGTGGGCCATTCGATCCCGACCCAAATCGATAGCTACGATGACG"
# -man must be handled before the 'no sequence given' fallback, otherwise
# 'primerdesign.pl --man' only prints the one-line usage.
pod2usage(-verbose => 2) if $flag_man;
$flag_help = 1 unless(defined $sequence);
pod2usage(1) if $flag_help;

my $obj = new Primer3();
my $primer3path = $obj->get_primer3_path();
unless(defined $primer3path && $primer3path ne '') {
    print STDERR "ERROR: primer3 is not installed (neither 'primer3_core' nor 'primer3' was found in PATH).\n";
    exit 1;
}
print STDERR "Primer3Path $primer3path\n" if($debug);

#
my %argHash = (
	id     => "exampleID",
    seq    => $sequence,
    maxn   => $param_maxn,
    product_range => $param_product_range
);

if($param_target) {
	# Anchored: '-target=chr1:10-20x' used to be silently accepted.
	if($param_target =~ /^\s*(\d+)-(\d+)\s*$/) {
		my $start = $1;
		my $end   = $2;
		if($end < $start) {
			print STDERR "Invalid range '$param_target': the end ($end) is before the start ($start)\n";
			exit 1;
		}
		$argHash{target_start} = $start;
		$argHash{target_end}   = $end;
		if($debug) {
			print STDERR "Range : ${start}-${end}\n";
		}
	} else {
		print STDERR "Cannot understand the range '$param_target'\n";
		print STDERR "Range must be in a format like '70-100'\n";
		exit 1;
	}
}

$argHash{primer_mintm}   = $param_mintm         if(defined $param_mintm);
$argHash{primer_maxtm}   = $param_maxtm         if(defined $param_maxtm);
$argHash{primer_opttm}   = $param_opttm         if(defined $param_opttm);
$argHash{primer_mingc}   = $param_mingc         if(defined $param_mingc);
$argHash{primer_maxgc}   = $param_maxgc         if(defined $param_maxgc);
$argHash{primer_optgc}   = $param_optgc         if(defined $param_optgc);
$argHash{primer_minsize} = $param_minsize       if(defined $param_minsize);
$argHash{primer_maxsize} = $param_maxsize       if(defined $param_maxsize);
$argHash{primer_optsize} = $param_optsize       if(defined $param_optsize);

my $primers = $obj->designprimer(%argHash);

sub nvl($)
{
	my $v = shift;
	return defined $v ? $v : '';
}

my $error = $primers->{error};
if(defined $error && $error ne '') {
	print STDERR "ERROR: $error\n";
	exit 1;
}
unless($flag_noheader) {
	print "#ForwardPrimer,ReversePrimer,ForwardTm,ReverseTm,ForwardPosition,ForwardLength,ReversePosition,ReverseLength\n";
}
for(@{$primers->{primer} || []}){
    print nvl($_->{leftsequence}), ",", nvl($_->{rightsequence}), ",";
    print nvl($_->{lefttm}), ",", nvl($_->{righttm}), ",";
    print nvl($_->{leftposition}), ",", nvl($_->{leftlength}), ",", nvl($_->{rightposition}), ",", nvl($_->{rightlength}), "\n";
}

=pod

=head1 NAME

primerdesign.pl - Designing a primer pair.

=head1 SYNOPSIS

primerdesign.pl [options] <sequence>

Options:
   -help            brief help message
   -man             full documentation

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<primerdesign.pl> will design a primer pair for the given sequence.
Various conditions can be specified at command line.

=cut

