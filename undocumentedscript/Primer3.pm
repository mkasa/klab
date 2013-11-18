use strict;

=pod

=head1 NAME

Primer3.pm - Wrapper for Primer3

=head1 SYNOPSIS

	use Primer3;

	my $primer3 = new Primer3();

		my $primers = $primer3->designprimer(
		id  => "Sequence1",
		seq => "ACATCAGGCGACGGAGCGACGAGCGAGCAGTCACACG",
		primer_mintm => 60.0,
		primer_maxtm => 70.0,
		primer_opttm => 75.5,
		primer_maxdifftm => 10.0,
		primer_mingc => 40.0,
		primer_maxgc => 60.0,
		primer_optgc => 45.0,
		primer_minsize => 16,
		primer_maxsize => 32,
		primer_optsize => 21,
		product_range => '75-100',
		target_start  => 37 # 0-origin, inclusive
		target_end    => 48 # 0-origin, exclusive
    );
    die if($primers->{error});
    for(@{$primers->{primer}){
   		print "$_->{leftsequence},$_->{rightsequence}\n";
   		print "$_->{lefttm},$_->{righttm}\n";
    }

=head1 DESCRIPTION

B<Primer3.pm> can execute Primer3 and get the results.

=head1 OPTION

	You can specify multiple targets like this,
	target => "1,20 300,20 400,20"

	You can also exclude specific regions from primer design
	excludedregion => "4,499, 1000,300"

	These region lists is
	startpos,len startpos,len ...
	where startpos is 0-origin index.

	You can see the number of the designed primer pairs
	$retval->{numprimerpairs}

=cut

package Primer3;

use Carp;
use Env::Path;
use FileHandle;
use File::Temp;

sub new($)
{
	my $class = shift;
	my $pathobj = Env::Path->PATH;
	my @primer3s = $pathobj->Whence('primer3*');
	my $primer3path;
	if(@primer3s) {
		$primer3path = shift @primer3s;
	} else {
		my @primer3cores = $pathobj->Whence('primer3_core*');
		$primer3path = shift @primer3cores;
	}
	my $self = {
		primer3path => $primer3path
	};
    bless($self, $class);
    return($self);
}

sub get_primer3_path()
{
	my $self = shift;
	return $self->{primer3path};
}

sub DESTROY{
    my $self = shift;
}

sub designprimer($$)
{
	my $self = shift;
	my %param = @_;
	my $sequenceid = $param{id};
	my $sequence   = $param{seq};
	my $primer3path = $self->{primer3path};

    croak "Sequence ID is not defined. " unless(defined $sequenceid);

	my $tmp = new File::Temp( UNLINK => 0, SUFFIX => '.dat' );
	print $tmp "PRIMER_SEQUENCE_ID=$sequenceid\n";
	print $tmp "SEQUENCE=$sequence\n";

    print $tmp "PRIMER_PICK_ANYWAY=1\n" if $param{pickanyway};

	my $target_start = $param{tstart};
	my $target_end   = $param{tend};
	if(defined $target_end && defined $target_start) {
		my $length = $target_end - $target_start;
		print $tmp "TARGET=$target_start,$length\n";
	} else {
		my $target = $param{target};
		print $tmp "TARGET=$target\n";
	}

	my $primer_opt = $param{primeroptlength};
	if(defined $primer_opt) {
		print $tmp "PRIMER_OPT_SIZE=$primer_opt\n";
	}

	my $excluded_region = $param{excludedregion};
	if(defined $excluded_region) {
		print $tmp "EXCLUDED_REGION=$excluded_region\n";
	}

	print $tmp "PRIMER_PRODUCT_SIZE_RANGE=$param{product_range}\n" if($param{product_range});
	print $tmp 	"PRIMER_NUM_NS_ACCEPTED=$param{maxn}\n" if($param{maxn});
#PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION=37,21
#PRIMER_MIN_SIZE=15
#PRIMER_MAX_SIZE=21
	print $tmp "PRIMER_FILE_FLAG=1\n";
	print $tmp "PRIMER_EXPLAIN_FLAG=1\n" if($param{explain});
    print $tmp "PRIMER_PRODUCT_MAX_TM=$param{product_maxtm}\n" if($param{product_maxtm});
    print $tmp "PRIMER_PRODUCT_MIN_TM=$param{product_mintm}\n" if($param{product_mintm});
    print $tmp "PRIMER_PRODUCT_OPT_TM=$param{product_opttm}\n" if($param{product_opttm});
    print $tmp "PRIMER_PRODUCT_OPT_SIZE=$param{product_optsize}\n" if($param{product_optsize});
    print $tmp "PRIMER_OPT_SIZE=$param{primer_optsize}\n" if($param{primer_optsize});
    print $tmp "PRIMER_MIN_SIZE=$param{primer_minsize}\n" if($param{primer_minsize});
    print $tmp "PRIMER_MAX_SIZE=$param{primer_maxsize}\n" if($param{primer_maxsize});
    print $tmp "PRIMER_OPT_TM=$param{primer_opttm}\n" if($param{primer_opttm});
    print $tmp "PRIMER_MIN_TM=$param{primer_mintm}\n" if($param{primer_mintm});
    print $tmp "PRIMER_MAX_TM=$param{primer_maxtm}\n" if($param{primer_maxtm});
    print $tmp "PRIMER_MAX_DIFF_TM=$param{primer_maxdifftm}\n" if($param{primer_maxdifftm});
    print $tmp "PRIMER_MIN_GC=$param{primer_mingc}\n" if($param{primer_mingc});
    print $tmp "PRIMER_MAX_GC=$param{primer_maxgc}\n" if($param{primer_maxgc});
    print $tmp "PRIMER_OPT_GCPERCENT=$param{primer_optgc}\n" if($param{primer_optgc});
    print $tmp "PRIMER_SALT_CONC=$param{primer_saltcon}\n" if($param{primer_saltcon});
    print $tmp "PRIMER_DNA_CONC=$param{primer_dnacon}\n" if($param{primer_dnacon});
	print $tmp "=\n";
	$tmp->close();

    sub tonumeric($) {
    	my $value = shift;
    	return 0 unless(defined $value);
    	return $value + 0;
    }
	my @elines = `$primer3path < $tmp`;
	# print "file : $tmp\n"; # for debug
	my $retval;
	$retval->{error} = "";
	$retval->{numprimerpairs} = 0;
 	for(@elines) {
 		# print;
 		chomp;
 		if(/^PRIMER_ERROR=(.*)/){
			$retval->{error} = $1;
 		} elsif(/^PRIMER_LEFT(_(\d+))?_PENALTY=([\d\.]*)/) {
			my $index = tonumeric($2);
			$retval->{primer}->[$index]->{leftpenalty} = $3;
 		} elsif(/^PRIMER_RIGHT(_(\d+))?_PENALTY=([\d\.]*)/) {
			my $index = tonumeric($2);
			$retval->{primer}->[$index]->{rightpenalty} = $3;
 		} elsif(/^PRIMER_LEFT(_(\d+))?_SEQUENCE=([A-Z]*)/) {
			my $index = tonumeric($2);
			$retval->{primer}->[$index]->{leftsequence} = $3;
			$retval->{numprimerpairs}++;
 		} elsif(/^PRIMER_RIGHT(_(\d+))?_SEQUENCE=([A-Z]*)/) {
			my $index = tonumeric($2);
			$retval->{primer}->[$index]->{rightsequence} = $3;
 		} elsif(/^PRIMER_LEFT(_(\d+))?=(\d*),(\d*)/) {
			my $index = tonumeric($2);
			$retval->{primer}->[$index]->{leftposition} = $3;
			$retval->{primer}->[$index]->{leftlength} = $4;
 		} elsif(/^PRIMER_RIGHT(_(\d+))?=(\d*),(\d*)/) {
			my $index = tonumeric($2);
			$retval->{primer}->[$index]->{rightposition} = $3;
			$retval->{primer}->[$index]->{rightlength} = $4;
 		} elsif(/^PRIMER_LEFT(_(\d+))?_TM=([\d\.]*)/) {
			my $index = tonumeric($2);
			$retval->{primer}->[$index]->{lefttm} = $3;
 		} elsif(/^PRIMER_RIGHT(_(\d+))?_TM=([\d\.]*)/) {
			my $index = tonumeric($2);
			$retval->{primer}->[$index]->{righttm} = $3;
 		} elsif(/^PRIMER_PRODUCT_SIZE(_(\d+))?=(\d*)/) {
			my $index = tonumeric($2);
			$retval->{primer}->[$index]->{productsize} = $3;
 		} elsif(/^PRIMER_LEFT(_(\d+))?_GC_PERCENT=([\d\.]*)/) {
			my $index = tonumeric($2);
			$retval->{primer}->[$index]->{leftgc} = $3;
		} elsif(/^PRIMER_RIGHT(_(\d+))?_GC_PERCENT=([\d\.]*)/) {
			my $index = tonumeric($2);
			$retval->{primer}->[$index]->{rightgc} = $3;
		}
	}
	unlink $tmp;
	return $retval;
}

1;
