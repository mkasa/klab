# BLASTM8Parse.pm - parse BLAST -m8 (a.k.a. -outfmt 6/7, tabular) output.
#
# Nothing is exported by default; import exactly the subs you use, or call
# them fully qualified as BLASTM8Parse::parseBLASTM8ResultLine(...).
#
# use BLASTM8Parse qw(parseBLASTM8ResultFile);
#
# my @matches = parseBLASTM8ResultFile("blastm8.txt");
# for(@matches) {
#   my $query_id   = $_->{qid};
#   my $subject_id = $_->{sid};
#   ... so on.
# }

package BLASTM8Parse;

use strict;
use warnings;

use Carp;
use Exporter 'import';

our @EXPORT_OK = qw(
	parseBLASTM8ResultFile
	parseBLASTM8ResultLine
	createBLASTM8Line
);
our %EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

sub parseBLASTM8ResultFile($) {
    my $m8file = shift;
    open(my $fh, '<', $m8file) or croak "Cannot open BLAST m8 file '$m8file': $!";
    my @retval;
    while(<$fh>) {
	next if(/^\s*$/);  # blank line
	next if(/^\s*#/);  # comment line (BLAST -outfmt 7)
	my $br = parseBLASTM8ResultLine($_);
	push(@retval, $br);
    }
    close $fh;
    return @retval;
}

sub parseBLASTM8ResultLine($) {
    my $line = shift;
    $line = '' unless(defined $line);
    $line =~ s/[\r\n]+\z//;
    my @args = split(/\t/, $line);
    my %rval;
    $rval{'qid'}      = $args[ 0]; # Query ID
    chop $rval{'qid'} if(defined $rval{'qid'} && $rval{'qid'} =~ /\|$/);
    $rval{'sid'}      = $args[ 1]; # Subject ID
    chop $rval{'sid'} if(defined $rval{'sid'} && $rval{'sid'} =~ /\|$/);
    $rval{'identity'} = $args[ 2]; # Identity
    $rval{'alignlen'} = $args[ 3]; # Alignment length
    $rval{'mismatch'} = $args[ 4]; # Number of mismatch
    $rval{'gapopen'}  = $args[ 5]; # Number of gap opening
    $rval{'qstart'}   = $args[ 6]; # Query start position   (1-origin)
    $rval{'qend'}     = $args[ 7]; # Query end position     (1-origin)
    $rval{'sstart'}   = $args[ 8]; # Subject start position (1-origin)
    $rval{'send'}     = $args[ 9]; # Subject end position   (1-origin)
    $rval{'eval'}     = $args[10]; # e-value
    $rval{'bitscore'} = $args[11]; # bit-score
    return \%rval;
}

sub createBLASTM8Line($) {
    my $parsedobj = shift;
    my @line = ( $parsedobj->{'qid'},
		 $parsedobj->{'sid'},
		 $parsedobj->{'identity'},
		 $parsedobj->{'alignlen'},
		 $parsedobj->{'mismatch'},
		 $parsedobj->{'gapopen'},
		 $parsedobj->{'qstart'},
		 $parsedobj->{'qend'},
		 $parsedobj->{'sstart'},
		 $parsedobj->{'send'},
		 $parsedobj->{'eval'},
		 $parsedobj->{'bitscore'}
		);
    return join("\t", map { defined $_ ? $_ : '' } @line);
}

1;
