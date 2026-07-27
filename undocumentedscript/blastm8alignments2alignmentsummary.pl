#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib $FindBin::Bin;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::Seq::SeqWithQuality;
use Bio::SeqIO;
use BLASTM8Parse;

my $flag_man     = 0;
my $flag_help    = 0;
my $debug        = 0;
my $flag_outputunique     = 0;
my $flag_outputnohit      = 0;
my $flag_outputrepetitive = 0;
my $flag_outputambiguous  = 0;
my $flag_outputmatches    = 0;
my $param_outputfilename       = undef;
my $param_leastBasesToDeclareUnique      = 60;
my $param_leastBasesToDeclareRepeat      = 300;
my $param_leastBasesRequiredForAlignemnt = 30;
my $param_maxdiagonal                    = 10;

GetOptions( 'help|?'    => \$flag_help,
	    'man'       => \$flag_man,
	    'lunique=i' => \$param_leastBasesToDeclareUnique,
            'lrepeat=i' => \$param_leastBasesToDeclareRepeat,
	    'lalign=i'  => \$param_leastBasesRequiredForAlignemnt,
	    'ldiag=i'   => \$param_maxdiagonal,
            'unique'    => \$flag_outputunique,
            'nohit'     => \$flag_outputnohit,
            'repeat'    => \$flag_outputrepetitive,
            'ambiguous' => \$flag_outputambiguous,
	    'match'     => \$flag_outputmatches,
	    'list=s'    => \$param_outputfilename,
	    'debug'     => \$debug
	    ) or pod2usage(2);

pod2usage(1) if($flag_help || @ARGV != 2);
pod2usage(-verbose => 2) if $flag_man; 

my $querySequenceFileName = shift;
my $blastM8ResultFileName = shift;

print STDERR "Least length required for unique span : $param_leastBasesToDeclareUnique bp\n";
print STDERR "Least length required for repeat span : $param_leastBasesToDeclareRepeat bp \n";
print STDERR "Least length required for alignment   : $param_leastBasesRequiredForAlignemnt bp\n";
print STDERR "Maximam width of diagonal             : $param_maxdiagonal bp\n";

my $ofh;
if(defined $param_outputfilename) {
    open($ofh, '>', $param_outputfilename) or die "Cannot open '$param_outputfilename': $!";
    print STDERR "File name to output list              : $param_outputfilename\n";
}
# -match writes whole BLAST -m8 lines. They must not be mixed into the plain
# read-name stream produced by -unique/-nohit/-repeat/-ambiguous, so they go to
# the -list file whenever one was given.
my $matchfh = defined $ofh ? $ofh : \*STDOUT;

my %readName2SequenceLength;
my @readNames;
{
    my $numSequences = 0;
    my $numTotalBases = 0;
    my $seqin  = Bio::SeqIO->new('-file'   => "< $querySequenceFileName",
				 '-format' => 'fasta'     );
    while(my $seqobj = $seqin->next_seq) {
	my $seqname = $seqobj->display_id;
	my $length = length($seqobj->seq);
	$numTotalBases += $readName2SequenceLength{$seqname} = $length;
	push(@readNames, $seqname);
	$numSequences++;
    }    
    print STDERR "$numSequences sequences loaded (${numTotalBases}bp in total)\n";
}

my $numNoHitSequences      = 0;
my $numAmbiguousSequences  = 0;
my $numUniqueSequences     = 0;
my $numRepetitiveSequences = 0;
my $numTotalSequences      = 0;

open(FH, '<', $blastM8ResultFileName) or die "Cannot open '$blastM8ResultFileName': $!";
my $currentQueryName = undef;
my @currentQueryMatches;
# Which query IDs the m8 file accounted for. The m8 file is NOT assumed to
# list the queries in FASTA order, nor exactly once each.
my %isQueryAccountedFor;

while(<FH>) {
    chomp;
    next if(/^\s*$/ or /^\s*#/);
    my $parsedBLASTM8Obj = parseBLASTM8ResultLine($_);
    my $queryID = $parsedBLASTM8Obj->{qid};
    my $alignmentLength = $parsedBLASTM8Obj->{alignlen};
    unless(defined $queryID && defined $alignmentLength) {
	print STDERR "WARNING: skipping malformed BLAST -m8 line: $_\n";
	next;
    }
    next unless($alignmentLength >= $param_leastBasesRequiredForAlignemnt);
    if(!defined $currentQueryName || $currentQueryName ne $queryID) {
	processOneQueryID($currentQueryName, \@currentQueryMatches, \%readName2SequenceLength) if(defined $currentQueryName);
	@currentQueryMatches = ();
        $currentQueryName = $queryID;
	if($isQueryAccountedFor{$queryID}) {
	    print STDERR "WARNING: query '$queryID' occurs in more than one block of '$blastM8ResultFileName'.\n";
	    print STDERR "         Each block is classified separately; sort the file by query ID to avoid this.\n";
	}
	$isQueryAccountedFor{$queryID} = 1;
    }
    push(@currentQueryMatches, $parsedBLASTM8Obj);
}
processOneQueryID($currentQueryName, \@currentQueryMatches, \%readName2SequenceLength) if(defined $currentQueryName);
close FH;

# Everything in the query FASTA that the m8 file never mentioned is a no-hit.
for my $readName (@readNames) {
    next if($isQueryAccountedFor{$readName});
    processNoHitQueryID($readName);
}

sub processNoHitQueryID {
    my $queryID = shift;
    $numTotalSequences++;
    print "$queryID\n" if($flag_outputnohit);
    print {$ofh} "$queryID\tNOHIT\n" if(defined $ofh);
    $numNoHitSequences++;
}

sub processOneQueryID {
    my $readName = shift;
    my $refArrayOfMatchesOfEachRead = shift;
    my $readName2SequenceLength = shift;

    my @readcoverage;
    unless(exists $readName2SequenceLength->{$readName}) {
	print STDERR "WARNING: not found '$readName'. Maybe you specified wrong file.\n";
	return;
    }
    $numTotalSequences++;
    if(@{$refArrayOfMatchesOfEachRead} == 0) {
        $numNoHitSequences++;
	print "$readName\n" if($flag_outputnohit);
	print {$ofh} "$readName\tNOHIT\n" if(defined $ofh);
        return;
    }
    my $readLength = $readName2SequenceLength->{$readName};
    for(my $i = 0; $i < $readLength; $i++) {
	$readcoverage[$i] = undef;
    }
    for my $refMatch (@{$refArrayOfMatchesOfEachRead}) {
	my $qst = $refMatch->{qstart} - 1; # to 0-origin
	my $qed = $refMatch->{qend} - 1;   # to 0-origin
	my $q1 = ($qst < $qed ? $qst : $qed);
	my $q2 = ($qst < $qed ? $qed : $qst);
	for(my $i = $q1; $i <= $q2; $i++) {
	    if(defined $readcoverage[$i]) {
		if(ref($readcoverage[$i]) eq 'HASH') {
		    $readcoverage[$i] = 2;
		} else {
		    $readcoverage[$i]++;
		}
	    } else {
		$readcoverage[$i] = $refMatch;
	    }
	}
    }
    my $maxs = undef;
    my $maxl = 0;
    {
        my $s = -1;
        my $l = 0;
        for(my $i = 0; $i < $readLength; $i++) {
            if(ref($readcoverage[$i]) eq 'HASH') {
                $s = $i if($l == 0);
                $l++;
            } else {
                if($maxl < $l) {
                    $maxl = $l;
                    $maxs = $s;
                }
                $s = -1; $l = 0;
            }
        }
	if($maxl < $l) {
	    $maxl = $l;
	    $maxs = $s;
	}
    }
    print STDERR "l = $maxl, s = ", (defined $maxs ? $maxs : 'n/a'), "\n" if ($debug);
    # -lunique (not -lalign) is the threshold that declares a read unique.
    my $isUnique     = $param_leastBasesToDeclareUnique <= $maxl;
    my $isRepetitive = 0;
    {
        my $s = -1;
        my $l = 0;
        for(my $i = 0; $i < $readLength; $i++) {
            if(defined $readcoverage[$i] && !(ref($readcoverage[$i]) eq 'HASH') && $readcoverage[$i] > 1) {
                $s = $i if($l == 0);
                $l++;
            } else {
                if($param_leastBasesToDeclareRepeat <= $l) {
                    $isRepetitive = 1; last;
                }
                $s = -1; $l = 0;
            }
        }
	if($param_leastBasesToDeclareRepeat <= $l) {
	    $isRepetitive = 1;
	}
    }
    my @fhashkeys;
    if($isUnique && !$isRepetitive) { # further examination
	my %numOccurence;
	{
	    # Collect the distinct matches that uniquely cover at least one base,
	    # in the order in which they first appear along the read. Hash
	    # references cannot be deduplicated with sort()/'!=', so use a
	    # 'seen' hash keyed on the stringified reference.
	    my %isAlreadySeen;
	    for(my $i = 0; $i < $readLength; $i++) {
		my $c = $readcoverage[$i];
		next unless(ref($c) eq 'HASH');
		$numOccurence{$c}++;
		push(@fhashkeys, $c) unless($isAlreadySeen{$c}++);
	    }
        }
	if(@fhashkeys > 1) {
	    # Only matches that uniquely cover at least -lunique bases are
	    # significant; a one base island must not demote the read.
	    my @significant = grep { $numOccurence{$_} >= $param_leastBasesToDeclareUnique } @fhashkeys;
	    if(@significant > 1) {
		my $fhash = $significant[0];
		my $ftarget_sequence = $fhash->{sid};
		my $fisplusstrand    = $fhash->{sstart} < $fhash->{send};
		my $fdiagonal        = $fhash->{sstart} + ($fisplusstrand ? - $fhash->{qstart} : + $fhash->{qstart});
		for(my $fpointer = 1; $fpointer < @significant; $fpointer++) {
		    my $h = $significant[$fpointer];
		    my $htarget_sequence = $h->{sid};
		    my $hisplusstarnd    = $h->{sstart} < $h->{send};
		    my $hdiagonal        = $h->{sstart} + ($hisplusstarnd ? - $h->{qstart} : + $h->{qstart});
		    if( ($htarget_sequence ne $ftarget_sequence) ||
			($hisplusstarnd != $fisplusstrand) ||
			(abs($fdiagonal - $hdiagonal) >= $param_maxdiagonal) ) {
			$isUnique = 0;
			last;
		    }
		}
	    }
	}
    }
    if(!$isUnique && !$isRepetitive) {
	$numNoHitSequences++;
	print "$readName\n" if($flag_outputnohit);
	print {$ofh} "$readName\tNOHIT\n" if(defined $ofh);
    } elsif($isUnique && !$isRepetitive) {
        $numUniqueSequences++;
	print "$readName\n" if($flag_outputunique);
	print {$ofh} "$readName\tUNIQUE\n" if(defined $ofh);
	if($flag_outputmatches) {
	    for(@fhashkeys){
		print {$matchfh} createBLASTM8Line($_), "\tUNIQUE\n";
	    }
	}
    } elsif(!$isUnique && $isRepetitive) {
        $numRepetitiveSequences++;
	print "$readName\n" if($flag_outputrepetitive);
	print {$ofh} "$readName\tREPEAT\n" if(defined $ofh);
	if($flag_outputmatches) {
	    print {$matchfh} createBLASTM8Line($_), "\tREPEAT\n" for(@{$refArrayOfMatchesOfEachRead});
	}
    } else {
        # print STDERR join(',', @readcoverage), "\n";
        $numAmbiguousSequences++;
	print "$readName\n" if($flag_outputambiguous);
	print {$ofh} "$readName\tAMBIGUOUS\n" if(defined $ofh);
	if($flag_outputmatches) {
	    print {$matchfh} createBLASTM8Line($_), "\tAMBIGUOUS\n" for(@{$refArrayOfMatchesOfEachRead});
	}
    }
}

if(defined $ofh) {
    close($ofh) or die "Cannot write '$param_outputfilename': $!";
}

print STDERR "\nOf all $numTotalSequences, \n";
print STDERR "\tNo hit     : $numNoHitSequences\n";
print STDERR "\tUnique     : $numUniqueSequences\n";
print STDERR "\tRepetitive : $numRepetitiveSequences\n";
print STDERR "\tAmbiguous  : $numAmbiguousSequences\n";

=pod

=head1 NAME

blastm8alignments2alignmentsummary - Create a summary of uniqueness of reads.

=head1 SYNOPSIS

blastm8alignments2alignmentsummary.pl [options] <query sequences> <blast -m8 result>

Options:
   -help            brief help message
   -man             full documentation

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-lalign>

Set the least number of bases for each (potentially partial) alignment.
If -align=100 is specified, alignments less than 100bp are thrown away, and will not be used for further analysis.

=item B<-lunique>

Set the least number of bases to classify a sequence as unique. For example,
when -unique=50 is set, a sequence is marked unique if and only if
at least consecutive 50bp of the given sequence is mapped onto one unique 
location of the target sequence.

=item B<-lrepeat>

Set the least number of bases to classify a sequence as repeat. For example,
when -lrepeat=100 is set, a sequence is marked repeat if and only if
at least consecutive 100bp of the given sequence is mapped onto multiple
locations of some target sequence(s).

=item B<-ldiag>

Set the maximum width of a diagonal of a hit chain. The default value is set to 10.

=item B<-unique>

If -unique option is given, only the names of unique sequences are output.

=item B<-nohit>

If -unique option is given, only the names of no-hit sequences are output.

=item B<-repeat>

If -repeat option is given, only the names of repetitive sequences are output.

=item B<-ambiguous>

If -ambiguous option is given, only the names of repetitive sequences are output.

=item B<-list>

When -list=filename option is given, a classification (nohit/repeat/ambiguous/list) is output for each sequence to the specified file.

=item B<-match>

When -match option is given, BLAST -m8 matches are output. Only difference between the input and the output is that if the sequence has unique match, non-unique matches are filtered out.
The matches are written to the file given by B<-list> if that option is used, and to the standard output otherwise, so that they are never interleaved with
the plain sequence names printed by B<-unique>/B<-nohit>/B<-repeat>/B<-ambiguous>.

=back

=head1 DESCRIPTION

B<This program> will read the given blast -m8 output and generates a
summary of uniquness of each sequence.

=cut
