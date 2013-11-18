#!/usr/bin/env perl

use strict;

use Getopt::Long;
use Pod::Usage;
use Bio::Perl;
use Bio::Seq;
use Bio::Tools::RestrictionEnzyme;

my $flag_man     = 0;
my $flag_help    = 0;
my $debug        = 0;

GetOptions( 'help|?'  => \$flag_help,
	    'man'     => \$flag_man,
	    'debug'   => \$debug
	    ) or pod2usage(2);
pod2usage(1) if $flag_help;
pod2usage(-verbose => 2) if $flag_man;

sub parseBLASTCSV($)
{
	my $line = shift;
	my ($qname, $qlen, $qst, $qed, $sname, $slen, $sst, $sed, $nummatch, $alignlen, $qastr, $sastr) = split(/,/);
	return {
		queryid      => $qname,
		querylen     => $qlen,
		querystart   => $qst,
		queryend     => $qed,
		subjectid    => $sname,
		subjectstart => $sst,
		subjectend   => $sed,
		matchlen     => $nummatch,
		alignlen     => $alignlen,
		queryalignmentstring   => $qastr,
		subjectalignmentstring => $sastr
	};
}

sub adjust_for_gap($$)
{
	my $gappedStr = shift;
	my @splitGappedStr = split(//, $gappedStr);

	my $positionListRef = shift;
	my $indexPositionListRef = 0;
	my @retval;

    return \@retval if(@$positionListRef <= 0);
    my $positionOnSeqStr    = 0;
	my $positionOnGappedStr = 0;
	for(@splitGappedStr) {
		if($$positionListRef[$indexPositionListRef] <= $positionOnSeqStr) {
			push(@retval, $positionOnGappedStr);
			last if(++$indexPositionListRef >= @$positionListRef);
		}
		$positionOnGappedStr++;
		unless($_ eq '-' || $_ eq ' ') {
			$positionOnSeqStr++;
		}
	}
	return \@retval;
}

sub cutstr($)
{
	my $listRef = shift;
	my $prevPos = 0;
	my $retval  = '';
	for(@$listRef) {
		my $diff = $_ - $prevPos;
		$retval .= ' 'x$diff if($diff > 0);
		$retval .= '|';
		$prevPos = $_ + 1;
	}
	return $retval;
}

sub diff_of_array($$)
{
	my $array1Ref = shift;
	my $array2Ref = shift;
	my @diff;
	my $i = 0;
	my $j = 0;
	while($i < @$array1Ref && $j < @$array2Ref) {
		if($$array1Ref[$i] < $$array2Ref[$j]) {
		    push(@diff, $$array1Ref[$i]);
			$i++;
		} elsif ($$array1Ref[$i] > $$array2Ref[$j]) {
		    push(@diff, $$array2Ref[$j]);
			$j++;
		} else {
		    $i++; $j++;
		}
	}
	push(@diff, $$array1Ref[$i++]) while($i < @$array1Ref);
	push(@diff, $$array2Ref[$j++]) while($j < @$array2Ref);
	return \@diff;
}

sub intersection_of_array($$)
{
	my $array1Ref = shift;
	my $array2Ref = shift;
	my @intersection;
	my $i = 0;
	my $j = 0;
	while($i < @$array1Ref && $j < @$array2Ref) {
		if($$array1Ref[$i] < $$array2Ref[$j]) {
			$i++;
		} elsif ($$array1Ref[$i] > $$array2Ref[$j]) {
			$j++;
		} else {
		    push(@intersection, $$array1Ref[$i]);
		    $i++; $j++;
		}
	}
	return \@intersection;
}

sub nmask_sequence($$)
{
	my $str1 = shift;
	my $str2 = shift;
	my $retval;
	my @str1arr = split(//, $str1);
	my @str2arr = split(//, $str2);
	for(my $i = 0; $i < @str1arr; $i++) {
		if($str2arr[$i] =~ /n/i && $str1arr[$i] ne '-') {
			$retval .= 'n';
		} else {
			$retval .= $str1arr[$i];
		}
	}
	return $retval;
}

my @restrictionEnzymeObjects;
{
	my @restrictionEnzymes = qw(
		HinfI MboI MseI MspI MvaI RsaI Sau96I ScrFI TaqI TspEI HaeIII DdeI EcoRI BamHI HindIII DraI EcoRV PstI SacI XhoI
	);
	for(@restrictionEnzymes) {
	    push(@restrictionEnzymeObjects,
	        new Bio::Tools::RestrictionEnzyme(-name => $_));
	}
	# TODO: ’Ç‰Á.
	# MnlI
}

print '<?xml version="1.0" encoding="US-ASCII" ?>', "\n";
print "<rflpmap>\n";
while(<>){
	chomp;
	next if(/^#/);
	my $bmatch = parseBLASTCSV($_);
	print "  <matchsegment queryid=\"$bmatch->{queryid}\" subjectid=\"$bmatch->{subjectid}\" ";
	print "querystart=\"$bmatch->{querystart}\" queryend=\"$bmatch->{queryend}\" ";
	print "subjectstart=\"$bmatch->{subjectstart}\" subjectend=\"$bmatch->{subjectend}\" >\n";
	print "    <queryseq>$bmatch->{queryalignmentstring}</queryseq>\n";
	print "    <subjectseq>$bmatch->{subjectalignmentstring}</subjectseq>\n";
	print "    <rflps>\n";
	for my $reobj (@restrictionEnzymeObjects) {
		my $rname = $reobj->name();
		my $qalignstr = $bmatch->{queryalignmentstring};
		my $salignstr = $bmatch->{subjectalignmentstring};
		my $qstr = $qalignstr;
		my $sstr = $salignstr;
		$qstr = nmask_sequence($qstr, $sstr);
		$sstr = nmask_sequence($sstr, $qstr);
		$qstr =~ s/-//g;
		$sstr =~ s/-//g;
		if($debug) {
			print "$qstr\n";
			print "$sstr\n";
		}
		my $qseqobj = new Bio::Seq(-seq => $qstr);
		my $sseqobj = new Bio::Seq(-seq => $sstr);
		my $qcutsites_qstr = $reobj->cut_locations($qseqobj);
		my $scutsites_sstr = $reobj->cut_locations($sseqobj);
		if($debug) {
	        print "For enzyme $rname\n";
	        print " Ungap:\n";
	        print "  ", join(',', @$qcutsites_qstr), "\n";
	        print "  ", join(',', @$scutsites_sstr), "\n";
	    }
		my $qcutsites_qalignstr = adjust_for_gap($qalignstr, $qcutsites_qstr);
		my $scutsites_salignstr = adjust_for_gap($salignstr, $scutsites_sstr);
		if($debug) {
	        print " Gapped:\n";
	        print "  ", join(',', @$qcutsites_qalignstr), "\n";
	        print "  ", join(',', @$scutsites_salignstr), "\n";
	    }
		my $diff         = diff_of_array($qcutsites_qalignstr, $scutsites_salignstr);
		my $intersection = intersection_of_array($qcutsites_qalignstr, $scutsites_salignstr);
		if(@$diff > 0) {
			my $recognitionstr = $reobj->seq->seq;
			my $cutafter       = $reobj->cuts_after;
			print "      <rflp enzymename=\"$rname\" recog=\"$recognitionstr\" cutafter=\"$cutafter\">\n";
			for(@$diff){
				print "        <site pos=\"$_\" />\n";
			}
			for(@$intersection){
				print "        <commonsite pos=\"$_\" />\n";
			}
			print "      </rflp>\n";
		}
	}
	print "    </rflps>\n";
	print "  </matchsegment>\n";
}
print "</rflpmap>\n";

=pod

=head1 NAME

blastcsv2rflpmap.pl - Making RFLP map from BLAST CSV

=head1 SYNOPSIS

blastcsv2rflpmap.pl [options...] [input BLAST CSV file(s)]

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

B<blastcsv2rflpmap.pl> will read the given input BLAST CSV file and
makes RFLP map in XML format.

=cut
