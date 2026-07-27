#!/usr/bin/env perl

use strict;
use warnings;

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

# Split one CSV record into fields, honouring RFC4180 style double quoting.
# NOTE: this is the same implementation as BLASTCSVParse.pm's splitCSVLine();
#       the two copies must be kept in sync.
sub splitCSVLine($)
{
	my $line = shift;
	return () unless(defined $line);
	$line =~ s/[\r\n]+\z//;
	my @fields;
	while($line =~ /\G(?:"((?:[^"]|"")*)"|([^,]*))(,|\z)/gc) {
		my $field;
		if(defined $1) {
			$field = $1;
			$field =~ s/""/"/g;
		} else {
			$field = $2;
		}
		push(@fields, $field);
		last if($3 eq '');
	}
	return @fields;
}

sub parseBLASTCSV($)
{
	my $line = shift;
	my ($qname, $qlen, $qst, $qed, $sname, $slen, $sst, $sed, $nummatch, $alignlen, $qastr, $sastr) = splitCSVLine($line);
	return {
		queryid      => $qname,
		querylen     => $qlen,
		querystart   => $qst,
		queryend     => $qed,
		subjectid    => $sname,
		subjectlen   => $slen,
		subjectstart => $sst,
		subjectend   => $sed,
		matchlen     => $nummatch,
		alignlen     => $alignlen,
		queryalignmentstring   => $qastr,
		subjectalignmentstring => $sastr
	};
}

# Escape a string so that it may be used as XML character data or inside a
# double quoted XML attribute value.
sub xml_escape($)
{
	my $s = shift;
	return '' unless(defined $s);
	$s =~ s/&/&amp;/g;
	$s =~ s/</&lt;/g;
	$s =~ s/>/&gt;/g;
	$s =~ s/"/&quot;/g;
	return $s;
}

# Translate positions given as 0-origin offsets into the UNGAPPED sequence
# into 0-origin column indices of the GAPPED alignment string.
#
# The returned array always has exactly as many elements as the input list
# (the previous implementation dropped the last site and mapped repeated
# positions to different columns, which made the diff/intersection below
# compare mismatched lists).
sub adjust_for_gap($$)
{
	my $gappedStr = shift;
	my $positionListRef = shift;

	my @retval;
	return \@retval unless(defined $positionListRef && @$positionListRef > 0);

	my @splitGappedStr = split(//, defined $gappedStr ? $gappedStr : '');
	# ungapped offset -> gapped column
	my @columnOfUngappedPosition;
	for(my $col = 0; $col < @splitGappedStr; $col++) {
		my $c = $splitGappedStr[$col];
		next if($c eq '-' || $c eq ' ');
		push(@columnOfUngappedPosition, $col);
	}
	my $lastColumn = @splitGappedStr > 0 ? $#splitGappedStr : 0;
	for my $pos (@$positionListRef) {
		if(!defined $pos || $pos < 0) {
			push(@retval, 0);
		} elsif($pos >= @columnOfUngappedPosition) {
			push(@retval, $lastColumn);
		} else {
			push(@retval, $columnOfUngappedPosition[$pos]);
		}
	}
	return \@retval;
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
	my $retval = '';
	my @str1arr = split(//, defined $str1 ? $str1 : '');
	my @str2arr = split(//, defined $str2 ? $str2 : '');
	for(my $i = 0; $i < @str1arr; $i++) {
		# $str2arr may be shorter than $str1arr on a malformed record;
		# treat the missing columns as unmasked instead of reading undef.
		if($i < @str2arr && $str2arr[$i] =~ /n/i && $str1arr[$i] ne '-') {
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
	# TODO: �ǉ�.
	# MnlI
}

# The input bytes are copied through verbatim, so the declared encoding must
# be one that accepts arbitrary bytes above 0x7f; US-ASCII made the
# declaration a lie (and the document unparsable) for any non-ASCII input.
print '<?xml version="1.0" encoding="UTF-8" ?>', "\n";
print "<rflpmap>\n";
while(<>){
	chomp;
	next if(/^#/);
	next if(/^\s*$/);
	my $bmatch = parseBLASTCSV($_);
	unless(defined $bmatch->{queryalignmentstring} && defined $bmatch->{subjectalignmentstring}) {
		print STDERR "WARNING: skipping malformed input line $.: not enough CSV fields\n";
		next;
	}
	print "  <matchsegment queryid=\"", xml_escape($bmatch->{queryid}), "\" subjectid=\"", xml_escape($bmatch->{subjectid}), "\" ";
	print "querystart=\"", xml_escape($bmatch->{querystart}), "\" queryend=\"", xml_escape($bmatch->{queryend}), "\" ";
	print "subjectstart=\"", xml_escape($bmatch->{subjectstart}), "\" subjectend=\"", xml_escape($bmatch->{subjectend}), "\" >\n";
	print "    <queryseq>", xml_escape($bmatch->{queryalignmentstring}), "</queryseq>\n";
	print "    <subjectseq>", xml_escape($bmatch->{subjectalignmentstring}), "</subjectseq>\n";
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
			# Diagnostics must never go to STDOUT: STDOUT carries the XML
			# document and any extra text would make it ill-formed.
			print STDERR "$qstr\n";
			print STDERR "$sstr\n";
		}
		my $qseqobj = new Bio::Seq(-seq => $qstr);
		my $sseqobj = new Bio::Seq(-seq => $sstr);
		# Bio::Tools::RestrictionEnzyme::cut_locations() numbers locations
		# starting with 1, while adjust_for_gap() works in 0-origin offsets.
		my $qcutsites_qstr = [ map { $_ - 1 } @{$reobj->cut_locations($qseqobj)} ];
		my $scutsites_sstr = [ map { $_ - 1 } @{$reobj->cut_locations($sseqobj)} ];
		if($debug) {
	        print STDERR "For enzyme $rname\n";
	        print STDERR " Ungap (0-origin):\n";
	        print STDERR "  ", join(',', @$qcutsites_qstr), "\n";
	        print STDERR "  ", join(',', @$scutsites_sstr), "\n";
	    }
		my $qcutsites_qalignstr = adjust_for_gap($qalignstr, $qcutsites_qstr);
		my $scutsites_salignstr = adjust_for_gap($salignstr, $scutsites_sstr);
		if($debug) {
	        print STDERR " Gapped (0-origin column):\n";
	        print STDERR "  ", join(',', @$qcutsites_qalignstr), "\n";
	        print STDERR "  ", join(',', @$scutsites_salignstr), "\n";
	    }
		my $diff         = diff_of_array($qcutsites_qalignstr, $scutsites_salignstr);
		my $intersection = intersection_of_array($qcutsites_qalignstr, $scutsites_salignstr);
		if(@$diff > 0) {
			my $recognitionstr = $reobj->seq->seq;
			my $cutafter       = $reobj->cuts_after;
			print "      <rflp enzymename=\"", xml_escape($rname), "\" recog=\"", xml_escape($recognitionstr), "\" cutafter=\"", xml_escape($cutafter), "\">\n";
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

The B<pos> attribute of E<lt>siteE<gt>/E<lt>commonsiteE<gt> is a 0-origin
column index into the gapped alignment strings reported in
E<lt>queryseqE<gt>/E<lt>subjectseqE<gt>.

With B<-debug>, diagnostics are written to standard error so that standard
output remains a well-formed XML document.

=cut
