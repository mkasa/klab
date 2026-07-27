# BLASTCSVParse.pm - parse one line of a BLAST CSV file.
#
#   use BLASTCSVParse qw(parseBLASTCSV);
#
#   my $m = parseBLASTCSV($line);
#   print $m->{queryid}, "\t", $m->{subjectid}, "\n";
#
# Nothing is exported by default; import exactly the subs you use, or call
# them fully qualified as BLASTCSVParse::parseBLASTCSV(...).
# This file is a module, so it carries no shebang.

package BLASTCSVParse;

use strict;
use warnings;

use Exporter 'import';

our @EXPORT_OK = qw(
	parseBLASTCSV
	splitCSVLine
);
our %EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

# Split one CSV record into fields, honouring RFC4180 style double quoting
# ("a,b" is one field, "" inside a quoted field is a literal quote).
# Text::CSV is not required (it is frequently absent on the machines these
# scripts run on), so the quoting is handled here.
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

1;
