# BLASTCSVParse.pm - parse one line of a BLAST CSV file.
#
#   my $m = parseBLASTCSV($line);
#   print $m->{queryid}, "\t", $m->{subjectid}, "\n";
#
# NOTE: this file intentionally does NOT declare a 'package'; parseBLASTCSV()
#       is defined directly in the caller's namespace (normally main::).
#       Adding a package would break existing callers, but beware of name
#       collisions. This file is also a module, so it carries no shebang.

use strict;
use warnings;

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
