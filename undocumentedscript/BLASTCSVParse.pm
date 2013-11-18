#!/usr/bin/env perl

use strict;

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

1;
