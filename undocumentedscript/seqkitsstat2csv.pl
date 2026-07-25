#!/usr/bin/env perl
# NOTE: '#!/usr/bin/perl' picks the deprecated system perl on current macOS
#       and the wrong interpreter on Homebrew setups; use env.

use strict;
use warnings;

my @header_columns;

sub split_header
{
    my $line = shift;
    $line =~ s/^\s+//;
    $line =~ s/\s+\z//;
    return () if($line eq '');
    return split(/\s+/, $line);
}

# The FIRST column of 'seqkit stat' output is a file name, which may itself
# contain whitespace. The remaining columns never do, so peel them off from
# the right and treat whatever is left as the file name.
sub split_row
{
    my ($line, $ncolumns) = @_;
    $line =~ s/^\s+//;
    $line =~ s/\s+\z//;
    return () if($line eq '');
    return ($line) if($ncolumns <= 1);
    my @tail;
    for(my $k = 0; $k < $ncolumns - 1; $k++) {
        if($line =~ s/\s+(\S+)\z//) {
            unshift(@tail, $1);
        } else {
            return ();   # fewer fields than the header: malformed row
        }
    }
    return ($line, @tail);
}

sub print_csv
{
    my @args = @_;
    my $first_flag = 0;
    for my $c (@args) {
        if($first_flag) {
            print ",";
        } else {
            $first_flag = 1;
        }
        $c = '' unless(defined $c);
        $c =~ s/"/""/g;   # a quote is escaped by DOUBLING it in CSV
        print "\"$c\"";
    }
    print "\n";
}

my $header_printed = 0;
while(<>) {
    chomp;
    s/\r\z//;
    # $. is reset by 'close ARGV if eof' below, so $. == 1 marks the header
    # line of EVERY input file, not just the first one.
    if($. == 1) {
        @header_columns = split_header($_);
        print_csv(@header_columns) unless($header_printed);
        $header_printed = 1;
        close ARGV if(eof);
        next;
    }
    if(/^\s*$/) {
        close ARGV if(eof);
        next;
    }
    my @columns = split_row($_, scalar(@header_columns));
    unless(@columns) {
        print STDERR "WARNING: skipping malformed line $. ('$_')\n";
        close ARGV if(eof);
        next;
    }
    # Strip the thousands separators seqkit puts in its numbers. Applied to
    # every numeric-looking column but the file name, so that the extra
    # columns of 'seqkit stat -a' are handled too (the old code hardcoded
    # columns 3..7 and autovivified @columns up to 8 entries).
    for(my $i = 1; $i < @columns; $i++) {
        $columns[$i] =~ s/,//g if($columns[$i] =~ /^[\d,]+(\.\d+)?$/);
    }
    print_csv(@columns);
    close ARGV if(eof);
}
unless($header_printed) {
    print STDERR "WARNING: no input (no header line was read)\n";
}
