#!/usr/bin/env perl
# NOTE: '#!/usr/bin/perl' picks the deprecated system perl on current macOS
#       and the wrong interpreter on Homebrew setups; use env.

use strict;
use warnings;

# Whether the current input file is a real TSV. Decided from its header line:
# a TSV must be split on /\t/ with a -1 limit, otherwise empty cells collapse
# and every later value shifts one column to the left.
my $is_tsv = 0;

sub split_row
{
    my $line = shift;
    return () unless(defined $line);
    return split(/\t/, $line, -1) if($is_tsv);
    # Whitespace separated (column aligned) table: leading whitespace must be
    # stripped, otherwise split(/\s+/) yields a phantom empty first field.
    $line =~ s/^\s+//;
    $line =~ s/\s+\z//;
    return () if($line eq '');
    return split(/\s+/, $line);
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
        $is_tsv = (index($_, "\t") >= 0) ? 1 : 0;
        my @header_columns = split_row($_);
        print_csv(@header_columns) unless($header_printed);
        $header_printed = 1;
        close ARGV if(eof);
        next;
    }
    my @columns = split_row($_);
    print_csv(@columns);
    close ARGV if(eof);
}
unless($header_printed) {
    print STDERR "WARNING: no input (no header line was read)\n";
}
