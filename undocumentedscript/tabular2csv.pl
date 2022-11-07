#!/usr/bin/perl

use strict;

my $header_line = <>; chomp $header_line;
my @header_columns = split(/\s+/, $header_line);

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
        $c =~ s/"/\"/g;
        print "\"$c\"";
    }
    print "\n";
}

print_csv(@header_columns);
while(<>) {
    chomp;
    my @columns = split(/\s+/);
    print_csv(@columns);
}

