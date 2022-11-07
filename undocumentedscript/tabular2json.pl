#!/usr/bin/perl

use strict;
use JSON;

my $header_line = <>; chomp $header_line;
my @header_columns = split(/\s+/, $header_line);

while(<>) {
    chomp;
    my @columns = split(/\s+/);
    my $json = {};
    for(my $i = 0; $i < @columns; $i++) {
        $json->{$header_columns[$i]} = $columns[$i];
    }
    print encode_json($json), "\n";
}

