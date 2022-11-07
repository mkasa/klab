#!/usr/bin/perl

use strict;
use JSON;

my $header_line = <>; chomp $header_line;
my @header_columns = split(/\s+/, $header_line);

print <<"EOF";
#{"name": "stats"}
#{"name": "file", "width": 20, "height": 1, "expr": "file", "keycol": 0, "fmtstr": "", "voffset": 0, "hoffset": 0, "aggstr": "", "type": "", "col": "Column"}
#{"name": "sum_len", "width": 12, "height": 1, "expr": "sum_len", "keycol": 0, "fmtstr": "{:,}", "voffset": 0, "hoffset": 0, "aggstr": "", "type": "int", "col": "Column"}
#{"name": "num_seqs", "width": 10, "height": 1, "expr": "num_seqs", "keycol": 0, "fmtstr": "{:,}", "voffset": 0, "hoffset": 0, "aggstr": "", "type": "int", "col": "Column"}
#{"name": "min_len", "width": 9, "height": 1, "expr": "min_len", "keycol": 0, "fmtstr": "{:,}", "voffset": 0, "hoffset": 0, "aggstr": "", "type": "int", "col": "Column"}
#{"name": "avg_len", "width": 9, "height": 1, "expr": "avg_len", "keycol": 0, "fmtstr": "{%.01f}", "voffset": 0, "hoffset": 0, "aggstr": "", "type": "float", "col": "Column"}
#{"name": "max_len", "width": 9, "height": 1, "expr": "max_len", "keycol": 0, "fmtstr": "{:,}", "voffset": 0, "hoffset": 0, "aggstr": "", "type": "int", "col": "Column"}
EOF
while(<>) {
    chomp;
    my @columns = split(/\s+/);
    my $json = {};
    for(my $i = 0; $i < @columns; $i++) {
        $columns[$i] =~ s/,//g if(3 <= $i && $i <= 7);
        $columns[$i] .= ".0" if($i eq 6 && $columns[$i] !~ /\./);
        $json->{$header_columns[$i]} = $columns[$i];
    }
    print encode_json($json), "\n";
}

