#!/usr/bin/env perl
# NOTE: '#!/usr/bin/perl' picks the deprecated system perl on current macOS
#       and the wrong interpreter on Homebrew setups; use env.

use strict;
use warnings;
use JSON::PP;

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
# the right and treat whatever is left as the file name. This also guarantees
# that a row never has more fields than the header, which used to make undef
# a hash key and collapse every surplus field into a single "" key.
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

my $header_line = <>;
unless(defined $header_line) {
    print STDERR "WARNING: no input (no header line was read)\n";
    exit 0;
}
chomp $header_line;
$header_line =~ s/\r\z//;
my @header_columns = split_header($header_line);

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
    s/\r\z//;
    next if(/^\s*$/);
    # Input bytes are UTF-8; decode them so that encode_json does not
    # re-encode them as if they were latin-1 (which produced mojibake).
    utf8::decode($_);
    my @columns = split_row($_, scalar(@header_columns));
    unless(@columns) {
        print STDERR "WARNING: skipping malformed line $.\n";
        next;
    }
    my $json = {};
    for(my $i = 0; $i < @columns; $i++) {
        my $value = $columns[$i];
        if($i > 0) {
            # Strip seqkit's thousands separators, then turn the result into
            # a real number so that encode_json emits 2000 rather than "2000"
            # (the VDS header above declares these columns as int/float).
            # Applies to every numeric column, so 'seqkit stat -a' works too.
            $value =~ s/,//g if($value =~ /^[\d,]+(\.\d+)?$/);
            $value += 0 if($value =~ /^-?(0|[1-9][0-9]*)(\.[0-9]+)?([eE][-+]?[0-9]+)?$/);
        }
        $json->{$header_columns[$i]} = $value;
    }
    print encode_json($json), "\n";
}
