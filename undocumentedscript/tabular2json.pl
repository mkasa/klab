#!/usr/bin/env perl
# NOTE: '#!/usr/bin/perl' picks the deprecated system perl on current macOS
#       and the wrong interpreter on Homebrew setups; use env.

use strict;
use warnings;
use JSON;

# Whether the current input file is a real TSV. Decided from its header line:
# a TSV must be split on /\t/ with a -1 limit, otherwise empty cells collapse
# and every later value ends up under the wrong key.
my $is_tsv = 0;

sub split_row
{
    my $line = shift;
    return () unless(defined $line);
    return split(/\t/, $line, -1) if($is_tsv);
    $line =~ s/^\s+//;
    $line =~ s/\s+\z//;
    return () if($line eq '');
    return split(/\s+/, $line);
}

my @header_columns;
my %warned;
my $header_read = 0;

while(<>) {
    chomp;
    s/\r\z//;
    # Input bytes are UTF-8; decode them so that encode_json does not
    # re-encode them as if they were latin-1 (which double-encoded them).
    utf8::decode($_);
    # $. is reset by 'close ARGV if eof' below, so $. == 1 marks the header
    # line of EVERY input file, not just the first one.
    if($. == 1) {
        $is_tsv = (index($_, "\t") >= 0) ? 1 : 0;
        my @columns = split_row($_);
        if($header_read) {
            # Second and later files: their header must be skipped, not
            # emitted as a data row.
            unless(join("\t", @columns) eq join("\t", @header_columns)) {
                print STDERR "WARNING: '$ARGV' has a different header; its columns are ignored\n";
            }
        } else {
            @header_columns = @columns;
            my %seen;
            for my $h (@header_columns) {
                if($seen{$h}++) {
                    print STDERR "WARNING: duplicated column name '$h': only the last one is kept in the JSON output\n";
                }
            }
            $header_read = 1;
        }
        close ARGV if(eof);
        next;
    }
    my @columns = split_row($_);
    my $json = {};
    for(my $i = 0; $i < @columns; $i++) {
        unless($i < @header_columns) {
            # More fields than the header: the key would be undef, which
            # stringifies to "" and silently collapses all surplus columns.
            unless($warned{surplus}++) {
                print STDERR "WARNING: line $. has more fields than the header; the surplus fields are ignored\n";
            }
            last;
        }
        my $value = $columns[$i];
        # Emit real JSON numbers for plain numeric values, so that downstream
        # comparisons are numeric. Values with a leading zero (e.g. '007' or a
        # zip code) are deliberately left as strings.
        $value += 0 if(defined $value && $value =~ /^-?(0|[1-9][0-9]*)(\.[0-9]+)?([eE][-+]?[0-9]+)?$/);
        $json->{$header_columns[$i]} = $value;
    }
    print encode_json($json), "\n";
    close ARGV if(eof);
}
unless($header_read) {
    print STDERR "WARNING: no input (no header line was read)\n";
}
