#!/usr/bin/env perl

use strict;

while(<>){
    chomp;
    s/\bstlp_std\:\://g;
    s/\:\:iterator/\:\:iter/g;
    s/\:\:const_iterator/\:\:c_iter/g;
    if(/^g(\+\+|cc)/) {
        print color(1), "$_", color(7), "\n";
    } elsif(/^[^:]+\:(\d+)\: warning\:/) {
        print color(3), "$_", color(7), "\n";
    } elsif(/^[^:]+\:(\d+)\: error\:/) {
        print color(2), "$_", color(7), "\n";
    } else {
        print "$_\n";
    }
}

sub color($) {
    my $cnum = shift;
    my $rv = "\x1b[";
    if($cnum == 0 || $cnum == 8) {
        $rv .= "30";
    } elsif($cnum == 1 || $cnum == 9) {
        $rv .= "34";
    } elsif($cnum == 2 || $cnum == 10) {
        $rv .= "32";
    } elsif($cnum == 3 || $cnum == 11) {
        $rv .= "36";
    } elsif($cnum == 4 || $cnum == 12) {
        $rv .= "31";
    } elsif($cnum == 5 || $cnum == 13) {
        $rv .= "35";
    } elsif($cnum == 6 || $cnum == 14) {
        $rv .= "33";
    } elsif($cnum == 7 || $cnum == 15) {
        $rv .= "37";
    }
    $rv .= "m";
    return $rv;
} 
