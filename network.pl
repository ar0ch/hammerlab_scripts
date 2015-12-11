#!/usr/bin/perl -w
# Aroon Chande
# Prepare files fo rnetwork analysis
use strict; use warnings;
use Algorithm::Combinatorics qw(combinations);
max
my $seq = [`seq -s , $max`];

my $pairs = combinations($seq, 2);

while (my $c = $pairs->next) {
    print "@$c\n";
}
