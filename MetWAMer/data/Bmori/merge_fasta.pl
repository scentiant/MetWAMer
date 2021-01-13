#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 22 December 2020

# merge_fasta.pl - a utility that reads in a Fasta-formatted
# file of nucleotide sequences, and reports the sequence data
# on one line.

my $seq="";
my $ct=0;

while(my $line=<>) {
  if($line=~m/^>/) { # header
    if($ct) {
      #$seq=toupper($seq);
      print $seq,"\n";
      $seq="";
    }
    ++$ct;

    print $line;
    next;
  }

  chomp($line);
  $line=~s/\s//g;
  $seq .= $line;
}

print $seq,"\n";

exit;
