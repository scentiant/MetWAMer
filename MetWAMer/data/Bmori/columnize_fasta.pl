#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 22 December 2020

# for pretty printing fasta sequences

while(my $line=<>) {
  if($line=~m/^>/) {
    print $line;
    next;
  }
  chomp($line);
  for(my $i=0;$i<length($line);$i += 60) {
    print substr($line,$i,60),"\n";
  }
}

exit 0;
