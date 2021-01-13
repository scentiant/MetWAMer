#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 20 December 2020

while(<>) {
  my $desc=$_;
  my $seq=<>;
  chomp($seq);
  if($seq!~/[^aAcCgGtT]/) {
    print "$desc$seq\n";
  }
  else {
    print STDERR "purged $desc";
  }
}

exit 0;
