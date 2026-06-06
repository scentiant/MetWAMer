#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (mespar1@gmail.com), 20 December 2020

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
