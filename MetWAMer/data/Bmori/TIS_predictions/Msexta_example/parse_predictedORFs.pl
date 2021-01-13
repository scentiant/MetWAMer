#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 12 January 2021

while(<>) {
  next if ($_ !~ /^>/);
  do {
    chomp($_);
    if($_!~/^>/){
      $_=~s/\d//g;
      $_=~s/\s//g;
    }
    print $_,"\n";
    $_=<>;
  } while($_=~/\S/);
}

exit(0);
