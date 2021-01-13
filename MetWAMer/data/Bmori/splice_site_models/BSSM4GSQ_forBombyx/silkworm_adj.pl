#!/usr/bin/perl -w

# Michael E. Sparks (michael.sparks2@usda.gov), 20 December 2020

while(<>) {
  if(/^>/) {
    chomp($_);
    /^(>.*?)\.1([-+]_PGL.*?)$/;
    print "$1$2\n";
  }
  else {
    print $_;
  }
}

exit 0;
