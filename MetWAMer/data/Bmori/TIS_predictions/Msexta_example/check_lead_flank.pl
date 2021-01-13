#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 13 January 2021

use constant UPSTREXTENT => 5;
use constant CONTENTSWATHLEN => 96;

my $fraglen=UPSTREXTENT + CONTENTSWATHLEN;

while(<>) {
  my $desc=$_;
  my $seq=<>;
  $seq=~/^(\w{$fraglen})/;
  my $frag=$1;
  my $ambig=0;
  ++$ambig while $frag =~ m/[^acgt]/ig;
  if($ambig < int($fraglen / 10)) {
    print $desc,$seq;
  }
}

exit 0;
