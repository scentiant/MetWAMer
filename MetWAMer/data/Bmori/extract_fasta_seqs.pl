#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 20 December 2020

my $FILE = shift or die; # sequences should be merged onto one line in advance

my %seqs=();
open(FNA,"<$FILE") or die "$!\n";
while(<FNA>) {
  chomp($_);
  $_=~/^>(\S+)/; # use first token as key
  my $key=$1;

  $_=<FNA>;
  chomp($_);
  my $seq=$_;
  if(! exists($seqs{$key}) ) {
    $seqs{$key}=$seq;
  }
  else {
    die "$key seen >= 2 times?\n";
  }
}
close FNA;

while(<>) {
  chomp($_);
  if(exists($seqs{$_})) {
    print ">$_\n",$seqs{$_},"\n";
  }
  else {
    print STDERR "Not seen: $_\n";
  }
}

exit 0;
