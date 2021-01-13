#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 13 January 2021

my $FILE = shift or die; # sequences should be merged onto one line in advance

my %seqs=();
open(FNA,"<$FILE") or die "$!\n";
while(<FNA>) {
  chomp($_);
  $_=~/^>(\S+.*?)$/;
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
  my $key=$_;
  if(! exists($seqs{$key}) ) { die "where is $key?\n"; }
  print ">$key\n",$seqs{$key},"\n";
}

exit 0;
