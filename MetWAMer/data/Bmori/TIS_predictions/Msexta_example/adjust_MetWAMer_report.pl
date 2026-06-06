#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (mespar1@gmail.com), 13 January 2021

while(<>) {
  my $desc=$_;
  $_=<>;
  chomp($_);
  print $_,"\t",$desc;
}

exit 0;
