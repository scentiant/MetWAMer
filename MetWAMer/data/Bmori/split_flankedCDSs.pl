#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 22 December 2020

my %tismap=();
open(MAP,"<TIS_to_medoid_mapping.txt") or die "$!\n";
while(<MAP>) {
  chomp($_);
  $_=~/^(\w+)\t(\d)$/;
  $tismap{$1}=$2;
}
close MAP;

# assume sequences merged onto a single line in advance
open(CDS,"<Bmori.codseq.withFlanks.fna") or die "$!\n";
open(CDS1,">Bmori.codseq.withFlanks.1.fna") or die "$!\n";
open(CDS2,">Bmori.codseq.withFlanks.2.fna") or die "$!\n";
open(CDS3,">Bmori.codseq.withFlanks.3.fna") or die "$!\n";

while(<CDS>) {
  my $desc=$_;
  $_=<CDS>;
  my $seq=$_;
  chomp($seq);
  $seq=~/^\w{96}(\w{11})\w/;
  my $bin=$tismap{$1};
  if($bin == 1)    { print CDS1 "${desc}$seq\n"; }
  elsif($bin == 2) { print CDS2 "${desc}$seq\n"; }
  elsif($bin == 3) { print CDS3 "${desc}$seq\n"; }
  else { die "?\n"; }
}

close CDS;
close CDS1;
close CDS2;
close CDS3;

exit 0;
