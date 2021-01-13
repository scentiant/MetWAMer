#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 20 December 2020

my %tismap=();
open(MAP,"<TIS_to_medoid_mapping.txt") or die "$!\n";
while(<MAP>) {
  chomp($_);
  $_=~/^(\w+)\t(\d)$/;
  $tismap{$1}=$2;
}
close MAP;

# assume sequences merged onto a single line
open(TIS,"<TIS_for_PWM.fna") or die "$!\n";
open(TIS1,">TIS_for_PWM.1.fna") or die "$!\n";
open(TIS2,">TIS_for_PWM.2.fna") or die "$!\n";
open(TIS3,">TIS_for_PWM.3.fna") or die "$!\n";

while(<TIS>) {
  my $desc=$_;
  $_=<TIS>;
  my $seq=$_;
  chomp($seq);
  $seq=~/^(\w{11})$/;
  my $bin=$tismap{$1};
  if($bin == 1)    { print TIS1 "${desc}$seq\n"; }
  elsif($bin == 2) { print TIS2 "${desc}$seq\n"; }
  elsif($bin == 3) { print TIS3 "${desc}$seq\n"; }
  else { die "?\n"; }
}

close TIS;
close TIS1;
close TIS2;
close TIS3;

exit 0;
