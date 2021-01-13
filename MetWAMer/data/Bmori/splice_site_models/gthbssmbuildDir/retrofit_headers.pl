#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 20 December 2020

my %md5=();
open(MAP,"<bmori.md5map.dat") or die "$!\n";
while(<MAP>) {
  my $kv=$_;
  chomp($kv);
  #md5:6d023d2783c5d6cec91716cdba6ba908:NW_004581679.1
  $kv=~m/^.*:(\S+)$/;
  $md5{$1}=$kv;
}
close MAP;

while(<>) {
  my $desc=$_;
  #>  1855828 1856051 1 NW_004582009+_PGL-83_AGS-1_PPS_1FALSE 1
  $desc=~s/^>  />/;
  $desc=~m/^(>\d+\s+\d+\s+[012]\s+)(\S+)([+-])_\S+/;
  #>349718 349860 0 md5:642e90b99b2ab029fa401f876d2b67f6:NW_004581786.1+
  my ($tk1,$tk2,$tk3)=($1,$2,$3);
  $tk2 .= ".1";
  print $tk1,$md5{$tk2},$tk3,"\n";

  $_=<>;
  print $_;
}

exit 0;
