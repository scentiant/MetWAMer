#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 13 January 2021

while(<>) {
  if($_!~/^>/) {
    print $_;
    next;
  }

  chomp($_);
  $_=~/^(>\S+[+-])\s+(\(\d.*?\d\))\s+\(MetWAMer ATG offset:\s+(\d+)\s+bp\)$/;
  my($ctg,$coords,$offset)=($1,$2,$3);

  if($ctg=~/\+$/) {
    my $coords2="";
    while($coords=~/(\d+  \d+)/g) {
      my $exon=$1;
      $exon=~/^(\d+)\s+(\d+)$/;
      my($start,$stop)=($1,$2);
      if($offset < ($stop - $start + 1)) {
        $start += $offset;
        $offset=0;
        $coords2 .= "$start  $stop,";
      }
      else {
        $offset -= $stop - $start + 1;
      }
    }
    chop($coords2);
    print "$ctg ($coords2)\n";
  }
  else { # reverse strand
    my $coords2="";
    while($coords=~/(\d+  \d+)/g) {
      my $exon=$1;
      $exon=~/^(\d+)\s+(\d+)$/;
      my($start,$stop)=($1,$2);
      if($offset < ($start - $stop + 1)) {
        $start -= $offset;
        $offset=0;
        $coords2 .= "$start  $stop,";
      }
      else {
        $offset -= $start - $stop + 1;
      }
    }
    chop($coords2);
    print "$ctg ($coords2)\n";
  }
}

exit 0;
