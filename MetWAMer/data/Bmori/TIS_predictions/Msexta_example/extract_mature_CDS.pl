#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 13 January 2021

while(<>) {
  my $desc=$_;
  $desc=~/^>(\S+)([+-])\s+\((\d.*?\d)\)/;
  my($ctg,$orient,$coords)=($1,$2,$3);
  print $_;

  my @exons=split(/,/,$coords);
  for (@exons) {
    /^(\d+)\s+(\d+)$/;
    my($estart,$estop)=($1,$2);
    my $seq=`indexFasSeq ./gDNA_templates/${ctg}.fas.ind $estart $estop`;
    if( ($estart == $estop) && ($orient eq '-') ) {
      chomp($seq);
      if   ($seq eq 'A') { $seq="T\n"; }
      elsif($seq eq 'C') { $seq="G\n"; }
      elsif($seq eq 'G') { $seq="C\n"; }
      elsif($seq eq 'T') { $seq="A\n"; }
      else               { $seq.="\n"; } # ambiguity's handled elsewhere
    }
    print $seq;
  }
}

exit 0;
