#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 12 January 2021

my $MAXMASKRAT = shift or die;
my $desc="";

while(<>) {
  if($_=~ /^>/) { $desc=$_; next; }
  chomp($_);
  my $len=length($_);
  $_=uc($_);
  my @seq=split(//,$_);
  my $x=0;
  for (@seq) {
    if($_ eq 'X') { ++$x; }
  }

  if( ($x / $len) > $MAXMASKRAT ) {
    print STDERR "\t",$x / $len,"\n$desc$_\n";
  }
  else {
    print $desc,$_,"\n";
  }
}

exit 0;
