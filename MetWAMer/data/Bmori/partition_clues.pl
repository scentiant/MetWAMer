#!/usr/bin/perl -w
use strict;
use List::Util qw( reduce );

# Michael E. Sparks (michael.sparks2@usda.gov), 20 December 2020

my $MEDFILE=shift @ARGV;
my $CLUSTFILE=shift @ARGV or die "$0 medoids.xml cluster_report.dat\n";

my %medoids=();
open(MEDOIDS,"<$MEDFILE") or die "$!\n";
while(<MEDOIDS>) {
  if(/<medoid id=/) {
    $_=~/id=\"(\d+)\">(\w+)</;
    my($id,$med)=($1,$2);
    $medoids{$id}=$med;
  }
}
close MEDOIDS;

open(CLUSTERS,"<$CLUSTFILE") or die "$!\n";
while(<CLUSTERS>) {
  my $line=$_;
  chomp($line);
  $line=~/\s(\w+)$/;
  my $site=$1;

  # embellish clustering report
  my %dists=();
  print "$line";
  foreach my $key (sort keys %medoids) {
    # count null chars after XOR'ing together
    $dists{$key}=($medoids{$key}^$site)=~tr/\0//c;
    print "\t",$dists{$key};
  }
  my $keyOfMinValue = reduce { $dists{$a} <= $dists{$b} ? $a : $b } keys %dists;
  print "\t$keyOfMinValue\n";
}
close CLUSTERS;

exit 0;
