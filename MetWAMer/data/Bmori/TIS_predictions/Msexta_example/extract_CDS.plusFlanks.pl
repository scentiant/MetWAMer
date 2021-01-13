#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 12 January 2021

use constant UPSTREXTENT => 5;
use constant DOWNSTREXTENT => 3;
use constant CONTENTSWATHLEN => 96;

while(<>) {
  my $desc=$_;
  $desc=~/^>(.*?\.\d+)([+-]).*?\((\d.*?\d)\)/;
  my ($gdna,$orient,$coords)=($1,$2,$3);

  my @exons=split(/,/,$coords);

  #>NC_051126.1+ (14589914  14589969,14590237  14590504)
  #>NC_051140.1- (2191813  2191745,2190688  2190576,2188570  2188441)

  $exons[0]=~/^(\d+)\s+(\d+)$/;
  my ($start,$stop)=($1,$2);
  if   ($orient eq '+') { $start-=(CONTENTSWATHLEN + UPSTREXTENT); }
  elsif($orient eq '-') { $start+=(CONTENTSWATHLEN + UPSTREXTENT); }
  else { die "not +-"; }
  if ($start < 0) {
    print STDERR $desc;
    next;
  }
  $exons[0]="$start  $stop";

  $exons[$#exons]=~/^(\d+)\s+(\d+)$/;
  ($start,$stop)=($1,$2);
  if    ($orient eq '+') { $stop+=(DOWNSTREXTENT + CONTENTSWATHLEN); }
  elsif ($orient eq '-') { $stop-=(DOWNSTREXTENT + CONTENTSWATHLEN); }
  else { die "not +-"; }
  if ($stop < 0) {
    print STDERR $desc;
    next;
  }
  $exons[$#exons]="$start  $stop";

  my $error=0;
  my $cds="";
  for (@exons) {
    /^(\d+)\s+(\d+)$/;
    my($estart,$estop)=($1,$2);
    my $seq=`indexFasSeq ./gDNA_templates/${gdna}.fas.ind $estart $estop`;
    if( ($estart == $estop) && ($orient eq '-') ) {
      chomp($seq);
      if   ($seq eq 'A') { $seq="T\n"; }
      elsif($seq eq 'C') { $seq="G\n"; }
      elsif($seq eq 'G') { $seq="C\n"; }
      elsif($seq eq 'T') { $seq="A\n"; }
      else               { $seq.="\n"; } # ambiguity's handled elsewhere
    }
    if($seq ne "") { $cds .= $seq; }
    else { $error=1; last; }
  }
  if($error == 0) { print $desc, $cds; }
  else { print STDERR $desc; }
}

exit 0;
