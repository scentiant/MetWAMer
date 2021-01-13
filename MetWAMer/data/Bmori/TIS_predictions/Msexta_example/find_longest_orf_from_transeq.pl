#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 12 January 2021

my $PEPFILEIN=shift;
my $CDSFILEIN=shift;
my $PEPFILEOUT=shift;
my $CDSFILEOUT=shift or die "$0 pepfile.fpa cdsfile.fna outpep.fpa outcds.fna\n";

open(PEP,"<$PEPFILEIN") or die "$!\n";
open(CDS,"<$CDSFILEIN") or die "$!\n";
open(OUTPEP,">$PEPFILEOUT") or die "$!\n";
open(OUTCDS,">$CDSFILEOUT") or die "$!\n";

while(<PEP>) {
  my $descPEP=$_;
  my $seqPEP=<PEP>;
  chomp($seqPEP);

  my $descCDS=<CDS>;
  my $seqCDS=<CDS>;

  my $lenlongest=-1; #stores length of longest protein fragment encountered
  my $longest; #stores the longest protein fragment encountered

  my @tokens = split(/\*/,$seqPEP);
  foreach my $token (@tokens) {
    my $toklen=length($token);
    if($toklen > $lenlongest) {
      $lenlongest = $toklen;
      $longest = $token;
    }
  }

  if($longest eq $seqPEP) { # no adjustments needed!
    chomp($descCDS);
    $descCDS=~s/\s+$//;
    print OUTPEP "$descCDS\n$seqPEP\n";
    print OUTCDS "$descCDS\n$seqCDS";
    next;
  }

  my $offsetNterm = 3 * index($seqPEP, $longest);

  #>NW_023592558.1+ (7909  8280,8518  8571,9031  9087)
  $descCDS=~/^(>.*?\.\d+([+-])\s+\()(\d.*?\d)(\))/;
  my ($bite1,$orient,$coords,$bite2)=($1,$2,$3,$4); # $bite2 keeps the regex honest

  my @exons=split(/,/,$coords);

  if($orient eq '+') {
    #>NW_023592558.1+ (7909  8280,8518  8571,9031  9087)

    #adjust for clipped stuff at N terminus
    my $deplete=$offsetNterm;
    my $i;
    for($i=0;$i<=$#exons;++$i) {
      $exons[$i]=~/^(\d+)\s+(\d+)$/;
      my($start,$stop)=($1,$2);

      if(($stop-$start+1) > $deplete) {
        $start += $deplete;
        $exons[$i]="$start  $stop";
        last;
      }
      else {
        $deplete -= $stop-$start+1;
      }
    }
    # from here, $i stores the left-most exon

    #adjust for clipped stuff at C terminus
    $deplete=3 * $lenlongest;
    my $j;
    for($j=$i;$j<=$#exons;++$j) {
      $exons[$j]=~/^(\d+)\s+(\d+)$/;
      my($start,$stop)=($1,$2);
  
      if(($stop-$start+1) < $deplete) {
        $deplete -= $stop-$start+1;
      }
      else {
        $stop=$start+$deplete-1;
        $exons[$j]="$start  $stop";
        last;
      }
    }
    # from here, $j stores the right-most exon

    my $desc=$bite1;
    for(my $k=$i;$k<$j;++$k) {
       $desc .= $exons[$k];
       $desc .= ",";
    }
    $desc .= $exons[$j];
    $desc .= $bite2;

    print OUTPEP "$desc\n$longest\n";
    print OUTCDS "$desc\n",substr($seqCDS,$offsetNterm,3*$lenlongest),"\n";
  }
  elsif($orient eq '-') {
    #>NC_051141.1- (4900556  4899767,4880231  4879993)

    #adjust for clipped stuff at N terminus
    my $deplete=$offsetNterm;
    my $i;
    for($i=0;$i<=$#exons;++$i) {
      $exons[$i]=~/^(\d+)\s+(\d+)$/;
      my($start,$stop)=($1,$2);

      if(($start-$stop+1) > $deplete) {
        $start -= $deplete;
        $exons[$i]="$start  $stop";
        last;
      }
      else {
        $deplete -= $start-$stop+1;
      }
    }
    # from here, $i stores the left-most exon

    #adjust for clipped stuff at C terminus
    $deplete=3 * $lenlongest;
    my $j;
    for($j=$i;$j<=$#exons;++$j) {
      $exons[$j]=~/^(\d+)\s+(\d+)$/;
      my($start,$stop)=($1,$2);
  
      if(($start-$stop+1) < $deplete) {
        $deplete -= $start-$stop+1;
      }
      else {
        $stop=$start-$deplete+1;
        $exons[$j]="$start  $stop";
        last;
      }
    }
    # from here, $j stores the right-most exon

    my $desc=$bite1;
    for(my $k=$i;$k<$j;++$k) {
       $desc .= $exons[$k];
       $desc .= ",";
    }
    $desc .= $exons[$j];
    $desc .= $bite2;

    print OUTPEP "$desc\n$longest\n";
    print OUTCDS "$desc\n",substr($seqCDS,$offsetNterm,3*$lenlongest),"\n";
  }
  else { die "not +-"; }
}

close PEP;
close CDS;
close OUTPEP;
close OUTCDS;

exit 0;
