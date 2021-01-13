#!/usr/bin/perl
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 13 January 2021

my $PRIORITY=10;

my $BADINTRON="intron-length-out-of-bounds\n";

# parameters used in Augustus' scripts
my $MIN_INTRON_LEN=41;
my $MAX_INTRON_LEN=350000;
my $ADJ=15;

my $cds_cnt=0; # serial no.
while(<>) {
  chomp($_);
  my $overallRecord=$_;
  #>NC_051115.1+ (1335602  1336002,1339520  1339745,1341010  1342137)
  #>NC_051115.1- (1403512  1403301,1401908  1401797,1401423  1400575)
  $overallRecord=~/^>(\S+)([+-])\s+\((\d.*?\d)\)$/;
  my($ctg_id,$orient,$coords)=($1,$2,$3);
  ++$cds_cnt;

  # if present, record the intron line items first
  my @introns=();
  if($coords=~/,/) {
    my $intron_coords=$coords;
    while($intron_coords=~/(\d+),(\d+)/g) {
      my($istart,$istop)=($1,$2);
      if($orient eq '+') { # forward strand
        ++$istart; --$istop;
        if( (($istop - $istart + 1) >= $MIN_INTRON_LEN) &&
            (($istop - $istart + 1) <= $MAX_INTRON_LEN) ) {
          my $record="$ctg_id\tgth2h\tintron\t$istart\t$istop\t.\t+\t.\tsrc=P;grp=gene.${cds_cnt};pri=$PRIORITY\n";
          push(@introns,$record);
        }
        else { # we want to preserve exon-intron interleaving order during printout
          push(@introns,$BADINTRON);
        }
      }
      else { # reverse strand
      --$istart; ++$istop;
        if( (($istart - $istop + 1) >= $MIN_INTRON_LEN) &&
            (($istart - $istop + 1) <= $MAX_INTRON_LEN) ) {
          my $record="$ctg_id\tgth2h\tintron\t$istop\t$istart\t.\t-\t.\tsrc=P;grp=gene.${cds_cnt};pri=$PRIORITY\n";
          push(@introns,$record);
        }
        else {
          push(@introns,$BADINTRON);
        }
      }
    }
  }

  # >= 1 exons per CDS
  my @exons=();
  while($coords=~/(\d+)  (\d+)/g) {
    my($estart,$estop)=($1,$2);
    if($orient eq '+') { # forward strand
      $estart += $ADJ; $estop -= $ADJ; # AUGUSTUS wants incorrect CDS ends
      if($estart > $estop) { # set to midpoint
        $estart = $estop = int(($estart + $estop)/2);
      }
      my $record="$ctg_id\tgth2h\tCDSpart\t$estart\t$estop\t.\t+\t.\tsrc=P;grp=gene.${cds_cnt};pri=$PRIORITY\n";
      push(@exons,$record);
    }
    else { # reverse strand
      $estart -= $ADJ; $estop += $ADJ;
      if($estart < $estop) { # set to midpoint
        $estart = $estop = int(($estart + $estop)/2);
      }
      my $record="$ctg_id\tgth2h\tCDSpart\t$estop\t$estart\t.\t-\t.\tsrc=P;grp=gene.${cds_cnt};pri=$PRIORITY\n";
      push(@exons,$record);
    }
  }

  if( ($#exons > 0) && ($#exons != $#introns + 1) ) { # sanity check
    print STDERR "exon-intron cnts off in $overallRecord ??\n";
  }
  # print to stdout (delete intron error lines w/ `grep -v intron-length-out-of-bounds`)
  for(my $i=0;$i<=$#exons;++$i) {
    print $exons[$i];
    if(($#exons > 0) && ($i < $#exons)) {
      print $introns[$i];
    }
  }
}

exit 0;
