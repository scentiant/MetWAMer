#!/usr/bin/perl
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 13 January 2021

my $PRIORITY=10;

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

        my $record="$ctg_id\tanchor\tintronpart\t$istart\t$istop\t.\t+\t.\tsrc=M;grp=gene.${cds_cnt};pri=$PRIORITY\n";

        my($donBeg,$donEnd)=($istart,$istart+1);
        my $don=`indexFasSeq ./gDNA_templates/${ctg_id}.fas.ind $donBeg $donEnd`;
        chomp($don);
        if($don eq 'GT') {
          $record.="$ctg_id\tanchor\tdss\t$donBeg\t$donBeg\t.\t+\t.\tsrc=M;grp=gene.${cds_cnt};pri=$PRIORITY\n";
        }

        my($accBeg,$accEnd)=($istop-1,$istop);
        my $acc=`indexFasSeq ./gDNA_templates/${ctg_id}.fas.ind $accBeg $accEnd`;
        chomp($acc);
        if($acc eq 'AG') {
          $record.="$ctg_id\tanchor\tass\t$accEnd\t$accEnd\t.\t+\t.\tsrc=M;grp=gene.${cds_cnt};pri=$PRIORITY\n";
        }

        push(@introns,$record);
      }
      else { # reverse strand
        --$istart; ++$istop;

        my $record="$ctg_id\tanchor\tintronpart\t$istop\t$istart\t.\t-\t.\tsrc=M;grp=gene.${cds_cnt};pri=$PRIORITY\n";

        my($donBeg,$donEnd)=($istart,$istart-1);
        my $don=`indexFasSeq ./gDNA_templates/${ctg_id}.fas.ind $donBeg $donEnd`;
        chomp($don);
        if($don eq 'GT') {
          $record.="$ctg_id\tanchor\tdss\t$donBeg\t$donBeg\t.\t-\t.\tsrc=M;grp=gene.${cds_cnt};pri=$PRIORITY\n";
        }

        my($accBeg,$accEnd)=($istop+1,$istop);
        my $acc=`indexFasSeq ./gDNA_templates/${ctg_id}.fas.ind $accBeg $accEnd`;
        chomp($acc);
        if($acc eq 'AG') {
          $record.="$ctg_id\tanchor\tass\t$accEnd\t$accEnd\t.\t-\t.\tsrc=M;grp=gene.${cds_cnt};pri=$PRIORITY\n";
        }
        
        push(@introns,$record);
      }
    }
  }

  # >= 1 exons per CDS
  my @exons=();
  my $firstp=1;
  while($coords=~/(\d+)  (\d+)/g) {
    my($estart,$estop)=($1,$2);
    if($orient eq '+') { # forward strand
      my $record="$ctg_id\tanchor\texonpart\t$estart\t$estop\t.\t+\t.\tsrc=M;grp=gene.${cds_cnt};pri=$PRIORITY\n";

      if($firstp) {
        my($metBeg,$metEnd)=($estart,$estart+2);
        my $seq=`indexFasSeq ./gDNA_templates/${ctg_id}.fas.ind $metBeg $metEnd`;
        chomp($seq);
        if($seq eq 'ATG') {
          $record.="$ctg_id\tanchor\tstart\t$metBeg\t$metEnd\t.\t+\t.\tsrc=M;grp=gene.${cds_cnt};pri=$PRIORITY\n";
        }
        $firstp=0;
      }

      push(@exons,$record);
    }
    else { # reverse strand
      my $record="$ctg_id\tanchor\texonpart\t$estop\t$estart\t.\t-\t.\tsrc=M;grp=gene.${cds_cnt};pri=$PRIORITY\n";

      if($firstp) {
        my($metBeg,$metEnd)=($estart,$estart-2);
        my $seq=`indexFasSeq ./gDNA_templates/${ctg_id}.fas.ind $metBeg $metEnd`;
        chomp($seq);
        if($seq eq 'ATG') {
          $record.="$ctg_id\tanchor\tstart\t$metEnd\t$metBeg\t.\t-\t.\tsrc=M;grp=gene.${cds_cnt};pri=$PRIORITY\n";
        }
        $firstp=0;
      }

      push(@exons,$record);
    }
  }

  $exons[$#exons]=~/exonpart\t(\d+)\t(\d+)\t/;
  my($start,$stop)=($1,$2);
  if($orient eq '+') { # forward strand
    my($stopBeg,$stopEnd)=($stop-2,$stop);
    my $seq=`indexFasSeq ./gDNA_templates/${ctg_id}.fas.ind $stopBeg $stopEnd`;
    chomp($seq);
    if( ($seq eq 'TAA') || ($seq eq 'TAG') || ($seq eq 'TGA') ) {
      $exons[$#exons].="$ctg_id\tanchor\tstop\t$stopBeg\t$stopEnd\t.\t+\t.\tsrc=M;grp=gene.${cds_cnt};pri=$PRIORITY\n";
    }
  }
  else {
    my($stopBeg,$stopEnd)=($start+2,$start);
    my $seq=`indexFasSeq ./gDNA_templates/${ctg_id}.fas.ind $stopBeg $stopEnd`;
    chomp($seq);
    if( ($seq eq 'TAA') || ($seq eq 'TAG') || ($seq eq 'TGA') ) {
      $exons[$#exons].="$ctg_id\tanchor\tstop\t$stopEnd\t$stopBeg\t.\t-\t.\tsrc=M;grp=gene.${cds_cnt};pri=$PRIORITY\n";
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
