#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 13 January 2021

# if the gene/group conforms to the ideal gene grammar,
# then report it as an anchor; otherwise, as a hint.

# An instance conforming to such a grammar:
#NC_051115.1     anchor  exonpart        1137245 1139350 .       -       .       src=M;grp=gene.22;pri=10
#NC_051115.1     anchor  start   1139348 1139350 .       -       .       src=M;grp=gene.22;pri=10
#NC_051115.1     anchor  intronpart      1135000 1137244 .       -       .       src=M;grp=gene.22;pri=10
#NC_051115.1     anchor  dss     1137244 1137244 .       -       .       src=M;grp=gene.22;pri=10
#NC_051115.1     anchor  ass     1135000 1135000 .       -       .       src=M;grp=gene.22;pri=10
#NC_051115.1     anchor  exonpart        1134967 1134999 .       -       .       src=M;grp=gene.22;pri=10
#NC_051115.1     anchor  stop    1134967 1134969 .       -       .       src=M;grp=gene.22;pri=10
#NC_051115.1     anchor  exonpart        1137245 1139350 .       -       .       src=M;grp=gene.23;pri=10
#NC_051115.1     anchor  start   1139348 1139350 .       -       .       src=M;grp=gene.23;pri=10
#NC_051115.1     anchor  intronpart      1135927 1137244 .       -       .       src=M;grp=gene.23;pri=10
#NC_051115.1     anchor  dss     1137244 1137244 .       -       .       src=M;grp=gene.23;pri=10
#NC_051115.1     anchor  ass     1135927 1135927 .       -       .       src=M;grp=gene.23;pri=10
#NC_051115.1     anchor  exonpart        1135819 1135926 .       -       .       src=M;grp=gene.23;pri=10
#NC_051115.1     anchor  stop    1135819 1135821 .       -       .       src=M;grp=gene.23;pri=10

sub check_grammar {
  my @states=@_;

  # cds should always begin as an exonpart
  # if last element of parse is stop codon, then preceding element should be either
  #   a) exonpart (if multi-exon or single-exon w/out canonical start) or
  #   b) start (if single-exon gene)
  if(($states[0] ne "exonpart") ||
     (($states[$#states] eq "stop") &&
      (($states[$#states-1] ne "exonpart") &&
       ($states[$#states-1] ne "start")))) {
    die "serious grammatical exception encountered\n";
  }

  # single-exon gene without a canonical start and stop site
  # CDS lacks canonical start site
  # CDS lacks canonical stop site
  if(($#states == 0) ||
     ($states[1] ne "start") ||
     ($states[$#states] ne "stop")) {
    return 0;
  }

  for(my $i=2;$i<=$#states;$i+=4) {
    if($states[$i] eq "stop") { return 1; }

    if($#states < ($i+3)) { return 0; }

    if(($states[$i] ne "intronpart") ||
       ($states[$i+1] ne "dss") ||
       ($states[$i+2] ne "ass") ||
       ($states[$i+3] ne "exonpart")) {
      return 0;
    }
  }

  die "should have exited from for-loop\n";
}

my $anchorsf=shift;
my $hintsf=shift or die "$0 anchor.gff hints.gff\n";

my $curr_grp="nOtHiNg.0";

open(ANCHORS,"<$anchorsf") or die "$!\n";
open(HINTS,"<$hintsf") or die "$!\n";

my $linea="";
my $lineh=<HINTS>;

my @anchors=();
my @parse=();

while($linea=<ANCHORS>) {
  $linea=~/anchor\t(\w+)\t.*?src=M;grp=(gene\.\d+);/;
  my($feature,$grp)=($1,$2);

  if($grp eq $curr_grp) {
    push(@anchors,$linea);
    push(@parse,$feature);
  }
  else {
    if($curr_grp ne "nOtHiNg.0") {
      my $grammar_good_p=&check_grammar(@parse);

      if($grammar_good_p) {
        for(my $i=0;$i<=$#anchors;++$i) { print $anchors[$i]; }
      }
      else {
        do {
          print $lineh;
          $lineh=<HINTS>;
        } while($lineh=~/P;grp=${curr_grp};pri/);
      }
    }

    # update for next gene/group
    $curr_grp=$grp;
    @anchors=();
    push(@anchors,$linea);
    @parse=();
    push(@parse,$feature);
    until($lineh=~/P;grp=${curr_grp};pri/) { $lineh=<HINTS>; }
  }
}

# process last gene/group
my $grammar_good_p=&check_grammar(@parse);

if($grammar_good_p) {
  for(my $i=0;$i<=$#anchors;++$i) { print $anchors[$i]; }
}
else {
  do {
    print $lineh;
    $lineh=<HINTS>;
  } while( (defined($lineh)) && ($lineh=~/P;grp=${curr_grp};pri/) );
}

close ANCHORS;
close HINTS;

exit 0;
