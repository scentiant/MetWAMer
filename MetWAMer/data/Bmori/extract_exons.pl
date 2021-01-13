#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 22 December 2020

# assumes fasta sequences were merged onto single lines in advance

# Relevant statements from MetWAMer C sources:
#
##define UPSTREXTENT    5 /* Region to consider upstream of ATG   */
##define DOWNSTREXTENT  3 /* Region to consider downstream of ATG */
##define CONTENTSWATHLEN 96 /* Length of fragment to apply a content *
#                            * sensor, i.e., a Markov chain, over.   */
#
#  prependlen=(int)(UPSTREXTENT);
#  appendlen=(int)(DOWNSTREXTENT);
#  if(content_sensors_p==TRUE) {
#    prependlen+=(int)(CONTENTSWATHLEN);
#    appendlen+=(int)(CONTENTSWATHLEN);
#  }

my $PREPENDLEN=5+96;
my $APPENDLEN=3+96;

if(scalar @ARGV) { # pass an argument at CLI to add no flanking content
  $PREPENDLEN=$APPENDLEN=0;
}

my $INDDIR="./gdna_inds/";
my $PATH2IFS="../../bin/";

MAINLOOP: while(my $line=<STDIN>) {
  $line=~/^>lcl\|(\S+)_cds.*?\[location=(.*?)\]/;
  my ($gdna,$coords)=($1,$2);

  # disallow genes with questionable exon boundaries
  if( ($coords=~m/[<>-]/) ||
      ($coords=~m/\d+,\d+\)/) ||
      ($coords=~m/\(\d+,\d+/) ||
      ($coords=~m/\d+,\d+,\d+/) ) {
    $line=<STDIN>;
    next MAINLOOP;
  }

  my @exons=();
  while($coords=~m/(\d+\.\.\d+)/g) {
    push(@exons,$1);
  }

  # ensure strictly increasing/decreasing order of boundaries
  if($coords!~m/^complement/) { #forward strand
    $exons[0]=~m/^(\S+)\.\.(\S+)$/;
    my($prevstart,$prevstop)=($1,$2);
    if($prevstart > $prevstop) {
      print STDERR "Problem detected at line $. (ignoring gene) : $line";
      $line=<STDIN>;
      next MAINLOOP;
    }
    for(my $i=1;$i<=$#exons;++$i) {
      $exons[$i]=~m/^(\S+)\.\.(\S+)$/;
      my ($start,$stop)=($1,$2);
      if( ($prevstop >= $start) || ($start > $stop) ) {
        print STDERR "Problem detected at line $. (ignoring gene) : $line";
        $line=<STDIN>;
        next MAINLOOP;
      }
      else {
        $prevstop=$stop;
      }
    }
  }
  else {
    $exons[$#exons]=~m/^(\S+)\.\.(\S+)$/;
    my($prevstop,$prevstart)=($1,$2);
    if($prevstart < $prevstop) {
      print STDERR "Problem detected at line $. (ignoring gene) : $line";
      $line=<STDIN>;
      next MAINLOOP;
    }
    for(my $i=$#exons-1;$i>=0;--$i) {
      $exons[$i]=~m/^(\S+)\.\.(\S+)$/;
      my($stop,$start)=($1,$2);
      if( ($prevstop <= $start) || ($start < $stop) ) {
        print STDERR "Problem detected at line $. (ignoring gene) : $line";
        $line=<STDIN>;
        next MAINLOOP;
      }
      else {
        $prevstop=$stop;
      }
    }
  }

  my $seq="";
  if($coords!~m/^complement/) { #forward strand
    # extend terminal exon(s)
    $exons[0]=~m/^(\S+)\.\.(\S+)$/;
    my($start,$stop)=($1,$2);
    $start -= $PREPENDLEN;
    if($start < 0) {
      print STDERR "Extended beyond gDNA start at line $. (parsing $start..$stop) : $line";
      $line=<STDIN>;
      next MAINLOOP;
    }
    $exons[0]="${start}..${stop}";

    $exons[$#exons]=~m/^(\S+)\.\.(\S+)$/;
    ($start,$stop)=($1,$2);
    $stop += $APPENDLEN;
    # we don't know the gDNA length at this point, but indexFasSeq
    # will signal if we've tried to read beyond the 3' end of it
    $exons[$#exons]="${start}..${stop}";

    # extract sequence data
    for(my $i=0;$i<=$#exons;++$i) {
      $exons[$i]=~m/^(\S+)\.\.(\S+)$/;
      ($start,$stop)=($1,$2);
      my $frag=`${PATH2IFS}indexFasSeq ${INDDIR}/${gdna}.ind $start $stop`;
      # infer the $stop value exceeded $gdna's length (most likely reason)
      if($frag eq "") {
        print STDERR "Extended beyond gDNA end (?) at line $. (parsing $start..$stop) : $line";
        $line=<STDIN>;
        next MAINLOOP;
      }
      chomp($frag);
      $seq .= $frag;
    }
  }
  else { #reverse
    # extend terminal exon(s)
    $exons[0]=~m/^(\S+)\.\.(\S+)$/;
    my($stop,$start)=($1,$2);
    $stop -= $APPENDLEN;
    $exons[0]="${stop}..${start}";
    if($stop < 0) {
      print STDERR "Extended beyond gDNA start at line $. (parsing $start..$stop) : $line";
      $line=<STDIN>;
      next MAINLOOP;
    }

    $exons[$#exons]=~m/^(\S+)\.\.(\S+)$/;
    ($stop,$start)=($1,$2);
    $start += $PREPENDLEN; # rely on indexFasSeq to signal out-of-bounds error
    $exons[$#exons]="${stop}..${start}";

    # extract sequence data
    @exons = reverse @exons;
    for(my $i=0;$i<=$#exons;++$i) {
      $exons[$i]=~m/^(\S+)\.\.(\S+)$/;
      ($stop,$start)=($1,$2);
      my $frag=`${PATH2IFS}indexFasSeq ${INDDIR}/${gdna}.ind $start $stop`;
      if($frag eq "") { # if extensions go wrong (barring I/O errors, etc.)
        print STDERR "Extended beyond gDNA end (?) at line $. (parsing $start..$stop) : $line";
        $line=<STDIN>;
        next MAINLOOP;
      }

      if($start == $stop) {
        chomp($frag);
        if   ($frag eq 'A') { $frag="T\n"; }
        elsif($frag eq 'C') { $frag="G\n"; }
        elsif($frag eq 'G') { $frag="C\n"; }
        elsif($frag eq 'T') { $frag="A\n"; }
        else                { $frag.="\n"; } # handle ambiguity elsewhere
      }

      chomp($frag);
      $seq .= $frag;
    }
  }

  if($seq!~m/[xn]/i) { # no ambiguous residues in CDS
    print "${line}${seq}\n";
  }

  $line=<STDIN>; # discard original (unextended) fasta sequence
}

exit 0;
