#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 22 December 2020

# assumes fasta sequences were merged onto single lines in advance

my $INDDIR="./gdna_inds/";
my $PATH2IFS="../../bin/";

MAINLOOP: while(my $line=<>) {
  $line=~/^>lcl\|(\S+)_cds.*?\[location=(.*?)\]/;
  my ($gdna,$coords)=($1,$2);

  # disallow genes with questionable intron boundaries
  if( ($coords=~m/[<>-]/) ||
      ($coords=~m/\d+,\d+\)/) ||
      ($coords=~m/\(\d+,\d+/) ||
      ($coords=~m/\d+,\d+,\d+/) ||
      ($coords!~m/,/) ) {
    $line=<>;
    next MAINLOOP;
  }

  my @introns=();
  while($coords=~m/(\.\.\d+,\d+)/g) {
    push(@introns,$1);
  }

  # ensure strictly increasing/decreasing order of boundaries
  if($coords!~m/^complement/) { #forward strand
    $introns[0]=~m/^\.\.(\S+),(\S+)$/;
    my($prevstart,$prevstop)=($1,$2);
    ++$prevstart; --$prevstop;
    if($prevstart > $prevstop) {
      print STDERR "Problem detected at line $. (ignoring gene) : $line";
      $line=<>;
      next MAINLOOP;
    }
    for(my $i=1;$i<=$#introns;++$i) {
      $introns[$i]=~m/^\.\.(\S+),(\S+)$/;
      my ($start,$stop)=($1,$2);
      ++$start; --$stop;
      if( ($prevstop >= $start) || ($start > $stop) ) {
        print STDERR "Problem detected at line $. (ignoring gene) : $line";
        $line=<>;
        next MAINLOOP;
      }
      else {
        $prevstop=$stop;
      }
    }
  }
  else {
    $introns[$#introns]=~m/^\.\.(\S+),(\S+)$/;
    my($prevstop,$prevstart)=($1,$2);
    --$prevstart; ++$prevstop;
    if($prevstart < $prevstop) {
      print STDERR "Problem detected at line $. (ignoring gene) : $line";
      $line=<>;
      next MAINLOOP;
    }
    for(my $i=$#introns-1;$i>=0;--$i) {
      $introns[$i]=~m/^\.\.(\S+),(\S+)$/;
      my($stop,$start)=($1,$2);
      --$start; ++$stop;
      if( ($prevstop <= $start) || ($start < $stop) ) {
        print STDERR "Problem detected at line $. (ignoring gene) : $line";
        $line=<>;
        next MAINLOOP;
      }
      else {
        $prevstop=$stop;
      }
    }
  }

  if($coords!~m/^complement/) { #forward strand
    # extract introns
    for(my $i=0;$i<=$#introns;++$i) {
      $introns[$i]=~m/^\.\.(\S+),(\S+)$/;
      my ($start,$stop)=($1,$2);
      ++$start; --$stop;
      my $frag=`${PATH2IFS}indexFasSeq ${INDDIR}/${gdna}.ind $start $stop`;
      # the following check is added for safety, but is unexpected:
      # intron coordinates are entirely contained in the gene under
      # consideration, and this gene is itself entirely contained in
      # the host genomic DNA sequence its coordinates are expressed in.
      if($frag eq "") {
        print STDERR "Problem detected at line $. (parsing $start..$stop) : $line";
        $line=<>;
        next MAINLOOP;
      }

      if($frag!~m/[xn]/i) { # no ambiguous residues in intron
        chomp($line);
        $line .= " (intron $start..$stop)";
        print "$line\n$frag";
      }
    }
  }
  else { #reverse
    for(my $i=0;$i<=$#introns;++$i) {
      $introns[$i]=~m/^\.\.(\S+),(\S+)$/;
      my ($stop,$start)=($1,$2);
      --$start; ++$stop;
      my $frag=`${PATH2IFS}indexFasSeq ${INDDIR}/${gdna}.ind $start $stop`;
      if($frag eq "") { # same reasoning as given above...
        print STDERR "Problem detected at line $. (parsing $start..$stop) : $line";
        $line=<>;
        next MAINLOOP;
      }

      # not that a single-base intron makes biological sense (AFAIK),
      # but we'll adjust $frag if it does somehow occur in NCBI's annotation
      if($start == $stop) {
        chomp($frag);
        if   ($frag eq 'A') { $frag="T\n"; }
        elsif($frag eq 'C') { $frag="G\n"; }
        elsif($frag eq 'G') { $frag="C\n"; }
        elsif($frag eq 'T') { $frag="A\n"; }
        else                { $frag.="\n"; } # handle ambiguity elsewhere
      }

      if($frag!~m/[xn]/i) { # permit no ambiguous residues in intron
        chomp($line);
        $line .= " (intron $start..$stop)";
        print "$line\n$frag";
      }
    }
  }

  $line=<>; # discard original fasta sequence
}

exit 0;
