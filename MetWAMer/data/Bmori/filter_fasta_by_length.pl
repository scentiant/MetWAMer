#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 22 December 2020

# filter_fasta_by_length.pl - utility that reads in a Fasta-formatted
# file of nucleotide sequences, and prints those sequences having length
# .ge. MINLEN to stdout.  The length dist could be written to stderr,
# if so desired.

my $MINLEN=shift or die "$0 min_len\n";
my $seqlen=-1;

my $desc="";
my $seq="";

while(my $line=<>) {
  if($line=~m/^>/) { # header
    if($seqlen >= 0) { # nothing to report for 1st instance
      if($seqlen >= $MINLEN) {
        print $desc,$seq,"\n";
      }
      #print STDERR $seqlen,"\n";
    }
    $seqlen=0;
    $desc=$line;
    $seq="";
    next;
  }

  chomp($line);
  $line=~s/\s//g;

  $seqlen += length($line);
  $seq .= $line;
}

# Report length of last instance...if any.
if($seqlen >= 0) { # nothing to report if were 0 instances
  if($seqlen >= $MINLEN) {
    print $desc,$seq,"\n";
  }
  #print STDERR $seqlen,"\n";
}

exit;
