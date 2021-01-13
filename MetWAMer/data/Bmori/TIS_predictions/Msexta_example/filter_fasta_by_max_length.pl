#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 13 January 2021

# filter_fasta_by_max_length.pl - utility that reads in a Fasta-formatted
# file of nucleotide sequences, and prints those sequences having length
# .le. to MAXLEN to stdout.  The length dist could be written to stderr,
# if so desired.

my $MAXLEN=shift or die "$0 max_len\n";
my $seqlen=-1;

my $desc="";
my $seq="";

while(my $line=<>) {
  if($line=~m/^>/) { # header
    if($seqlen >= 0) { # nothing to report for 1st instance
      if($seqlen <= $MAXLEN) {
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
  if($seqlen <= $MAXLEN) {
    print $desc,$seq,"\n";
  }
  #print STDERR $seqlen,"\n";
}

exit;
