#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 22 December 2020

# users may wish to confirm they won't bump into any
# number-of-files-in-a-folder limits on their local filesystem.
# also, ensure the script can access indexFasSeq at the rel. path below

while(my $line=<>) {
  $line=~m/^>(\S+)\s/;
  my $fh=$1;

  open(FOUT,">${fh}.fas") or die "$!\n";
  print FOUT $line;
  # assumes sequence was merged onto one line in advance
  $line=<>;
  print FOUT $line;
  close FOUT;

  `../../../bin/indexFasSeq ${fh}.fas ; mv ${fh}.fas.ind ${fh}.ind`;
}

exit 0;
