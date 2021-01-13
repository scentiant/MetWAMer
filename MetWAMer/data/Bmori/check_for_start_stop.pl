#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 22 December 2020

# ensure only canonical TIS structures are used for training

# assumes fasta sequences were merged onto single lines in advance!

# Relevant statements copied from MetWAMer C sources:
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

while(my $desc=<>) {
  if($desc!~m/^>/) {
    die "fasta\n";
  }

  my $seq=<>;
  chomp($seq);
  my $seqr = reverse $seq;

  # recall we're looking at stop codons in reverse
  # orientation (not complemented, just reversed)
  if( ($seq=~m/^\w{101}ATG/) &&
      (($seqr=~m/^\w{99}AAT/) ||
       ($seqr=~m/^\w{99}GAT/) ||
       ($seqr=~m/^\w{99}AGT/))) {
    print "${desc}${seq}\n";
  }
  else {
    print STDERR "Ignoring ${desc}${seq}\n";
  }
}

exit 0;
