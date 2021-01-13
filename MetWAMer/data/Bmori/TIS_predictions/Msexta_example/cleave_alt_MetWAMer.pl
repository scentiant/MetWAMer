#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 13 January 2021

use constant UPSTREXTENT => 5;
use constant CONTENTSWATHLEN => 96;

#if Met found in first 45 bases of CDS, accept unconditionally
#else
#  if Met in first 2/10ths of overall CDS length (get this from header) -> keep
#  else -> ignore
#
#153     >projectID.1182+ (121581  121592,121850  122987,124014  124072)
#153     >projectID.1182+ (168199  168510)
#
#>>> 153 - 101
#52
#>>> 52.0 / 1209
#0.043010752688172046
#>>> 52.0 / 309
#0.16828478964401294

use constant ACCEPTRANGE => 45; # if Met in first 15 aa's, accept it
# technically, this could be a <= 43 test, but 45 works fine, too:
# 1   2   3      4        5           14       15
#123 456 789 10_11_12 13_14_15 ... 40_41_42 43_44_45
use constant COVERTHRESH => 0.2;

while(<>) {
  $_=~/^(\S+).*?([+-])\s+\((\d.*?\d)\)/;
  my ($pred,$orient,$coords)=($1,$2,$3);

  die if ($pred < 105); # -1 & 102 entries should not be in this file!
  $pred-=(CONTENTSWATHLEN+UPSTREXTENT);
  die if ($pred < 4); # sanity check on constants

  my $cdslen=0;
  while($coords=~/(\d+)  (\d+)/g) {
    my($start,$stop)=($1,$2);
    if($orient eq '+') {
      $cdslen += $stop - $start + 1;
    }
    else {
      $cdslen += $start - $stop + 1;
    }
  }

  if( ($pred < ACCEPTRANGE) || ( ($pred / $cdslen) < COVERTHRESH ) ) { print $_; }
  else { print STDERR $_; }
}

exit 0;
