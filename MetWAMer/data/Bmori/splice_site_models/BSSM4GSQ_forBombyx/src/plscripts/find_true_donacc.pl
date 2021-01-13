#!/usr/bin/perl -w

# find_true_donacc.pl - Gathers flanking dinucleotide data.

# Michael E Sparks (mespar1@iastate.edu)
# Last modified: 6-15-03

# Copyright (c) 2003 Michael E Sparks
# All rights reserved.

# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

use strict;

my ($infile, $don, $acc, $line,
    $don_start, $don_stop,
    $acc_start, $acc_stop,
    $bin, $index
   ) = "";

if (@ARGV != 5) {
  die "Usage: find_true_donacc.pl infile bin index don_dinuc acc_dinuc\n";
}
else {
  $infile = shift(@ARGV); # file to process
  $bin = shift(@ARGV); # path to indexFasSeq executable
  $index = shift(@ARGV); # path to indexFasSeq-derived indices
  $don = shift(@ARGV); $don = uc($don); # donor dinucleotide
  $acc = shift(@ARGV); $acc = uc($acc); # acceptor dinucleotide
}

open (FIN, "<$infile")
  or die "Can't open $infile";
open (FOUTd, ">${infile}.${don}_${acc}.truedons")
  or die "Can't open ${infile}.${don}_${acc}.truedons";
open (FOUTa, ">${infile}.${don}_${acc}.trueaccs")
  or die "Can't open ${infile}.${don}_${acc}.trueaccs";

while ($line = <FIN>) {
  if ($line !~ m/^>/) { next; }
  else {
    chomp($line);
    $line =~ m/^>\s+(\d+) (\d+)\s+\d\s+(\S+)[+-]_/;
    my ($i_start,$i_stop,$temp_id)=($1,$2,$3);

    my $test_seq = `$bin/indexFasSeq $index/${temp_id}.fas.ind $i_start $i_stop`;
    $test_seq = uc($test_seq);

    # Determine if intron has the correct don/acc pair we want data for
    if ($test_seq !~ m/^$don.*?$acc$/) {
      next;
    }
    else {
      if ($i_start < $i_stop) {
        $don_start = $i_start - 50;
        $don_stop  = $i_start + 51;
        $acc_start = $i_stop  - 51;
        $acc_stop  = $i_stop  + 50;
      }
      else {
        $don_start = $i_start + 50;
        $don_stop  = $i_start - 51;
        $acc_start = $i_stop  + 51;
        $acc_stop  = $i_stop  - 50;
      }

      # Get donor/acceptor data
      my $don_seq = `$bin/indexFasSeq $index/${temp_id}.fas.ind $don_start $don_stop`;
      my $acc_seq = `$bin/indexFasSeq $index/${temp_id}.fas.ind $acc_start $acc_stop`;

      if (($don_seq eq "") || ($acc_seq eq "")) {
        print STDERR "Note: Flanks extend beyond genomic template at
$bin/indexFasSeq $index/${temp_id}.fas.ind $don_start $don_stop
  .OR.
$bin/indexFasSeq $index/${temp_id}.fas.ind $acc_start $acc_stop
(All's well, though!)\n";
        next;
      }
      else {
        $don_seq = uc($don_seq);
        $acc_seq = uc($acc_seq);
      }

      print FOUTd "$line\n$don_seq\n";
      print FOUTa "$line\n$acc_seq\n";
    }
  }
}

close(FIN);
close(FOUTd);
close(FOUTa);

exit;
