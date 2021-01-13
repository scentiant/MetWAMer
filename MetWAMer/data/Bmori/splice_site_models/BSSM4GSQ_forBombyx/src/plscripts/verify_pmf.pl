#!/usr/bin/perl -w

# verify_pmf.pl - Verify PMFs
#   Apart from the all-0.0000 entries near don/acc dinucleotides, we want
#   to guarantee that no row deviates from a probability mass function
#   (sums to 1.0000) by any more than 0.0005.  We will test each line
#   individually and uniquely sort the results.  Your output should
#   only be "Problematic case! : 0.0000 0.0000 0.0000 0.0000"

# Michael E Sparks (mespar1@iastate.edu)
# Last modified: 5 April 2007

# Copyright (c) 2005,2007 Michael E Sparks
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
use POSIX;

my %problems=();
while (my $line=<>) {
  if ($line =~ m/^Term/) {
    for (my $i=0;$i<6;++$i) {
      $line=<>;
    }
  }
  if ($line =~ m/(\d\.\d{4}) (\d\.\d{4}) (\d\.\d{4}) (\d\.\d{4})/) {
    if ( fabs(($1 + $2 + $3 + $4) - 1.0000) > 0.0005) {
      chomp($line);
      $problems{$line}=1;
    }
  }
}

foreach my $problem (keys(%problems)) {
  print $problem,"\n";
}

exit 0;
