#!/usr/bin/perl -w

# verify_pmf.pl - Verify PMFs
#   Apart from the all-0.0000 entries near ATGs, and the
#   first slot's equilibrium frequencies, we want to
#   guarantee that no row in our (right) stochastic matrices
#   deviates from a probability mass function by any more
#   than 0.0005.

# Michael E Sparks (mespar1@gmail.com)
# Last modified: 11 July 2007

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

my $i=0;
my %problems=();
while(my $line=<>) {
  ++$i;
  if($line=~m/(\d\.\d{4}) (\d\.\d{4}) (\d\.\d{4}) (\d\.\d{4})/) {
    if($1==$2&&$1==$3&&$1==$4) {
      next;
    }
    elsif(fabs(($1+$2+$3+$4)-1.0000)>0.0005) {
      chomp($line);
      $problems{$line}=$i;
    }
  }
}

foreach my $problem (keys(%problems)) {
  print "Line ",$problems{$problem},": $problem\n";
}

exit 0;
