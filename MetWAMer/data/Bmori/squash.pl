#!/usr/bin/perl -w

# squash.pl

# Michael E Sparks (mespar1@iastate.edu)
# Last modified: 2 Feb 2003

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

my ($line1,$seq1,$desc1,
    $line2,$seq2,$desc2) = "";
my $counter=2;

$line1 = <>;
chomp($line1);
while ($line2 = <>) {
  chomp($line2);

  $line1 =~ m/^(.*?) (.*?)$/;
  ($seq1,$desc1)=($1,$2);
  $line2 =~ m/^(.*?) (.*?)$/;
  ($seq2,$desc2)=($1,$2);

  if ($seq1 ne $seq2) { # case where redundant entry doesn't present
    print STDOUT "$line1\n";
  }
  else {
    #print STDERR "Line $counter is redundant\n";
    $line2 .= " AND $desc1";
  }
    
  $line1=$line2;
  $counter++;
}

print STDOUT "$line1\n";

exit 0;
