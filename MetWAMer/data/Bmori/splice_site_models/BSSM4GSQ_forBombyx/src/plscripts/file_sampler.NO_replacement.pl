#!/usr/bin/perl -w
use strict;

# file_sampler.NO_replacement.pl
# Michael E Sparks (mespar1@iastate.edu)
# Last modified: 31 October 2006

# This is a script to randomly sample some specified number of FASTA
# formatted sequences from a file WITHOUT REPLACEMENT.

# Copyright (c) 2006 Michael E Sparks
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

my $NUMSAMPLES = shift or die
"file_sampler.NO_replacement.pl num_samples < file_to_sample_from\n";
if ($NUMSAMPLES < 1) {
  die "Must sample 1 or more FASTA sequences!\n";
}

my @seqs=();
$/ = '>';
while(<>) {
  next unless $_ =~ /\w/;
  my $item = $_;
  $item=~s/^/>/;
  $item=~s/\n>/\n/;
  push(@seqs,$item);
}
if($NUMSAMPLES > ($#seqs + 1)) {
  print STDERR "Trying to sample more elements than in file!?? ";
  print STDERR "Shuffling file instead.\n";
  $NUMSAMPLES = $#seqs + 1;
}

my %keys=();
for(my $i = 0; $i < $NUMSAMPLES; ) {
  my $randnum=int(rand ($#seqs + 1));
  if(!exists($keys{$randnum})) {
    $keys{$randnum}=1;
    ++$i;
  }
}

foreach my $index (keys(%keys)) {
  print $seqs[$index];
}

exit 0;
