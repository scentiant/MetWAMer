#!/usr/bin/perl -w
# Michael E Sparks (mespar1@gmail.com)

# This script can be used to eliminate ambiguous and too-short training
# sequence entries that may be present in a multi-FASTA data file.

# Copyright (C) 2005,2006 Michael E Sparks
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

# Import packages
use strict;

# Global "Constant" Parameters taken from ``immpractical.h"
# This variable is free to change from [0..5]
my $MAXORDER=5;

# Variable declarations
my($line,$comment,$sequence)="";
my $firstseq=1;

# Main Application
MAIN : {
  # Reminder to use redirect operators
  print STDERR "  Reading data from STDIN and writing to STDOUT...\n";

  # Process the data
  while($line=<>) {
    chomp($line);

    # Store comment line
    if ($line =~ m/^>/) {
      if ( (!$firstseq) &&
           ($sequence !~ m/[^AaCcGgTt]/) &&
           (length($sequence) > ($MAXORDER+1))
         ) {
        print STDOUT "$comment\n$sequence\n";
        $comment=$line;
        $sequence="";
      }
      elsif ($firstseq) {
        $comment=$line;
        $firstseq=0;
      }
      else {
        print STDERR "    Eliminated: $comment\n";
        $comment=$line;
        $sequence="";
      }
    }
    else { # accrue sequence
      $sequence .= $line;
    }
  } # end while

  if ( ($sequence !~ m/[^AaCcGgTt]/) &&
       (length($sequence) > ($MAXORDER+1))
     ) {
    print STDOUT "$comment\n$sequence\n";
  }
  else {
    print STDERR "    Eliminated: $comment\n";
  }

  print STDERR "  Finished cleaning the data!\n";
} # end MAIN

exit 0;
