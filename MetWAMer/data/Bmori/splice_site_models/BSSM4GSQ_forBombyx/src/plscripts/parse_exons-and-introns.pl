#!/usr/bin/perl -w

# parse_exons-and-introns.pl - Does what it says. Says what it does.

# Michael E Sparks (mespar1@iastate.edu)
# Last modified: 22 May 2007

# Copyright (c) 2003,2007 Michael E Sparks
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

# Import Packages
use strict;

# Variable declarations
my ($in_file, $index, $substring_path, $out0, $out1, $out2, $outI0, $outI1, $outI2,
    $line, $lineINTRONS, $phase, $e_start, $e_stop, $i_start, $i_stop, $reference_source) = "";
my @phases = ();

# Main application
MAIN: {
  # Verify command line parameters
  if (@ARGV != 9) {
    die
  "
    Usage: ./parse_exons.pl in_file index substring_path out0 out1 out2 outI0 outI1 outI2

     Note: in_file is an output file as produced by get_support_pps.pl
           index is an indexFasSeq index'ed file containing an appropriate template
           substring_path is the path to the indexFasSeq executable
           out0 stores C O D | (phase 0 exons)
           out1 stores C | O D (phase 1 exons)
           out2 stores C O | D (phase 2 exons)
           outI0 stores phase 0 introns
           outI1 stores phase 1 introns
           outI2 stores phase 2 introns

  ";
  }
  else { # Parse in the parameters
    $in_file = shift(@ARGV);
    $index = shift(@ARGV);
    $substring_path = shift(@ARGV);
    $out0 = shift(@ARGV);
    $out1 = shift(@ARGV);
    $out2 = shift(@ARGV);
    $outI0 = shift(@ARGV);
    $outI1 = shift(@ARGV);
    $outI2 = shift(@ARGV);
  }

  # Open the necessary files 
  open(FIN,"<$in_file") or die "Can't open $in_file $!\n";
  open(OUT0,">$out0") or die "Can't open $out0 $!\n";
  open(OUT1,">$out1") or die "Can't open $out1 $!\n";
  open(OUT2,">$out2") or die "Can't open $out2 $!\n";
  open(OUTI0,">$outI0") or die "Can't open $outI0 $!\n";
  open(OUTI1,">$outI1") or die "Can't open $outI1 $!\n";
  open(OUTI2,">$outI2") or die "Can't open $outI2 $!\n";

  # Process the input data
  while ($line = <FIN>) {
    chomp($line);
    if ($line !~ m/^>/) { next; }
    else {
      $lineINTRONS=$line; # Parallel string so I can parse out introns as well; see below
      $line =~ m/^>(.*?)\s+\(/;
      $reference_source=$1;

      $phase=0;

      # Initialize the phase array
      @phases=();
    }

    # Process the exons that comprise the PPS ORF
    while ($line =~ m/(\d+)  (\d+)/g) {
      ($e_start, $e_stop) = ($1, $2);

      my $exon_seq = `$substring_path/indexFasSeq $index $e_start $e_stop`;
      if ($exon_seq eq "") {
        print STDERR "Note: Flanks extend beyond genomic template at
$substring_path/indexFasSeq $index $e_start $e_stop
(All's well, though!)\n";
        next;
      }
      else {
        $exon_seq = uc($exon_seq);
      }

      # Print the exon to the appropriate file...if it is unambiguous!
      if ($exon_seq !~ m/[^ACGT]/) {
        if ($phase == 0) {
          print OUT0 ">$e_start $e_stop $phase $reference_source\n$exon_seq\n";
        }
        elsif ($phase == 1) {
          print OUT1 ">$e_start $e_stop $phase $reference_source\n$exon_seq\n";
        }
        elsif ($phase == 2) {
          print OUT2 ">$e_start $e_stop $phase $reference_source\n$exon_seq\n";
        }
        else {
          die "Error in phase assignment!\n";
        }
      }

      # Determine the phase of the (potentially) next exon
      my $exon_length = length($exon_seq);
      $phase = ($phase + $exon_length) % 3;

      push(@phases,$phase);
    } # end while loop for exon processing

    # Now, we process the succeeding intron, if it exists
    while ($lineINTRONS =~ m/(\d+),(\d+)/g) {
      ($i_start, $i_stop) = ($1, $2);
      if ($i_start < $i_stop) {
        ++$i_start;
        --$i_stop;
      }
      elsif ($i_start > $i_stop) {
        --$i_start;
        ++$i_stop;
      }
      else {
        die "Encountered a 1nt intron!??\n";
      }
      
      my $intron_seq = `$substring_path/indexFasSeq $index $i_start $i_stop`;
      if ($intron_seq eq "") {
        print STDERR "Note: Flanks extend beyond genomic template at
$substring_path/indexFasSeq $index $i_start $i_stop
(All's well, though!)\n";
        next;
      }
      else {
        $intron_seq = uc($intron_seq);
      }
     
      # Print the intron to the appropriate file...if it is unambiguous!
      if ($intron_seq !~ m/[^ACGT]/) {
        $phase=shift(@phases);
        if ($phase == 0) {
          print OUTI0 ">$i_start $i_stop $phase $reference_source\n$intron_seq\n";
        }
        elsif ($phase == 1) {
          print OUTI1 ">$i_start $i_stop $phase $reference_source\n$intron_seq\n";
        }
        elsif ($phase == 2) {
          print OUTI2 ">$i_start $i_stop $phase $reference_source\n$intron_seq\n";
        }
        else {
          die "Error in phase assignment!\n";
        }
      }
    } # end while, intron parsing
  } # end while

  # Clean up
  close(FIN);
  close(OUT0); close(OUT1); close(OUT2);
  close(OUTI0); close(OUTI1); close(OUTI2);
} # end main

exit 0;
