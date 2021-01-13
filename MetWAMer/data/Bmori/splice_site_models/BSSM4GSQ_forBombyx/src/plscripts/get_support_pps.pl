#!/usr/bin/perl -w

# Michael E Sparks (mespar1@iastate.edu)
# Last modified: 6 April 2007

# get_support_pps.pl - Gets supporting PPS information
#                      from spliced alignment data.

# Copyright (c) 2004,2007 Michael E Sparks
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

my $USAGE = "$0 [plain|gthxml] cutoff goodexonct infile\n";
my ($intype,$cutoff,$infile,$count,$line,$reference,$score,
    $test,$ex_count,$ex_coors,$estart,$estop) = "";
my %references;

if($#ARGV != 3) {
  die "$USAGE";
}
else { # can't shift a zero from the command line
  $intype=$ARGV[0]; # plain text or xml
  $cutoff=$ARGV[1]; # min score of qualifying exon
  $count=$ARGV[2]; # Min number adequate scoring exons per structure
  $infile=$ARGV[3]; # Input file to process
}

if (($intype ne "plain") && ($intype ne "gthxml")) { die "$USAGE"; }
open(FIN,"<$infile") or die "Can't open $infile! $!\n";

# This loop will populate %reference, which stores the names of all
# transcripts such that at least $count many exons in their alignment
# exhibit an alignment score >= $cutoff.
# Names of qualifying aligned transcript sequences are stored as keys,
# and a list of exons as values, the latter of which allow us to distinguish
# among alignments at multiple loci on a given template.
while ($line=<FIN>) {
  if ($intype eq "plain") {
    if ($line !~ m/^Predicted gene structure/) { next; }
    else {
      $line=<FIN>;
      for($line=<FIN>,$ex_count=0,$ex_coors="";
          ($line =~ m/\S/);
          $line=<FIN>) {
        if ($line =~ m/Intron/) { next; }
        else {
          if ($line !~ m/Exon/)  {
            if ($line !~ m/PPA/) { die "Unexpected error occurred.\n"; }
            else { last; }
          }
          $line =~ m/Exon\s+\d+\s+(\d+)\s+(\d+).*?score:\s+(\S+)$/;
          ($estart,$estop,$score)=($1,$2,$3);
          $ex_coors .= "$estart $estop\t";
          if ($score >= $cutoff) { # testing score of each exon
            ++$ex_count;
          }
        }
      }
      chop($ex_coors);
      until ($line =~ m/^MATCH/) { $line=<FIN>; }
      $line =~ m/^MATCH\s+\S+\s+(\S+)[+-]/; $reference=$1;
      if($ex_count >= $count) {
        if(! defined $references{$reference} ) {
          $references{$reference}=();
        }
        push @{$references{$reference}}, $ex_coors;
      }
    }
  }
  elsif ($intype eq "gthxml") {
    if ($line !~ m/<predicted_gene_structure/) { next; }
    else {
      $line=<FIN>; # discard opening <exon-intron_info> tag
      for($line=<FIN>,$ex_count=0,$ex_coors="";
          ($line !~ m/<\/exon-intron_info>/);
          $line=<FIN>) {
        if ($line =~ m/<intron/) {
          until($line =~ m/<\/intron/) { $line=<FIN>; }
        }
        else {
          if ($line !~ m/<exon/) { die "Unexpected error occurred.\n"; }
          $line=<FIN>;
          $line =~ m/g_start=\"(\d+)\" g_stop=\"(\d+)\"/;
          ($estart,$estop)=($1,$2);
          $ex_coors .= "$estart $estop\t";
          $line=<FIN>;
          $line =~ m/r_score=\"(\S+)\"/; $score=$1;
          if ($score >= $cutoff) { ++$ex_count; }
          $line=<FIN>; # exon element closing tag
        }
      }
      chop($ex_coors);
      until ($line =~ m/<MATCH_line/) { $line=<FIN>; }
      $line =~ m/ref_id=\"(\S+)\"/; $reference=$1;
      if($ex_count >= $count) {
        if(! defined $references{$reference} ) {
          $references{$reference}=();
        }
        push @{$references{$reference}}, $ex_coors;
      }
    }
  }
  else { die "Unexpected error\n"; }
}

# rewind file
seek(FIN,0,0);

# Report PPS entries that have adequate EST supporting evidence
PPSLOOP: for ($line=<FIN>;defined($line); ) {
  if ($intype eq "plain") {
    if ($line !~ m/^PGL/) { $line=<FIN>; next; }
    until ($line =~ m/^\s{2}PGS.*?\S+[+-]$/) { $line=<FIN>; }
    while ($line =~ m/^\s{2}PGS.*?\S+[+-]$/) {
      $ex_coors="";
      while($line =~ m/(\d+)  (\d+)/g) {
        $ex_coors .= "$1 $2\t";
      }
      chop($ex_coors);
      $line =~ m/^\s{2}PGS.*?(\S+)[+-]$/; $test=$1;
      my $testflag=0;
      if (defined($references{$test})) {
        foreach my $ex_string (@{$references{$test}}) {
          if($ex_coors eq $ex_string) { # adjust code if > 1 good
                                        # EST alignment is required
                                        # to support a valid PPS
            $testflag=1;
            last;
          }
        }
      }
      if ($testflag) { # Print PPS
        until ($line =~ m/^Maximal non-overlapping open reading frames/) {
          $line=<FIN>;
        }
        $line=<FIN>;
        if ($line !~ m/^none/) {
          $line=<FIN>;
          # We only want to consider the first PPS,
          # as it is (usually) the longest, and therefore,
          # the most reliable.
          if ($line !~ m/frame/) { die "Inconsistent file format.\n"; }
          if($line !~ m/\(.*?,.*?\)\s+\(frame/) { # Must contain >= 1 intron
            next PPSLOOP;
          }
          print $line;

          # parse protein sequence
          for($line=<FIN>;($line !~ m/^$/);$line=<FIN>) {
            $line =~ s/[\d\s]//g;
            print "$line\n";
          }
        }
        last;
      }
      $line=<FIN>;
    } # end while
  }
  elsif ($intype eq "gthxml") {
    if ($line !~ m/<PGL_line/) { $line=<FIN>; next; }
    until ($line =~ m/<PGS_line/) { $line=<FIN>; }
    while ($line =~ m/<PGS_line/) {
      $line=<FIN>;
      $line=<FIN>; # first <exon element
      $ex_coors="";
      until($line =~ m/<\/gDNA_exon_coordinates/) {
        $line=~m/<exon start=\"(\d+)\" stop=\"(\d+)\"/;
        $ex_coors .= "$1 $2\t";
        $line=<FIN>;
      }
      chop($ex_coors);
      $line=<FIN>;
      $line =~ m/id=\"(\S+)\"/; $test=$1;
      my $testflag=0;
      if (defined($references{$test})) {
        foreach my $ex_string (@{$references{$test}}) {
          if($ex_coors eq $ex_string) { # adjust code if > 1 good
                                        # EST alignment is required
                                        # to support a valid PPS
            $testflag=1;
            last;
          }
        }
      }
      if ($testflag) { # Print PPS
        until ($line =~ m/<probable_ORFs/) { $line=<FIN>; }
        $line=<FIN>;
        if ($line !~ m/<none/) {
          # We only want to consider the first PPS,
          # as it is (usually) the longest, and therefore,
          # the most reliable.
          if($line !~ m/<orf_entry/) { die "File format error.\n"; }
          my $printable="";
          $line=<FIN>;
          $line=<FIN>;
          $line=~m/id=\"(\S+)\" strand=\"([+-])\"/;
          $printable = ">$1$2_";
          $line=<FIN>;
          $line=~m/PGL_serial=\"(\d+)\" AGS_serial=\"(\d+)\" PPS_serial=\"(\d+)\"/;
          $printable .= "PGL-$1_AGS-$2_PPS_$3 (";
          until ($line =~ m/<exon_boundaries/) { $line=<FIN>; }
          my $passfirst=0;
          for($line=<FIN>;($line=~m/<exon/);$line=<FIN>) {
            if ($passfirst) { $printable .= ","; }
            $line=~m/<exon start=\"(\d+)\" stop=\"(\d+)\"/;
            $printable .= "$1  $2";
            $passfirst+=1;
          }
          if ($passfirst < 2) { # case with no introns
            last;
          }
          print $printable,")\t(frame '";
          $line=<FIN>;
          $line=~m/<frame>(\d)/;
          print "$1';";
          $line=<FIN>;
          $line=~m/<number_coding_nucleotides>(\d+)/;
          printf("%6i bp,",$1);
          $line=<FIN>;
          $line=~m/<number_encoded_amino_acids>(\d+)/;
          printf("%5i residues) \n",$1);

          # print protein sequence
          until ($line =~ m/<predicted_protein_sequence/) { $line=<FIN>; }
          $line=~m/>(\S+)</;
          my $sequence=$1;
          while(length($sequence)) {
            $sequence =~s/^(.{1,60})//;
            print "$1\n";
          }
        } # end if not none
        last;
      }
      else {
        $line=<FIN>; # discard </PGS_line> and try again
        $line=<FIN>;
      }
    } # end PGS_line while
  }
  else { die "Unexpected error\n"; }
}

close FIN;

exit 0;
