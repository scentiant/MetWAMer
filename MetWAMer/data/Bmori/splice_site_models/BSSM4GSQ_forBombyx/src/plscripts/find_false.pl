#!/usr/bin/perl -w

# find_false.pl - This script parses all false GT or AG sites
#   in either exons or introns

# Michael E Sparks (mespar1@iastate.edu)
# Last modified: 4-11-04

# Copyright (c) 2004 Michael E Sparks
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
use Getopt::Std;

my %types = (
  "e" => "exon",
  "i" => "intron",
);
my %pole = (
  "d" => "donor",
  "a" => "acceptor",
);
my $USAGE="./find_false.pl [eE|iI] in_file bin ind_dir don_dinuc acc_dinuc out_dir\n";
my $FLANK=50; # symmetric extension

my ($don,$acc,$start,$stop,$temp_id,$bin,$type,$line,$linelen,
    $infile,$index_dir,$desc,$out_dir,$in_dir,$phase,
    $faldon_start,$faldon_stop,$falacc_start,$falacc_stop,
    $faldon_phase,$falacc_phase) = "";

# get parms
if (@ARGV != 7) { die "$USAGE"; }
else {
  use vars qw($opt_e $opt_E $opt_i $opt_I);
  getopts('eEiI');

  if    ($opt_e || $opt_E) { $type=$types{e} }
  elsif ($opt_i || $opt_I) { $type=$types{i} }
  else                     { die "$USAGE\n"; }

  $infile=shift(@ARGV); # file to process
  $bin=shift(@ARGV); # path to indexFasSeq executable
  $index_dir=shift(@ARGV); # path to indexFasSeq'ed indices
  $don=shift(@ARGV); # donor dinucleotide
  $acc=shift(@ARGV); # acceptor dinucleotide
  $out_dir=shift(@ARGV); # directory to write output to
}

# open output files
if ($type eq $types{i}) {
  open (FIN, "<$infile")
    or die "Can't open $infile";
  $infile =~ s/^(.*\/)//g;
  open (FOUTd, ">${out_dir}/${infile}.${don}_${acc}.falsedons")
    or die "Can't open ${out_dir}/${infile}.${don}_${acc}.falsedons";
  open (FOUTa, ">${out_dir}/${infile}.${don}_${acc}.falseaccs")
    or die "Can't open ${out_dir}/${infile}.${don}_${acc}.falseaccs";
}
else { # ($type eq $types{e})
  open (FIN, "<$infile")
    or die "Can't open $infile";
  $infile =~ s/^(.*\/)//g;
  open (FOUTd0, ">>${out_dir}/${infile}.${don}_${acc}.falsedons0")
    or die "Can't open ${out_dir}/${infile}.${don}_${acc}.falsedons0";
  open (FOUTa0, ">>${out_dir}/${infile}.${don}_${acc}.falseaccs0")
    or die "Can't open ${out_dir}/${infile}.${don}_${acc}.falseaccs0";
  open (FOUTd1, ">>${out_dir}/${infile}.${don}_${acc}.falsedons1")
    or die "Can't open ${out_dir}/${infile}.${don}_${acc}.falsedons1";
  open (FOUTa1, ">>${out_dir}/${infile}.${don}_${acc}.falseaccs1")
    or die "Can't open ${out_dir}/${infile}.${don}_${acc}.falseaccs1";
  open (FOUTd2, ">>${out_dir}/${infile}.${don}_${acc}.falsedons2")
    or die "Can't open ${out_dir}/${infile}.${don}_${acc}.falsedons2";
  open (FOUTa2, ">>${out_dir}/${infile}.${don}_${acc}.falseaccs2")
    or die "Can't open ${out_dir}/${infile}.${don}_${acc}.falseaccs2";
}

# process input
while ($line=<FIN>) {
  chomp($line);
  $desc=$line;
  if ($type eq $types{i}) {
    $desc =~ m/^>\s+(\d+) (\d+)\s+\d\s+(\S+)[+-]_/;
    ($start,$stop,$temp_id)=($1,$2,$3);
  }
  else { # ($type eq $types{e})
    $desc =~ m/^>\s+(\d+) (\d+) (\d)\s+(\S+)[+-]_/;
    ($start,$stop,$phase,$temp_id)=($1,$2,$3,$4);
  }

  $line=<FIN>; chomp($line);
  $linelen=length($line);

  # Introns
  if ($type eq $types{i}) {
    while ($line =~ m/$don/gi) {
      $faldon_start=pos($line) - 1;
      if ( (pos($line) - 1) != 1 ) {
        $faldon_start += $start;
        --$faldon_start;
        $faldon_stop = $faldon_start + ($FLANK + 1);
        $faldon_start -= $FLANK;
        &print_seq($pole{d},$temp_id,$faldon_start,$faldon_stop,$desc);
      }
    }
    while ($line =~ m/$acc/gi) {
      $falacc_start=pos($line) - 1;
      if ( (pos($line) - 1) != ($linelen - 1) ) {
        $falacc_start += $start;
        --$falacc_start;
        $falacc_stop = $falacc_start + ($FLANK + 1);
        $falacc_start -= $FLANK;
        &print_seq($pole{a},$temp_id,$falacc_start,$falacc_stop,$desc);
      }
    }
  }

  # Exons
  elsif ($type eq $types{e}) {
    while ($line =~ m/$don/gi) {
      $faldon_phase=(pos($line) - 1) % 3;
      $faldon_phase += $phase;
      $faldon_phase %= 3;
      $faldon_start=pos($line) - 1;
      $faldon_start += $start;
      --$faldon_start;
      $faldon_stop = $faldon_start + ($FLANK + 1);
      $faldon_start -= $FLANK;
      &print_seq_CDNA($pole{d},$faldon_phase,$temp_id,$faldon_start,$faldon_stop,$desc);
    }
    while ($line =~ m/$acc/gi) {
      $falacc_phase=(pos($line) - 1) % 3;
      $falacc_phase += $phase;
      $falacc_phase %= 3;
      $falacc_start=pos($line) - 1;
      $falacc_start += $start;
      --$falacc_start;
      $falacc_stop = $falacc_start + ($FLANK + 1);
      $falacc_start -= $FLANK;
      &print_seq_CDNA($pole{a},$falacc_phase,$temp_id,$falacc_start,$falacc_stop,$desc);
    }
  }
  else {
    die "Unexpected error\n";
  }
} # end while

exit 0;

sub print_seq {
  my $var=shift(@_); # donor dinucleotide
  my $temp_id=shift(@_); # genomic template to extract from
  my $beg=shift(@_); # start
  my $end=shift(@_); # stop
  my $line=shift(@_) or die "Improper print_seq usage\n";
  $line .= "FALSE";

  if ($var eq $pole{d}) {
    # Get donor data
    my $don_seq = `$bin/indexFasSeq $index_dir/${temp_id}.fas.ind $beg $end`;
    if ($don_seq eq "") {
      print STDERR "Note: Flanks extend beyond genomic template at
$bin/indexFasSeq $index_dir/${temp_id}.fas.ind $beg $end
(All's well, though!)\n";
      return;
    }
    else { $don_seq = uc($don_seq); }
    if ($don_seq =~ m/^.{50}$don.{50}$/) {
      print FOUTd "$line\n$don_seq\n";
    }
  }
  elsif ($var eq $pole{a}) {
    # Get acceptor data
    my $acc_seq = `$bin/indexFasSeq $index_dir/${temp_id}.fas.ind $beg $end`;
    if ($acc_seq eq "") {
      print STDERR "Note: Flanks extend beyond genomic template at
$bin/indexFasSeq $index_dir/${temp_id}.fas.ind $beg $end
(All's well, though!)\n";
      return;
    }
    else { $acc_seq = uc($acc_seq); }
    if ($acc_seq =~ m/^.{50}$acc.{50}$/) {
      print FOUTa "$line\n$acc_seq\n";
    }
  }
  else {
    die "Unexpected var in print_seq\n";
  }

  return;
} # end print_seq

# Function to print the false sites from the cdna.  Will mark the phase of the false
# site; these can be sorted in a post-processing step from the automation script.
sub print_seq_CDNA {
  my $var=shift(@_); # acceptor dinucleotide
  my $phase=shift(@_); # 0,1,2
  my $temp_id=shift(@_); # genomic template to extract from
  my $beg=shift(@_); # start
  my $end=shift(@_); # stop
  my $line=shift(@_) or die "Improper print_seq_CDNA usage\n";
  $line .= "FALSE";

  if ($var eq $pole{d}) {
    # Get donor data
    my $don_seq = `$bin/indexFasSeq $index_dir/${temp_id}.fas.ind $beg $end`;
    if ($don_seq eq "") {
      print STDERR "Note: Flanks extend beyond genomic template at
$bin/indexFasSeq $index_dir/${temp_id}.fas.ind $beg $end
(All's well, though!)\n";
      return;
    }
    else { $don_seq = uc($don_seq); }
    if ($don_seq =~ m/^.{50}$don.{50}$/) {
      if    ($phase == 0) { print FOUTd0 "$line $phase\n$don_seq\n"; }
      elsif ($phase == 1) { print FOUTd1 "$line $phase\n$don_seq\n"; }
      elsif ($phase == 2) { print FOUTd2 "$line $phase\n$don_seq\n"; }
      else { die "Bad phase passed\n"; }
    }
  }
  elsif ($var eq $pole{a}) {
    # Get acceptor data
    my $acc_seq = `$bin/indexFasSeq $index_dir/${temp_id}.fas.ind $beg $end`;
    if ($acc_seq eq "") {
      print STDERR "Note: Flanks extend beyond genomic template at
$bin/indexFasSeq $index_dir/${temp_id}.fas.ind $beg $end
(All's well, though!)\n";
      return;
    }
    else { $acc_seq = uc($acc_seq); }
    if ($acc_seq =~ m/^.{50}$acc.{50}$/) {
      if    ($phase == 0) { print FOUTa0 "$line $phase\n$acc_seq\n"; }
      elsif ($phase == 1) { print FOUTa1 "$line $phase\n$acc_seq\n"; }
      elsif ($phase == 2) { print FOUTa2 "$line $phase\n$acc_seq\n"; }
      else { die "Bad phase passed\n"; }
    }
  }
  else {
    die "Unexpected var\n";
  }

  return;
} # end print_seq_CDNA
