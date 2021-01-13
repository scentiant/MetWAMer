#!/usr/bin/perl -w
# Michael E Sparks (mespar1@gmail.com)
# 14 Dec 2005

# This script will take input, presumed to be in frame 0
# and length of an exact multiple of 3, and it will write
# output to 6 files, prefixed with $HANDLE, and with
# extensions as follows:
# Model 0 (C g t) (forward)
#       1 (c G t) (forward)
#       2 (c g T) (forward)
#       3 (A c g) (reverse)
#       4 (a C g) (reverse)
#       5 (a c G) (reverse)

# Copyright (C) 2005  Michael E Sparks
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

use strict;

my $HANDLE = shift or die "script.pl HANDLE < data.fas\n";
my @file = ();

for(my $i=0;$i<6;++$i) {
  $file[$i] = $HANDLE.'.'.$i;
  `cat /dev/null > $file[$i]`;
}

my $sequence = "";
my $desc=<>;
if ($desc !~ m/^>/) { die "Something's amiss!\n"; }
else { chomp($desc); }

while(my $line=<>) {
  chomp($line);
  if ($line =~ m/^>/) {
    for(my $i=0;$i<6;++$i) {
      my $result = &framemaker($i,$sequence);
      my $newlen = length($result);
      `echo "$desc (revised length: $newlen)\n$result" >> $file[$i]`;
    }
    $sequence = "";
    $desc=$line;
  }
  else { $sequence .= $line; }
}
for(my $i=0;$i<6;++$i) {
  my $result = &framemaker($i,$sequence);
  my $newlen = length($result);
  `echo "$desc (revised length: $newlen)\n$result" >> $file[$i]`;
}

exit 0;

sub framemaker {
  my $frame = shift(@_);
  my $string = shift(@_) or die;

  if ($frame == 0) {
    ;
  }
  elsif ($frame == 1) {
    $string = reverse $string;
    chop $string;
    $string = reverse $string;
  }
  elsif ($frame == 2) {
    $string = reverse $string;
    chop $string;
    chop $string;
    $string = reverse $string;
  }
  elsif ($frame == 3) {
    $string = &revcmp($string);
    if ((length($string) % 3) eq 0) {
      ;
    }
    elsif ((length($string) % 3) == 1) {
      $string = reverse $string;
      chop $string;
      $string = reverse $string;
    }
    else { # (length($string) % 3 == 2)
      $string = reverse $string;
      chop $string;
      chop $string;
      $string = reverse $string;
    }
  }
  elsif ($frame == 4) {
    $string = &revcmp($string);
    if (length($string) % 3 == 0) {
      $string = reverse $string;
      chop $string;
      $string = reverse $string;
    }
    elsif (length($string) % 3 == 1) {
      $string = reverse $string;
      chop $string;
      chop $string;
      $string = reverse $string;
    }
    else { # (length($string) % 3 == 2)
      ;
    }
  }
  elsif ($frame == 5) {
    $string = &revcmp($string);
    if (length($string) % 3 == 0) {
      $string = reverse $string;
      chop $string;
      chop $string;
      $string = reverse $string;
    }
    elsif (length($string) % 3 == 1) {
      ;
    }
    else { # (length($string) % 3 == 2)
      $string = reverse $string;
      chop $string;
      $string = reverse $string;
    }
  }

  return uc $string;
}

sub revcmp {
  my $string = shift(@_) or die;
  my @oldstring = split(//,$string);

  my @newstring = ();
  for(my $i=$#oldstring,my $j=0;$i>=0;--$i,++$j) {
    if    (($oldstring[$i] eq 'A') || ($oldstring[$i] eq 'a')) {
      $newstring[$j] = 'T';
    }
    elsif (($oldstring[$i] eq 'C') || ($oldstring[$i] eq 'c')) {
      $newstring[$j] = 'G';
    }
    elsif (($oldstring[$i] eq 'G') || ($oldstring[$i] eq 'g')) {
      $newstring[$j] = 'C';
    }
    elsif (($oldstring[$i] eq 'T') || ($oldstring[$i] eq 't')) {
      $newstring[$j] = 'A';
    }
    else {
      $newstring[$j] = uc $oldstring[$i];
    }
  }

  my $newstring = join("",@newstring);
  return $newstring;
}
