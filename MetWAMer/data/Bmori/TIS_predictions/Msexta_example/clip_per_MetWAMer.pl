#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 13 January 2021

use constant UPSTREXTENT => 5;
use constant DOWNSTREXTENT => 3;
use constant CONTENTSWATHLEN => 96;

my $FILE = shift or die; # sequences should be merged onto one line in advance

my %seqs=();
open(FNA,"<$FILE") or die "$!\n";
while(<FNA>) {
  chomp($_);
  $_=~/^>(\S+.*?)$/;
  my $key=$1;

  $_=<FNA>;
  chomp($_);
  my $seq=$_;
  if(! exists($seqs{$key}) ) {
    $seqs{$key}=$seq;
  }
  else {
    die "$key seen >= 2 times?\n";
  }
}
close FNA;

while(<>) {
  chomp($_);
  $_=~/^(\d+)\s+(>\S+\s.*?)$/;
  my ($atg,$desc)=($1,$2);
  my $ctgid=$desc;
  $ctgid=~s/^>//;
  --$atg;
  $desc.=" (MetWAMer ATG offset: ";
  $desc.=$atg-(CONTENTSWATHLEN+UPSTREXTENT);
  $desc.=" bp)\n";

  if(!exists($seqs{$ctgid})) {
    print STDERR "Not seen: $ctgid\n";
  }

  my $threePrimeClip=DOWNSTREXTENT+CONTENTSWATHLEN;
  my $seq=$seqs{$ctgid};
  $seq=~s/^\S{$atg}//;
  $seq=~s/\S{$threePrimeClip}$//;
  print $desc,$seq,"\n";
}

exit 0;
