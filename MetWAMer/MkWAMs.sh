#!/usr/bin/env bash

# MkWAMs.sh - Script to guide derivation of start-Met WAMs using
#             gthXML data sets as training materials.

# Michael E Sparks (mespar1@gmail.com)
# Last modified: 20 July 2013

# Copyright (c) 2007,2013 Michael E Sparks
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

BIN="bin/"
PLBIN="src/plscripts/"
PREREQS="${BIN}/indexFasSeq \
         ${BIN}/parse_pps_codseqs \
         ${BIN}/print_MetWAM \
         ${BIN}/train_MetWAM"
RNGDIR="specs/"
TEMPLATES="data/gthxml_sample_input/fas/"
THREADINGS="data/gthxml_sample_input/xml/"

echo -e "\nRunning script ... `date`\n"
start_time=`date +%s`

# Verify that necessary executables are present
for prereq in $PREREQS; do
  if [ ! -e $prereq ]; then
    # build IMMpractical library
    cd IMMpractical/src; make redo; cd -
    # build MetWAMer executables
    cd src; make redo; make clean; cd -
    echo ""
  fi
done

# Process each gthXML file in turn
cat /dev/null > codseqs.fas
for file in `ls ${THREADINGS}`; do
  handle=`echo "$file" | cut -f 1 -d'.'`
  echo -n "Processing $handle ... "

  $BIN/indexFasSeq $TEMPLATES/$handle.fas
  $BIN/parse_pps_codseqs \
    -f $THREADINGS/$file \
    -g $TEMPLATES/$handle.fas.ind \
    -r $RNGDIR/GenomeThreader.rng \
    >> codseqs.fas

  echo "done"
done
rm -f $TEMPLATES/*.fas.ind
perl -e '$i=1; while($line=<>){ if($line=~m/^>/){ $line=~s/^> \d+/> $i/; ++$i; } print $line; }' -i codseqs.fas

# Build the WAM for translation start Met's and
# assert that it is probabilistically sound.
$BIN/train_MetWAM -h 1 -f codseqs.fas -r parmsT -o 6 -e 3
$BIN/train_MetWAM -h 0 -f codseqs.fas -r parmsF -o 6 -e 3
rm -f codseqs.fas
$BIN/print_MetWAM parmsT.MetWAM 2>&1 | $PLBIN/verify_pmf.pl
$BIN/print_MetWAM parmsF.MetWAM 2>&1 | $PLBIN/verify_pmf.pl

end_time=`date +%s`
et=`echo "$end_time - $start_time" | bc`
let "minutes=et / 60"
echo -e "\nScript finished!  Took ~$minutes minutes\a\n"

exit 0
