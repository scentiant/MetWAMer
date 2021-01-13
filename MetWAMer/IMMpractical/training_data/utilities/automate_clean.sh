#!/bin/bash
# Script to automate use of clean_data.pl.  Just copy it (and
# clean_data.pl) into the directory with your training data
# and issue as is.

FILES=`ls ./file[01]`

for file in $FILES; do
  echo "Processing $file"
  cat $file | ./clean_data.pl > $file.tmp
  mv $file.tmp $file
done

exit 0
