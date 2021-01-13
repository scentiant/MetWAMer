#!/bin/bash
# Michael E Sparks (michael.sparks2@usda.gov)
# 2 June 2006, updated 4 January 2021 by MES

# Script to randomly sample into $BUCKETS equally-sized
# pools for ${BUCKETS}-fold cross-validation. To produce
# a 1/10, 1/10, 4/5 split for typical deleted interpolation, 
# set $BUCKETS (below) to 10, run, and concatenate 8 of the
# resultant files.

# Copyright (C) 2005,2006  Michael E Sparks
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

USAGE="
$0 num_buckets lines_per_sample file(s)_to_split

  Example: $0 5 2 file0
           -> Partitions file0, a Fasta file with each
              sequence entirely on one line, into five
              roughly equal-sized buckets
"

# Get parameters
if [ "$#" -lt 3 ]; then
  echo -e "\a$USAGE" ; exit 0
else
  BUCKETS=$1      # How many ways to equally split?
  LINESPERSAM=$2  # How many lines constitute an atomic,
                  # sampleable element.
  file=$3         # File to split
fi 

# Update (4 January 2021): this script was developed at a time when
# I was transitioning from tcsh to bash and still learning the ropes
# of the latter, and when datasets were generally much smaller than
# they are nowadays. Although the code performs correctly, it is
# prohibitively slow. I've decided to reimplement it in the if
# clause below, letting standard linux utilities do the heavy
# lifting rather than the shell interpreter language.
# If GNU awk is not available on your system, the script will fall
# back to the prior version.
if command -v gawk > /dev/null 2>&1; then # start of reimplementation

if [ -n "$(find . -maxdepth 1 -name "${file}.shuffled.*" \
-type f -print0 | xargs --null -r ls)" ]; then
  echo "Please remove or rename ${file}.shuffled.* and restart."
  exit 1
fi

cat "$file" | \
  gawk -v lps=$LINESPERSAM '{printf("%s%s",$0,(NR%lps==0)?"\n":"\0")}' | \
  sed 's/\x0$/\n/' | shuf | \
  split -n r/${BUCKETS} --numeric-suffixes=1 - "${file}.shuffled.tmp."

find . -maxdepth 1 -name "${file}.shuffled.tmp.*" -type f -print0 | \
while IFS= read -r -d '' f; do
  cat "$f" | tr '\0' '\n' \
    > "$file".shuffled.`echo "$f" | rev | cut -f 1 -d '.' | rev`
  rm "$f"
done

else # prior code begins here (bypassed by default as of 4 January 2021)

echo "*WARNING* : GNU awk (gawk) not detected - using old partitioner"

# Don't modify these.
TRUE=1
FALSE=0
FLOOR=0
SEDSCRIPT="sedcom.tmp"
MAXSEDLEN=250

# The following array will assist with quickly determining which elements
# of the input file have not been sampled so as to print only those lines.
# Keyed by primary key (no relation to anything else, just arbitrary),
# values are LINE NUMBERS OF START OF SAMPLED ELEMENTS
declare -a Previous

# Function to heap sort the Previous array
heapsort () {
  if [[ -z $1 ]]; then
    echo "No parameter passed to heapsort!"
    exit 1
  else
    limit=$1
  fi

  for (( i = `expr \( $limit / 2 \) - 1` ; i >= 0 ; --i )); do
    sift $i $limit
  done

  for (( i = `expr $limit - 1` ; i >= 1 ; --i )); do
    temp=${Previous[0]}
    Previous[0]=${Previous[$i]}
    Previous[$i]=$temp
    sift 0 `expr $i - 1` 
  done
}

# Function that heapsort depends on
sift () {
  if [[ -z $2 ]]; then
    echo "Must pass two parameters to sift function!"
    exit 1
  else
   root=$1
   lowlim=$2    
  fi

  done=$FALSE
  while [[ `expr $root \* 2` -le $lowlim && $done -eq $FALSE ]]; do
    if [[ `expr $root \* 2` -eq $lowlim ]]; then
      childmax=`expr $root \* 2`
    elif [[ ${Previous[`expr $root \* 2`]} -gt \
            ${Previous[`expr \( $root \* 2 \) + 1`]} ]]; then
      childmax=`expr $root \* 2`
    else
      childmax=`expr \( $root \* 2 \) + 1`
    fi

    if [[ ${Previous[$root]} -lt ${Previous[$childmax]} ]]; then
      temp=${Previous[$root]}
      Previous[$root]=${Previous[$childmax]}
      Previous[$childmax]=$temp
      root=$childmax
    else
      done=$TRUE
    fi
  done
}

# process file --------------------------------------------------------------

# Determine how many units to sample
totct=`grep -c '^>' $file`
ctN=`expr $totct / $BUCKETS`

# Announce the sample sizes
echo "$file ($totct): $ctN (~1/$BUCKETS of total)"

FILESIZE=`wc -l $file | awk '{print $1}'`
let "test=$FILESIZE % $LINESPERSAM"
if [[ $test -ne 0 ]]; then
  echo -e "\aInconsistency with $file!"
  exit 1
fi

# Array of file output names
declare -a file
# Specify output filenames
i=1
while [[ $i -le $BUCKETS ]]; do
  file[$i]="${file}.shuffled.${i}"
  cat /dev/null > ${file[$i]}
  let "i+=1"
done

# We want to randomly sample from the input file without replacement.
# To track what input items have been previously sampled, we use an
# array whose elements are 0 if the line has not been sampled, 1 if so.
# Keyed from 0 to $FILESIZE / $LINESPERSAM (total number of elements
# in input file, key corresponds .NOT. to line number of input where a
# particular data unit is sampled from, BUT any given unit's serial
# occurrence in the input, i.e., first, second, third, ... , n-th).
# Note that Sampledstat[0] is never used!!
declare -a Sampledstat

# This array will allow us to
# consider only sample-able data points, refreshed every so often
# to reasonably ensure that we are sampling productively, i.e., not
# consider only previously sampled elements as their proportion in
# the urn begins to increase.
# Index (key) is just primary key (no relation, 0..size-1), and value
# is index of sample-able elements from Sampledstat
# Note that the zeroth element IS used, in the C array sense
declare -a Samplefromme

# init all to 0; Note that Sampledstat[0] is never used;
for (( i = 0 ; i <= `expr $FILESIZE / $LINESPERSAM` ; ++i )); do
  Sampledstat[$i]=$FALSE
done
for (( i = 0 ; i < `expr $FILESIZE / $LINESPERSAM` ; ++i )); do
  Samplefromme[$i]=`expr $i + 1`
done
# CURRCEIL gives the upper limit, exclusive, for the Samplefromme array
CURRCEIL=`expr $FILESIZE / $LINESPERSAM`

# Produce the random sample for H
j=1 #indexes first $BUCKETS-1 files
i=0 #tracks how many units have been sampled
randloopct=0 # helpful in determining if we are sampling unproductively
# A seed trick I borrowed from the BASH advanced scripting guide
RANDOM=$(head -1 /dev/urandom | od -N 1 | awk '{ print $2 }')
while [[ $j -le `expr $BUCKETS - 1` ]]; do
  SEDCOM="sed -n \";"
  echo "Now sampling for the ${j} file"

  while [[ $i -lt `expr $ctN \* $j` ]]; do
    # Obtain random number
    random=$RANDOM
    let "random %= $CURRCEIL" # Since sampleable indices from Samplefromme
                              # are < $CURRCEIL
    if [[ $random -lt $FLOOR ]]; then
      random=0
    fi

    # Has this unit been previously sampled?
    if [[ ${Sampledstat[${Samplefromme[$random]}]} -eq $TRUE ]]; then
      # No good, try again
      if [ "$randloopct" -eq 100 ]; then
        CURRCEIL=0
        #Reset Samplefromme to makes things a bit more productive
        for (( k = 1 ; k <= `expr $FILESIZE / $LINESPERSAM` ; ++k )); do
          if [[ ${Sampledstat[$k]} -eq $FALSE ]]; then
            Samplefromme[$CURRCEIL]=$k
            let "CURRCEIL+=1"
          fi
        done
        let "CURRCEIL+=1" # valid indices in Samplefromme < CURRCEIL

        randloopct=0
        RANDOM=$(head -1 /dev/urandom | od -N 1 | awk '{ print $2 }')
      else
        let "randloopct+=1"
      fi
      continue 1
    else
      Sampledstat[${Samplefromme[$random]}]=$TRUE
    fi

    # Spit out the data
    let "START=(${Samplefromme[$random]} * $LINESPERSAM) - \
               ($LINESPERSAM - 1)"
    Previous[$i]=$START
    let "STOP=$START + $LINESPERSAM - 1"

    if [ "$START" -eq "$STOP" ]; then
      echo "Error: start == stop!! ($START)"
      exit 1
    else
      SEDCOM="${SEDCOM} ${START},${STOP} p;"
    fi

    let "i+=1"

    if [[ `expr $i % $MAXSEDLEN` -eq 0 ]]; then
      # sed commands can only grow so long, so this
      # will keep it trimmed back a bit.
      SEDCOM="${SEDCOM} \""
      (echo $SEDCOM > $SEDSCRIPT && \
      chmod 700 $SEDSCRIPT && \
      cat $file | ./$SEDSCRIPT >> ${file[$j]} && \
      rm -f $SEDSCRIPT) || \
      (echo "Error associated with ${SEDSCRIPT}!" ; exit 1)
      SEDCOM="sed -n \";"
    fi
  done
  if [[ `expr $i % $MAXSEDLEN` -ne 0 ]]; then
    # Caught a straggler!
    SEDCOM="${SEDCOM} \""
    (echo $SEDCOM > $SEDSCRIPT && \
    chmod 700 $SEDSCRIPT && \
    cat $file | ./$SEDSCRIPT >> ${file[$j]} && \
    rm -f $SEDSCRIPT) || \
    (echo "Error associated with ${SEDSCRIPT}!" ; exit 1)
  fi
  let "j+=1"
done

echo "Heap sorting after first `expr $BUCKETS - 1` files \
      (`expr $ctN \* \( $BUCKETS - 1 \)` elements)"
heapsort `expr $ctN \* \( $BUCKETS - 1 \)`

# Now populate the last file with entries not in the others
# That is, print those parts of $file not listed in the
# Previous array
echo "Now sampling for the ${BUCKETS} file"

SEDCOM="sed -n \";"
begin=1
for (( i = 0 ; i < `expr $ctN \* \( $BUCKETS - 1 \)` ; ++i )); do
  if [[ ${Previous[$i]} -eq $begin ]]; then
    let "begin+=$LINESPERSAM"
    continue
  fi
  end=`expr ${Previous[$i]} - 1`

  if [ "$begin" -eq "$end" ]; then
    echo "Error: begin == end!! ($begin)"
    exit 1
  else
    SEDCOM="${SEDCOM} ${begin},${end} p;"
  fi

  begin=`expr ${Previous[$i]} + $LINESPERSAM`

  if [[ `expr $i % $MAXSEDLEN` -eq 0 ]]; then
    # sed commands can only grow so long, so this
    # will keep it trimmed back a bit.
    SEDCOM="${SEDCOM} \""
    (echo $SEDCOM > $SEDSCRIPT && \
    chmod 700 $SEDSCRIPT && \
    cat $file | ./$SEDSCRIPT >> ${file[$BUCKETS]} && \
    rm -f $SEDSCRIPT) || \
    (echo "Error associated with ${SEDSCRIPT}!" ; exit 1)
    SEDCOM="sed -n \";"
  fi
done

# Handle end game, if need be.
let "i-=1"
if [[ $begin -lt $FILESIZE ]]; then
  SEDCOM="${SEDCOM} ${begin},${FILESIZE} p; \""
  (echo $SEDCOM > $SEDSCRIPT && \
  chmod 700 $SEDSCRIPT && \
  cat $file | ./$SEDSCRIPT >> ${file[$BUCKETS]} && \
  rm -f $SEDSCRIPT) || \
  (echo "Error associated with ${SEDSCRIPT}!" ; exit 1)
elif [[ `expr $i % $MAXSEDLEN` -ne 0 ]]; then
  SEDCOM="${SEDCOM} \""
  (echo $SEDCOM > $SEDSCRIPT && \
  chmod 700 $SEDSCRIPT && \
  cat $file | ./$SEDSCRIPT >> ${file[$BUCKETS]} && \
  rm -f $SEDSCRIPT) || \
  (echo "Error associated with ${SEDSCRIPT}!" ; exit 1)
fi

# Finished for this input file
echo "Finished with random sampling"
j=1
totct2=0
while [[ $j -le $BUCKETS ]]; do
  tmpct=`grep -c '^>' ${file[$j]}`
  echo "  ${file[$j]} : $tmpct"
  let "totct2+=$tmpct"
  let "j+=1"
done
if [ "$totct2" -ne "$totct" ]; then
  echo "Sed hiccoughed unexpectedly! (Note: That's bad.)"
  exit 1
fi

echo "Done"

fi

exit 0
