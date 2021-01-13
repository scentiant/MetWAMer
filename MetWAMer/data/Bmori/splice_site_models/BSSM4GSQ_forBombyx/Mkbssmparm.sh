#!/usr/bin/env bash

# Mkbssmparm.sh - Manages weight array matrix development
#                 given appropriate training data.

# Michael E Sparks (mespar1@iastate.edu)
# Last modified: 8 August 2004

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

EXECPATH="bin/"
PREREQS="${EXECPATH}BSSM_build ${EXECPATH}BSSM_print"
DATAROOT="${PWD}/output/training_data/"
acceptor="AG"

# verify that necessary executable is present
for prereq in $PREREQS; do
  if [ ! -e $prereq ]; then
    cd src; make; cd -
  fi
done

echo -en "\nPlease name your BSSM parm file, e.g., \"foo.bssm\" : "
read response
PARMFILE=$response

echo -e "\nExpected splice site training data root directory is :"
echo "$DATAROOT"
while true; do
  echo -n "Is this correct? (y/n) : "
  read response

  case "$response" in
  "y") break
       ;;
  "n") echo "Please indicate the root directory"
       echo -n "  of your splice site training data: "
       read response
       DATAROOT=$response
       break
       ;;
  *)
       continue
  esac
done

if [ ! -e $DATAROOT ]; then
  echo "Can't find that directory!"
  exit 0
fi

DONOR="GT GC"
for donor in $DONOR; do
  # Does the user wish to train for a given terminus type?
  while true; do
    echo -en "\nDo you want to train the $donor model? (y/n) : "
    read response

    case "$response" in
    "n") echo "  Ignoring the $donor model..."
         break
         ;;
    "y") response="y"
         DONINPUT="$DATAROOT/${donor}_${acceptor}/don/"
         ACCINPUT="$DATAROOT/${donor}_${acceptor}/acc/"

         echo "  Training the $donor model..."

         echo "Processing donor data"
         FILES="T1_don T2_don T0_don F1_don F2_don F0_don Fi_don"
         for file in $FILES; do
           echo "  $file"
           length=`wc -l $DONINPUT/$file | awk '{print $1}'`
           let "length /= 2"
           $EXECPATH/BSSM_build $PARMFILE $length $donor $DONINPUT/$file
         done

         echo "Processing acceptor data"
         FILES="T1_acc T2_acc T0_acc F1_acc F2_acc F0_acc Fi_acc"
         for file in $FILES; do
           echo "  $file"
           length=`wc -l $ACCINPUT/$file | awk '{print $1}'`
           let "length /= 2"
           $EXECPATH/BSSM_build $PARMFILE $length $donor $ACCINPUT/$file
         done

         break
         ;;
    *)
         continue
         ;;
    esac
  done
done

echo -en "\nPlease name your BSSM ascii output file, e.g., \"foo.bssm.ascii\" : "
read response
ASCIIFILE=$response
$EXECPATH/BSSM_print $PARMFILE > $ASCIIFILE

echo "Finished"

exit 0
