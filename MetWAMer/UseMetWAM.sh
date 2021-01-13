#!/usr/bin/env bash

# UseMetWAM.sh - Script to demonstrate use of the MetWAMer.gthXML
#                program on gthXML and PASIF xml documents.

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
PREREQS="$BIN/indexFasSeq \
         $BIN/MetWAMer.gthXML"
DATADIR="data/"
PARMDIR="$DATADIR/Athal/"
RNGDIR="specs/"

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

echo -en "\nAbout to demo the MetWAMer program ( press Return to proceed ) : "
read response
./$BIN/indexFasSeq $DATADIR/rice.fas
./$BIN/MetWAMer.gthXML \
  -g $DATADIR/rice.fas.ind \
  -f $DATADIR/rice.gthxml \
  -r $RNGDIR/GenomeThreader.rng \
  -a 4 \
  -k 3 \
  -s \
  -d 1 \
  -l $PARMDIR/meds.medoids.xml \
  -w $PARMDIR/parms.MetPWM \
  -m $PARMDIR/parmsT.MetWAM \
  -x $PARMDIR/parmsF.MetWAM \
  -t 1 \
  -n $PARMDIR/neur.CHI2.sigunit \
  -c 3 $PARMDIR/MCparms.CHI2 \
  > rice.plus_MetWAMer.gthxml
while true; do
  echo -en "\nView results? (y/n) : "
  read response

  case "$response" in
  "y") less rice.plus_MetWAMer.gthxml
       break
       ;;
  "n") break
       ;;
  *)   echo -en "\a"
       continue
  esac
done
# Validate results, if jing (www.thaiopensource.com/relaxng/jing.html)
# is installed on the user's system.
if [[ -e `which jing 2> /dev/null` ]]; then
  while true; do
    echo -en "\nValidate results with jing? (y/n) : "
    read response
  
    case "$response" in
    "y") jing $RNGDIR/GenomeThreader.rng rice.plus_MetWAMer.gthxml
         break
         ;;
    "n") break
         ;;
    *)   echo -en "\a"
         continue
    esac
  done
fi

# Demo MetWAMer on PASIF output, if interested
while true; do
  echo -en "\nDemo MetWAMer on PASIF output, also? (y/n) : "
  read response

  case "$response" in
  "y") ./$BIN/MetWAMer.gthXML \
         -p \
         -a 1 \
         -m $PARMDIR/parmsT.MetWAM \
         -x $PARMDIR/parmsF.MetWAM \
         -f $DATADIR/sample.PASIFxml \
         -g $DATADIR/sample.fas \
         -r $RNGDIR/PASIF.rng > \
         sample.plus_MetWAMer.PASIFxml
       while true; do
         echo -en "\nView results? (y/n) : "
         read response
       
         case "$response" in
         "y") less sample.plus_MetWAMer.PASIFxml
              break
              ;;
         "n") break
              ;;
         *)   echo -en "\a"
              continue
         esac
       done
       # Validate results, if jing (www.thaiopensource.com/relaxng/jing.html)
       # is installed on the user's system.
       if [[ -e `which jing 2> /dev/null` ]]; then
         while true; do
           echo -en "\nValidate results with jing? (y/n) : "
           read response
         
           case "$response" in
           "y") jing $RNGDIR/PASIF.rng sample.plus_MetWAMer.PASIFxml
                break
                ;;
           "n") break
                ;;
           *)   echo -en "\a"
                continue
           esac
         done
       fi
       break
       ;;
  "n") echo -en "\n"
       break
       ;;
  *)   echo -en "\a"
       continue
  esac
done

# Cleanup lint
rm rice.plus_MetWAMer.gthxml $DATADIR/rice.fas.ind
if [ -e sample.plus_MetWAMer.PASIFxml ]; then
  rm -f sample.plus_MetWAMer.PASIFxml
fi

exit 0
