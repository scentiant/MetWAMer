#!/usr/bin/env bash

# Mktraindata.sh - Script to create BSSM training data

# Michael E Sparks (michael.sparks2@usda.gov)
# Recent modifications:
# 3 Sep 2021 (made key on which to sort explicit)
# 14 Dec 2020 (made a few idiosyncratic adjustments
#   to accommodate Bombyx mori training data)

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

ROOT=$PWD
BIN="bin/"
PLBIN="src/plscripts/"
NONREDUNDANT="src/plscripts/nonredundant_tools/"
PREREQS="${BIN}indexFasSeq"
OUT="output/training_data/"
EIOUT="output/exons_introns/"
TEMPLATES="input/fas/"
THREADINGS="input/gsq/"

# If your THREADINGS files are in GeneSeqer's plain text
# format, change the FORMAT variable to "plain".
# If you are using GenomeThreader -xmlout output or
# GSQ2XML.pl-converted GeneSeqer output
# (see http://www.genomethreader.org),
# set the following variable to "gthxml".
#FORMAT="plain"
FORMAT="gthxml"

CUTOFF=1.0
GOODEXONCT=3

DON="GT"
#DON="GC"
ACC="AG"

echo -e "\nRunning script ... `date`\n"
start_time=`date +%s`

# Phase 1 -------

# verify that necessary executables are present
for prereq in $PREREQS; do
  if [ ! -e $prereq ]; then
    cd src; make; cd -
  fi
done

# clear out any pre-existing result data
if [ -e $EIOUT/${DON}_${ACC} ]; then
  rm -rf $EIOUT/${DON}_${ACC}
fi
if [ -e $OUT/${DON}_${ACC} ]; then
  rm -rf $OUT/${DON}_${ACC}
fi

# process each GeneSeqer output file in turn
for file in `ls ${THREADINGS}`; do
  handle=`echo "$file" | cut -f 1 -d'.'`
  echo -n "Processing $handle ... "
  $PLBIN/get_support_pps.pl $FORMAT $CUTOFF $GOODEXONCT \
    $THREADINGS/$file > $handle.tmp

  # build the models
  if [ -s $handle.tmp ]; then
    echo -n "Useful stuff in $handle, getting it ... "

    $BIN/indexFasSeq $TEMPLATES/$handle.fas

    $PLBIN/parse_exons-and-introns.pl $handle.tmp \
      $TEMPLATES/$handle.fas.ind $BIN \
      $OUT/${handle}_0  $OUT/${handle}_1  $OUT/${handle}_2 \
      $OUT/${handle}_I0 $OUT/${handle}_I1 $OUT/${handle}_I2

    echo "done"
  else
    echo "Nothing useful in $handle"
  fi
  rm -f $handle.tmp
done

# concatenate the model files and eliminate redundant/empty entries
cat $OUT/*_I[012] | $NONREDUNDANT/fasta2convertable.pl | \
  sed '/^ \*\*\*/d' | sort -k1,1 | $NONREDUNDANT/squash.pl   | \
  $NONREDUNDANT/convert2fasta.pl | \
  ./silkworm_adj.pl > $EIOUT/fileI
HYPOS="0 1 2 I0 I1 I2"
for i in $HYPOS; do
  cat $OUT/*_${i} | $NONREDUNDANT/fasta2convertable.pl | \
    sed '/^ \*\*\*/d' | sort -k1,1 | $NONREDUNDANT/squash.pl | \
    $NONREDUNDANT/convert2fasta.pl | \
    ./silkworm_adj.pl > $EIOUT/file${i}
  rm -f $OUT/*_${i}
done

# Phase 2 -------

# Collect true splice sites
echo -n "Processing true intron sites ... "
i=0
while [ "$i" -le 2 ]; do
  $PLBIN/find_true_donacc.pl $EIOUT/fileI$i $BIN $TEMPLATES $DON $ACC
  mv $EIOUT/fileI$i.* $OUT
  let "i+=1"
done
echo "done"

# Collect false within-intron sites
echo -n "Processing false intron sites ... "
$PLBIN/find_false.pl -i $EIOUT/fileI $BIN $TEMPLATES $DON $ACC $OUT 2> /dev/null
echo "done"

# Produce false within-exon sites
echo -n "Processing false exon sites ... "
i=0
while [ "$i" -le 2 ]; do
  $PLBIN/find_false.pl -e $EIOUT/file${i} $BIN $TEMPLATES $DON $ACC $OUT 2> /dev/null
  let "i+=1"
done
echo "done"

i=0
while [ "$i" -le 2 ]; do
  cat $OUT/file[012]*.falsedons$i > exon.${DON}_${ACC}.falsedons$i
  rm -f $OUT/file[012]*.falsedons$i
  mv exon.${DON}_${ACC}.falsedons$i $OUT/exon.${DON}_${ACC}.falsedons$i

  cat $OUT/file[012]*.falseaccs$i > exon.${DON}_${ACC}.falseaccs$i
  rm -f $OUT/file[012]*.falseaccs$i
  mv exon.${DON}_${ACC}.falseaccs$i $OUT/exon.${DON}_${ACC}.falseaccs$i

  let "i+=1"
done

# Clear empty entries
echo -n "Cleaning output ... "
for file in `find $OUT -name "*(dons|accs)*" -type f -print`; do
  cat $OUT/$file | $NONREDUNDANT/fasta2convertable.pl | sed '/^ \*\*\*/d' | \
    $NONREDUNDANT/convert2fasta.pl $PLBIN/clean_data.pl > $OUT/$file.tmp 2> /dev/null
  mv $OUT/$file.tmp $OUT/$file
done
echo "done"

# No more substring parsing is necessary
rm -f $TEMPLATES/*.ind

# Phase 3 -------

# Random sampling for the model
echo -n "Sampling false site data to appropriate sizes ... "
cd $OUT
i=0
while [ "$i" -le 2 ]; do
  length=`wc -l fileI$i.${DON}_${ACC}.truedons | awk '{print $1}'`
  let "length /= 2"
  DONLEN[$i]=$length
  ../../$PLBIN/file_sampler.NO_replacement.pl $length \
    < exon.${DON}_${ACC}.falsedons$i > exon.${DON}_${ACC}.falsedons$i.tmp
  mv exon.${DON}_${ACC}.falsedons$i.tmp exon.${DON}_${ACC}.falsedons$i

  length=`wc -l fileI$i.${DON}_${ACC}.trueaccs | awk '{print $1}'`
  let "length /= 2"
  ACCLEN[$i]=$length
  ../../$PLBIN/file_sampler.NO_replacement.pl $length \
    < exon.${DON}_${ACC}.falseaccs$i > exon.${DON}_${ACC}.falseaccs$i.tmp
  mv exon.${DON}_${ACC}.falseaccs$i.tmp exon.${DON}_${ACC}.falseaccs$i

  let "i+=1"
done

# Determine max size of the 3
if [ ${DONLEN[0]} -ge ${DONLEN[1]} ] && [ ${DONLEN[0]} -ge ${DONLEN[2]} ]; then
  MAX=${DONLEN[0]}
elif [ ${DONLEN[1]} -ge ${DONLEN[0]} ] && [ ${DONLEN[1]} -ge ${DONLEN[2]} ]; then
  MAX=${DONLEN[1]}
else
  MAX=${DONLEN[2]}
fi
../../$PLBIN/file_sampler.NO_replacement.pl $MAX \
  < fileI.${DON}_${ACC}.falsedons > fileI.${DON}_${ACC}.falsedons.tmp
mv fileI.${DON}_${ACC}.falsedons.tmp fileI.${DON}_${ACC}.falsedons

if [ ${ACCLEN[0]} -ge ${ACCLEN[1]} ] && [ ${ACCLEN[0]} -ge ${ACCLEN[2]} ]; then
  MAX=${ACCLEN[0]}
elif [ ${ACCLEN[1]} -ge ${ACCLEN[0]} ] && [ ${ACCLEN[1]} -ge ${ACCLEN[2]} ]; then
  MAX=${ACCLEN[1]}
else
  MAX=${ACCLEN[2]}
fi
../../$PLBIN/file_sampler.NO_replacement.pl $MAX \
  < fileI.${DON}_${ACC}.falseaccs > fileI.${DON}_${ACC}.falseaccs.tmp
mv fileI.${DON}_${ACC}.falseaccs.tmp fileI.${DON}_${ACC}.falseaccs
cd $ROOT
echo "done"

# Adjust phase nomenclature to match Volker's
mkdir $EIOUT/${DON}_${ACC}
mv ${EIOUT}/file*                  ${EIOUT}/${DON}_${ACC}
mv $EIOUT/${DON}_${ACC}/file1      $EIOUT/${DON}_${ACC}/file1.tmp
mv $EIOUT/${DON}_${ACC}/file2      $EIOUT/${DON}_${ACC}/file2.tmp
mv $EIOUT/${DON}_${ACC}/file0      $EIOUT/${DON}_${ACC}/file1
mv $EIOUT/${DON}_${ACC}/file1.tmp  $EIOUT/${DON}_${ACC}/file2
mv $EIOUT/${DON}_${ACC}/file2.tmp  $EIOUT/${DON}_${ACC}/file0
mv $EIOUT/${DON}_${ACC}/fileI1     $EIOUT/${DON}_${ACC}/fileI1.tmp
mv $EIOUT/${DON}_${ACC}/fileI2     $EIOUT/${DON}_${ACC}/fileI2.tmp
mv $EIOUT/${DON}_${ACC}/fileI0     $EIOUT/${DON}_${ACC}/fileI1
mv $EIOUT/${DON}_${ACC}/fileI1.tmp $EIOUT/${DON}_${ACC}/fileI2
mv $EIOUT/${DON}_${ACC}/fileI2.tmp $EIOUT/${DON}_${ACC}/fileI0

mkdir $OUT/${DON}_${ACC} && \
mkdir $OUT/${DON}_${ACC}/don && \
mkdir $OUT/${DON}_${ACC}/acc

mv $OUT/*${DON}_${ACC}*dons* $OUT/${DON}_${ACC}/don/
cd $OUT/${DON}_${ACC}/don/
mv exon.${DON}_${ACC}.falsedons0 F1_don.ambig
mv exon.${DON}_${ACC}.falsedons1 F2_don.ambig
mv exon.${DON}_${ACC}.falsedons2 F0_don.ambig
mv fileI.${DON}_${ACC}.falsedons Fi_don.ambig
mv fileI0.${DON}_${ACC}.truedons T1_don.ambig
mv fileI1.${DON}_${ACC}.truedons T2_don.ambig
mv fileI2.${DON}_${ACC}.truedons T0_don.ambig
../../../../filter_ambiguous.pl < F1_don.ambig 1> F1_don 2>/dev/null
../../../../filter_ambiguous.pl < F2_don.ambig 1> F2_don 2>/dev/null
../../../../filter_ambiguous.pl < F0_don.ambig 1> F0_don 2>/dev/null
../../../../filter_ambiguous.pl < Fi_don.ambig 1> Fi_don 2>/dev/null
../../../../filter_ambiguous.pl < T1_don.ambig 1> T1_don 2>/dev/null
../../../../filter_ambiguous.pl < T2_don.ambig 1> T2_don 2>/dev/null
../../../../filter_ambiguous.pl < T0_don.ambig 1> T0_don 2>/dev/null
cd $ROOT

mv $OUT/*${DON}_${ACC}*accs* $OUT/${DON}_${ACC}/acc/
cd $OUT/${DON}_${ACC}/acc/
mv exon.${DON}_${ACC}.falseaccs0 F1_acc.ambig
mv exon.${DON}_${ACC}.falseaccs1 F2_acc.ambig
mv exon.${DON}_${ACC}.falseaccs2 F0_acc.ambig
mv fileI.${DON}_${ACC}.falseaccs Fi_acc.ambig
mv fileI0.${DON}_${ACC}.trueaccs T1_acc.ambig
mv fileI1.${DON}_${ACC}.trueaccs T2_acc.ambig
mv fileI2.${DON}_${ACC}.trueaccs T0_acc.ambig
../../../../filter_ambiguous.pl < F1_acc.ambig 1> F1_acc 2>/dev/null
../../../../filter_ambiguous.pl < F2_acc.ambig 1> F2_acc 2>/dev/null
../../../../filter_ambiguous.pl < F0_acc.ambig 1> F0_acc 2>/dev/null
../../../../filter_ambiguous.pl < Fi_acc.ambig 1> Fi_acc 2>/dev/null
../../../../filter_ambiguous.pl < T1_acc.ambig 1> T1_acc 2>/dev/null
../../../../filter_ambiguous.pl < T2_acc.ambig 1> T2_acc 2>/dev/null
../../../../filter_ambiguous.pl < T0_acc.ambig 1> T0_acc 2>/dev/null
cd $ROOT

end_time=`date +%s`
et=`echo "$end_time - $start_time" | bc`
let "minutes=et / 60"
echo -e "\nScript finished!  Took ~$minutes minutes\a\n"

exit 0
