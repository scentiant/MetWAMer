# Michael E Sparks (mespar1@gmail.com)
# Last modified : 20 July 2013

# Copyright (c) 2004,2007,2010,2013 Michael E Sparks
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

BIN=../bin/
CC=gcc
# THREEPERIODIC is defined to invoke the proper definitions
#   from the IMMpractical Markov chain model library.
#
# The following parameter settings work cleanly on GCC 4.8.1,
# but a less stringent combination is used for distribution.
CFLAGS=-c -std=c99 -pedantic-errors -O2 -Wall -Werror -I../include \
-I../IMMpractical/include/ \
-DTHREEPERIODIC
LINK=-L../IMMpractical/lib/ -limmpractical -lm -lc
LIBXMLCFLAGS=`xml2-config --cflags`
LIBXMLLIBS=`xml2-config --libs`

executables = \
  build_featvecs \
  calc_medoids \
  eval_MetWAM \
  indexFasSeq \
  MetWAMer.gthXML \
  MetWAMer.CDS \
  parse_tigrxml_cds \
  parse_pps_codseqs \
  print_MetWAM \
  train_MC \
  train_MetPWM \
  train_MetWAM \
  train_perceptron

build_featvecs_objs = \
  build_featvecs.o \
  classifiers.o \
  index_utils.o \
  MetWAM_utils.o \
  sequence_parse.o
calc_medoids_objs = \
  calc_medoids.o \
  classifiers.o \
  cluster.o \
  index_utils.o \
  MetWAM_utils.o \
  sequence_parse.o
eval_MetWAM_objs = \
  classifiers.o \
  eval_MetWAM.o \
  index_utils.o \
  MetWAM_utils.o \
  sequence_parse.o
indexFasSeq_objs = \
  indexFasSeq.o \
  index_utils.o \
  sequence_parse.o
MetWAMer_gthXML_objs = \
  classifiers.o \
  index_utils.o \
  MetWAM_utils.o \
  MetWAMer.gthXML.o \
  sequence_parse.o
MetWAMer_mRNA_objs = \
  classifiers.o \
  index_utils.o \
  MetWAM_utils.o \
  MetWAMer.CDS.o \
  sequence_parse.o
parse_tigrxml_cds_objs = \
  classifiers.o \
  index_utils.o \
  MetWAM_utils.o \
  parse_tigrxml_cds.o \
  sequence_parse.o
parse_pps_codseqs_objs = \
  classifiers.o \
  index_utils.o \
  MetWAM_utils.o \
  parse_pps_codseqs.o \
  sequence_parse.o
print_MetWAM_objs = \
  print_MetWAM.o
train_MetPWM_objs = \
  sequence_parse.o \
  train_MetPWM.o
train_MetWAM_objs = \
  sequence_parse.o \
  train_MetWAM.o
train_perceptron_objs = \
  classifiers.o \
  train_perceptron.o

ALL : $(executables)
	mv $(executables) $(BIN)

build_featvecs : $(build_featvecs_objs)
	$(CC) $(build_featvecs_objs) -o $@ $(LINK) $(LIBXMLLIBS)
calc_medoids : $(calc_medoids_objs)
	$(CC) $(calc_medoids_objs) -o $@ $(LINK) $(LIBXMLLIBS)
eval_MetWAM : $(eval_MetWAM_objs)
	$(CC) $(eval_MetWAM_objs) -o $@ $(LINK) $(LIBXMLLIBS)
indexFasSeq : $(indexFasSeq_objs)
	$(CC) $(indexFasSeq_objs) -o $@
MetWAMer.gthXML : $(MetWAMer_gthXML_objs)
	$(CC) $(MetWAMer_gthXML_objs) -o $@ $(LINK) $(LIBXMLLIBS)
MetWAMer.CDS : $(MetWAMer_mRNA_objs)
	$(CC) $(MetWAMer_mRNA_objs) -o $@ $(LINK) $(LIBXMLLIBS)
parse_tigrxml_cds : $(parse_tigrxml_cds_objs)
	$(CC) $(parse_tigrxml_cds_objs) -o $@ $(LINK) $(LIBXMLLIBS)
parse_pps_codseqs : $(parse_pps_codseqs_objs)
	$(CC) $(parse_pps_codseqs_objs) -o $@ $(LINK) $(LIBXMLLIBS)
print_MetWAM : $(print_MetWAM_objs)
	$(CC) $(print_MetWAM_objs) -o $@
train_MC : train_MC.o
	$(CC) train_MC.o -o $@ $(LINK)
train_MetPWM : $(train_MetPWM_objs)
	$(CC) $(train_MetPWM_objs) -o $@
train_MetWAM : $(train_MetWAM_objs)
	$(CC) $(train_MetWAM_objs) -o $@
train_perceptron : $(train_perceptron_objs)
	$(CC) $(train_perceptron_objs) -o $@ $(LINK)

build_featvecs.o : build_featvecs.c
	$(CC) $(CFLAGS) $(LIBXMLCFLAGS) build_featvecs.c -o $@
calc_medoids.o : calc_medoids.c
	$(CC) $(CFLAGS) $(LIBXMLCFLAGS) -DSHOWCLUSTERING calc_medoids.c -o $@
classifiers.o : classifiers.c
	$(CC) $(CFLAGS) classifiers.c -o $@
cluster.o : cluster.c
	$(CC) $(CFLAGS) cluster.c -o $@
eval_MetWAM.o : eval_MetWAM.c
	$(CC) $(CFLAGS) $(LIBXMLCFLAGS) eval_MetWAM.c -o $@
indexFasSeq.o : indexFasSeq.c
	$(CC) $(CFLAGS) indexFasSeq.c -o $@
index_utils.o : index_utils.c
	$(CC) $(CFLAGS) index_utils.c -o $@
MetWAMer.gthXML.o : MetWAMer.gthXML.c
	$(CC) $(CFLAGS) $(LIBXMLCFLAGS) MetWAMer.gthXML.c -o $@
MetWAMer.CDS.o : MetWAMer.CDS.c
	$(CC) $(CFLAGS) $(LIBXMLCFLAGS) MetWAMer.CDS.c -o $@
MetWAM_utils.o : MetWAM_utils.c
	$(CC) $(CFLAGS) $(LIBXMLCFLAGS) MetWAM_utils.c -o $@
parse_tigrxml_cds.o : parse_tigrxml_cds.c
	$(CC) $(CFLAGS) $(LIBXMLCFLAGS) parse_tigrxml_cds.c -o $@
parse_pps_codseqs.o : parse_pps_codseqs.c
	$(CC) $(CFLAGS) $(LIBXMLCFLAGS) parse_pps_codseqs.c -o $@
print_MetWAM.o : print_MetWAM.c
	$(CC) $(CFLAGS) $(LIBXMLCFLAGS) print_MetWAM.c -o $@
sequence_parse.o : sequence_parse.c
	$(CC) $(CFLAGS) sequence_parse.c -o $@
train_MC.o : train_MC.c
	$(CC) $(CFLAGS) train_MC.c -o $@
train_MetPWM.o : train_MetPWM.c
	$(CC) $(CFLAGS) $(LIBXMLCFLAGS) train_MetPWM.c -o $@
train_MetWAM.o : train_MetWAM.c
	$(CC) $(CFLAGS) $(LIBXMLCFLAGS) train_MetWAM.c -o $@
train_perceptron.o : train_perceptron.c
	$(CC) $(CFLAGS) -DTRACKERROR train_perceptron.c -o $@

clean :
	rm -f *.o
cleaner : clean
	rm -f $(BIN)/*
redo : cleaner
	make
