# Makefile for train_MM.x and use_MM.x and

# Michael E Sparks (mespar1@gmail.com)

# Copyright (C) 2005,2006,2013  Michael E Sparks
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

# Relevant directories ########################################################

BIN=../bin
DATADIR=../training_data
INCLUDE=../include
LIB=../lib
PARMDIR=../parm_files
SRC=.

# C compiler of choice ########################################################

CC=gcc

# Compiler flag descriptions ##################################################

#These are set to seemingly reasonable defaults.
# ADJUSTNULLS If any particular prob mass function
#   contains a null probability (can happen if a
#   history is observed in the training data, but
#   that history + some 3'-terminal nucleotide is
#   not), this code will adjust the distribution
#   based on the simple smoothing rule described
#   in Sparks and Brendel (2005) Bioinformatics
#   21Suppl3:iii20-iii30.
# BUCKBOUNDCHECK reports bucket boundaries as
#   determined in the history tying phase of
#   training.
# DMMCREPORT reports the PMF deduced by the DMMM
#   for a given test substrate instance.
# DMPARMPRINT prints the min, max, mu-sigma,
#   and mu+sigma results for each oligomer
#   considered in the DMMM training process
# EMPTYEQUIV will alternatively compile code
#   for propagating pseudocounts in response
#   to non-occurring histories, an experimental
#   alternative on Jelinek and Salzberg's suggestions
#   for handling non-occurring oligomers.
# PROBREPORT reports PMF status of final, smoothed
#   probability estimates for each (occurring) history.
#   (Note that, if a history does not occur and
#   EMPTYEQUIV is not defined, then the PMF that will
#   be ascertained will correspond to that of the
#   largest of the smaller pretexts/histories that
#   did occur, which will thus be reported.)
# QUANTSPEC compiles a library version that will
#   calculate quantiles based on G+C composition
#   of overlapping windows in the training data.
#   Quantile-specific transition probability estimates
#   will be produced for each training method by
#   binning oligomer counts based on window content.
#   Both THREEPERIODIC and SHADOWSTRAND and, optionally,
#   WINDOWGC, must be defined if QUANTSPEC is.
# ROOTREPORT will compile in code that reports
#   max likelihood estimates of lambdahat values
#   as the training progresses.
# SHADOWSTRAND develops a seven-class model for
#   the six coding phases and noncoding sequences.
#   Defining this option REQUIRES for THREEPERIODIC
#   to be set, also (see below)!
# SIMPLEXDBG prints to stdout progress of the simplex
#   algorithm when solving some instance of an LP problem.
# STRICTTEST mandates that there be no ambiguous
#   symbols encountered in an input test sequence
#   to compute the smoothed log-probability of.
#   If not set, then a residue is randomly chosen
#   in place of the ambiguous symbol, and that is
#   instead used in computing the estimate.
# STRICTTRAIN mandates that no ambiguous symbols
#   are permitted in the training data.
# TESTWINDOWGC compiles a QUANTSPEC-versioned
#   library that, for a given test sequence, will
#   compute the G+C composition of a window of
#   2*FLANK+MAXORDER+1 bases, and assay the
#   MAXORDER+1-mer centered in it with the appropriate
#   quantile-specific parameters. If not defined,
#   these will be selected on the basis of the
#   overall G+C composition of the test sequence.
#   If the test sequence is < 2*FLANK+MAXORDER+1 bases
#   long, overall G+C composition will be used.
# THREEPERIODIC compiles the library to accommodate
#   models taking the periodicity of nucleic acid
#   sequences into account.
CFLAGS=-c -std=c99 -pedantic-errors -Wall -Werror -O2 -I${INCLUDE} \
       -DADJUSTNULLS \
       -DBUCKBOUNDCHECK \
       -DROOTREPORT \
       -DSTRICTTRAIN \
       -DPROBREPORT \
       -DTHREEPERIODIC \
       #-DSHADOWSTRAND \
       #-DQUANTSPEC \
       #-DTESTWINDOWGC \
       #-DDMMCREPORT \
       #-DDMPARMPRINT \
       #-DEMPTYEQUIV \
       #-DSIMPLEXDBG \
       #-DSTRICTTEST \

# Default target ##############################################################

DEFAULT : ${BIN}/train_MM.x \
          ${BIN}/use_MM.x

# Objects #####################################################################

STDOBJ=immpractical.o \
       dimm_utils.o \
       chisquare_utils.o \
       fo_utils.o \
       dm_utils.o \
       sequence_parse.o \
       simplex.o \
       sorting.o

immpractical.o : ${SRC}/immpractical.c \
                 ${INCLUDE}/immpractical.h \
                 ${INCLUDE}/dimm_utils.h \
                 ${INCLUDE}/chisquare_utils.h \
                 ${INCLUDE}/fo_utils.h \
                 ${INCLUDE}/dm_utils.h \
                 ${INCLUDE}/sequence_parse.h
	${CC} ${CFLAGS} ${SRC}/immpractical.c -o $@
#	${CC} ${CFLAGS} -fPIC ${SRC}/immpractical.c -o $@

fo_utils.o : ${SRC}/fo_utils.c \
             ${INCLUDE}/fo_utils.h \
             ${INCLUDE}/immpractical.h \
             ${INCLUDE}/sequence_parse.h
	${CC} ${CFLAGS} ${SRC}/fo_utils.c -o $@
#	${CC} ${CFLAGS} -fPIC ${SRC}/fo_utils.c -o $@

dm_utils.o : ${SRC}/dm_utils.c \
             ${INCLUDE}/fo_utils.h \
             ${INCLUDE}/dm_utils.h \
             ${INCLUDE}/immpractical.h \
             ${INCLUDE}/sequence_parse.h \
             ${INCLUDE}/simplex.h 
	${CC} ${CFLAGS} ${SRC}/dm_utils.c -o $@
#	${CC} ${CFLAGS} -fPIC ${SRC}/dm_utils.c -o $@

chisquare_utils.o : ${SRC}/chisquare_utils.c \
                    ${INCLUDE}/chisquare_utils.h \
                    ${INCLUDE}/immpractical.h
	${CC} ${CFLAGS} ${SRC}/chisquare_utils.c -o $@
#	${CC} ${CFLAGS} -fPIC ${SRC}/chisquare_utils.c -o $@

dimm_utils.o : ${SRC}/dimm_utils.c \
               ${INCLUDE}/dimm_utils.h \
               ${INCLUDE}/immpractical.h \
               ${INCLUDE}/sorting.h
	${CC} ${CFLAGS} ${SRC}/dimm_utils.c -o $@
#	${CC} ${CFLAGS} -fPIC ${SRC}/dimm_utils.c -o $@

sequence_parse.o : ${SRC}/sequence_parse.c \
                   ${INCLUDE}/sequence_parse.h \
                   ${INCLUDE}/immpractical.h
	${CC} ${CFLAGS} ${SRC}/sequence_parse.c -o $@
#	${CC} ${CFLAGS} -fPIC ${SRC}/sequence_parse.c -o $@

sorting.o : ${SRC}/sorting.c \
            ${INCLUDE}/sorting.h \
            ${INCLUDE}/immpractical.h
	${CC} ${CFLAGS} ${SRC}/sorting.c -o $@
#	${CC} ${CFLAGS} -fPIC ${SRC}/sorting.c -o $@

simplex.o : ${SRC}/simplex.c \
            ${INCLUDE}/simplex.h
	${CC} ${CFLAGS} ${SRC}/simplex.c -o $@
#	${CC} ${CFLAGS} -fPIC ${SRC}/simplex.c -o $@

train_MM.o : ${SRC}/train_MM.c \
              ${INCLUDE}/immpractical.h
	${CC} ${CFLAGS} ${SRC}/train_MM.c -o $@

use_MM.o : ${SRC}/use_MM.c \
            ${INCLUDE}/immpractical.h
	${CC} ${CFLAGS} ${SRC}/use_MM.c -o $@

# Library #####################################################################

# Libraries for driver programs to link with
LINK=-L${LIB} -limmpractical -lm

# commands for building static library
STDLIB=${LIB}/libimmpractical.a
${LIB}/libimmpractical.a : ${STDOBJ}
	ar -rc $@ ${STDOBJ}
	ranlib $@

# commands for building shared object library
#STDLIB=${LIB}/libimmpractical.so
#NOTE THAT THIS NOMENCLATURE ASSUMES YOU'RE RUNNING A BASH INTERPRETER!
#${LIB}/libimmpractical.so : ${STDOBJ}
#	ld -shared -soname ${LIB}/libimmpractical.so.1 \
#          -o ${LIB}/libimmpractical.so.1.0 -lc ${STDOBJ}
#	/sbin/ldconfig -v -n ${LIB}
#	ln -sf ${LIB}/libimmpractical.so.1 ${LIB}/libimmpractical.so
#	export LD_LIBRARY_PATH=.:${LD_LIBRARY_PATH}

# Driver programs #############################################################

${BIN}/train_MM.x : ${STDLIB} train_MM.o
	${CC} train_MM.o -o $@ ${LINK}
${BIN}/use_MM.x : ${STDLIB} use_MM.o
	${CC} use_MM.o -o $@ ${LINK}

# Training ####################################################################

#NOTE THAT THIS NOMENCLATURE ASSUMES YOU'RE RUNNING A BASH INTERPRETER!
PAGER=2>&1 | less

# Various training algorithms are described in immpractical.h (see enum type)
ALGOCODE=1

#parameter files to build
TRAIN=Arabidopsis

Arabidopsis : ${BIN}/train_MM.x
	${BIN}/train_MM.x ${ALGOCODE} \
                           ${DATADIR}/Arabidopsis/D_data/file0.D \
                           ${DATADIR}/Arabidopsis/D_data/file1.D \
                           ${DATADIR}/Arabidopsis/H_data/file0.H \
                           ${DATADIR}/Arabidopsis/H_data/file1.H \
                           ${PARMDIR}/Arabidopsis ${PAGER}

# Testing #####################################################################

#demo for testing parameterization
test : ${BIN}/use_MM.x
	cat ${DATADIR}/Arabidopsis/T_data/file0.T | \
          ${BIN}/use_MM.x ${ALGOCODE} ${PARMDIR}/Arabidopsis.* ${PAGER}

# Additional targets ##########################################################
all :
	make clean; make; make train test
clean :
	rm -f ${SRC}/*.o ${LIB}/*.a ${BIN}/*.x ${PARMDIR}/Arabidopsis.*
dataunzip :
	find ${DATADIR} -name "file*" | xargs bunzip2
datazip :
	find ${DATADIR} -name "file*" | xargs bzip2
default : ${DEFAULT}
redo :
	make clean ; make
train : ${TRAIN}
