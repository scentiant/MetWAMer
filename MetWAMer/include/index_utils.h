/* index_utils.h
 * Michael E Sparks (mespar1@gmail.com)
 * Last modified: 17 August 2007
 *
 * Functions used by indexFasSeq and related utilities.
 *
 * Copyright (c) 2007 Michael E Sparks
 * All rights reserved.
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#ifndef INDEX_UTILS_H
#define INDEX_UTILS_H

#include <ctype.h>
#include "sequence_parse.h"

/* Standard genetic code.  This lookup table should not be initialized as *
 * a table, because we don't necessarily know which of [0-4] Ade, Cyt,    *
 * Gua, and Thy map to: Call the INIT_STD_GEN_CODE macro, instead.        */
char stdgencode[NTALFSIZE][NTALFSIZE][NTALFSIZE];
#define INIT_STD_GEN_CODE { \
  stdgencode[Thy][Thy][Thy]='F'; \
  stdgencode[Thy][Thy][Cyt]='F'; \
  stdgencode[Thy][Thy][Ade]='L'; \
  stdgencode[Thy][Thy][Gua]='L'; \
  stdgencode[Cyt][Thy][Thy]='L'; \
  stdgencode[Cyt][Thy][Cyt]='L'; \
  stdgencode[Cyt][Thy][Ade]='L'; \
  stdgencode[Cyt][Thy][Gua]='L'; \
  stdgencode[Ade][Thy][Thy]='I'; \
  stdgencode[Ade][Thy][Cyt]='I'; \
  stdgencode[Ade][Thy][Ade]='I'; \
  stdgencode[Ade][Thy][Gua]='M'; \
  stdgencode[Gua][Thy][Thy]='V'; \
  stdgencode[Gua][Thy][Cyt]='V'; \
  stdgencode[Gua][Thy][Ade]='V'; \
  stdgencode[Gua][Thy][Gua]='V'; \
  stdgencode[Thy][Cyt][Thy]='S'; \
  stdgencode[Thy][Cyt][Cyt]='S'; \
  stdgencode[Thy][Cyt][Ade]='S'; \
  stdgencode[Thy][Cyt][Gua]='S'; \
  stdgencode[Cyt][Cyt][Thy]='P'; \
  stdgencode[Cyt][Cyt][Cyt]='P'; \
  stdgencode[Cyt][Cyt][Ade]='P'; \
  stdgencode[Cyt][Cyt][Gua]='P'; \
  stdgencode[Ade][Cyt][Thy]='T'; \
  stdgencode[Ade][Cyt][Cyt]='T'; \
  stdgencode[Ade][Cyt][Ade]='T'; \
  stdgencode[Ade][Cyt][Gua]='T'; \
  stdgencode[Gua][Cyt][Thy]='A'; \
  stdgencode[Gua][Cyt][Cyt]='A'; \
  stdgencode[Gua][Cyt][Ade]='A'; \
  stdgencode[Gua][Cyt][Gua]='A'; \
  stdgencode[Thy][Ade][Thy]='Y'; \
  stdgencode[Thy][Ade][Cyt]='Y'; \
  stdgencode[Thy][Ade][Ade]='-'; \
  stdgencode[Thy][Ade][Gua]='-'; \
  stdgencode[Cyt][Ade][Thy]='H'; \
  stdgencode[Cyt][Ade][Cyt]='H'; \
  stdgencode[Cyt][Ade][Ade]='Q'; \
  stdgencode[Cyt][Ade][Gua]='Q'; \
  stdgencode[Ade][Ade][Thy]='N'; \
  stdgencode[Ade][Ade][Cyt]='N'; \
  stdgencode[Ade][Ade][Ade]='K'; \
  stdgencode[Ade][Ade][Gua]='K'; \
  stdgencode[Gua][Ade][Thy]='D'; \
  stdgencode[Gua][Ade][Cyt]='D'; \
  stdgencode[Gua][Ade][Ade]='E'; \
  stdgencode[Gua][Ade][Gua]='E'; \
  stdgencode[Thy][Gua][Thy]='C'; \
  stdgencode[Thy][Gua][Cyt]='C'; \
  stdgencode[Thy][Gua][Ade]='-'; \
  stdgencode[Thy][Gua][Gua]='W'; \
  stdgencode[Cyt][Gua][Thy]='R'; \
  stdgencode[Cyt][Gua][Cyt]='R'; \
  stdgencode[Cyt][Gua][Ade]='R'; \
  stdgencode[Cyt][Gua][Gua]='R'; \
  stdgencode[Ade][Gua][Thy]='S'; \
  stdgencode[Ade][Gua][Cyt]='S'; \
  stdgencode[Ade][Gua][Ade]='R'; \
  stdgencode[Ade][Gua][Gua]='R'; \
  stdgencode[Gua][Gua][Thy]='G'; \
  stdgencode[Gua][Gua][Cyt]='G'; \
  stdgencode[Gua][Gua][Ade]='G'; \
  stdgencode[Gua][Gua][Gua]='G'; \
}

/* Macro to write the complement of base nt into ntcomp. *
 * Both arguments should be of type char.                */
#define TAKE_COMPLEMENT(nt,ntcomp) { \
  char c=(char)toupper((int)nt); \
  switch(c) { \
    case 'A' : \
      ntcomp='T'; \
      break; \
    case 'C' : \
      ntcomp='G'; \
      break; \
    case 'G' : \
      ntcomp='C'; \
      break; \
    case 'T' : \
      ntcomp='A'; \
      break; \
    default : /* Leave ambiguous, non-nucleotide symbol as-is */ \
      ntcomp=c; \
  } \
}

/* Function to import an indexFasSeq'ed sequence from the filesystem */
char *slurpindex(
  const char *infile,
  char *seq
);

/* Function to efficiently parse substrings *
 * of an index in a one-shot manner.  start *
 * and stop should be supplied using zero-  *
 * based coordinates.                       */
char *oneshotparse(
  const char *infile,
  char *seq,
  int start,
  int stop
);

/* Function to parse a sequence fragment from hostseq, *
 * using a ZERO-BASED coordinate system.  While for a  *
 * one-shot fragment parse, an fseek() operation would *
 * be more efficient, this code was designed with      *
 * repeated extractions in mind, for which importing   *
 * the data into memory is quicker.                    */
char *parseseqfrag(
  char *hostseq,
  int start,
  int stop
);

/* Returns a char string of the translated cds */
char *translate_cds(
  char *cds,
  char gencode[NTALFSIZE][NTALFSIZE][NTALFSIZE]
);

#endif
