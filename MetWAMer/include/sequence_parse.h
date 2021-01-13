/* sequence_parse.h
 * Michael Sparks (mespar1@gmail.com)
 * Last modified : 19 August 2007
 *
 * Code for parsing FASTA formatted sequences.
 *
 * Copyright (c) 2005,2006,2007 Michael E Sparks
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

#ifndef SEQUENCE_PARSE_H
#define SEQUENCE_PARSE_H

#define MAXDESCLINE 16384 /* Max length of FASTA description line */

#define NTALFSIZE 4 /* Cardinality of nucleotide alphabet */
#define Ade 0
#define Cyt 1
#define Gua 2
#define Thy 3

/* Macro to recover the character form of a *
 * base when given its integer translation. */
#define RECOVER_BASE(ntint,ntchar) { \
  switch(ntint) { \
    case Ade : \
      ntchar='A'; \
      break; \
    case Cyt : \
      ntchar='C'; \
      break; \
    case Gua : \
      ntchar='G'; \
      break; \
    case Thy : \
      ntchar='T'; \
      break; \
    default : /* Leave ambiguous, non-nucleotide symbol as-is */ \
      FATALERROR("Err (): int form of base not recognized.")\
  } \
}

/* Returns TRUE if seq points to ATG, else returns FALSE */
#define IS_MET_P(seq,i) \
  (seq[(i)]==Ade&&seq[(i)+1]==Thy&&seq[(i)+2]==Gua?TRUE:FALSE)

/* Returns TRUE if seq points to a stop codon, FALSE otherwise */
#define IS_STOP_P(seq,i) \
  ((seq[(i)]==Thy&&seq[(i)+1]==Ade&&seq[(i)+2]==Ade)||\
   (seq[(i)]==Thy&&seq[(i)+1]==Ade&&seq[(i)+2]==Gua)||\
   (seq[(i)]==Thy&&seq[(i)+1]==Gua&&seq[(i)+2]==Ade)?TRUE:FALSE)

/* Returns a random value in the interval *
 * [0,NTALFSIZE-1]. Requires <time.h>.    */
int random_residue(void);

/* get_fasta:  A general function to parse FASTA entries  *
 * Input :                                                *
 *   A file pointer: This must be connected to a stream   *
 *                   containing Fasta-formatted sequences *
 *   A char pointer: This must be set to NULL.  The       *
 *                   function will attempt to allocate    *
 *                   sufficient memory to store the seq   *
 *   An int pointer: This is a single integer.  The       *
 *                   function will store the length of    *
 *                   the current sequence in it.          *
 *   A char pointer: Should have sufficient storage when  *
 *                   called for recording the sequence's  *
 *                   description data (MAXDESCLINE+1).    *
 * Output :                                               *
 *   If a sequence was remaining in file, it is returned. *
 *   Otherwise, NULL is returned.                         */
char *get_fasta(FILE *file,char *seq,int *seqlength,char *desc);

/* Function returns the integer translation of a character    *
 * in the alphabet [AaCcGgTtUu].  Depending on preprocessing  *
 * conditions, any other non-whitespace character encountered *
 * will either throw an error or a random residue is returned *
 * via random_residue.                                        *
 *                                                            *
 * Transliteration scheme                                     *
 *     A   -> 0                                               *
 *     C   -> 1                                               *
 *     G   -> 2                                               *
 *     T/U -> 3                                               */
int trans(char c);

/* Return complement of argument (base) */
int basecomp(int x);

#endif
