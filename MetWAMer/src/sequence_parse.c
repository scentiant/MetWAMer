/* sequence_parse.c
 * Michael E Sparks (mespar1@gmail.com)
 * Last modified: 20 July 2013
 *
 * Code for parsing FASTA formatted sequences.
 *
 * Copyright (c) 2005,2006,2007,2013 Michael E Sparks
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

/* Includes */
#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "errors.h"
#include "sequence_parse.h"

/* Function definitions ******************************************************/

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
char *get_fasta(FILE *file,char *seq,int *seqlength,char *desc)
{
  int c,          /* For reading in letters   */
      length,     /* Total length of sequence */
      desclength; /* Total length of sequence *
                   * description              */
  fpos_t fpos;

  /* Verify and parse comment line */
  if((c=fgetc(file))!=(int)'>'||(bool)(feof(file)))
    if((bool)feof(file))
      return(NULL);
    else /* Improper format */
      FATALERROR("Err (get_fasta): This requires FASTA formatted input!\n")
  else {
    desclength=1;
    desc[0]=c;
    /* Scroll past the comment line and record it */
    while((c=fgetc(file))!=(int)'\n'&&c!=EOF)
      if(++desclength<=MAXDESCLINE)
        desc[desclength-1]=c;
      /* Else the description line exceeds permissable lengths */

    if(c==EOF) /* We should always expect sequence after the comment line */
      FATALERROR("Err (get_fasta): A comment, but...no sequence!??\n")
    else if(desclength<=MAXDESCLINE)
      desc[desclength]='\0';
    else
      desc[MAXDESCLINE]='\0';
  }

  /* Determine sequence length */
  (void)fgetpos(file,&fpos);
  for(length=0;(c=fgetc(file))!=(int)'>'&&c!=EOF;)
    if(!isspace(c))
      ++length;

  /* Copy data into sequence pointer */
  if((seq=(char*)malloc(sizeof(char)*(length+1)))!=NULL) {
    for(fsetpos(file,&fpos),length=0;
        (c=fgetc(file))!=(int)'>'&&c!=EOF;)
      if(!isspace(c))
        seq[length++]=(char)c;
    seq[length]='\0';
    *seqlength=length;

    /* Re-position for the next invocation */
    if((bool)ungetc(c,file))
      return(seq);
    else
      FATALERROR("Err (get_fasta): \
Error putting character back to input stream.\n")
  }
  else
    FATALERROR("Err (get_fasta): Out of memory!\n")

  FATALERROR("Err (get_fasta): Unexpected condition!\n")
  return(NULL);
} /* end get_fasta */

/* Returns a random value in the interval [0,NTALFSIZE-1]. */
int random_residue(void)
{
  return(rand()%NTALFSIZE);
} /* end random_residue */

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
int trans(char c)
{
  int cint=-1; /* int translation of c */

  if(!isspace((int)c)) {
    switch(c) {
      case 'A' : case 'a' :
        cint=Ade;
        break;
      case 'C' : case 'c' :
        cint=Cyt;
        break;
      case 'G' : case 'g' :
        cint=Gua;
        break;
      case 'T' : case 't' : case 'U' : case 'u' :
        cint=Thy;
        break;
      default :
        #ifdef STRICTTEST
        FATALERROR("Err (trans): Ambiguous test char encountered \
and strictness requested!\n")
        #endif
        cint=random_residue();
    } /* end switch */
  } /* end if */
  else
    FATALERROR("Err (trans): Whitespace passed to trans fxn!\n")

  return(cint);
} /* end trans */

/* Return complement of argument (base) */
int basecomp(int x) {
  int comp=-1;

  if(x==Ade)
    comp=Thy;
  else if(x==Cyt)
    comp=Gua;
  else if(x==Gua)
    comp=Cyt;
  else if(x==Thy)
    comp=Ade;
  else {
    (void)snprintf(errbuff,MAXERRLEN,"Err (basecomp): \
Bad base int representation: %i\n",x);
    FATALERROR(errbuff)
  }
  return(comp);
}
