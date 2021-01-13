/* sequence_parse.c
 * Michael E Sparks (mespar1@gmail.com)
 *
 * Contains code for generic functions meant to manipulate
 * biological sequences stored in common-format files.
 *
 * Copyright (C) 2005,2006 Michael E Sparks
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/* Includes */
#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "immpractical.h"
#include "sequence_parse.h"

/* Local function prototypes *************************************************/

/* Returns a random value in the interval *
 * [0,ALFSIZE-1]. Requires <time.h>.      */
static int random_residue(void);

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
 *                   the current sequence in it.          */
char *get_fasta(FILE *file,char *seq,int *seqlength)
{
  int c,      /* For reading in letters   */
      length; /* Total length of sequence */
  fpos_t fpos;

  /* Verify comment line */
  if ( ((c=fgetc(file)) != (int)'>') || (bool)(feof(file)) ) {
    if ((bool)feof(file)) {
      *seqlength=0;
      return(NULL);
    }
    else { /* Improper format */
      fprintf(stderr,
        "Err (get_fasta): This requires FASTA formatted input!\n");
      *seqlength=0;
      return(NULL);
    }
  }
  else {
    /* Scroll past the comment line and ignore it */
    while ( ((c=fgetc(file)) != (int)'\n') && (c != EOF) )
      ;
    /* We should always expect sequence after the comment line */
    if (c == EOF) {
      fprintf(stderr,"Err (get_fasta): A comment, but...no sequence!??\n");
      *seqlength=0;
      return(NULL);
    }
  }
  /* Determine sequence length */
  (void)fgetpos(file,&fpos);
  for(length=0;((c=fgetc(file)) != (int)'>') && (bool)(c != EOF); )
    if(!isspace(c))
      ++length;

  /* Copy data into sequence pointer */
  if ((seq=(char*)malloc(sizeof(char)*(length+1)))!=NULL) {
    for((void)fsetpos(file,&fpos),length=0;
        ((c=fgetc(file)) != (int)'>') && (bool)(c != EOF);
       )
      if(!isspace(c))
        seq[length++] = (char)c;
    seq[length]='\0';
    *seqlength=length;

    /* Re-position for the next invocation */
    if ((bool)ungetc(c,file))
      return(seq);
    else {
      fprintf(stderr,
        "Err (get_fasta): Error putting character back to input stream.\n");
      exit(EXIT_FAILURE);
    }
  }
  else {
    fprintf(stderr,"Err (get_fasta): Out of memory!\n");
    exit(EXIT_FAILURE);
  }
} /* end get_fasta */

/* Returns a random value in the interval *
 * [0,ALFSIZE-1]. Requires <time.h>.      */
static int random_residue(void)
{
  srand((unsigned int)time(NULL));
  return(rand() % (ALFSIZE-1));
} /* end random_residue */

/* Function returns the integer translation of a character     *
 * in the alphabet [AaCcGgTtUu].  Depending on preprocessing   *
 * conditions, any other non-whitespace character encountered  *
 * will either throw and error or a random residue is returned *
 * via random_residue.                                         *
 *                                                             *
 * Transliteration scheme                                      *
 *     A   -> 0                                                *
 *     C   -> 1                                                *
 *     G   -> 2                                                *
 *     T/U -> 3                                                */
int trans(char c)
{
  int cint; /* int translation of c */

  if(!isspace((int)c)) {
    switch (c) {
      case 'A' : case 'a' :
        cint=0;
        break;
      case 'C' : case 'c' :
        cint=1;
        break;
      case 'G' : case 'g' :
        cint=2;
        break;
      case 'T' : case 't' : case 'U' : case 'u' :
        cint=3;
        break;
      default :
        #ifdef STRICTTRAIN
        fprintf(stderr,"Err (trans): Error in input file sequence! ");
        fprintf(stderr,
          "Err (trans): (no ambiguity -- encountered the symbol %c)\n",c);
        exit(EXIT_FAILURE);
        #endif
        #ifdef STRICTTEST
        fprintf(stderr,"Err (trans): Ambiguous test char encountered \
and strictness requested!\n");
        exit(EXIT_FAILURE);
        #endif
        cint=random_residue();
    } /* end switch */
  } /* end if */
  else {
    fprintf(stderr,"Err (trans): Whitespace passed to trans fxn!\n");
    exit(EXIT_FAILURE);
  }

  return(cint);
} /* end trans */

/* Return complement of argument (base) */
int basecomp(int x) {
  int comp;
  if (x == 0)
    comp = 3;
  else if (x == 1)
    comp = 2;
  else if (x == 2)
    comp = 1;
  else if (x == 3)
    comp = 0;
  else {
    fprintf(stderr,"Err (basecomp): Bad base int representation: %i\n",x);
    exit(EXIT_FAILURE);
  }
  return(comp);
}
