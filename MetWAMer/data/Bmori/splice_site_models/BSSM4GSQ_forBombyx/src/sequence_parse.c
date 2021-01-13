/* sequence_parse.c
 * Michael Sparks (mespar1@iastate.edu)
 * Last modified : 3 April 2007
 *
 * Code for parsing FASTA formatted sequences.
 *
 * Please Note: Transliteration scheme
 *   A   -> 0
 *   C   -> 1
 *   G   -> 2
 *   T/U -> 3
 *
 * Copyright (c) 2004,2007 Michael E Sparks
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

#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "sequence_parse.h"

/* get_fasta:  A general function to parse FASTA entries  *
 * Input :                                                *
 *   A file pointer: This must be connected to a stream   *
 *                   containing Fasta-formatted sequences *
 *   An int pointer: This must have sufficient space      *
 *                   allocated to store the integer       *
 *                   translation of the sequence in the   *
 *                   file.                                *
 *   An int pointer: This is a single integer.  The       *
 *                   function will store the length of    *
 *                   the current sequence in it.          *
 * Output : Returns TRUE if a sequence has been           *
 *            encountered.                                *
 *          Returns FALSE if the end of file has been     *
 *            reached.                                    */

int get_fasta(FILE *file,int *seq,int *seqlength,const int alphsize) {
  int c,      /* For reading in letters   */
      length, /* Total length of sequence */
      i;      /* Iterator variable        */

  /* Verify comment line */
  if ( ((c=fgetc(file)) != '>') || (feof(file)) ) {
    if (feof(file))
      return(FALSE);
    else { /* Improper format */
      fprintf(stderr,"Error: This requires FASTA formatted input!\n");
      return(FALSE);
    }
  }

  /* Scroll past the comment line and ignore it */
  while ( ((c=fgetc(file)) != '\n') && (c != EOF) )
    ;
  /* We should always expect sequence after the comment line */
  if (c == EOF) {
    fprintf(stderr,"Error: A comment, but...no sequence!??\n");
    exit(EXIT_FAILURE);
  }

  /* Copy data into sequence pointer */
  for(i=length=0;((c=fgetc(file)) != '>') && (c != EOF); ) {
    if(!isspace(c)) {
      switch (c) {
        case 'A' :
        case 'a' :
          c=0;
          break;
        case 'C' :
        case 'c' :
          c=1;
          break;
        case 'G' :
        case 'g' :
          c=2;
          break;
        case 'T' :
        case 't' :
        case 'U' :
        case 'u' :
          c=3;
          break;
        default :
          #ifdef STRICTTRAIN
          fprintf(stderr,"Error in input file sequence! \
(no ambiguity -- encountered the symbol %c)\n",c);
          exit(EXIT_FAILURE);
          #endif
          c=random_residue(alphsize);
      } /* end switch */
      seq[i++]=c;
      ++length;
    }
  } /* end for */
  *seqlength=length;

  /* Re-position for the next invocation */
  ungetc(c,file); 
  return(TRUE);
} /* end get_fasta */

/* Returns a random value in the interval [0,ALFSIZE-1] */
int random_residue(const int alphsize) {
  srand(time(NULL)); /* Requires <time.h> */
  return(rand()%(alphsize-1));
} /* end random_residue */
