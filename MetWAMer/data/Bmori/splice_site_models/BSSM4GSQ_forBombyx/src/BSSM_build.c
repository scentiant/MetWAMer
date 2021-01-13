/* BSSM_build.c
 * Michael Sparks (mespar1@iastate.edu)
 * Last modified : 3 April 2007
 *
 * Main Implementation file for "BSSM_build.x", a program to 
 * build a 1st-order Markov model (position weight array) for
 * splice sites.
 *
 * Please Note: Transliteration scheme
 *   A   -> 0
 *   C   -> 1
 *   G   -> 2
 *   T/U -> 3
 *
 * Input data files are--strictly!--named as follows, using
 * an obvious schema:
 *   F0_don  F1_don  F2_don  Fi_don  T0_don  T1_don  T2_don
 *   F0_acc  F1_acc  F2_acc  Fi_acc  T0_acc  T1_acc  T2_acc
 *
 * Phase is denoted as follows (Brendel conventions):
 *   1 -> C O D |
 *   2 -> C | O D
 *   0 -> C O | D
 *
 * Copyright (c) 2003,2004,2006,2007 Michael E Sparks
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

/* Include Statements */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "probdef.h"
#include "sequence_parse.h"
#include "bssm_utils.h"

/* String literals */
#define ARGCT 5
#define USAGE "Usage: %s BSSM_parmfile (grep -c '^>' in_file) \
G[T|C] in_file\n"

/* Main Application */
int main(int argc, char *argv[]) {
  Bssmparm my_bssmparm;     /* Stores model parameterization  */
  FILE *f_ptr1=NULL,        /* File pointers                  */
       *f_ptr2=NULL;
  int sequence[STRINGSIZE], /* For storing int equivalent of  *
                             * sequence                       */
      length_token,         /* Stores length of string        */
      number_seqs,          /* Row size, num seqs             */
      **seq_matrix=NULL,    /* Stores sequence data           */
      dimension1,           /* donor or acceptor?             */
      dimension2,           /* T1, T2, T0, F1, F2, F0, or Fi? */
      i,j;                  /* Iterator variables             */

  /* open input file */
  if (argc != ARGCT) {
    fprintf(stderr,USAGE,argv[0]);
    exit(EXIT_FAILURE);
  }

  /* Determine if the BSSM parameter file already exists *
   * in the filesystem.  If so, update. Else, create it. */
  if ((f_ptr1=fopen(argv[1],"rb"))!=NULL) {
    if (fread(&my_bssmparm,sizeof(Bssmparm),1,f_ptr1) != 1) {
      fprintf(stderr,"Error importing Bssmparm data from %s!\n",argv[1]);
      fclose(f_ptr1);
      exit(EXIT_FAILURE);
    }
  }
  else {
    /* Initialize our parameter object */
    init_bssm(&my_bssmparm);

    /* Prepare a file for writing Bssmparm object to */
    if ((f_ptr1=fopen(argv[1],"wb"))==NULL) {
      fprintf(stderr,"Can't open %s for binary writing!\n",argv[1]);
      exit(EXIT_FAILURE);
    }
  }

  /* Open input file */
  if ((f_ptr2=fopen(argv[4],"rt"))==NULL) {
    fprintf(stderr,"Can't open %s for ascii reading!\n",argv[4]);
    fclose(f_ptr1);
    exit(EXIT_FAILURE);
  }

  /* Allocate memory for the sequence storage matrix */
  number_seqs=(int)atoi(argv[2]);
  seq_matrix=(int**)malloc(number_seqs*sizeof(int*));
  if (seq_matrix == NULL) {
    fprintf(stderr,"Insufficient Memory!\n");
    fclose(f_ptr1);
    fclose(f_ptr2);
    exit(EXIT_FAILURE);
  }
  else
    for(i=0;i<number_seqs;++i) {
      seq_matrix[i]=(int*)malloc(STRINGSIZE*sizeof(int));
      if (seq_matrix[i] == NULL) {
        fprintf(stderr,"Insufficient Memory!\n");
        while (--i >= 0)
          free(seq_matrix[i]);
        free(seq_matrix);
        fclose(f_ptr1);
        fclose(f_ptr2);
        exit(EXIT_FAILURE);
      }
    }

  /* Process each input sequence string */
  for(i=0;get_fasta(f_ptr2,sequence,&length_token,ALFSIZE);++i) {
    if (length_token != STRINGSIZE) {
      fprintf(stderr,"Inconsistency with input file detected!\n");
      fclose(f_ptr1);
      fclose(f_ptr2);
      exit(EXIT_FAILURE);
    }
    else
      /* read in the sequence */
      for (j=0;j<STRINGSIZE;++j)
        seq_matrix[i][j]=sequence[j];
  }

  /* Since we opened the input file, we don't need its name anymore, but *
   * we need the tokens that comprise the name (per a rigid input file   *
   * naming scheme described in the header of this file) as information  *
   * to key into appropriate dimensions of the Bssmparm hypo7table       */
  parse_dimensions(argv[4],&dimension1,&dimension2);

  /* Record the parameterization */
  build_bssm(seq_matrix,number_seqs,&my_bssmparm,
             argv[3],dimension1,dimension2);

  /* Close files */
  fclose(f_ptr1);
  fclose(f_ptr2);

  /* Overwrite the binary file with our new object. */
  if ((f_ptr1=fopen(argv[1],"wb"))==NULL) {
    fprintf(stderr,"Unable to open %s for binary writing!\n",argv[1]);
    exit(EXIT_FAILURE);
  }
  else {
    /* Write parameterization to output file */
    if (fwrite(&my_bssmparm,sizeof(Bssmparm),1,f_ptr1) != 1) {
      fprintf(stderr,"Error exporting Bssmparm data to %s!\n",argv[1]);
      fclose(f_ptr1);
      exit(EXIT_FAILURE);
    }
    else
      fclose(f_ptr1);
  }

  /* Free memory */
  for(i=0;i<number_seqs;++i)
    free(seq_matrix[i]);
  free(seq_matrix);

  return(EXIT_SUCCESS);
} /* end main */
