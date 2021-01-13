/* bssm_utils.c
 *
 * Michael Sparks (mespar1@iastate.edu)
 * Last modified : 6 April 2007
 *
 * This is a collection of functions associated with
 * manipulating Bssmparm objects for training purposes.
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

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bssm_utils.h"
#include "probdef.h"

void init_bssm(Bssmparm *bssm) {
  int i,j,k,l,m; /* iterator variables */

  bssm->gtmodelset=FALSE;
  bssm->gcmodelset=FALSE;

  /* Init GT|AG and GC|AG models */
  for(i=0;i<TERMCOUNT;++i)
    for(j=0;j<HYPOTHESIS7;++j)
      for(k=0;k<WINSIZE+2;++k)
        for(l=0;l<ALFSIZE;++l)
          for(m=0;m<ALFSIZE;++m) {
            bssm->gtmodel.hypotables.hypo7table[i][j][k][l][m]=INITVAL_FLT;
            bssm->gcmodel.hypotables.hypo7table[i][j][k][l][m]=INITVAL_FLT;
          }

  return;
}

void build_bssm(int **seq_matrix, int num_entries, Bssmparm *bssm,
                char *type, int dim1, int dim2) {
  int mono_ct[STRINGSIZE-1][ALFSIZE],        /* Mononuc freq          */
      di_ct[STRINGSIZE-1][ALFSIZE][ALFSIZE], /* Dinuc freq            */
      i,j,k;                                 /* Iterator variables    */

  /* Verify dinucleotide terminus type */
  if(!strcmp(type,"GT"))
    bssm->gtmodelset=TRUE;
  else if (!strcmp(type,"GC"))
    bssm->gcmodelset=TRUE;
  else {
    fprintf(stderr,"Error: Unknown donor dinucleotide type specified!\n");
    exit(EXIT_FAILURE);
  }

  /* Record the ``obvious" data into the Bssmparm structure */
  if(!strcmp(type,"GT")) {
    /* We prefer the seven-hypothesis model to that with two. */
    bssm->gtmodel.hypothesisnum=HYPOTHESIS7;
    /* Just set maximum extent of splice signal *
     * to its maximum possible value.           */
    bssm->gtmodel.wsizedonorleft=MAXSPLICESIG;
    bssm->gtmodel.wsizedonorright=MAXSPLICESIG;
    bssm->gtmodel.wsizeacceptorleft=MAXSPLICESIG;
    bssm->gtmodel.wsizeacceptorright=MAXSPLICESIG;
  }
  else {
    /* We prefer the seven-hypothesis model to that with two. */
    bssm->gcmodel.hypothesisnum=HYPOTHESIS7;
    /* Just set maximum extent of splice signal *
     * to its maximum possible value.           */
    bssm->gcmodel.wsizedonorleft=MAXSPLICESIG;
    bssm->gcmodel.wsizedonorright=MAXSPLICESIG;
    bssm->gcmodel.wsizeacceptorleft=MAXSPLICESIG;
    bssm->gcmodel.wsizeacceptorright=MAXSPLICESIG;
  }

  /* Inits of local variables */
  for(i=0;i<(STRINGSIZE-1);++i)
    for(j=0;j<ALFSIZE;++j) {
      mono_ct[i][j]=INITVAL_INT;
      for(k=0;k<ALFSIZE;++k)
        di_ct[i][j][k]=INITVAL_INT; 
    }
  
  /* tabulate mono/dinucleotide frequencies */
  for(i=0;i<(STRINGSIZE-1);++i)
    for(j=0;j<num_entries;++j) {
      /* Note that the first dimension of mono_ct indexes the  *
       * position (only the first STRINGSIZE-1 are considered) *
       * that nucleotide seq_matrix[j][i] occurs in.           */
      ++mono_ct[i][seq_matrix[j][i]];        

      /* Note that the first dimension of di_ct indexes the    *
       * position that the first base of dinucleotide          *
       * [seq_matrix[j][i]][seq_matrix[j][i+1]] occurs in.     *
       * Only the first STRINGSIZE-1 positions are considered. */
      ++di_ct[i][seq_matrix[j][i]][seq_matrix[j][i+1]];
    }

  /* Record equilibrium frequencies (1st ``slot" in transition freqs). *
   * Note that only base i is ``important", but we populate all bases  *
   * indexed by j, for consistency.                                    */
  for(i=0;i<ALFSIZE;++i)
    for(j=0;j<ALFSIZE;++j) {
      if(!strcmp(type,"GT"))
        bssm->gtmodel.hypotables.hypo7table[dim1][dim2][0][i][j] =
          ((double)mono_ct[0][i]) / num_entries;
      else
        bssm->gcmodel.hypotables.hypo7table[dim1][dim2][0][i][j] =
          ((double)mono_ct[0][i]) / num_entries;
    }

  /* Populate the remaining transition frequencies.                         *
   * Note that bssm->g[t|c]model.hypotables.hypo7table[dim1][dim2][k][i][j] *
   * records the frequency that dinucleotide [i][j] occurs with base j      *
   * occurring in position k, relative to dinucs [i][z], z != j, also       *
   * occurring in that position in the training data.                       */
  for(k=1;k<STRINGSIZE;++k) {
    for(i=0;i<ALFSIZE;++i) { /* rows: give probability of transition *
                              * to a different base j predicated on  *
                              * the first base being i.              */
      for(j=0;j<ALFSIZE;++j) { /* columns: base being transitioned into */
        if(!strcmp(type,"GT")) {
          if(mono_ct[k-1][i] != 0)
            bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
              ((double)di_ct[k-1][i][j]) / mono_ct[k-1][i];
          else
            bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] = 0.0;
        }
        else {
          if(mono_ct[k-1][i] != 0)
            bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
              ((double)di_ct[k-1][i][j]) / mono_ct[k-1][i];
          else
            bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] = 0.0;
        }
      } /* end col for */

      /* Smooth zero-probability transition probabilities (empty counts) *
       * on a row-by-row basis.                                          *
       * Note that, since the ascii parameterization output is good only *
       * up to a precision of 1e-4, we will consider any input less than *
       * PROBMIN to be in need of serious smoothing.                     */
      if ( k < (STRINGSIZE - 52) || k > (STRINGSIZE - 50) ) {
        if (!strcmp(type,"GT")) {
          /* check if all four columns are < MAXFAULT */
          if ((bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][0]
                < MAXFAULT) &&
              (bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][1]
                < MAXFAULT) &&
              (bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][2]
                < MAXFAULT) &&
              (bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][3]
                < MAXFAULT))
            /* Adjust all probabilities */
            for(j=0;j<ALFSIZE;++j)
              bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                EQUIPROB;
          /* check if any of the four columns are < MAXFAULT */
          else if ((bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][0]
                     < MAXFAULT) ||
                   (bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][1]
                     < MAXFAULT) ||
                   (bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][2]
                     < MAXFAULT) ||
                   (bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][3]
                     < MAXFAULT))
            for (j=0;j<ALFSIZE;++j) {
              if (bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j]
                   < MAXFAULT)
                bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  PROBMIN;  
              else
                bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  (bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j]
                  * (1 - (4 * PROBMIN))) + PROBMIN; 
            }
          else /* All's well. */
            ;
        }
        else { /* GC donor */
          /* check if all four columns are < MAXFAULT */
          if ((bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][0]
                < MAXFAULT) &&
              (bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][1]
                < MAXFAULT) &&
              (bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][2]
                < MAXFAULT) &&
              (bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][3]
                < MAXFAULT))
            /* Adjust all probabilities */
            for(j=0;j<ALFSIZE;++j)
              bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                EQUIPROB;
          /* check if any of the four columns are < MAXFAULT */
          else if ((bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][0]
                     < MAXFAULT) ||
                   (bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][1]
                     < MAXFAULT) ||
                   (bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][2]
                     < MAXFAULT) ||
                   (bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][3]
                     < MAXFAULT))
            for (j=0;j<ALFSIZE;++j) {
              if (bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j]
                    < MAXFAULT)
                bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  PROBMIN;  
              else
                bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  (bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j]
                  * (1 - (4 * PROBMIN))) + PROBMIN; 
            }
          else /* All's well. */
            ;
        } /* end GT/GC if/else */
      } /* end k if */
    } /* end row for */

    /* Now, handle probabilities proximal to don/acc dinucleotides */
    if (k == 50) { /* one column of 1's, all else 0's */
      if (dim1 == 0) { /* donor */
        for(i=0;i<ALFSIZE;++i) {
          for(j=0;j<ALFSIZE;++j) {
            if (j == 2) {
              if (!strcmp(type,"GT"))
                bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  MAXPROB;
              else
                bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  MAXPROB;
            }
            else {
              if (!strcmp(type,"GT"))
                bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  NULLPROB;
              else
                bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  NULLPROB;
            }
          }
        }
      }
      else { /* acceptor */
        for(i=0;i<ALFSIZE;++i) {
          for(j=0;j<ALFSIZE;++j) {
            if (j == 0) {
              if (!strcmp(type,"GT"))
                bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  MAXPROB;
              else
                bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  MAXPROB;
            }
            else {
              if (!strcmp(type,"GT"))
                bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  NULLPROB;
              else
                bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  NULLPROB;
            }
          }
        }
      }
    } /* end k == 50 if */

    if (k == 51) { /* only one 1, all else 0's */
      if (dim1 == 0) { /* donor */
        for(i=0;i<ALFSIZE;++i) {
          for(j=0;j<ALFSIZE;++j) {
            if (!strcmp(type,"GT")) {
              if ( i == 2 && j == 3 )
                bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  MAXPROB;
              else
                bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  NULLPROB;
            }
            else {
              if ( i == 2 && j == 1 )
                bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  MAXPROB;
              else
                bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  NULLPROB;
            }
          }
        }
      }
      else { /* acceptor */
        for(i=0;i<ALFSIZE;++i) {
          for(j=0;j<ALFSIZE;++j) {
            if (!strcmp(type,"GT")) {
              if ( i == 0 && j == 2 )
                bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  MAXPROB;
              else
                bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  NULLPROB;
            }
            else {
              if ( i == 0 && j == 2 )
                bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  MAXPROB;
              else
                bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  NULLPROB;
            }
          }
        }
      } /* end don/acc if */
    } /* end k == 51 if */

    if (k == 52) { /* one row of non-0's, all else 0 */
      if (dim1 == 0) { /* donor */
        for(i=0;i<ALFSIZE;++i) {
          if (!strcmp(type,"GT")) {
            if (i == 3) { 
              /* check if all four columns are < MAXFAULT */
              if((bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][0]
                   < MAXFAULT) &&
                 (bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][1]
                   < MAXFAULT) &&
                 (bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][2]
                   < MAXFAULT) &&
                 (bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][3]
                   < MAXFAULT))
                /* Adjust all probabilities */
                for (j=0;j<ALFSIZE;++j)
                  bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                    EQUIPROB;
              /* check if any of the four columns are < MAXFAULT */
              else if((bssm->gtmodel.hypotables.hypo7table[dim1][dim2]
                        [k][i][0] < MAXFAULT) ||
                      (bssm->gtmodel.hypotables.hypo7table[dim1][dim2]
                        [k][i][1] < MAXFAULT) ||
                      (bssm->gtmodel.hypotables.hypo7table[dim1][dim2]
                        [k][i][2] < MAXFAULT) ||
                      (bssm->gtmodel.hypotables.hypo7table[dim1][dim2]
                        [k][i][3] < MAXFAULT))
                for(j=0;j<ALFSIZE;++j) {
                  if(bssm->gtmodel.hypotables.hypo7table[dim1][dim2]
                      [k][i][j] < MAXFAULT)
                    bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                      PROBMIN;  
                  else
                    bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                      (bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j]
                      * (1 - (4 * PROBMIN))) + PROBMIN; 
                }
              else /* All's well. */
                ;
            }
            else
              for(j=0;j<ALFSIZE;++j)
                bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  NULLPROB;
          }
          else {
            if (i == 1) { 
              /* check if all four columns are < MAXFAULT */
              if((bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][0]
                   < MAXFAULT) &&
                 (bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][1]
                   < MAXFAULT) &&
                 (bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][2]
                   < MAXFAULT) &&
                 (bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][3]
                   < MAXFAULT))
                /* Adjust all probabilities */
                for(j=0;j<ALFSIZE;++j)
                  bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                    EQUIPROB;
              /* check if any of the four columns are < MAXFAULT */
              else if((bssm->gcmodel.hypotables.hypo7table[dim1][dim2]
                        [k][i][0] < MAXFAULT) ||
                      (bssm->gcmodel.hypotables.hypo7table[dim1][dim2]
                        [k][i][1] < MAXFAULT) ||
                      (bssm->gcmodel.hypotables.hypo7table[dim1][dim2]
                        [k][i][2] < MAXFAULT) ||
                      (bssm->gcmodel.hypotables.hypo7table[dim1][dim2]
                        [k][i][3] < MAXFAULT))
                for(j=0;j<ALFSIZE;++j) {
                  if(bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j]
                      < MAXFAULT)
                    bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                      PROBMIN;  
                  else
                    bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                      (bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j]
                      * (1 - (4 * PROBMIN))) + PROBMIN; 
                }
              else /* All's well. */
                ;
            }
            else
              for(j=0;j<ALFSIZE;++j)
                bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  NULLPROB;
          } /* end GT/GC if */
        } /* end i for */
      }
      else { /* acceptor */
        for(i=0;i<ALFSIZE;++i) {
          if (i == 2) { 
            if (!strcmp(type,"GT")) {
              /* check if all four columns are < MAXFAULT */
              if((bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][0]
                   < MAXFAULT) &&
                 (bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][1]
                   < MAXFAULT) &&
                 (bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][2]
                   < MAXFAULT) &&
                 (bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][3]
                   < MAXFAULT))
                /* Adjust all probabilities */
                for (j=0;j<ALFSIZE;++j)
                  bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                    EQUIPROB;
              /* check if any of the four columns are < MAXFAULT */
              else if((bssm->gtmodel.hypotables.hypo7table[dim1][dim2]
                        [k][i][0] < MAXFAULT) ||
                      (bssm->gtmodel.hypotables.hypo7table[dim1][dim2]
                        [k][i][1] < MAXFAULT) ||
                      (bssm->gtmodel.hypotables.hypo7table[dim1][dim2]
                        [k][i][2] < MAXFAULT) ||
                      (bssm->gtmodel.hypotables.hypo7table[dim1][dim2]
                        [k][i][3] < MAXFAULT))
                for(j=0;j<ALFSIZE;++j) {
                  if(bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j]
                      < MAXFAULT)
                    bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                      PROBMIN;  
                  else
                    bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                      (bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j]
                      * (1 - (4 * PROBMIN))) + PROBMIN; 
                }
              else /* All's well. */
                ;
            }
            else {
              /* check if all four columns are < MAXFAULT */
              if((bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][0]
                   < MAXFAULT) &&
                 (bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][1]
                   < MAXFAULT) &&
                 (bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][2]
                   < MAXFAULT) &&
                 (bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][3]
                   < MAXFAULT))
                /* Adjust all probabilities */
                for(j=0;j<ALFSIZE;++j)
                  bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                    EQUIPROB;
              /* check if any of the four columns are < MAXFAULT */
              else if((bssm->gcmodel.hypotables.hypo7table[dim1][dim2]
                        [k][i][0] < MAXFAULT) ||
                      (bssm->gcmodel.hypotables.hypo7table[dim1][dim2]
                        [k][i][1] < MAXFAULT) ||
                      (bssm->gcmodel.hypotables.hypo7table[dim1][dim2]
                        [k][i][2] < MAXFAULT) ||
                      (bssm->gcmodel.hypotables.hypo7table[dim1][dim2]
                        [k][i][3] < MAXFAULT))
                for(j=0;j<ALFSIZE;++j) {
                  if(bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j]
                      < MAXFAULT)
                    bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                      PROBMIN;  
                  else
                    bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                      (bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j]
                      * (1 - (4 * PROBMIN))) + PROBMIN; 
                }
              else /* All's well. */
                ;
            } /* end GT/GC if */
          } /* end i==2 if */
          else
            for(j=0;j<ALFSIZE;++j) {
              if (!strcmp(type,"GT"))
                bssm->gtmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  NULLPROB;
              else
                bssm->gcmodel.hypotables.hypo7table[dim1][dim2][k][i][j] =
                  NULLPROB;
            }
        } /* end i for */
      } /* end don/acc if */
    } /* end k == 52 if */
  } /* end k for */

  return;
} /* end build_bssm */

/* Function to delineate the dimensions of our BSSM parameterization *
 * object to be updated, as per a naming convention described in     *
 * the header of bssm_utils.h                                        */
void parse_dimensions(char *filename,int *dim1, int *dim2) {
  char *chp1, 
       *chp2;

  chp1=strrchr(filename,'/');
  if (chp1 != NULL)
    ++chp1;
  else
    chp1=filename;

  chp2=strtok(chp1,"_");
  /* Set per Volker Brendel's splice site phase conventions */
  if (!strcmp(chp2,"T1"))
    *dim2=0;
  else if (!strcmp(chp2,"T2"))
    *dim2=1;
  else if (!strcmp(chp2,"T0"))
    *dim2=2;
  else if (!strcmp(chp2,"F1"))
    *dim2=3;
  else if (!strcmp(chp2,"F2"))
    *dim2=4;
  else if (!strcmp(chp2,"F0"))
    *dim2=5;
  else if (!strcmp(chp2,"Fi"))
    *dim2=6;
  else {
    fprintf(stderr,"Improper filename passed!\n");
    exit(EXIT_FAILURE);
  }

  chp2=strtok(NULL,"_");
  if (!strcmp(chp2,"don"))
    *dim1=0;
  else if (!strcmp(chp2,"acc"))
    *dim1=1;
  else {
      fprintf(stderr,"Improper filename passed!\n");
      exit(EXIT_FAILURE);
  }

  return;
} /* end parse_dimensions */

/* Function to print the contents of the bssm parameterization *
 * to stdout for debugging purposes.                           */
void echo_bssm(Bssmparm *bssm) {
  int i,j,k,l,m;

  fprintf(stderr,"Is the GT model set? -> %i\n",bssm->gtmodelset);
  fprintf(stderr,"Is the GC model set? -> %i\n\n",bssm->gcmodelset);

  if (bssm->gtmodelset == TRUE) {
    fprintf(stdout,"REPORTING GT|AG MODEL PARAMETERIZATION\n");
    for(i=0;i<TERMCOUNT;++i) {
      for(j=0;j<HYPOTHESIS7;++j) {
        fprintf(stdout,"\n\nTerminal: %i, Hypothesis: %i",i,j);
        for(k=0;k<STRINGSIZE;++k) {
          fprintf(stdout,"\n");
          for(l=0;l<ALFSIZE;++l) {
            fprintf(stdout,"\n");
            for(m=0;m<ALFSIZE;++m) {
              fprintf(stdout,"%.4f ",
                bssm->gtmodel.hypotables.hypo7table[i][j][k][l][m]);
    }}}}}
    fprintf(stdout,"\n\n");
  }

  if (bssm->gcmodelset == TRUE) {
    fprintf(stdout,"REPORTING GC|AG MODEL PARAMETERIZATION\n");
    for(i=0;i<TERMCOUNT;++i) {
      for(j=0;j<HYPOTHESIS7;++j) {
        fprintf(stdout,"\n\nTerminal: %i, Hypothesis: %i",i,j);
        for(k=0;k<STRINGSIZE;++k) {
          fprintf(stdout,"\n");
          for(l=0;l<ALFSIZE;++l) {
            fprintf(stdout,"\n");
            for(m=0;m<ALFSIZE;++m) {
              fprintf(stdout,"%.4f ",
                bssm->gcmodel.hypotables.hypo7table[i][j][k][l][m]);
    }}}}}
    fprintf(stdout,"\n\n");
  }

  return;
} /* end echo_bssm */
