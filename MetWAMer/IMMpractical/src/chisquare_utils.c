/* chisquare_utils.c
 * Michael E Sparks (mespar1@gmail.com)
 *
 * This file implements routines specific to computing chi-square-
 * based interpolated Markov model oligomer probabilities.
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

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "chisquare_utils.h"
#include "immpractical.h"
#ifdef PROBREPORT
#include <stdio.h>
#endif
#ifdef ADJUSTNULLS
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

/* Local struct defs *******************************************************/

typedef struct { /* struct for temporarily storing weights of pretexts */
  double lamhis0,
         lamhis1[ALFSIZE],
         lamhis2[ALFSIZE][ALFSIZE],
         lamhis3[ALFSIZE][ALFSIZE][ALFSIZE],
         lamhis4[ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE],
         lamhis5[ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE];
} lambdatmpT;

/* Local function/subroutine prototypes ************************************/

/* Given a populated contingency table, this function returns the  *
 * appropriate pretext weighting for the history it was called on. */
static double calc_chisquareweight(double obs_matrix[][ALFSIZE],
  int historyct);

/* This function computes the chi**2 confidence value (complement of *
 * P-value, i.e., area under left tail of the curve) to assist with  *
 * derivation of pretext weightings per 4th eqn in Salzberg et al.   *
 * The formulae are based on the classic Peizer and Pratt 1968 paper *
 * on approximating P-values of various statistical distributions:   *
 * Journal of the American Statistical Association, 63:1416-56       */
static double calc_chisq_conf(double chisq_stat);

/* Function definitions ******************************************************/

/* Subroutine to coordinate development of final, smoothed oligo *
 * likelihoods by the chi**2 interpolated training method.       */
void CHISQUARE_final_probs(immprobT *finalprob,relfreqsT *rfreqs,
  countarrayT *tal,int model)
{
  double obs_matrix[HYPONUM][ALFSIZE]; /* contingency table. I expect this *
           * to be 2X4! The first row corresponds to relative frequencies  *
           * of the current order + 3'-terminal nucleotide (which is       *
           * specified using the 2nd dimension), whereas the second row    *
           * corresponds to smoothed IMM probabilities of the next-        *
           * shorter history + 3'-terminal nucleotide. (truncate from 5')  *
           * Though this data structure is storing counts, the type is a   *
           * float. This is because expected values might not be integral  *
           * (which is legitimate for chi**2 independence tests), and I    *
           * prefer not to mix types when taking differences.              */
  int i,j,k,l,m,z,   /* base iterator variables                            */
      base;          /* base iterator specifically used to denote          *
                      * the reducing terminus of the sugar ;)              */
  lambdatmpT lambda; /* for storing weights of pretexts                    */
  #ifndef EMPTYEQUIV
  int order, /* iterator variable                                          */
      tmpshortcts[MAXORDER];  /* A temporary array to store counts of      *
        * successive 5'-truncated histories, to address training for       *
        * oligos whose histories did not occur in the training dataset     */
  double tmpfprobs[MAXORDER]; /* A temporary array to store final          *
        * smoothed probabilities of oligomers based on shorter pretext     *
        * lengths.                                                         */
  /* Note that in the previous two arrays, the zeroth-element is not used; *
   * this is just to make the indexing more intuitive.                     */
  #endif
  #ifdef PROBREPORT
  double pmf_test, /* Used to ascertain that a valid PMF was formed */
         pfix;     /* Statistically based adjustment for null probs */
  #else
  #ifdef ADJUSTNULLS
  double pmf_test, /* Used to ascertain that a valid PMF was formed */
         pfix;     /* Statistically based adjustment for null probs */
  #endif
  #endif

  /* Sequentially develop the model for oligos    *
   * based on histories of order 0 up to MAXORDER */
  #ifdef PROBREPORT
  fprintf(stdout,"\nNow reporting PMFs for model %u:\n\n",model);
  #endif

  /* This for-loop + switch facilitate MAXORDER ranging from [0..5] */
  for(z=0;z<=MAXORDER;++z) {
    switch(z) {
      case(0) :
        /* Smoothed probabilities of mononucleotides are just ML estimates, *
         * assuming the weight of an empty history is certainty (I can't    *
         * imagine that it's not!).                                         */
        lambda.lamhis0=CERTPROB;
        for(base=0;base<ALFSIZE;++base)
          /* Note that mult by lambda.lamhis0 would just yield the identity */
          finalprob->prob1[model][base] = rfreqs->MLhis0[base];
        #ifdef PROBREPORT
        /* Have we computed a valid PMF? */
        for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
          fprintf(stdout,"%u=%.4f,",base,finalprob->prob1[model][base]);
          pmf_test += finalprob->prob1[model][base];
        }
        fprintf(stdout," Sum of prob's for null history is %.2f ",pmf_test);
        if( fabs(pmf_test - CERTPROB) < MAXDIFF )
          fprintf(stdout,"valid\n");
        else
          fprintf(stdout,"problematic!\n");
        #endif
        break;

      case(1) :
        /* Set weights of mononucleotide pretexts */
        for(i=0;i<ALFSIZE;++i) {
          if (tal->count1[model][i] >= MINRELIABLECHISQCT)
            lambda.lamhis1[i]=CERTPROB;
          else if (tal->count1[model][i] == 0) {
            /* If this ever happened, something is terribly amiss! */
            fprintf(stdout,"Err (CHISQUARE_final_probs): \
Training data unrealistically sparse!\n");
            exit(EXIT_FAILURE);
          }
          else {
            /* Populate contingency table */
            for(base=0;base<ALFSIZE;++base) {
              obs_matrix[0][base]=(double)tal->count2[model][i][base];
              obs_matrix[1][base]=finalprob->prob1[model][base] *
                (double)tal->count1[model][i];
            }
            /* determine weight */
            lambda.lamhis1[i] = calc_chisquareweight(obs_matrix,
              tal->count1[model][i]);
          } 
        }
        /* Now record final, smoothed dinucleotide probabilities */
        for(i=0;i<ALFSIZE;++i) {
          for(base=0;base<ALFSIZE;++base)
            finalprob->prob2[model][i][base] =
              (lambda.lamhis1[i] * rfreqs->MLhis1[i][base]) +
              ((CERTPROB - lambda.lamhis1[i]) *
               finalprob->prob1[model][base]);
          #ifdef PROBREPORT
          for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
            fprintf(stdout,"%u=%.4f,",base,
              finalprob->prob2[model][i][base]);
            pmf_test += finalprob->prob2[model][i][base];
          }
          fprintf(stdout,
            " Sum of prob's for history %u is %.2f ",i,pmf_test);
          if( fabs(pmf_test - CERTPROB) < MAXDIFF )
            fprintf(stdout,"valid\n");
          else
            fprintf(stdout,"problematic!\n");
          #endif
          #ifdef ADJUSTNULLS
          /* First, see if any null probability adjustments need be made */
          if ( (fabs(finalprob->prob2[model][i][0]) < MAXDIFF) ||
               (fabs(finalprob->prob2[model][i][1]) < MAXDIFF) ||
               (fabs(finalprob->prob2[model][i][2]) < MAXDIFF) ||
               (fabs(finalprob->prob2[model][i][3]) < MAXDIFF)
             ) {
            /* This can occur if (provided EMPTYEQUIV is not defined),  *
             * though an oligo's history did occur .GE.                 *
             * MINRELIABLECHISQCT times, the history + some 3'-terminal *
             * nucleotide was not observed in the training data.        */
            fprintf(stdout,"Making null probability adjustments...\n");
            pfix=(double)MAX(PfixOMEGA,1.0 / tal->count1[model][i]);
            for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
              if( fabs(finalprob->prob2[model][i][base]) < MAXDIFF )
                finalprob->prob2[model][i][base] = pfix;
              else
                finalprob->prob2[model][i][base] =
                  (finalprob->prob2[model][i][base] * (1 - (4 * pfix))) + pfix;

              fprintf(stdout,"%u=%.4f,",base,
                finalprob->prob2[model][i][base]);
              pmf_test += finalprob->prob2[model][i][base];
            } /* end for */
            fprintf(stdout,
              " Sum of prob's (revised) for history %u is %.2f ",i,pmf_test);
            if( fabs(pmf_test - CERTPROB) < MAXDIFF )
              fprintf(stdout,"valid\n");
            else
              fprintf(stdout,"problematic!\n");
          } /* end if (nulls exist in current pmf) */
          #endif
        } /* end pretext for loop(s) */
        break;

      case(2) :
        /* Set weights of dinucleotide pretexts */
        for(i=0;i<ALFSIZE;++i) {
        for(j=0;j<ALFSIZE;++j) {
          if (tal->count2[model][i][j] >= MINRELIABLECHISQCT)
            lambda.lamhis2[i][j]=CERTPROB;
          else if (tal->count2[model][i][j] == 0)
            /* Here, we're in the business of setting weights for pretexts. *
             * If one doesn't occur, then the problem is moot, as the final *
             * smoothed probability for the unobserved oligomer will be set *
             * based on the longest of the 5'-truncated subhistories that   *
             * was seen.  We'll tackle this problem in the code below.      */
            continue;
          else {
            /* Populate contingency table */
            for(base=0;base<ALFSIZE;++base) {
              obs_matrix[0][base]=(double)tal->count3[model][i][j][base];
              obs_matrix[1][base]=finalprob->prob2[model][j][base] *
                (double)tal->count2[model][i][j];
            }
            /* determine weight */
            lambda.lamhis2[i][j] = calc_chisquareweight(obs_matrix,
              tal->count2[model][i][j]);
          }
        }}
        /* Now record final, smoothed trinucleotide probabilities */
        for(i=0;i<ALFSIZE;++i) {
        for(j=0;j<ALFSIZE;++j) {
          for(base=0;base<ALFSIZE;++base) {
            if (tal->count2[model][i][j] == 0) {
              #ifdef EMPTYEQUIV
              /* Note that induce_pseudocounts would have   *
               * been called in the imm_probmaker function. */
              fprintf(stdout,"Err (CHISQUARE_final_probs): \
induce_pseudocounts left empty devel history!\n");
              exit(EXIT_FAILURE);
              #else
              /* In this case, there is only one possible shorter history, *
               * a single nucleotide (which *must* occur), based on code   *
               * earlier in this function.  So, we would set the final     *
               * probability for the current oligo (based on the non-      *
               * occurring history) to the only possible compensatory one  *
               * that can be assigned to it.                               */
              finalprob->prob3[model][i][j][base]=
                finalprob->prob2[model][j][base];
              #endif
            }
            else
              finalprob->prob3[model][i][j][base] =
                (lambda.lamhis2[i][j] * rfreqs->MLhis2[i][j][base]) +
                ((CERTPROB - lambda.lamhis2[i][j]) *
                 finalprob->prob2[model][j][base]);
          } /* end terminal 3'-nt for loop */
          #ifdef PROBREPORT
          for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
            fprintf(stdout,"%u=%.4f,",base,
              finalprob->prob3[model][i][j][base]);
            pmf_test += finalprob->prob3[model][i][j][base];
          }
          fprintf(stdout,
            " Sum of prob's for history %u%u is %.2f ",i,j,pmf_test);
          if( fabs(pmf_test - CERTPROB) < MAXDIFF )
            fprintf(stdout,"valid\n");
          else
            fprintf(stdout,"problematic!\n");
          #endif
          #ifdef ADJUSTNULLS
          if ( (fabs(finalprob->prob3[model][i][j][0]) < MAXDIFF) ||
               (fabs(finalprob->prob3[model][i][j][1]) < MAXDIFF) ||
               (fabs(finalprob->prob3[model][i][j][2]) < MAXDIFF) ||
               (fabs(finalprob->prob3[model][i][j][3]) < MAXDIFF)
             ) {
            fprintf(stdout,"Making null probability adjustments...\n");
            pfix=MAX(PfixOMEGA,1.0 / tal->count2[model][i][j]);
            for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
              if( fabs(finalprob->prob3[model][i][j][base]) < MAXDIFF )
                finalprob->prob3[model][i][j][base] = pfix;
              else
                finalprob->prob3[model][i][j][base] =
                  (finalprob->prob3[model][i][j][base] *
                  (1 - (4 * pfix))) + pfix;

              fprintf(stdout,"%u=%.4f,",base,
                finalprob->prob3[model][i][j][base]);
              pmf_test += finalprob->prob3[model][i][j][base];
            } /* end for */
            fprintf(stdout,
              " Sum of prob's (revised) for history %u%u is %.2f ",
              i,j,pmf_test);
            if( fabs(pmf_test - CERTPROB) < MAXDIFF )
              fprintf(stdout,"valid\n");
            else
              fprintf(stdout,"problematic!\n");
          } /* end if (nulls exist in current pmf) */
          #endif
        }} /* end pretext for loop(s) */
        break;

      case(3) :
        /* Set weights of trinucleotide pretexts */
        for(i=0;i<ALFSIZE;++i) {
        for(j=0;j<ALFSIZE;++j) {
        for(k=0;k<ALFSIZE;++k) {
          if (tal->count3[model][i][j][k] >= MINRELIABLECHISQCT)
            lambda.lamhis3[i][j][k]=CERTPROB;
          else if (tal->count3[model][i][j][k] == 0)
            continue;
          else {
            /* Populate contingency table */
            for(base=0;base<ALFSIZE;++base) {
              obs_matrix[0][base]=(double)tal->count4[model][i][j][k][base];
              obs_matrix[1][base]=finalprob->prob3[model][j][k][base] *
                (double)tal->count3[model][i][j][k];
            }
            /* determine weight */
            lambda.lamhis3[i][j][k] = calc_chisquareweight(obs_matrix,
              tal->count3[model][i][j][k]);
          }
        }}}
        /* Now record final, smoothed tetranucleotide probabilities */
        for(i=0;i<ALFSIZE;++i) {
        for(j=0;j<ALFSIZE;++j) {
        for(k=0;k<ALFSIZE;++k) {
          for(base=0;base<ALFSIZE;++base) {
            if (tal->count3[model][i][j][k] == 0) {
              #ifdef EMPTYEQUIV
              /* Note that induce_pseudocounts would have   *
               * been called in the imm_probmaker function. */
              fprintf(stdout,"Err (CHISQUARE_final_probs): \
induce_pseudocounts left empty devel history!\n");
              exit(EXIT_FAILURE);
              #else
              tmpshortcts[2]=tal->count2[model][j][k];
              tmpshortcts[1]=tal->count1[model][k];
              tmpfprobs[2]=finalprob->prob3[model][j][k][base];
              tmpfprobs[1]=finalprob->prob2[model][k][base];
              for(order=2;order>0;--order) {
                if (tmpshortcts[order] > 0) {
                  finalprob->prob4[model][i][j][k][base]=
                    tmpfprobs[order];
                  break;
                }
              }
              #endif
            }
            else
              finalprob->prob4[model][i][j][k][base] =
                (lambda.lamhis3[i][j][k] * rfreqs->MLhis3[i][j][k][base]) +
                ((CERTPROB - lambda.lamhis3[i][j][k]) *
                 finalprob->prob3[model][j][k][base]);
          } /* end terminal 3'-nt for loop */
          #ifdef PROBREPORT
          for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
            fprintf(stdout,"%u=%.4f,",base,
              finalprob->prob4[model][i][j][k][base]);
            pmf_test += finalprob->prob4[model][i][j][k][base];
          }
          fprintf(stdout,
            " Sum of prob's for history %u%u%u is %.2f ",i,j,k,pmf_test);
          if( fabs(pmf_test - CERTPROB) < MAXDIFF )
            fprintf(stdout,"valid\n");
          else
            fprintf(stdout,"problematic!\n");
          #endif
          #ifdef ADJUSTNULLS
          if ( (fabs(finalprob->prob4[model][i][j][k][0]) < MAXDIFF) ||
               (fabs(finalprob->prob4[model][i][j][k][1]) < MAXDIFF) ||
               (fabs(finalprob->prob4[model][i][j][k][2]) < MAXDIFF) ||
               (fabs(finalprob->prob4[model][i][j][k][3]) < MAXDIFF)
             ) {
            fprintf(stdout,"Making null probability adjustments...\n");
            pfix=MAX(PfixOMEGA,1.0 / tal->count3[model][i][j][k]);
            for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
              if( fabs(finalprob->prob4[model][i][j][k][base]) < MAXDIFF )
                finalprob->prob4[model][i][j][k][base] = pfix;
              else
                finalprob->prob4[model][i][j][k][base] =
                  (finalprob->prob4[model][i][j][k][base] *
                  (1 - (4 * pfix))) + pfix;

              fprintf(stdout,"%u=%.4f,",base,
                finalprob->prob4[model][i][j][k][base]);
              pmf_test += finalprob->prob4[model][i][j][k][base];
            } /* end for */
            fprintf(stdout,
              " Sum of prob's (revised) for history %u%u%u is %.2f ",
              i,j,k,pmf_test);
            if( fabs(pmf_test - CERTPROB) < MAXDIFF )
              fprintf(stdout,"valid\n");
            else
              fprintf(stdout,"problematic!\n");
          } /* end if (nulls exist in current pmf) */
          #endif
        }}} /* end pretext for loop(s) */
        break;

      case(4) :
        /* Set weights of tetranucleotide pretexts */
        for(i=0;i<ALFSIZE;++i) {
        for(j=0;j<ALFSIZE;++j) {
        for(k=0;k<ALFSIZE;++k) {
        for(l=0;l<ALFSIZE;++l) {
          if (tal->count4[model][i][j][k][l] >= MINRELIABLECHISQCT)
            lambda.lamhis4[i][j][k][l]=CERTPROB;
          else if (tal->count4[model][i][j][k][l] == 0)
            continue;
          else {
            /* Populate contingency table */
            for(base=0;base<ALFSIZE;++base) {
              obs_matrix[0][base]=(double)tal->count5[model][i][j][k][l][base];
              obs_matrix[1][base]=finalprob->prob4[model][j][k][l][base] *
                (double)tal->count4[model][i][j][k][l];
            }
            /* determine weight */
            lambda.lamhis4[i][j][k][l] = calc_chisquareweight(obs_matrix,
              tal->count4[model][i][j][k][l]);
          }
        }}}}
        /* Now record final, smoothed pentanucleotide probabilities */
        for(i=0;i<ALFSIZE;++i) {
        for(j=0;j<ALFSIZE;++j) {
        for(k=0;k<ALFSIZE;++k) {
        for(l=0;l<ALFSIZE;++l) {
          for(base=0;base<ALFSIZE;++base) {
            if (tal->count4[model][i][j][k][l] == 0) {
              #ifdef EMPTYEQUIV
              /* Note that induce_pseudocounts would have   *
               * been called in the imm_probmaker function. */
              fprintf(stdout,"Err (CHISQUARE_final_probs): \
induce_pseudocounts left empty devel history!\n");
              exit(EXIT_FAILURE);
              #else
              tmpshortcts[3]=tal->count3[model][j][k][l];
              tmpshortcts[2]=tal->count2[model][k][l];
              tmpshortcts[1]=tal->count1[model][l];
              tmpfprobs[3]=finalprob->prob4[model][j][k][l][base];
              tmpfprobs[2]=finalprob->prob3[model][k][l][base];
              tmpfprobs[1]=finalprob->prob2[model][l][base];
              for(order=3;order>0;--order) {
                if (tmpshortcts[order] > 0) {
                  finalprob->prob5[model][i][j][k][l][base]=
                    tmpfprobs[order];
                  break;
                }
              }
              #endif
            }
            else
              finalprob->prob5[model][i][j][k][l][base] =
                (lambda.lamhis4[i][j][k][l] *
                 rfreqs->MLhis4[i][j][k][l][base]) +
                ((CERTPROB - lambda.lamhis4[i][j][k][l]) *
                 finalprob->prob4[model][j][k][l][base]);
          } /* end terminal 3'-nt for loop */
          #ifdef PROBREPORT
          for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
            fprintf(stdout,"%u=%.4f,",base,
              finalprob->prob5[model][i][j][k][l][base]);
            pmf_test += finalprob->prob5[model][i][j][k][l][base];
          }
          fprintf(stdout,
            " Sum of prob's for history %u%u%u%u is %.2f ",i,j,k,l,pmf_test);
          if( fabs(pmf_test - CERTPROB) < MAXDIFF )
            fprintf(stdout,"valid\n");
          else
            fprintf(stdout,"problematic!\n");
          #endif
          #ifdef ADJUSTNULLS
          /* First, see if any null probability adjustments need be made */
          if ( (fabs(finalprob->prob5[model][i][j][k][l][0]) < MAXDIFF) ||
               (fabs(finalprob->prob5[model][i][j][k][l][1]) < MAXDIFF) ||
               (fabs(finalprob->prob5[model][i][j][k][l][2]) < MAXDIFF) ||
               (fabs(finalprob->prob5[model][i][j][k][l][3]) < MAXDIFF)
             ) {
            fprintf(stdout,"Making null probability adjustments...\n");
            pfix=MAX(PfixOMEGA,1.0 / tal->count4[model][i][j][k][l]);
            for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
              if( fabs(finalprob->prob5[model][i][j][k][l][base]) < MAXDIFF )
                finalprob->prob5[model][i][j][k][l][base] = pfix;
              else
                finalprob->prob5[model][i][j][k][l][base] =
                  (finalprob->prob5[model][i][j][k][l][base] *
                  (1 - (4 * pfix))) + pfix;

              fprintf(stdout,"%u=%.4f,",base,
                finalprob->prob5[model][i][j][k][l][base]);
              pmf_test += finalprob->prob5[model][i][j][k][l][base];
            } /* end for */
            fprintf(stdout,
              " Sum of prob's (revised) for history %u%u%u%u is %.2f ",
              i,j,k,l,pmf_test);
            if( fabs(pmf_test - CERTPROB) < MAXDIFF )
              fprintf(stdout,"valid\n");
            else
              fprintf(stdout,"problematic!\n");
          } /* end if (nulls exist in current pmf) */
          #endif
        }}}} /* end pretext for loop(s) */
        break;

      case(5) :
        /* Set weights of pentanucleotide pretexts */
        for(i=0;i<ALFSIZE;++i) {
        for(j=0;j<ALFSIZE;++j) {
        for(k=0;k<ALFSIZE;++k) {
        for(l=0;l<ALFSIZE;++l) {
        for(m=0;m<ALFSIZE;++m) {
          if (tal->count5[model][i][j][k][l][m] >= MINRELIABLECHISQCT)
            lambda.lamhis5[i][j][k][l][m]=CERTPROB;
          else if (tal->count5[model][i][j][k][l][m] == 0)
            continue;
          else {
            /* Populate contingency table */
            for(base=0;base<ALFSIZE;++base) {
              obs_matrix[0][base]=
                (double)tal->count6[model][i][j][k][l][m][base];
              obs_matrix[1][base]=finalprob->prob5[model][j][k][l][m][base] *
                (double)tal->count5[model][i][j][k][l][m];
            }
            /* determine weight */
            lambda.lamhis5[i][j][k][l][m]=
              calc_chisquareweight(obs_matrix,
              tal->count5[model][i][j][k][l][m]);
          }
        }}}}}
        /* Now record final, smoothed hexanucleotide probabilities */
        for(i=0;i<ALFSIZE;++i) {
        for(j=0;j<ALFSIZE;++j) {
        for(k=0;k<ALFSIZE;++k) {
        for(l=0;l<ALFSIZE;++l) {
        for(m=0;m<ALFSIZE;++m) {
          for(base=0;base<ALFSIZE;++base) {
            if (tal->count5[model][i][j][k][l][m] == 0) {
              #ifdef EMPTYEQUIV
              /* Note that induce_pseudocounts would have   *
               * been called in the imm_probmaker function. */
              fprintf(stdout,"Err (CHISQUARE_final_probs): \
induce_pseudocounts left empty devel history!\n");
              exit(EXIT_FAILURE);
              #else
              tmpshortcts[4]=tal->count4[model][j][k][l][m];
              tmpshortcts[3]=tal->count3[model][k][l][m];
              tmpshortcts[2]=tal->count2[model][l][m];
              tmpshortcts[1]=tal->count1[model][m];
              tmpfprobs[4]=finalprob->prob5[model][j][k][l][m][base];
              tmpfprobs[3]=finalprob->prob4[model][k][l][m][base];
              tmpfprobs[2]=finalprob->prob3[model][l][m][base];
              tmpfprobs[1]=finalprob->prob2[model][m][base];
              for(order=4;order>0;--order) {
                if (tmpshortcts[order] > 0) {
                  finalprob->prob6[model][i][j][k][l][m][base]=
                    tmpfprobs[order];
                  break;
                }
              }
              #endif
            }
            else
              finalprob->prob6[model][i][j][k][l][m][base] =
                (lambda.lamhis5[i][j][k][l][m] *
                 rfreqs->MLhis5[i][j][k][l][m][base]) +
                ((CERTPROB - lambda.lamhis5[i][j][k][l][m]) *
                 finalprob->prob5[model][j][k][l][m][base]);
          } /* end terminal 3'-nt for loop */
          #ifdef PROBREPORT
          for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
            fprintf(stdout,"%u=%.4f,",base,
              finalprob->prob6[model][i][j][k][l][m][base]);
            pmf_test += finalprob->prob6[model][i][j][k][l][m][base];
          }
          fprintf(stdout,
            " Sum of prob's for history %u%u%u%u%u is %.2f ",
            i,j,k,l,m,pmf_test);
          if( fabs(pmf_test - CERTPROB) < MAXDIFF )
            fprintf(stdout,"valid\n");
          else
            fprintf(stdout,"problematic!\n");
          #endif
          #ifdef ADJUSTNULLS
          if ( (fabs(finalprob->prob6[model][i][j][k][l][m][0]) < MAXDIFF) ||
               (fabs(finalprob->prob6[model][i][j][k][l][m][1]) < MAXDIFF) ||
               (fabs(finalprob->prob6[model][i][j][k][l][m][2]) < MAXDIFF) ||
               (fabs(finalprob->prob6[model][i][j][k][l][m][3]) < MAXDIFF)
             ) {
            fprintf(stdout,"Making null probability adjustments...\n");
            pfix=MAX(PfixOMEGA,1.0 / tal->count5[model][i][j][k][l][m]);
            for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
              if( fabs(finalprob->prob6[model][i][j][k][l][m][base]) < MAXDIFF)
                finalprob->prob6[model][i][j][k][l][m][base] = pfix;
              else
                finalprob->prob6[model][i][j][k][l][m][base] =
                  (finalprob->prob6[model][i][j][k][l][m][base] *
                  (1 - (4 * pfix))) + pfix;

              fprintf(stdout,"%u=%.4f,",base,
                finalprob->prob6[model][i][j][k][l][m][base]);
              pmf_test += finalprob->prob6[model][i][j][k][l][m][base];
            } /* end for */
            fprintf(stdout,
              " Sum of prob's (revised) for history %u%u%u%u%u is %.2f ",
              i,j,k,l,m,pmf_test);
            if( fabs(pmf_test - CERTPROB) < MAXDIFF )
              fprintf(stdout,"valid\n");
            else
              fprintf(stdout,"problematic!\n");
          } /* end if (nulls exist in current pmf) */
          #endif
        }}}}} /* end pretext for loop(s) */
        break;

      default:
        fprintf(stderr,"Err (CHISQUARE_final_probs): \
Invalid MAXORDER encountered!\n");
        exit(EXIT_FAILURE);
    }
  }

  return;
} /* end CHISQUARE_final_probs */

/* Given a populated contingency table, this function returns the  *
 * appropriate pretext weighting for the history it was called on. */
static double calc_chisquareweight(double obs_matrix[][ALFSIZE],
  int historyct)
{
  double rowsums[HYPONUM],             /* row totals                  */
         colsums[ALFSIZE],             /* column totals               */
         tablesum,                     /* table total                 */
         exp_matrix[HYPONUM][ALFSIZE], /* evil twin of obs_matrix,    *
                                        * stores expected frequencies */
         chisquare_stat,               /* chi-square statistic        */
         chi2confidence,               /* Complement of stat's P-val  */
         result;                       /* return to calling function  */
  int a,b;                             /* row/col iterator variables  */

  /* Compute rowsums */
  for(a=0;a<HYPONUM;++a)
    for(b=0,rowsums[a]=0.0;b<ALFSIZE;++b)
      rowsums[a] += obs_matrix[a][b];

  /* Compute colsums */
  for(b=0;b<ALFSIZE;++b)
    for(a=0,colsums[b]=0.0;a<HYPONUM;++a)
      colsums[b] += obs_matrix[a][b];

  /* Compute tablesum */
  for(a=0,tablesum=0.0;a<HYPONUM;++a)
    tablesum += rowsums[a];

  /* Populate expected contingency table */
  for(a=0;a<HYPONUM;++a) {
    for(b=0;b<ALFSIZE;++b) {
      exp_matrix[a][b]=(rowsums[a] * colsums[b]) / tablesum;
      if (fabs(exp_matrix[a][b]-5.0) > DBL_EPSILON &&
          exp_matrix[a][b] < 5.0)
        /* The chi**2 test is only valid if all cells in the expected *
         * matrix are .GT. 5.  If there are sufficiently few counts   *
         * that this is a real concern, then just return null.        */
        return(NULLPROB);
  }}

  /* Calculate chisquare statistic */
  chisquare_stat=0.0;
  for(a=0;a<HYPONUM;++a)
    for(b=0;b<ALFSIZE;++b)
      chisquare_stat +=
        pow((obs_matrix[a][b]-exp_matrix[a][b]),2) / exp_matrix[a][b];

  /* Compute associated chi**2 confidence value */
  chi2confidence=calc_chisq_conf(chisquare_stat);

  /* Are the distributions significantly different? */
  if(fabs(chi2confidence-CONFCUTOFF) > DBL_EPSILON &&
     chi2confidence < CONFCUTOFF)
    result=NULLPROB;
  else
    result=(chi2confidence*(double)historyct)/(double)MINRELIABLECHISQCT;

  return(result);
} /* end calc_chisquareweight */

/* This function computes the chi**2 confidence value (complement of *
 * P-value, i.e., area under left tail of the curve) to assist with  *
 * derivation of pretext weightings per 4th eqn in Salzberg et al.   *
 * The formulae are based on the classic Peizer and Pratt 1968 paper *
 * on approximating P-values of various statistical distributions:   *
 * Journal of the American Statistical Association, 63:1416-56       */
static double calc_chisq_conf(double chisq_stat)
{
  int DEGFREE = 3; /* 2x4 contingency table to compare frequency dists */
  double RATIOLOWBD = 0.5E-6,
                  z = chisq_stat / 2.0,
                  c = 1.0,
                  g = 1.0,
                  a = DEGFREE / 2.0,
                  d = a,
                 d3 = d + 2.0,
               d3sq = pow(d3,2);

  while(fabs((c/g)-RATIOLOWBD) > DBL_EPSILON &&
        (c/g) > RATIOLOWBD) {
    a += 1.0;
    c *= z/a;
    g += c;
  }

  g = g*exp((d*log(z))-
       (d3*(1.0/(12.0*d3sq)*(1.0-1.0/d3sq*(1.0/30.0-1.0/d3sq*
       (1.0/105.0-1.0/(140.0*d3sq))))))-(d3-0.5)*log(d3)+d3-z)*(d+1.0);

  return(g/sqrt(2.0*PI));
} /* end calc_chisq_conf */
