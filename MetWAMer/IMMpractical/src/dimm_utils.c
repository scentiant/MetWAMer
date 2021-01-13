/* dimm_utils.c
 * Michael E Sparks (mespar1@gmail.com)
 *
 * This file contains code specific to training deleted interpolated
 * Markov models.  Top-down and bottom-up deleted interpolated
 * Markov models are implemented.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dimm_utils.h"
#include "immpractical.h"
#include "sorting.h"
#ifdef ADJUSTNULLS
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

/* Local function/subroutine prototypes ************************************/

/* This function handles the bisection method; it's an  *
 * abstracted wrapper for the fxns4(TDDI|BUDI) fxns.    *
 * Set topdown to TRUE for top-down training, FALSE for *
 * bottom-up.                                           */
static double calc_root(int *iteration,double HI,double LO,countarrayT *talC,
  relfreqsT *rfreqs,immprobT* fprob,PsupiT *incsmooth,int model,int maxord,
  int order,ind_listT *sorted_listD,int startIND,int stopIND,int topdown);

/* Function to compute result of using candidate root given   *
 * the top-down training algorithm.                           *
 * Note: call with derivative == TRUE to compute second-order *
 *       fxn, derivative == FALSE for first-order fxn         */
static double fxns4TDDI(countarrayT *talC,relfreqsT* rfreqs,
  immprobT* fprob,int model,int order,ind_listT *sorted_listD,
  int startIND,int stopIND,double rootcandidate,int derivative);

/* Function to compute result of using candidate root given   *
 * the bottom-up training algorithm.                          *
 * Note: call with derivative == TRUE to compute second-order *
 *       fxn, derivative == FALSE for first-order fxn         */
static double fxns4BUDI(countarrayT *talC,relfreqsT *rfreqs,
  PsupiT *incsmooth,int model,int currmaxord,int currorder,
  ind_listT *sorted_listD,int startIND,int stopIND,
  double rootcandidate,int derivative);

/* Function definitions ******************************************************/

/* Subroutine to coordinate development of final, smoothed oligo *
 * likelihoods by the top-down or bottom-up training methods.    *
 * Set topdown to TRUE for top-down, FASLE, for bottom-up.       *
 * Note that ind_listT structures are defined in sorting.h       */
void DIMM_final_probs(immprobT *finalprob,relfreqsT *rfreqs,
  countarrayT *talD,countarrayT *talC,int model,int algo)
{
  ind_listT sorted_listD[(int)pow(ALFSIZE,MAXORDER)]; /* buffer-based   *
                         * space for for sorting oligos by counts       */
  int list_len,         /* length of seqment of sorted_listD to process */
      buckstart_IND=-1, /* Store index of start of the bucket           */
      buckstop_IND=-1,  /* Store index of stop of the bucket            */
      Buck_start_ct,    /* Bj-1 (variable for a given order)            */
      Buck_stop_ct,     /* Bj   (variable for a given order)            */
      BuckL,            /* BL   (fixed for a given order)               */
      tot_ctC,          /* Total number of C counts in a current bucket */
      iteration,        /* Keeps track of recursion depth               */
      order,            /* tracks 1..MAXORDER (k in Jelinek paper)      */
      innerorder,       /* bottom-up training inner loop (i in Jelinek) */
      topdown,          /* top-down or bottom-up truth value            */
      base,             /* 3'-most nt for an oligo, or 'v' in Jelinek   */
      i,j;              /* iterator variables                           */
  double lambdahat[MAXORDER+1]
                  [MAXORDER+1], /* store lambdahat estimates. The     *
           * first dimension corresponds to order (or 'k'), the       *
           * second to innerorder (or 'i'). Note that lambdahats      *
           * need not be associated with specific oligos, but rather  *
           * specific buckets. Also note that the second order is set *
           * to a constant index (0, arbitrarily chosen) in the case  *
           * of top-down training since there is no innerorder loop   *
           * to contend with there.                                   */
         tmprfreq; /* tmp var used to copy unsmooted ML estimates     */
  #ifndef EMPTYEQUIV
  int k,                     /* iterator variable                     */
      tmpshortcts[MAXORDER]; /* A temporary array to store counts of  *
           * successive 5'-truncated histories, to address training   *
           * for oligos whose histories did not occur in the training *
           * dataset                                                  */
  double tmpfprob[MAXORDER]; /* A temporary array to store final      *
           * smoothed probabilities of oligomers based on shorter     *
           * pretext lengths. This is related to the interpolation    *
           * that Jelinek describes in Section 2 for non-occurring    *
           * histories.                                               */
  /* Note that in the previous two arrays, the zeroth-element is not  *
   * used.                                                            */
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
  PsupiT incsmooth; /* store incrementally smoothed ML estimates */

  /* Determine requested training algorithm */
  if(algo == TDDIX)
    topdown=TRUE;
  else if(algo == BUDIX)
    topdown=FALSE;
  else {
    fprintf(stderr,
      "Err (DIMM_final_probs): Improper algorithm requested!\n");
    exit(EXIT_FAILURE);
  }
    
  /* history depth==0 requires no bucketing: "clearly, C0 is a singleton". *
   * Thus, we do not iterate over elements (histories) of the bucket.      */
  if(topdown == FALSE) /* (BUDI context) initialize P^0[v|h_0] per eqn 12a */
    for(base=0;base<ALFSIZE;++base)
      incsmooth.Psuphis0[0][base] =
        ((CERTPROB - XI) * rfreqs->MLhis0[base]) + (XI * EQUIPROB);
  iteration=1;
  lambdahat[0][0]=calc_root(&iteration,HBOUND,LBOUND,talC,rfreqs,
    finalprob,&incsmooth,model,0,0,sorted_listD,
    buckstart_IND,buckstop_IND,topdown);
  #ifdef ROOTREPORT
  fprintf(stdout,"  Root for this bucket is %.12f\n",lambdahat[0][0]);
  #endif

  /* Record these smoothed, final mononucleotide probabilities */
  for(base=0;base<ALFSIZE;++base) {
    if(topdown == TRUE) /* top-down (11c) */
      finalprob->prob1[model][base] =
        (lambdahat[0][0] * EQUIPROB) +
        ((CERTPROB - lambdahat[0][0]) * rfreqs->MLhis0[base]);
    else /* bottom-up (12d) */
      finalprob->prob1[model][base] =
        (lambdahat[0][0] * incsmooth.Psuphis0[0][base]) +
        ((CERTPROB - lambdahat[0][0]) * EQUIPROB);
  } /* end base for */

  #ifdef PROBREPORT
  /* Have we computed a valid PMF? */
  for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
    fprintf(stdout,"%u=%.4f,",base,finalprob->prob1[model][base]);
    pmf_test += finalprob->prob1[model][base];
  }
  fprintf(stdout," Sum of prob's for null history (model %u) is %.2f ",
    model,pmf_test);
  if( fabs(pmf_test - CERTPROB) < MAXDIFF )
    fprintf(stdout,"valid\n");
  else
    fprintf(stdout,"problematic!\n");
  #endif

  /* Other orders can be handled more generally */
  for(order=1;order<=MAXORDER;++order) {
    /* Input the relevant data into the array of ind_listT elements. */
    copy_data(order,model,talD,sorted_listD);

    /* Calculate how many distinct elements are used in the tally array *
     * for a particular order; this is needed for sorting purposes.     */
    list_len=(int)pow(ALFSIZE,(double)order);

    /* Sort the array of ind_listT elements *
     * Function definition in sorting.c     */
    quicksort(sorted_listD,0,list_len-1);

    /* Scroll past the non-occurring histories. (0 counts)     *
     * Note in proof: If using the induce_pseudocounts         *
     * function, there should never be a non-occurring history *
     * in the development data counts (would have been called  *
     * in the imm_probmaker function)                               */
    for(i=0;sorted_listD[i].ct_pret==0;++i)
      ;
    #ifdef EMPTYEQUIV
    if(i!=0) {
      /* Note that induce_pseudocounts would have   *
       * been called in the imm_probmaker function. */
      fprintf(stdout,"Err (DIMM_final_probs): \
induce_pseudocounts left empty devel history!\n");
      exit(EXIT_FAILURE);
    }
    #else
    /* In section 2 of Jelinek's report: "In practice, at least some    *
     * such histories never occur in the development set.  In these     *
     * cases, we set"...yadda, yadda, yadda.  This means we set the     *
     * final smoothed probability for those oligomers whose historical  *
     * pretexts were not encountered in the training data to the final  *
     * probability of one of its 5'-truncated suboligos, namely the one *
     * whose prefix (history, suboligo less 3'-terminal nt) did occur   *
     * in the training data.                                            */
    if(i!=0) {
      for(j=0;j<i;++j) { /* consider empty histories for current order   */
        for(base=0;base<ALFSIZE;++base) { /* If the history didn't occur *
          * in the first place, then it's impossible for the history +   *
          * 3'-terminal nucleotide oligo to occur.                       */
          switch(order) {
            case 1 :
              /* Let's be realistic here.  For this branch to execute,   *
               * one of the ALFSIZE possible mononucleotide histories    *
               * was not encountered in training.  That is unacceptable. */
              fprintf(stdout,"Err (DIMM_final_probs): \
Training data unrealistically sparse!\n");
              exit(EXIT_FAILURE);
              break;
            case 2 :
              /* In this case, there is only one possible shorter history,   *
               * a single nucleotide (which *must* occur).  This is also     *
               * quite unacceptable, but...we will set the final probability *
               * for the current oligo (based on the non-occurring history)  *
               * to the only possible compensatory one that can be assigned  *
               * to it.                                                      */
              finalprob->prob3[model][sorted_listD[j].key[0]]
                                     [sorted_listD[j].key[1]][base] =
                finalprob->prob2[model][sorted_listD[j].key[1]][base];
              break;
            case 3 :
              /* Here, we are dealing with an empty history, so we must, *
               * for each of all four possible oligos with this history, *
               * set their final smoothed probabilities to that based on *
               * the first in a series of 5'-truncated oligos whose      *
               * pretext count is .GT. zero in the training data.        */
              tmpshortcts[2]=talD->count2[model][sorted_listD[j].key[1]]
                                                [sorted_listD[j].key[2]];
              tmpshortcts[1]=talD->count1[model][sorted_listD[j].key[2]];
              tmpfprob[2]=finalprob->prob3[model][sorted_listD[j].key[1]]
                                                 [sorted_listD[j].key[2]]
                                                 [base];
              tmpfprob[1]=finalprob->prob2[model][sorted_listD[j].key[2]]
                                                 [base];
              for(k=order-1;k>0;--k) {
                if (tmpshortcts[k] > 0) {
                  finalprob->prob4[model][sorted_listD[j].key[0]]
                                         [sorted_listD[j].key[1]]
                                         [sorted_listD[j].key[2]][base]=
                    tmpfprob[k];
                  break;
                }
              }
              break;
            case 4 :
              tmpshortcts[3]=talD->count3[model][sorted_listD[j].key[1]]
                                                [sorted_listD[j].key[2]]
                                                [sorted_listD[j].key[3]];
              tmpshortcts[2]=talD->count2[model][sorted_listD[j].key[2]]
                                                [sorted_listD[j].key[3]];
              tmpshortcts[1]=talD->count1[model][sorted_listD[j].key[3]];
              tmpfprob[3]=finalprob->prob4[model][sorted_listD[j].key[1]]
                                                 [sorted_listD[j].key[2]]
                                                 [sorted_listD[j].key[3]]
                                                 [base];
              tmpfprob[2]=finalprob->prob3[model][sorted_listD[j].key[2]]
                                                 [sorted_listD[j].key[3]]
                                                 [base];
              tmpfprob[1]=finalprob->prob2[model][sorted_listD[j].key[3]]
                                                 [base];
              for(k=order-1;k>0;--k) {
                if (tmpshortcts[k] > 0) {
                  finalprob->prob5[model][sorted_listD[j].key[0]]
                                         [sorted_listD[j].key[1]]
                                         [sorted_listD[j].key[2]]
                                         [sorted_listD[j].key[3]][base]=
                    tmpfprob[k];
                  break;
                }
              }
              break;
            case 5 :
              tmpshortcts[4]=talD->count4[model][sorted_listD[j].key[1]]
                                                [sorted_listD[j].key[2]]
                                                [sorted_listD[j].key[3]]
                                                [sorted_listD[j].key[4]];
              tmpshortcts[3]=talD->count3[model][sorted_listD[j].key[2]]
                                                [sorted_listD[j].key[3]]
                                                [sorted_listD[j].key[4]];
              tmpshortcts[2]=talD->count2[model][sorted_listD[j].key[3]]
                                                [sorted_listD[j].key[4]];
              tmpshortcts[1]=talD->count1[model][sorted_listD[j].key[4]];
              tmpfprob[4]=finalprob->prob5[model][sorted_listD[j].key[1]]
                                                 [sorted_listD[j].key[2]]
                                                 [sorted_listD[j].key[3]]
                                                 [sorted_listD[j].key[4]]
                                                 [base];
              tmpfprob[3]=finalprob->prob4[model][sorted_listD[j].key[2]]
                                                 [sorted_listD[j].key[3]]
                                                 [sorted_listD[j].key[4]]
                                                 [base];
              tmpfprob[2]=finalprob->prob3[model][sorted_listD[j].key[3]]
                                                 [sorted_listD[j].key[4]]
                                                 [base];
              tmpfprob[1]=finalprob->prob2[model][sorted_listD[j].key[4]]
                                                 [base];
              for(k=order-1;k>0;--k) {
                if (tmpshortcts[k] > 0) {
                  finalprob->prob6[model][sorted_listD[j].key[0]]
                                         [sorted_listD[j].key[1]]
                                         [sorted_listD[j].key[2]]
                                         [sorted_listD[j].key[3]]
                                         [sorted_listD[j].key[4]][base]=
                    tmpfprob[k];
                  break;
                }
              }
              break;
            default :
              fprintf(stderr,
                "Err (DIMM_final_probs): Invalid case encountered!\n");
              exit(EXIT_FAILURE);
          } /* end switch */
        } /* end base for */
      } /* end empty history index for */
    } /* end if i!=0 */
    #endif
    buckstart_IND=i; /* we will need this later */

    /* Partition histories into buckets. */
    Buck_start_ct=sorted_listD[buckstart_IND].ct_pret;
    BuckL=sorted_listD[list_len-1].ct_pret+1;

    #ifdef BUCKBOUNDCHECK
    fprintf(stdout,"\nOverall bucket bounds (model %u, order %u) \
start_ct %u (index %u) stop_ct %u (index %u)\n",
      model,order,sorted_listD[buckstart_IND].ct_pret,
      buckstart_IND,BuckL-1,list_len-1);
    #endif

    /* process all histories of interest, one bucket at a time */
    for(tot_ctC=0;i<list_len;tot_ctC=0) {
      Buck_stop_ct=(int)ceil(BUCKRATIO*(double)Buck_start_ct);

      /* Pool the last two buckets if the final one is too narrow. */
      if( ((int)ceil(BUCKRATIO*(double)Buck_stop_ct)) >= BuckL )
        Buck_stop_ct=BuckL;

      /* If not processing last bucket, verify that sufficient counts occur *
       * in the held-out data or extend Buck_stop_ct (and index i) if not.  */
      if(Buck_stop_ct != BuckL) {
        while( i < list_len &&
          (sorted_listD[i].ct_pret <= Buck_stop_ct || tot_ctC < RELIABLECT)) {

          /* Update tot_ctC (counts of histories from H set per bucket) */
          switch (order) {
            case 1 :
              tot_ctC += talC->count1[model][sorted_listD[i].key[0]];
              break;
            case 2 :
              tot_ctC += talC->count2[model][sorted_listD[i].key[0]]
                                            [sorted_listD[i].key[1]];
              break;
            case 3 :
              tot_ctC += talC->count3[model][sorted_listD[i].key[0]]
                                            [sorted_listD[i].key[1]]
                                            [sorted_listD[i].key[2]];
              break;
            case 4 :
              tot_ctC += talC->count4[model][sorted_listD[i].key[0]]
                                            [sorted_listD[i].key[1]]
                                            [sorted_listD[i].key[2]]
                                            [sorted_listD[i].key[3]];
              break;
            case 5 :
              tot_ctC += talC->count5[model][sorted_listD[i].key[0]]
                                            [sorted_listD[i].key[1]]
                                            [sorted_listD[i].key[2]]
                                            [sorted_listD[i].key[3]]
                                            [sorted_listD[i].key[4]];
              break;
            default :
              fprintf(stderr,
                "Err (DIMM_final_probs): Invalid case encountered!\n");
              exit(EXIT_FAILURE);
          }
          ++i;
        } /* end while */
      } /* end if */
      else {
        i=list_len;
        /* For the record, compute how many counts occur in the *
         * held-out data for the last bucket of the order.      */
        for(j=buckstart_IND,tot_ctC=0;j<i;++j) {
          switch(order) {
            case 1 :
              tot_ctC += talC->count1[model][sorted_listD[j].key[0]];
              break;
            case 2 :
              tot_ctC += talC->count2[model][sorted_listD[j].key[0]]
                                            [sorted_listD[j].key[1]];
              break;
            case 3 :
              tot_ctC += talC->count3[model][sorted_listD[j].key[0]]
                                            [sorted_listD[j].key[1]]
                                            [sorted_listD[j].key[2]];
              break;
            case 4 :
              tot_ctC += talC->count4[model][sorted_listD[j].key[0]]
                                            [sorted_listD[j].key[1]]
                                            [sorted_listD[j].key[2]]
                                            [sorted_listD[j].key[3]];
              break;
            case 5 :
              tot_ctC += talC->count5[model][sorted_listD[j].key[0]]
                                            [sorted_listD[j].key[1]]
                                            [sorted_listD[j].key[2]]
                                            [sorted_listD[j].key[3]]
                                            [sorted_listD[j].key[4]];
              break;
            default :
              fprintf(stderr,
                "Err (DIMM_final_probs): Invalid case encountered!\n");
              exit(EXIT_FAILURE);
          } /* end switch */
        } /* end for */
      } /* end else */

      /* At this point, we have a bucket to process.    *
       * Note: buckstart_IND is prepared at this point. */
      buckstop_IND=i-1;
      #ifdef BUCKBOUNDCHECK
      fprintf(stdout,"Bucket start_ct %u (index %u) stop_ct %u (index %u)\n",
        sorted_listD[buckstart_IND].ct_pret,buckstart_IND,
        sorted_listD[buckstop_IND].ct_pret,buckstop_IND);
      if(tot_ctC < RELIABLECT)
        fprintf(stdout,"  Warning: Insufficient held-out counts! (%i < %i)\n",
          tot_ctC,RELIABLECT);
      else
        fprintf(stdout,"  held-out counts: (%i >= %i)\n",tot_ctC,RELIABLECT);
      #endif

      /* So, process the bucket. */
      iteration=1;
      /* The top-down method doesn't have an *
       * internal order loop to fuss with.   */
      if(topdown == TRUE) {
        /* First, compute the lambdahat estimate for this bucket */
        iteration=1;
        lambdahat[order][0]=calc_root(&iteration,HBOUND,LBOUND,talC,rfreqs,
          finalprob,&incsmooth,model,order,0,sorted_listD,
          buckstart_IND,buckstop_IND,topdown);
        #ifdef ROOTREPORT
        fprintf(stdout,"  Root for this bucket is %.12f\n",
          lambdahat[order][0]);
        #endif

        /* Then, store final, smoothed oligo probabilities, *
         * for all histories in the bucket, per 11c.        */
        for(j=buckstart_IND;j<=buckstop_IND;++j) {
          for(base=0;base<ALFSIZE;++base) {
            /* This switch statement is well-defined, as order ranges  *
             * from 1 upto MAXORDER in the enclosing loop.  Therefore, *
             * the DI.TD value for a shorter pretext length will have  *
             * been computed prior to it's being needed in the         *
             * calculations below.                                     */
            switch(order) {
              case 1 :
                finalprob->prob2[model][sorted_listD[j].key[0]][base]=
                  (lambdahat[order][0] *
                  finalprob->prob1[model][base]) +
                  ((CERTPROB - lambdahat[order][0]) *
                  rfreqs->MLhis1[sorted_listD[j].key[0]][base]);
                break;
              case 2 :
                finalprob->prob3[model][sorted_listD[j].key[0]]
                                       [sorted_listD[j].key[1]][base]=
                  (lambdahat[order][0] *
                  finalprob->prob2[model][sorted_listD[j].key[1]][base]) +
                  ((CERTPROB - lambdahat[order][0]) *
                  rfreqs->MLhis2[sorted_listD[j].key[0]]
                                [sorted_listD[j].key[1]][base]);
                break;
              case 3 :
                finalprob->prob4[model][sorted_listD[j].key[0]]
                                       [sorted_listD[j].key[1]]
                                       [sorted_listD[j].key[2]][base]=
                  (lambdahat[order][0] *
                  finalprob->prob3[model][sorted_listD[j].key[1]]
                                         [sorted_listD[j].key[2]][base]) +
                  ((CERTPROB - lambdahat[order][0]) *
                  rfreqs->MLhis3[sorted_listD[j].key[0]]
                                [sorted_listD[j].key[1]]
                                [sorted_listD[j].key[2]][base]);
                break;
              case 4 :
                finalprob->prob5[model][sorted_listD[j].key[0]]
                                       [sorted_listD[j].key[1]]
                                       [sorted_listD[j].key[2]]
                                       [sorted_listD[j].key[3]][base]=
                  (lambdahat[order][0] *
                  finalprob->prob4[model][sorted_listD[j].key[1]]
                                         [sorted_listD[j].key[2]]
                                         [sorted_listD[j].key[3]][base]) +
                  ((CERTPROB - lambdahat[order][0]) *
                  rfreqs->MLhis4[sorted_listD[j].key[0]]
                                [sorted_listD[j].key[1]]
                                [sorted_listD[j].key[2]]
                                [sorted_listD[j].key[3]][base]);
                break;
              case 5 :
                finalprob->prob6[model][sorted_listD[j].key[0]]
                                       [sorted_listD[j].key[1]]
                                       [sorted_listD[j].key[2]]
                                       [sorted_listD[j].key[3]]
                                       [sorted_listD[j].key[4]][base]=
                  (lambdahat[order][0] *
                  finalprob->prob5[model][sorted_listD[j].key[1]]
                                         [sorted_listD[j].key[2]]
                                         [sorted_listD[j].key[3]]
                                         [sorted_listD[j].key[4]][base]) +
                  ((CERTPROB - lambdahat[order][0]) *
                  rfreqs->MLhis5[sorted_listD[j].key[0]]
                                [sorted_listD[j].key[1]]
                                [sorted_listD[j].key[2]]
                                [sorted_listD[j].key[3]]
                                [sorted_listD[j].key[4]][base]);
                break;
              default :
                fprintf(stderr,
                  "Err (DIMM_final_probs): Invalid case encountered!\n");
                exit(EXIT_FAILURE);
            } /* end switch */
          } /* end base for */
        } /* end bucket index for */
        /* And that's all for TDDI-smoothed probs */
      } /* end if */
      /* The bottom-up algorithm requires an inner  *
       * loop that goes from current order downto 0 */
      else if(topdown == FALSE) {
        /* Init Psupi data structure for the current order [1..MAXORDER], *
         * for ALL histories in the current bucket, per eqn 12a           */
        for(j=buckstart_IND;j<=buckstop_IND;++j) {
          for(base=0;base<ALFSIZE;++base) {
            switch(order) {
              case 1 :
                incsmooth.Psuphis1[order]
                                  [sorted_listD[j].key[0]][base] =
                  ((CERTPROB - XI) *
                  rfreqs->MLhis1[sorted_listD[j].key[0]]
                                [base]) + (XI * EQUIPROB);
                break;
              case 2 :
                incsmooth.Psuphis2[order]
                                  [sorted_listD[j].key[0]]
                                  [sorted_listD[j].key[1]][base] =
                  ((CERTPROB - XI) *
                  rfreqs->MLhis2[sorted_listD[j].key[0]]
                                [sorted_listD[j].key[1]]
                                [base]) + (XI * EQUIPROB);
                break;
              case 3 :
                incsmooth.Psuphis3[order]
                                  [sorted_listD[j].key[0]]
                                  [sorted_listD[j].key[1]]
                                  [sorted_listD[j].key[2]][base] =
                  ((CERTPROB - XI) *
                  rfreqs->MLhis3[sorted_listD[j].key[0]]
                                [sorted_listD[j].key[1]]
                                [sorted_listD[j].key[2]]
                                [base]) + (XI * EQUIPROB);
                break;
              case 4 :
                incsmooth.Psuphis4[order]
                                  [sorted_listD[j].key[0]]
                                  [sorted_listD[j].key[1]]
                                  [sorted_listD[j].key[2]]
                                  [sorted_listD[j].key[3]][base] =
                  ((CERTPROB - XI) *
                  rfreqs->MLhis4[sorted_listD[j].key[0]]
                                [sorted_listD[j].key[1]]
                                [sorted_listD[j].key[2]]
                                [sorted_listD[j].key[3]]
                                [base]) + (XI * EQUIPROB);
                break;
              case 5 :
                incsmooth.Psuphis5[order]
                                  [sorted_listD[j].key[0]]
                                  [sorted_listD[j].key[1]]
                                  [sorted_listD[j].key[2]]
                                  [sorted_listD[j].key[3]]
                                  [sorted_listD[j].key[4]][base] =
                  ((CERTPROB - XI) *
                  rfreqs->MLhis5[sorted_listD[j].key[0]]
                                [sorted_listD[j].key[1]]
                                [sorted_listD[j].key[2]]
                                [sorted_listD[j].key[3]]
                                [sorted_listD[j].key[4]]
                                [base]) + (XI * EQUIPROB);
                break;
              default :
                fprintf(stderr,
                  "Err (DIMM_final_probs): Invalid order encountered!\n");
                exit(EXIT_FAILURE);
            } /* end switch */
          } /* end base for */
        } /* end bucket index for */

        /* Now, with all Psupi data points initialized for the max order *
         * (k) of this bucket, we can calculate roots, i.e., set         *
         * lambdahat, and update Psupi for next-lower innerorders (i's)  */
        for(innerorder=order;innerorder>=0;--innerorder) { /* "i" loop */

          /* Everything is ready to compute lambdahat for the *
           * current bucket and order,innerorder combination. */
          iteration=1;
          lambdahat[order][innerorder]=
            calc_root(&iteration,HBOUND,LBOUND,talC,rfreqs,
            finalprob,&incsmooth,model,order,innerorder,sorted_listD,
            buckstart_IND,buckstop_IND,topdown);
          #ifdef ROOTREPORT
          fprintf(stdout,"  Root for this bucket is %.12f\n",
            lambdahat[order][innerorder]);
          #endif

          /* Update Psupi data structure for the next-lower innerorder 12c */
          for(j=buckstart_IND;j<=buckstop_IND;++j) {
            for(base=0;base<ALFSIZE;++base) {
              /* First, record unsmoothed ML estimate (for oligos with *
               * histories of the current innerorder-1 length!) in a   *
               * local variable.  This prevents having to code nested  *
               * switches in the switch for *order* below. (key-in for *
               * appropriate ML estimates depends only on *innerorder* */
              switch(innerorder-1) { /* or 'i-1' */
                case -1 :
                  /* Note that, while order==0 (k==0) was handled as     *
                   * a special one-off case earlier in this function,    *
                   * innerorder==0 (i==0) is to be expected in this loop */
                  tmprfreq=EQUIPROB;
                  break;
                case 0 :
                  tmprfreq=rfreqs->MLhis0[base];
                  break;
                case 1 :
                  /* index the history by [order-X] because shorter   *
                   * histories, in order to maintain connect with the *
                   * 3'-terminal nucleotide (base), must be shortened *
                   * from the 5' end.                                 */
                  tmprfreq=rfreqs->MLhis1[sorted_listD[j].key[order-1]]
                                         [base];
                  break;
                case 2 :
                  tmprfreq=rfreqs->MLhis2[sorted_listD[j].key[order-2]]
                                         [sorted_listD[j].key[order-1]]
                                         [base];
                  break;
                case 3 :
                  tmprfreq=rfreqs->MLhis3[sorted_listD[j].key[order-3]]
                                         [sorted_listD[j].key[order-2]]
                                         [sorted_listD[j].key[order-1]]
                                         [base];
                  break;
                case 4 :
                  tmprfreq=rfreqs->MLhis4[sorted_listD[j].key[order-4]]
                                         [sorted_listD[j].key[order-3]]
                                         [sorted_listD[j].key[order-2]]
                                         [sorted_listD[j].key[order-1]]
                                         [base];
                  break;
                default :
                  fprintf(stderr,"Err (DIMM_final_probs): \
Invalid innerorder encountered!\n");
                  exit(EXIT_FAILURE);
              } /* end switch */

              /* Then, actually update incsmooth for the next-lower   *
               * historical context (remember, init for incsmooth was *
               * addressed earlier) according to eqn 12c. (key-in     *
               * depends exclusively on the *order* var).             */
              if (innerorder > 0) {
                switch(order) {
                  case 1 :
                    incsmooth.Psuphis1[innerorder-1]
                                      [sorted_listD[j].key[0]][base] =
                      (lambdahat[order][innerorder] *
                      incsmooth.Psuphis1[innerorder]
                                        [sorted_listD[j].key[0]]
                                        [base]) +
                      ((CERTPROB - lambdahat[order][innerorder]) *
                      tmprfreq);
                    break;
                  case 2 :
                    incsmooth.Psuphis2[innerorder-1]
                                      [sorted_listD[j].key[0]]
                                      [sorted_listD[j].key[1]][base] =
                      (lambdahat[order][innerorder] *
                      incsmooth.Psuphis2[innerorder]
                                        [sorted_listD[j].key[0]]
                                        [sorted_listD[j].key[1]]
                                        [base]) +
                      ((CERTPROB - lambdahat[order][innerorder]) *
                      tmprfreq);
                    break;
                  case 3 :
                    incsmooth.Psuphis3[innerorder-1]
                                      [sorted_listD[j].key[0]]
                                      [sorted_listD[j].key[1]]
                                      [sorted_listD[j].key[2]][base] =
                      (lambdahat[order][innerorder] *
                      incsmooth.Psuphis3[innerorder]
                                        [sorted_listD[j].key[0]]
                                        [sorted_listD[j].key[1]]
                                        [sorted_listD[j].key[2]]
                                        [base]) +
                      ((CERTPROB - lambdahat[order][innerorder]) *
                      tmprfreq);
                    break;
                  case 4 :
                    incsmooth.Psuphis4[innerorder-1]
                                      [sorted_listD[j].key[0]]
                                      [sorted_listD[j].key[1]]
                                      [sorted_listD[j].key[2]]
                                      [sorted_listD[j].key[3]][base] =
                      (lambdahat[order][innerorder] *
                      incsmooth.Psuphis4[innerorder]
                                        [sorted_listD[j].key[0]]
                                        [sorted_listD[j].key[1]]
                                        [sorted_listD[j].key[2]]
                                        [sorted_listD[j].key[3]]
                                        [base]) +
                      ((CERTPROB - lambdahat[order][innerorder]) *
                      tmprfreq);
                    break;
                  case 5 :
                    incsmooth.Psuphis5[innerorder-1]
                                      [sorted_listD[j].key[0]]
                                      [sorted_listD[j].key[1]]
                                      [sorted_listD[j].key[2]]
                                      [sorted_listD[j].key[3]]
                                      [sorted_listD[j].key[4]][base] =
                      (lambdahat[order][innerorder] *
                      incsmooth.Psuphis5[innerorder]
                                        [sorted_listD[j].key[0]]
                                        [sorted_listD[j].key[1]]
                                        [sorted_listD[j].key[2]]
                                        [sorted_listD[j].key[3]]
                                        [sorted_listD[j].key[4]]
                                        [base]) +
                      ((CERTPROB - lambdahat[order][innerorder]) *
                      tmprfreq);
                    break;
                  default :
                    fprintf(stderr,"Err (DIMM_final_probs): \
Invalid order encountered!\n");
                    exit(EXIT_FAILURE);
                } /* end switch */
              }
              else { /* innerorder == 0 */
                /* By way of eqn 12d, we can now compute the final, *
                 * smoothed probability for the oligomer, which we  *
                 * will go ahead and store in the immprobT data     *
                 * structure instance.                              */
                switch(order) {
                  case 1 :
                    finalprob->prob2[model][sorted_listD[j].key[0]][base]=
                      (lambdahat[order][0] *
                      incsmooth.Psuphis1[0][sorted_listD[j].key[0]][base]) +
                      ((CERTPROB - lambdahat[order][0]) * tmprfreq);
                    break;
                  case 2 :
                    finalprob->prob3[model][sorted_listD[j].key[0]]
                                           [sorted_listD[j].key[1]][base]=
                      (lambdahat[order][0] *
                      incsmooth.Psuphis2[0][sorted_listD[j].key[0]]
                                           [sorted_listD[j].key[1]][base]) +
                      ((CERTPROB - lambdahat[order][0]) * tmprfreq);
                    break;
                  case 3 :
                    finalprob->prob4[model][sorted_listD[j].key[0]]
                                           [sorted_listD[j].key[1]]
                                           [sorted_listD[j].key[2]][base]=
                      (lambdahat[order][0] *
                      incsmooth.Psuphis3[0][sorted_listD[j].key[0]]
                                           [sorted_listD[j].key[1]]
                                           [sorted_listD[j].key[2]][base]) +
                      ((CERTPROB - lambdahat[order][0]) * tmprfreq);
                    break;
                  case 4 :
                    finalprob->prob5[model][sorted_listD[j].key[0]]
                                           [sorted_listD[j].key[1]]
                                           [sorted_listD[j].key[2]]
                                           [sorted_listD[j].key[3]][base]=
                      (lambdahat[order][0] *
                      incsmooth.Psuphis4[0][sorted_listD[j].key[0]]
                                           [sorted_listD[j].key[1]]
                                           [sorted_listD[j].key[2]]
                                           [sorted_listD[j].key[3]][base]) +
                      ((CERTPROB - lambdahat[order][0]) * tmprfreq);
                    break;
                  case 5 :
                    finalprob->prob6[model][sorted_listD[j].key[0]]
                                           [sorted_listD[j].key[1]]
                                           [sorted_listD[j].key[2]]
                                           [sorted_listD[j].key[3]]
                                           [sorted_listD[j].key[4]][base]=
                      (lambdahat[order][0] *
                      incsmooth.Psuphis5[0][sorted_listD[j].key[0]]
                                           [sorted_listD[j].key[1]]
                                           [sorted_listD[j].key[2]]
                                           [sorted_listD[j].key[3]]
                                           [sorted_listD[j].key[4]][base]) +
                      ((CERTPROB - lambdahat[order][0]) * tmprfreq);
                    break;
                  default :
                    fprintf(stderr,"Err (DIMM_final_probs): \
Invalid order encountered!\n");
                    exit(EXIT_FAILURE);
                } /* end switch */
              } /* end innerorder == 0 else */
            } /* end base for */
          } /* end bucket index for */
        } /* end innerorder for */
      } /* end else */
      else {
        fprintf(stderr,"Err (DIMM_final_probs): \
Improper algorithm requested!\n");
        exit(EXIT_FAILURE);
      }

      /* That's all for this bucket. Optional testing block follows */
      #ifdef PROBREPORT
      /* scan through bucket to assure PMF status */
      for(j=buckstart_IND;j<=buckstop_IND;++j) { /* pretext/history */
        /* Have we computed a valid PMF? */
        switch(order) {
          case 1 :
            for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
              fprintf(stdout,"%u=%.4f,",base,
                finalprob->prob2[model][sorted_listD[j].key[0]]
                                       [base]);
              pmf_test += finalprob->prob2[model][sorted_listD[j].key[0]]
                                                 [base];
            }
            fprintf(stdout,
              " Sum of prob's for history %u is %.2f ",
              sorted_listD[j].key[0],pmf_test);
            if( fabs(pmf_test - CERTPROB) < MAXDIFF )
              fprintf(stdout,"valid\n");
            else
              fprintf(stdout,"problematic!\n");
            break;
          case 2 :
            for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
              fprintf(stdout,"%u=%.4f,",base,
                finalprob->prob3[model][sorted_listD[j].key[0]]
                                       [sorted_listD[j].key[1]]
                                       [base]);
              pmf_test += finalprob->prob3[model][sorted_listD[j].key[0]]
                                                 [sorted_listD[j].key[1]]
                                                 [base];
            }
            fprintf(stdout,
              " Sum of prob's for history %u%u is %.2f ",
              sorted_listD[j].key[0],
              sorted_listD[j].key[1],pmf_test);
            if( fabs(pmf_test - CERTPROB) < MAXDIFF )
              fprintf(stdout,"valid\n");
            else
              fprintf(stdout,"problematic!\n");
            break;
          case 3 :
            for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
              fprintf(stdout,"%u=%.4f,",base,
                finalprob->prob4[model][sorted_listD[j].key[0]]
                                       [sorted_listD[j].key[1]]
                                       [sorted_listD[j].key[2]]
                                       [base]);
              pmf_test += finalprob->prob4[model][sorted_listD[j].key[0]]
                                                 [sorted_listD[j].key[1]]
                                                 [sorted_listD[j].key[2]]
                                                 [base];
            }
            fprintf(stdout,
              " Sum of prob's for history %u%u%u is %.2f ",
              sorted_listD[j].key[0],
              sorted_listD[j].key[1],
              sorted_listD[j].key[2],pmf_test);
            if( fabs(pmf_test - CERTPROB) < MAXDIFF )
              fprintf(stdout,"valid\n");
            else
              fprintf(stdout,"problematic!\n");
            break;
          case 4 :
            for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
              fprintf(stdout,"%u=%.4f,",base,
                finalprob->prob5[model][sorted_listD[j].key[0]]
                                       [sorted_listD[j].key[1]]
                                       [sorted_listD[j].key[2]]
                                       [sorted_listD[j].key[3]]
                                       [base]);
              pmf_test += finalprob->prob5[model][sorted_listD[j].key[0]]
                                                 [sorted_listD[j].key[1]]
                                                 [sorted_listD[j].key[2]]
                                                 [sorted_listD[j].key[3]]
                                                 [base];
            }
            fprintf(stdout,
              " Sum of prob's for history %u%u%u%u is %.2f ",
              sorted_listD[j].key[0],
              sorted_listD[j].key[1],
              sorted_listD[j].key[2],
              sorted_listD[j].key[3],pmf_test);
            if( fabs(pmf_test - CERTPROB) < MAXDIFF )
              fprintf(stdout,"valid\n");
            else
              fprintf(stdout,"problematic!\n");
            break;
          case 5 :
            for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
              fprintf(stdout,"%u=%.4f,",base,
                finalprob->prob6[model][sorted_listD[j].key[0]]
                                       [sorted_listD[j].key[1]]
                                       [sorted_listD[j].key[2]]
                                       [sorted_listD[j].key[3]]
                                       [sorted_listD[j].key[4]]
                                       [base]);
              pmf_test += finalprob->prob6[model][sorted_listD[j].key[0]]
                                                 [sorted_listD[j].key[1]]
                                                 [sorted_listD[j].key[2]]
                                                 [sorted_listD[j].key[3]]
                                                 [sorted_listD[j].key[4]]
                                                 [base];
            }
            fprintf(stdout,
              " Sum of prob's for history %u%u%u%u%u is %.2f ",
              sorted_listD[j].key[0],
              sorted_listD[j].key[1],
              sorted_listD[j].key[2],
              sorted_listD[j].key[3],
              sorted_listD[j].key[4],pmf_test);
            if( fabs(pmf_test - CERTPROB) < MAXDIFF )
              fprintf(stdout,"valid\n");
            else
              fprintf(stdout,"problematic!\n");
            break;
          default :
            fprintf(stderr,
              "Err (DIMM_final_probs): Invalid order encountered!\n");
            exit(EXIT_FAILURE);
        } /* end switch */
      } /* end bucket index for */
      #endif
      #ifdef ADJUSTNULLS
      /* Potentially adjust prob mass functions *
       * containing null probabilities          */
      for(j=buckstart_IND;j<=buckstop_IND;++j) { /* pretext/history */
        switch(order) {
          case 1 :
            /* First, see if any null probability adjustments need be made */
            if ( (fabs(finalprob->prob2[model][sorted_listD[j].key[0]][0])
                   < MAXDIFF) ||
                 (fabs(finalprob->prob2[model][sorted_listD[j].key[0]][1])
                   < MAXDIFF) ||
                 (fabs(finalprob->prob2[model][sorted_listD[j].key[0]][2])
                   < MAXDIFF) ||
                 (fabs(finalprob->prob2[model][sorted_listD[j].key[0]][3])
                   < MAXDIFF)
               ) {
              /* This can occur if (provided EMPTYEQUIV is not defined),   *
               * though an oligo's history did occur .GT. zero times, the  *
               * history + some 3'-terminal nucleotide was not observed in *
               * the training data.                                        */
              fprintf(stdout,"Making null probability adjustments...\n");
              pfix=(double)MAX(PfixOMEGA,
                1.0 / talD->count1[model][sorted_listD[j].key[0]]);
              for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
                if( fabs(finalprob->prob2[model][sorted_listD[j].key[0]]
                                                [base]) < MAXDIFF )
                  finalprob->prob2[model][sorted_listD[j].key[0]][base] =
                    pfix;
                else
                  finalprob->prob2[model][sorted_listD[j].key[0]][base] =
                    (finalprob->prob2[model][sorted_listD[j].key[0]][base] *
                    (1 - (4 * pfix))) + pfix;

                fprintf(stdout,"%u=%.4f,",base,
                  finalprob->prob2[model][sorted_listD[j].key[0]][base]);
                pmf_test += finalprob->prob2[model][sorted_listD[j].key[0]]
                                                   [base];
              } /* end for */
              fprintf(stdout,
                " Sum of prob's (revised) for history %u is %.2f ",
                sorted_listD[j].key[0],pmf_test);
              if( fabs(pmf_test - CERTPROB) < MAXDIFF )
                fprintf(stdout,"valid\n");
              else
                fprintf(stdout,"problematic!\n");
            } /* end if (nulls exist in current pmf) */
            break;
          case 2 :
            if ( (fabs(finalprob->prob3[model][sorted_listD[j].key[0]]
                                              [sorted_listD[j].key[1]][0])
                   < MAXDIFF) ||
                 (fabs(finalprob->prob3[model][sorted_listD[j].key[0]]
                                              [sorted_listD[j].key[1]][1])
                   < MAXDIFF) ||
                 (fabs(finalprob->prob3[model][sorted_listD[j].key[0]]
                                              [sorted_listD[j].key[1]][2])
                   < MAXDIFF) ||
                 (fabs(finalprob->prob3[model][sorted_listD[j].key[0]]
                                              [sorted_listD[j].key[1]][3])
                   < MAXDIFF)
               ) {
              fprintf(stdout,"Making null probability adjustments...\n");
              pfix=(double)MAX(PfixOMEGA,
                1.0 / talD->count2[model][sorted_listD[j].key[0]]
                                         [sorted_listD[j].key[1]]);
              for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
                if( fabs(finalprob->prob3[model][sorted_listD[j].key[0]]
                                                [sorted_listD[j].key[1]]
                                                [base]) < MAXDIFF )
                  finalprob->prob3[model][sorted_listD[j].key[0]]
                                         [sorted_listD[j].key[1]][base] =
                    pfix;
                else
                  finalprob->prob3[model][sorted_listD[j].key[0]]
                                         [sorted_listD[j].key[1]][base] =
                    (finalprob->prob3[model][sorted_listD[j].key[0]]
                                            [sorted_listD[j].key[1]][base] *
                    (1 - (4 * pfix))) + pfix;

                fprintf(stdout,"%u=%.4f,",base,
                  finalprob->prob3[model][sorted_listD[j].key[0]]
                                         [sorted_listD[j].key[1]][base]);
                pmf_test += finalprob->prob3[model][sorted_listD[j].key[0]]
                                                   [sorted_listD[j].key[1]]
                                                   [base];
              } /* end for */
              fprintf(stdout,
                " Sum of prob's (revised) for history %u%u is %.2f ",
                sorted_listD[j].key[0],
                sorted_listD[j].key[1],pmf_test);
              if( fabs(pmf_test - CERTPROB) < MAXDIFF )
                fprintf(stdout,"valid\n");
              else
                fprintf(stdout,"problematic!\n");
            } /* end if (nulls exist in current pmf) */
            break;
          case 3 :
            if ( (fabs(finalprob->prob4[model][sorted_listD[j].key[0]]
                                              [sorted_listD[j].key[1]]
                                              [sorted_listD[j].key[2]][0])
                   < MAXDIFF) ||
                 (fabs(finalprob->prob4[model][sorted_listD[j].key[0]]
                                              [sorted_listD[j].key[1]]
                                              [sorted_listD[j].key[2]][1])
                   < MAXDIFF) ||
                 (fabs(finalprob->prob4[model][sorted_listD[j].key[0]]
                                              [sorted_listD[j].key[1]]
                                              [sorted_listD[j].key[2]][2])
                   < MAXDIFF) ||
                 (fabs(finalprob->prob4[model][sorted_listD[j].key[0]]
                                              [sorted_listD[j].key[1]]
                                              [sorted_listD[j].key[2]][3])
                   < MAXDIFF)
               ) {
              fprintf(stdout,"Making null probability adjustments...\n");
              pfix=(double)MAX(PfixOMEGA,
                1.0 / talD->count3[model][sorted_listD[j].key[0]]
                                         [sorted_listD[j].key[1]]
                                         [sorted_listD[j].key[2]]);
              for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
                if( fabs(finalprob->prob4[model][sorted_listD[j].key[0]]
                                                [sorted_listD[j].key[1]]
                                                [sorted_listD[j].key[2]]
                                                [base]) < MAXDIFF )
                  finalprob->prob4[model][sorted_listD[j].key[0]]
                                         [sorted_listD[j].key[1]]
                                         [sorted_listD[j].key[2]][base] =
                    pfix;
                else
                  finalprob->prob4[model][sorted_listD[j].key[0]]
                                         [sorted_listD[j].key[1]]
                                         [sorted_listD[j].key[2]][base] =
                    (finalprob->prob4[model][sorted_listD[j].key[0]]
                                            [sorted_listD[j].key[1]]
                                            [sorted_listD[j].key[2]][base] *
                    (1 - (4 * pfix))) + pfix;

                fprintf(stdout,"%u=%.4f,",base,
                  finalprob->prob4[model][sorted_listD[j].key[0]]
                                         [sorted_listD[j].key[1]]
                                         [sorted_listD[j].key[2]][base]);
                pmf_test += finalprob->prob4[model][sorted_listD[j].key[0]]
                                                   [sorted_listD[j].key[1]]
                                                   [sorted_listD[j].key[2]]
                                                   [base];
              } /* end for */
              fprintf(stdout,
                " Sum of prob's (revised) for history %u%u%u is %.2f ",
                sorted_listD[j].key[0],
                sorted_listD[j].key[1],
                sorted_listD[j].key[2],pmf_test);
              if( fabs(pmf_test - CERTPROB) < MAXDIFF )
                fprintf(stdout,"valid\n");
              else
                fprintf(stdout,"problematic!\n");
            } /* end if (nulls exist in current pmf) */
            break;
          case 4 :
            if ( (fabs(finalprob->prob5[model][sorted_listD[j].key[0]]
                                              [sorted_listD[j].key[1]]
                                              [sorted_listD[j].key[2]]
                                              [sorted_listD[j].key[3]][0])
                   < MAXDIFF) ||
                 (fabs(finalprob->prob5[model][sorted_listD[j].key[0]]
                                              [sorted_listD[j].key[1]]
                                              [sorted_listD[j].key[2]]
                                              [sorted_listD[j].key[3]][1])
                   < MAXDIFF) ||
                 (fabs(finalprob->prob5[model][sorted_listD[j].key[0]]
                                              [sorted_listD[j].key[1]]
                                              [sorted_listD[j].key[2]]
                                              [sorted_listD[j].key[3]][2])
                   < MAXDIFF) ||
                 (fabs(finalprob->prob5[model][sorted_listD[j].key[0]]
                                              [sorted_listD[j].key[1]]
                                              [sorted_listD[j].key[2]]
                                              [sorted_listD[j].key[3]][3])
                   < MAXDIFF)
               ) {
              fprintf(stdout,"Making null probability adjustments...\n");
              pfix=(double)MAX(PfixOMEGA,
                1.0 / talD->count4[model][sorted_listD[j].key[0]]
                                         [sorted_listD[j].key[1]]
                                         [sorted_listD[j].key[2]]
                                         [sorted_listD[j].key[3]]);
              for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
                if( fabs(finalprob->prob5[model][sorted_listD[j].key[0]]
                                                [sorted_listD[j].key[1]]
                                                [sorted_listD[j].key[2]]
                                                [sorted_listD[j].key[3]]
                                                [base]) < MAXDIFF )
                  finalprob->prob5[model][sorted_listD[j].key[0]]
                                         [sorted_listD[j].key[1]]
                                         [sorted_listD[j].key[2]]
                                         [sorted_listD[j].key[3]][base] =
                    pfix;
                else
                  finalprob->prob5[model][sorted_listD[j].key[0]]
                                         [sorted_listD[j].key[1]]
                                         [sorted_listD[j].key[2]]
                                         [sorted_listD[j].key[3]][base] =
                    (finalprob->prob5[model][sorted_listD[j].key[0]]
                                            [sorted_listD[j].key[1]]
                                            [sorted_listD[j].key[2]]
                                            [sorted_listD[j].key[3]][base] *
                    (1 - (4 * pfix))) + pfix;

                fprintf(stdout,"%u=%.4f,",base,
                  finalprob->prob5[model][sorted_listD[j].key[0]]
                                         [sorted_listD[j].key[1]]
                                         [sorted_listD[j].key[2]]
                                         [sorted_listD[j].key[3]][base]);
                pmf_test += finalprob->prob5[model][sorted_listD[j].key[0]]
                                                   [sorted_listD[j].key[1]]
                                                   [sorted_listD[j].key[2]]
                                                   [sorted_listD[j].key[3]]
                                                   [base];
              } /* end for */
              fprintf(stdout,
                " Sum of prob's (revised) for history %u%u%u%u is %.2f ",
                sorted_listD[j].key[0],
                sorted_listD[j].key[1],
                sorted_listD[j].key[2],
                sorted_listD[j].key[3],pmf_test);
              if( fabs(pmf_test - CERTPROB) < MAXDIFF )
                fprintf(stdout,"valid\n");
              else
                fprintf(stdout,"problematic!\n");
            } /* end if (nulls exist in current pmf) */
            break;
          case 5 :
            if ( (fabs(finalprob->prob6[model][sorted_listD[j].key[0]]
                                              [sorted_listD[j].key[1]]
                                              [sorted_listD[j].key[2]]
                                              [sorted_listD[j].key[3]]
                                              [sorted_listD[j].key[4]][0])
                   < MAXDIFF) ||
                 (fabs(finalprob->prob6[model][sorted_listD[j].key[0]]
                                              [sorted_listD[j].key[1]]
                                              [sorted_listD[j].key[2]]
                                              [sorted_listD[j].key[3]]
                                              [sorted_listD[j].key[4]][1])
                   < MAXDIFF) ||
                 (fabs(finalprob->prob6[model][sorted_listD[j].key[0]]
                                              [sorted_listD[j].key[1]]
                                              [sorted_listD[j].key[2]]
                                              [sorted_listD[j].key[3]]
                                              [sorted_listD[j].key[4]][2])
                   < MAXDIFF) ||
                 (fabs(finalprob->prob6[model][sorted_listD[j].key[0]]
                                              [sorted_listD[j].key[1]]
                                              [sorted_listD[j].key[2]]
                                              [sorted_listD[j].key[3]]
                                              [sorted_listD[j].key[4]][3])
                   < MAXDIFF)
               ) {
              fprintf(stdout,"Making null probability adjustments...\n");
              pfix=(double)MAX(PfixOMEGA,
                1.0 / talD->count5[model][sorted_listD[j].key[0]]
                                         [sorted_listD[j].key[1]]
                                         [sorted_listD[j].key[2]]
                                         [sorted_listD[j].key[3]]
                                         [sorted_listD[j].key[4]]);
              for(base=0,pmf_test=0.0;base<ALFSIZE;++base) {
                if( fabs(finalprob->prob6[model][sorted_listD[j].key[0]]
                                                [sorted_listD[j].key[1]]
                                                [sorted_listD[j].key[2]]
                                                [sorted_listD[j].key[3]]
                                                [sorted_listD[j].key[4]]
                                                [base]) < MAXDIFF )
                  finalprob->prob6[model][sorted_listD[j].key[0]]
                                         [sorted_listD[j].key[1]]
                                         [sorted_listD[j].key[2]]
                                         [sorted_listD[j].key[3]]
                                         [sorted_listD[j].key[4]][base] =
                    pfix;
                else
                  finalprob->prob6[model][sorted_listD[j].key[0]]
                                         [sorted_listD[j].key[1]]
                                         [sorted_listD[j].key[2]]
                                         [sorted_listD[j].key[3]]
                                         [sorted_listD[j].key[4]][base] =
                    (finalprob->prob6[model][sorted_listD[j].key[0]]
                                            [sorted_listD[j].key[1]]
                                            [sorted_listD[j].key[2]]
                                            [sorted_listD[j].key[3]]
                                            [sorted_listD[j].key[4]][base] *
                    (1 - (4 * pfix))) + pfix;

                fprintf(stdout,"%u=%.4f,",base,
                  finalprob->prob6[model][sorted_listD[j].key[0]]
                                         [sorted_listD[j].key[1]]
                                         [sorted_listD[j].key[2]]
                                         [sorted_listD[j].key[3]]
                                         [sorted_listD[j].key[4]][base]);
                pmf_test += finalprob->prob6[model][sorted_listD[j].key[0]]
                                                   [sorted_listD[j].key[1]]
                                                   [sorted_listD[j].key[2]]
                                                   [sorted_listD[j].key[3]]
                                                   [sorted_listD[j].key[4]]
                                                   [base];
              } /* end for */
              fprintf(stdout,
                " Sum of prob's (revised) for history %u%u%u%u%u is %.2f ",
                sorted_listD[j].key[0],
                sorted_listD[j].key[1],
                sorted_listD[j].key[2],
                sorted_listD[j].key[3],
                sorted_listD[j].key[4],pmf_test);
              if( fabs(pmf_test - CERTPROB) < MAXDIFF )
                fprintf(stdout,"valid\n");
              else
                fprintf(stdout,"problematic!\n");
            } /* end if (nulls exist in current pmf) */
            break;
          default :
            fprintf(stderr,
              "Err (DIMM_final_probs): Invalid order encountered!\n");
            exit(EXIT_FAILURE);
        } /* end switch */
      } /* end bucket index for */
      #endif

      /* Update indices for the next bucket to be processed. */
      buckstart_IND=buckstop_IND + 1; /* = i */
      if(buckstart_IND < list_len)
        Buck_start_ct=sorted_listD[buckstart_IND].ct_pret;
      /* else, loop will not execute again */
    } /* end bucketing for */
  } /* end order for */

  return;
} /* end DIMM_final_probs */

/* This function handles the bisection method; it's an  *
 * abstracted wrapper for the fxns4(TDDI|BUDI) fxns.    *
 * Set topdown to TRUE for top-down training, FALSE for *
 * bottom-up.                                           */
static double calc_root(int *iteration,double HI,double LO,countarrayT *talC,
  relfreqsT *rfreqs,immprobT* fprob,PsupiT *incsmooth,int model,int maxord,
  int order,ind_listT *sorted_listD,int startIND,int stopIND,int topdown)
{
  double currhi=HI,
         currlo=LO,
         currmid=((currhi-currlo)/2.0)+currlo,
         resulthi,resultlo,resultmid;
  
  /* compute results for lo,hi,mid possible root solutions */
  if (topdown == TRUE) { /* top-down training algorithm */
    resultlo=fxns4TDDI(talC,rfreqs,fprob,model,maxord,sorted_listD,
      startIND,stopIND,currlo,TRUE);
    resulthi=fxns4TDDI(talC,rfreqs,fprob,model,maxord,sorted_listD,
      startIND,stopIND,currhi,TRUE);
    resultmid=fxns4TDDI(talC,rfreqs,fprob,model,maxord,sorted_listD,
      startIND,stopIND,currmid,TRUE);
  }
  else if (topdown == FALSE) { /* bottom-up training algorithm */
    resultlo=fxns4BUDI(talC,rfreqs,incsmooth,model,maxord,order,sorted_listD,
      startIND,stopIND,currlo,TRUE);
    resulthi=fxns4BUDI(talC,rfreqs,incsmooth,model,maxord,order,sorted_listD,
      startIND,stopIND,currhi,TRUE);
    resultmid=fxns4BUDI(talC,rfreqs,incsmooth,model,maxord,order,sorted_listD,
      startIND,stopIND,currmid,TRUE);
  }
  else {
    fprintf(stderr,"Err (calc_root): Improper algorithm requested!\n");
    exit(EXIT_FAILURE);
  }

  /* If the root doesn't exist in the (0,1) interval of f'x, we *
   * will set lambda to argmax(fx) chosen from either 0 or 1.   */
  if( (fabs(resulthi) > DBL_EPSILON && fabs(resultlo) > DBL_EPSILON) &&
      ((resulthi < NULLPROB && resultlo < NULLPROB) ||
       (resulthi > NULLPROB && resultlo > NULLPROB))
    ) {
    #ifdef ROOTREPORT
    fprintf(stdout,"  Root DNE in (0<x<1) for f'x,\
 model %u (resulthi %.4f,resultlo %.4f,",model,resulthi,resultlo);
    if (topdown == TRUE) /* only maxord is relevant */
      fprintf(stdout,"maxorder %u) ",maxord);
    else /* both maxord and order are relevant for BUDI */
      fprintf(stdout,"maxorder %u,innerorder %u) ",maxord,order);
    #endif
    currhi=CERTPROB;
    currlo=NULLPROB;

    if (topdown == TRUE) { /* top-down training algorithm */
      resultlo=fxns4TDDI(talC,rfreqs,fprob,model,maxord,sorted_listD,
        startIND,stopIND,currlo,FALSE);
      resulthi=fxns4TDDI(talC,rfreqs,fprob,model,maxord,sorted_listD,
        startIND,stopIND,currhi,FALSE);
    }
    else if (topdown == FALSE) { /* bottom-up training algorithm */
      resultlo=fxns4BUDI(talC,rfreqs,incsmooth,model,maxord,order,sorted_listD,
        startIND,stopIND,currlo,FALSE);
      resulthi=fxns4BUDI(talC,rfreqs,incsmooth,model,maxord,order,sorted_listD,
        startIND,stopIND,currhi,FALSE);
    }
    else {
      fprintf(stderr,"Err (calc_root): Improper algorithm requested!\n");
      exit(EXIT_FAILURE);
    }

    /* Now, we can return a result to the calling   *
     * function--no need to proceed with bisection. */
    if(fabs(resulthi-resultlo) > DBL_EPSILON && resulthi > resultlo) {
      #ifdef ROOTREPORT
      fprintf(stdout,"(argmax of fx == 1.0)\n");
      #endif
      return(currhi);
    }
    else {
      #ifdef ROOTREPORT
      fprintf(stdout,"(argmax of fx == 0.0)\n");
      #endif
      return(currlo);
    }
  }

  /* So, I do these same conditional checks in the if-else if     *
   * clause below, but I'll do so here, too, so I don't multiply  *
   * bury this (preprocessor dependent) print stuff in the former */
  #ifdef ROOTREPORT
  if( fabs(resulthi) < MAXDIFF ||
      fabs(resultlo) < MAXDIFF ||
      fabs(resultmid) < MAXDIFF
    ) {
    fprintf(stdout,"  Root exists in (0<x<1) for f'x, model %u,",model);
    if (topdown == TRUE) /* only maxord is relevant */
      fprintf(stdout,"maxorder %u\n",maxord);
    else /* both maxord and order are relevant for BUDI */
      fprintf(stdout,"maxorder %u,innerorder %u\n",maxord,order);
  }
  #endif

  /* Are any of these reasonable roots to return to calling function? */
  if(fabs(resulthi) < MAXDIFF)
    return(currhi);
  else if(fabs(resultlo) < MAXDIFF)
    return(currlo);
  else if(fabs(resultmid) < MAXDIFF)
    return(currmid);
  else
    ;

  /* As a precaution, assert that none of the results might really be zero */
  if(fabs(resulthi) < DBL_EPSILON ||
     fabs(resultlo) < DBL_EPSILON ||
     fabs(resultmid) < DBL_EPSILON
    ) {
    fprintf(stderr,"Err (calc_root): \
Ignored valid root somehow! You set MAXDIFF < DBL_EPSILON perhaps?\n");
    exit(EXIT_FAILURE);
  }

  /* If we couldn't find reasonable root in MAXDEPTH recursions, *
   * we'll just say something's amiss and bail out.              */
  if(*iteration >= MAXITER) { /* .GT. for if the impossible happens! */
    fprintf(stderr,"Err (calc_root): Couldn't find root at \
reasonable recursion depth!\n");
    exit(EXIT_FAILURE);
  }
  else
    *iteration += 1;

  /* Determine whether currmid will become currhi *
   * or currlo when closing in on the root.       */
  if(resulthi < NULLPROB) {
    if(resultmid < NULLPROB)
      currhi=currmid;
    else
      currlo=currmid;
  }
  else if(resulthi > NULLPROB) {
    if(resultmid > NULLPROB)
      currhi=currmid;
    else
      currlo=currmid;
  }
  else {
    fprintf(stderr,"Err (calc_root): An unexpected error occurred!\n");
    exit(EXIT_FAILURE);
  }

  /* recurse */
  return(
    calc_root(iteration,currhi,currlo,talC,rfreqs,fprob,incsmooth,
      model,maxord,order,sorted_listD,startIND,stopIND,topdown)
  );
} /* end calc_root */

/* Function to compute result of using candidate root.        *
 * Note: call with derivative == TRUE to compute second-order *
 *       fxn, derivative == FALSE for first-order fxn         */
static double fxns4TDDI(countarrayT *talC,relfreqsT* rfreqs,
  immprobT* fprob,int model,int order,ind_listT *sorted_listD,
  int startIND,int stopIND,double rootcandidate,int derivative)
{
  int Ccount,         /* per history in bucket, its count in C       *
                       * dataset (Cc(v,hi))                          */
      hist[MAXORDER], /* used for keying into data structures by     *
                       * history                                     */
      base,           /* iterate through nt alphabet                 */
      i,j;            /* iterator variables                          */
  double MLcurr,      /* ML estimate of current history (PrML[v|hi]) */
         DIprev,      /* TDDI likelihood for next shorter history    *
                       * (PrDI[v|h_{i-1})                            */
         result=0.0;  /* return to calling function                  */

  /* Order 0 is a special case. Note in this fxn, *
   * order ranges from 0..MAXORDER (k)            */
  if (order == 0) {
    if(derivative == TRUE) { /* compute using the derivative */
      for(base=0;base<ALFSIZE;++base) {
        result += talC->count1[model][base] *
          ((EQUIPROB - rfreqs->MLhis0[base]) /
           ((rootcandidate * (EQUIPROB - rfreqs->MLhis0[base])) +
             rfreqs->MLhis0[base]));
      }
    }
    else if(derivative == FALSE) { /* compute first-order fxn */
      for(base=0;base<ALFSIZE;++base) {
        result += talC->count1[model][base] *
          log((rootcandidate * EQUIPROB) + ((CERTPROB - rootcandidate) *
            rfreqs->MLhis0[base]));
      }
    }
    else {
      fprintf(stderr,"Err (fxns4TDDI): Which order fxn do you want?\n");
      exit(EXIT_FAILURE);
    }
  } /* end order 0 if */
  else { /* order .GE. 1 */
    for(i=startIND;i<=stopIND;++i) {
      for(base=0;base<ALFSIZE;++base) {
        /* Populate variables used in formula */
        for(j=0;j<order;++j)
          hist[j]=sorted_listD[i].key[j];

        switch(order) {
          case 1 :
            Ccount = talC->count2[model][hist[0]][base];
            MLcurr = rfreqs->MLhis1[hist[0]][base];
            DIprev = fprob->prob1[model][base];
            break;
          case 2 :
            Ccount = talC->count3[model][hist[0]]
                                        [hist[1]][base];
            MLcurr = rfreqs->MLhis2[hist[0]]
                                   [hist[1]][base];
            DIprev = fprob->prob2[model][hist[1]][base];
            break;
          case 3 :
            Ccount = talC->count4[model][hist[0]]
                                        [hist[1]]
                                        [hist[2]][base];
            MLcurr = rfreqs->MLhis3[hist[0]]
                                   [hist[1]]
                                   [hist[2]][base];
            DIprev = fprob->prob3[model][hist[1]]
                                        [hist[2]][base];
            break;
          case 4 :
            Ccount = talC->count5[model][hist[0]]
                                        [hist[1]]
                                        [hist[2]]
                                        [hist[3]][base];
            MLcurr = rfreqs->MLhis4[hist[0]]
                                   [hist[1]]
                                   [hist[2]]
                                   [hist[3]][base];
            DIprev = fprob->prob4[model][hist[1]]
                                        [hist[2]]
                                        [hist[3]][base];
            break;
          case 5 :
            Ccount = talC->count6[model][hist[0]]
                                        [hist[1]]
                                        [hist[2]]
                                        [hist[3]]
                                        [hist[4]][base];
            MLcurr = rfreqs->MLhis5[hist[0]]
                                   [hist[1]]
                                   [hist[2]]
                                   [hist[3]]
                                   [hist[4]][base];
            DIprev = fprob->prob5[model][hist[1]]
                                        [hist[2]]
                                        [hist[3]]
                                        [hist[4]][base];
            break;
          default :
            fprintf(stderr,"Err (fxns4TDDI): Invalid case encountered!\n");
            exit(EXIT_FAILURE);
        }

        if(derivative == TRUE) /* compute using the derivative */
          result += Ccount * ((DIprev - MLcurr) /
            ((rootcandidate * (DIprev - MLcurr)) + MLcurr));
        else if(derivative == FALSE) /* compute first-order fxn */
          result += Ccount * log((rootcandidate * DIprev) +
              ((CERTPROB - rootcandidate) * MLcurr));
        else {
          fprintf(stderr,"Err (fxns4TDDI): Which order fxn do you want?\n");
          exit(EXIT_FAILURE);
        }
      } /* end base for */
    } /* end history for */
  } /* end else for order > 0 */

  return(result);
} /* end fxns4TDDI */

/* Function to compute result of using candidate root given   *
 * the bottom-up training algorithm.                          *
 * Note: call with derivative == TRUE to compute second-order *
 *       fxn, derivative == FALSE for first-order fxn         */
static double fxns4BUDI(countarrayT *talC,relfreqsT *rfreqs,
  PsupiT *incsmooth,int model,int currmaxord,int currorder,
  ind_listT *sorted_listD,int startIND,int stopIND,
  double rootcandidate,int derivative)
{
  int Ccount,         /* per history in bucket, its count in  *
                       * C dataset (Cc(v,hi))                 */
      hist[MAXORDER], /* used for keying into data structures *
                       * by history                           */
      base,           /* iterate through nt alphabet          */
      i,j;            /* general-purpose iterator variables   */
  double MLstor,      /* store relfreq for oligo with next-   *
                       * shorter history/pretext length.      */
         Psupistor,   /* Store current Psupi value based on   *
                       * current order,innerorder combination */
         result=0.0;  /* return to calling function           */

  /* Nest history loop within base one, b/c if currmaxord==0,     *
   * there are no bucket elements to iterate over, only 3' bases. */
  for(base=0;base<ALFSIZE;++base) {
    if(currmaxord==0) {
      /* Life's a little simpler with this condition...remember that  *
       * zeroth-order histories comprise a singleton bucket -> no     *
       * nested loop over histories in the bucket to worry about.     *
       * Note also that currorder must be 0 when this branch executes */
      if(currorder!=0) { /* just a sanity check */
        fprintf(stderr,
          "Err (fxns4BUDI): currmaxord==0, but currorder!=0 !??\n");
        exit(EXIT_FAILURE);
      }

      /* Compute result.  Note that incsmooth->Psuphis0[0][base]    *
       * was set in the final_prob function which called calc_root, *
       * which called this function.                                */
      if(derivative == TRUE) /* compute using the derivative (12b)  */
        result += talC->count1[model][base] *
          ((incsmooth->Psuphis0[0][base] - EQUIPROB) /
          ((rootcandidate * (incsmooth->Psuphis0[0][base] -
          EQUIPROB)) + EQUIPROB));
      else if(derivative == FALSE) /* compute first-order fxn (12b) */
        result += talC->count1[model][base] *
          log((rootcandidate * incsmooth->Psuphis0[0][base]) +
          ((CERTPROB - rootcandidate) * EQUIPROB));
      else {
        fprintf(stderr,"Err (fxns4BUDI): Which order fxn do you want?\n");
        exit(EXIT_FAILURE);
      }
    }
    else { /* more general case; currmaxord .GE. 1 */
      for(i=startIND;i<=stopIND;++i) {
        /* Go ahead and copy currmaxord history into hist array *
         * This will be useful for keying into various data     *
         * structures based on histories in a general manner.   */
        for(j=0;j<currmaxord;++j)
          hist[j]=sorted_listD[i].key[j];

        /* The take-home point of these next two switch statements *
         * is that they set variables that will allow us to deal   *
         * with certain data in a generic way relative to the      *
         * currmaxord variable.                                    */

        /* Generically assign Pr_{ML}[v|h_{i-1}] for eqn 12b */
        switch(currorder-1) {
          case -1 :
            MLstor = EQUIPROB;
            break;
          case 0 :
            MLstor = rfreqs->MLhis0[base];
            break;
          case 1 :
            MLstor = rfreqs->MLhis1[hist[currmaxord-1]][base];
            break;
          case 2 :
            MLstor = rfreqs->MLhis2[hist[currmaxord-2]]
                                   [hist[currmaxord-1]][base];
            break;
          case 3 :
            MLstor = rfreqs->MLhis3[hist[currmaxord-3]]
                                   [hist[currmaxord-2]]
                                   [hist[currmaxord-1]][base];
            break;
          case 4 :
            MLstor = rfreqs->MLhis4[hist[currmaxord-4]]
                                   [hist[currmaxord-3]]
                                   [hist[currmaxord-2]]
                                   [hist[currmaxord-1]][base];
            break;
          case 5 :
            /* Fooled you! Provided MAXORDER .LE. 5 .AND. MAXORDER .GE. 0,   *
             * (see immpractical.h), currorder-1 would only range from -1..4 */
          default :
            fprintf(stderr,"Err (fxns4BUDI): Invalid case encountered!\n");
            exit(EXIT_FAILURE);
        } /* end currorder switch */
        /* Generically assign Cc(v,h_k) and P^i[v|h_k] for eqn 12b */
        switch(currmaxord) {
          case 1 :
            Ccount = talC->count2[model][hist[0]][base];
            Psupistor = incsmooth->Psuphis1[currorder][hist[0]][base];
            break;
          case 2 :
            Ccount = talC->count3[model][hist[0]]
                                        [hist[1]][base];
            Psupistor = incsmooth->Psuphis2[currorder][hist[0]]
                                                      [hist[1]][base];
            break;
          case 3 :
            Ccount = talC->count4[model][hist[0]]
                                        [hist[1]]
                                        [hist[2]][base];
            Psupistor = incsmooth->Psuphis3[currorder][hist[0]]
                                                      [hist[1]]
                                                      [hist[2]][base];
            break;
          case 4 :
            Ccount = talC->count5[model][hist[0]]
                                        [hist[1]]
                                        [hist[2]]
                                        [hist[3]][base];
            Psupistor = incsmooth->Psuphis4[currorder][hist[0]]
                                                      [hist[1]]
                                                      [hist[2]]
                                                      [hist[3]][base];
            break;
          case 5 :
            Ccount = talC->count6[model][hist[0]]
                                        [hist[1]]
                                        [hist[2]]
                                        [hist[3]]
                                        [hist[4]][base];
            Psupistor = incsmooth->Psuphis5[currorder][hist[0]]
                                                      [hist[1]]
                                                      [hist[2]]
                                                      [hist[3]]
                                                      [hist[4]][base];
            break;
          default :
            fprintf(stderr,"Err (fxns4BUDI): Invalid case encountered!\n");
            exit(EXIT_FAILURE);
        } /* end currmaxord switch */

        /* Now, we can update our result */
        if(derivative == TRUE) /* compute using the derivative */
          result += Ccount * ((Psupistor - MLstor) /
            ((rootcandidate * (Psupistor - MLstor)) + MLstor));
        else if(derivative == FALSE) /* compute first-order fxn */
          result += Ccount * log((rootcandidate * Psupistor) +
            ((CERTPROB - rootcandidate) * MLstor));
        else {
          fprintf(stderr,"Err (fxns4BUDI): Which order fxn do you want?\n");
          exit(EXIT_FAILURE);
        }
      } /* end history index for */
    } /* end (currmaxord != 0) else */
  } /* end base for */

  return(result);
} /* end fxns4BUDI */
