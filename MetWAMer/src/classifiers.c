/* classifiers.c
 * Michael Sparks (mespar1@gmail.com)
 * Last modified : 20 July 2013
 *
 * Contains classifier-related functions for use in MetWAMer.
 *
 * Copyright (c) 2007,2013 Michael E Sparks
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "classifiers.h"
#include "errors.h"

/* Neuron-related classifier methods *****************************************/

/* Wraps the process of adjusting an LTU's weight array */
void update_weight_vector(
  featvecT features,  /* Feature vector of instance     */
  LTUparmT *weights,  /* Weight vector to modify        */
  double target,      /* Records true class of instance */
  double output,      /* Output of activation unit      */
  int classifiercode, /* see classifiers.h              */
  double learnrate    /* Learning rate                  */
) {
  double cachedfactor=0; /* Cache factor common to all dimensions */

  switch(classifiercode) {
    case LINEARUNITX : /* T4.1 of Mitchell 1997 */
      cachedfactor=learnrate*(target-output);
      break;
    case SIGMOIDUNITX : /* See Fig. 20.21 of Russell/Norvig AIMA 2ed */
      cachedfactor=learnrate*(target-output)*(output*(1.0-output));
      break;
    default :
      FATALERROR("\nErr (update_for_single_weight): \
Improper classifier method requested!\n\n")
  } /* end switch */

  weights->wzero+=cachedfactor; /* Imagine a feature of 1 */
  weights->wmetwam+=cachedfactor*features.metwam;
  weights->wmcratio+=cachedfactor*features.mcratio;

  return;
} /* end update_weight_vector */

/* Add components of vector arguments */
LTUparmT add_weight_vectors(
  LTUparmT ltu1,
  LTUparmT ltu2
) {
  LTUparmT ltuSUM;

  ltuSUM.wzero=ltu1.wzero+ltu2.wzero;
  ltuSUM.wmetwam=ltu1.wmetwam+ltu2.wmetwam;
  ltuSUM.wmcratio=ltu1.wmcratio+ltu2.wmcratio;

  return(ltuSUM);
} /* end add_weight_vectors */

/* Assay an input instance using a trained neuron */
double eval_activation_unit(
  featvecT features,     /* Feature vector of instance   */
  LTUparmT trainedmodel, /* Learned weights for the LTU  */
  int *classification,   /* Records (Boolean) result of  *
                          * classification               */
  double *fctofdotprod,  /* Stores result of activation  *
                          * fct on dot product of weight *
                          * and feature vectors          */
  int classifiercode     /* see classifiers.h            */
) {
  /* Classification boundaries */
  const double linear_threshold=0.0,  /* (+) = +1, (-) = -1 */
               sigmoid_threshold=0.5; /* (+) = +1, (-) =  0 */

  /* compute dot product of weight and feature arrays */
  *fctofdotprod=trainedmodel.wzero+
    trainedmodel.wmetwam*features.metwam+
    trainedmodel.wmcratio*features.mcratio;

  *classification=FALSE;
  switch(classifiercode) {
    case LINEARUNITX : /* Unthresholded neuron */
      if(*fctofdotprod>linear_threshold)
        *classification=TRUE;
      break;
    case SIGMOIDUNITX : /* Soft-thresholded neuron */
      /* Transform dot product using logistic fct */
      *fctofdotprod=1.0/(1.0+exp(-1.0**fctofdotprod));
      if(*fctofdotprod>sigmoid_threshold)
        *classification=TRUE;
      break;
    default :
      FATALERROR("\nErr (eval_activation_unit): \
Improper classifier method requested!\n\n")
  }

  return(*fctofdotprod);
} /* end eval_activation_unit */

/* Function to import a trained LTUparmT object */
void import_trained_perceptron(char *filename,LTUparmT *parms) {
  FILE *fptr=NULL; /* For connecting to binary input file */

  /* Open parameter file */
  if((fptr=fopen(filename,"rb"))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (import_trained_perceptron): \
Can't open %s for binary reading!\n",filename);
    FATALERROR(errbuff)
  }
  /* Read in the LTUparmT object */
  else if(fread(parms,sizeof(LTUparmT),1,fptr)!=1) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (import_trained_perceptron): \
Error reading from %s!\n",filename);
    FATALERROR(errbuff)
  }
  else
    (void)fclose(fptr);

  return;
} /* end import_trained_perceptron */
