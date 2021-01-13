/* classifiers.h
 * Michael Sparks (mespar1@gmail.com)
 * Last modified : 14 August 2007
 *
 * Contains classifier-related functions for use in MetWAMer.
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

#ifndef CLASSIFIERS_H
#define CLASSIFIERS_H

#include <immpractical.h>

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/* Data structures for handling machine learning classifiers. ****************/

#define MAXERROFWEIGHTS    10E-4 /* Target sum-of-squares error for  *
                                  * learned neuron weights.          */
#define AVAILCLASSIFIERCT      2 /* Number of classification methods *
                                  * available to PASIF.              */
enum { /* Which dimension corresponds to which classifier? *
        * It is critical that ORGOBJFCTX be the last elt   */
  LINEARUNITX=0,
  SIGMOIDUNITX
};
/* English descriptions of the available classifiers */
static const char classdesc[AVAILCLASSIFIERCT][MAXLOCSTRING] = {
  "Unthresholded Perceptron (Linear unit)",
  "Sigmoid Unit"
};
static const char classifiertags[AVAILCLASSIFIERCT][MAXLOCSTRING] = {
  "linunit",
  "sigunit"
};

/* Structure representing descriptive feature vector of a spliced isoform */
typedef struct {
  double
    metwam,  /* WAM-based, signal-sensing score */
    mcratio; /* MC-based, content-sensing score */
} featvecT;

/* Data structure for a trained perceptron classifier */
typedef struct {
  /* Weights for descriptive attributes */
  double
    wzero, /* Offset weight */
    wmetwam,
    wmcratio;
} LTUparmT;

/* Function prototypes *******************************************************/

/* Wraps the process of adjusting an LTU's weight array */
void update_weight_vector(
  featvecT features,  /* Feature vector of instance     */
  LTUparmT *weights,  /* Weight vector to modify        */
  double target,      /* Records true class of instance */
  double output,      /* Output of activation unit      */
  int classifiercode, /* see classifiers.h              */
  double learnrate    /* Learning rate                  */
);

/* Add components of vector arguments */
LTUparmT add_weight_vectors(
  LTUparmT ltu1,
  LTUparmT ltu2
);

/* Assay an input instance using a trained neuron */
double eval_activation_unit(
  featvecT features,     /* Feature vector of instance         */
  LTUparmT trainedmodel, /* Learned weights for the LTU        */
  int *classification,   /* Records (Boolean) result of        *
                          * classification                     */
  double *fctofdotprod,  /* Stores result of activation fct    *
                          * on dot product of weight/feat vecs */
  int classifiercode     /* see classifiers.h                  */
);

/* Function to import a trained LTUparmT object */
void import_trained_perceptron(char *filename,LTUparmT *parms);

#endif
