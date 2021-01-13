/* MetWAM_utils.h
 * Michael Sparks (mespar1@gmail.com)
 * Last modified : 24 February 2007
 *
 * Functions germane to calculations derived from a trained MetWAM
 *
 * Copyright (c) 2007,2008 Michael E Sparks
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

#ifndef METWAM_UTILS_H
#define METWAM_UTILS_H

#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <limits.h>
#include "classifiers.h"
#include "sequence_parse.h"

#define UPDATEDXMLOUTSTREAM stdout /* File stream to direct MetWAMer- *
                                    * updated gthXML content to.      */

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

#define UPSTREXTENT    5 /* Region to consider upstream of ATG   */
#define DOWNSTREXTENT  3 /* Region to consider downstream of ATG */
#define SITELEN        3 /* length(ATG)                          */
/* Length of complete fragment to apply start-Met WAM over */
#define STRINGSIZE (UPSTREXTENT+SITELEN+DOWNSTREXTENT)

#define CONTENTSWATHLEN 96 /* Length of fragment to apply a content *
                            * sensor, i.e., a Markov chain, over.   */

#define MetTFGAPnt   105 /* Gap length (in nts) between the first *
                          * (true) Met and downstream ones to use *
                          * for modeling false sites.             */

#define HYPOTHESISCT   2 /* True or False initiation sites        */
#define TRUEMODIND     0 /* Index under the true start Met model  */
#define FALSEMODIND    1 /* Index under the false start Met model */

#define NOTISPRED     -1 /* Code returned by MetWAMer.CDS when *
                          * no TIS was predicted for a CDS.    */
#define NOTMULT3     -20 /* Code returned by searchforbestinitincds *
                          * when CDS length is not a multiple of 3  */

#define CODINGPRIOR  (1/6.0) /* Coding prior probability    */
#define NONCODPRIOR  (1/2.0) /* Noncoding prior probability */

#define MFCTHRESH    0.0  /* Minimal threshold for labeling an ATG *
                           * as a TIS under the MFCWLLKR method    */

#define METHODCOUNT 5 /* Number of methods MetWAMer uses *
                       * to identify trans init codons   */
enum { /* Allow no value to be .LT. 0 here!! */
  LLKRX=0,
  WLLKRX,
  BAYESX,
  MFCWLLKRX,
  NFCWLLKRX,
};
static const char methoddesc[METHODCOUNT][128]={
  "Log-likelihood ratios",
  "Weighted log-likelihood ratios",
  "Bayesian method",
  "WLLKR with flank contrasting, multiplicative-based",
  "WLLKR with flank contrasting, perceptron-based"
};

/* Data structures for determining which cluster-specific *
 * parameter deployment strategy the user wises to use.   */
#define DEPLOYCOUNT 3
static const char deploydesc[DEPLOYCOUNT][128]={
  "Hamming-distance based",
  "Position weight matrix-based",
  "Weight array matrix-based"
};
enum { /* Allow no value to be .LT. 0 here!! */
  HAMMODX=0,
  PWMMODX,
  WAMMODX,
  HAMSTATX,
  PWMSTATX,
  WAMSTATX,
};
int cluster_parm_selection; /* global variable */

/* Minimum length of an open reading frame *
 * for which we predict a start-Methionine */
#define MIN_PEPTIDE_LENGTH       50
/* Maximum length of an open reading frame *
 * for which we predict a start-Methionine */
#define MAX_PEPTIDE_LENGTH  INT_MAX
/* Macro to confirm that the length of an open reading frame (peplen), *
 * for which MetWAMer has predicted a start-Met, satisfies basic       *
 * constraints on peptide sequence length.                             */
#define PEPTIDE_LENGTH_CONSTRAINTS_MET_P(cdslen) \
  ((cdslen)/3>=(MIN_PEPTIDE_LENGTH)&& \
   (cdslen)/3<=(MAX_PEPTIDE_LENGTH)?TRUE:FALSE)

/* Structure to store start-Met PWM */
typedef struct {
  float PWMtable[STRINGSIZE][NTALFSIZE];
} ATGparmPWM;
    
/* Structure to store start-Met WAM */
typedef struct {
  float WAMtable[STRINGSIZE][NTALFSIZE][NTALFSIZE];
} ATGparm;
    
/* Macro that returns TRUE if pos, the first base of a (start) *
 * codon, has adequate flanking sequence to allow for signal-  *
 * based scoring, and FALSE otherwise. Note that up and down   *
 * denote upstream and downstream, respectively; typically,    *
 * uplim will be .LT. downlim. :)                              */
#define MET_VALID_TO_SCORE_P(pos,uplim,downlim) \
  ((pos)-(UPSTREXTENT)>=(uplim)&& \
   (pos)+2+(DOWNSTREXTENT)<=(downlim)?TRUE:FALSE)

/* Macro that returns TRUE if pos, the first base of a (start) *
 * codon, has adequate flanking sequence to allow for content- *
 * based scoring, and FALSE otherwise. Note that up and down   *
 * denote upstream and downstream, respectively; typically,    *
 * uplim will be .LT. downlim. :)                              */
#define MET_FLANKS_VALID_TO_SCORE_P(pos,uplim,downlim) \
  ((pos)-(UPSTREXTENT)-(CONTENTSWATHLEN)>=(uplim)&& \
   (pos)+2+(DOWNSTREXTENT)+(CONTENTSWATHLEN)<=(downlim)?TRUE:FALSE)

/* Macro to make the expression of translation     *
 * initiation prediction method recording clearer. */
#define BUFFSIZE  2048
#define METHOD_PRINTER(node,k,meth) { \
  char locmethodbuff[BUFFSIZE<<1]; \
  if(k==1) \
    SET_STRING_AS_ATTR_VAL(node,"method",meth) \
  else { \
    (void)snprintf(locmethodbuff,BUFFSIZE, \
      "Clustered (using %i groups)",k); \
    (void)strncat(locmethodbuff,meth,BUFFSIZE); \
    SET_STRING_AS_ATTR_VAL(node,"method",locmethodbuff) \
  } \
}

/* Structure for representing individual exons in a genomic sequences */
typedef struct exoncoors {
  int start, /* exon start position (1st exon base) */
      stop;  /* exon stop position (final exon base) */
  struct exoncoors *nextelt; /* Pointer to next element in list */
} exoncoorsT;

/* Structure to encapsulate data for use with the *
 * flanking content log-likelihood ration method. */
typedef struct fcwllkraux {
  int markovchainmeth,    /* markov chain to deploy */
      dmcloudindex,       /* cloud type for DMMM's  */
      neurontype;         /* perceptron type        */
  LTUparmT *localneuron;  /* learned perceptron     */
  void *markovchainparms; /* trained markov chain   */
} fcwllkrauxT;

/* Calculate log-likelihood ratio of a putative translation *
 * initiation codon, whose first base (A) corresponds to    *
 * *(seq+(UPSTREXTENT)).                                    */
long double calc_llkr(
  ATGparm parms_true,
  ATGparm parms_false,
  int *seq
);

/* Function to coordinate prediction of translation *
 * initiation sites based on <orf_entry> data in    *
 * gthXML documents.                                */
void predicttransinitsites(
  int transinitmeth,
  int strata_ct,
  char **medoids,
  ATGparm *ATGparms,
  ATGparmPWM *ATGparmsPWM,
  fcwllkrauxT *fcwllkrparms,
  char *template,
  const char *xmldoc,
  const char *relaxngdoc
);

/* Function to coordinate prediction of translation   *
 * initiation sites based on <isoform> data, produced *
 * specifically by the PASIF system.                  */
void predicttransinitsites4PASIF(
  int transinitmeth,
  int strata_ct,
  char **medoids,
  ATGparm *ATGparms,
  ATGparmPWM *ATGparmsPWM,
  fcwllkrauxT *fcwllkrparms,
  FILE *loci,
  const char *xmldoc,
  const char *relaxngdoc
);

/* Finds the most likely translation initiation site in a CDS */
int searchforbestinitincds(
  int transinitmeth,
  int strata_ct,
  char **medoids,
  ATGparm *ATGparms,
  ATGparmPWM *ATGparmsPWM,
  fcwllkrauxT *fcwllkrparms,
  int *sequence,
  int seqlen
);

/* Function to modify xpcontext such that registerednamespaces *
 * will be registered, and thus available to XPath queries.    *
 * registerednamespaces format is                              *
 * abbrev1=verbatim1 abbrev2=verbatim2 ... abbrev3=verbatim3   */
int reg_ns(
  xmlXPathContextPtr xpcontext,
  const xmlChar *registerednamespaces
);

/* Returns TRUE if gthXML doc validates against RELAX NG schema, FALSE o.w. */
int validateagainstRNG(
  const char *relaxngdoc,
  xmlDocPtr gthxmldoc
);

/* Free memory in the linked list representing a structure */
void freestructure(exoncoorsT *structure);

/* Computes the distance between two data items, *
 * namely the number of nucleotides differing    *
 * between the pair, i.e., edit distance.        */
double eval_dist_kernel(
  char *dat1,
  char *dat2
);

#endif
