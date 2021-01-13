/* immpractical.h
 * Michael E Sparks (mespar1@gmail.com)
 *
 * Header file contains various definitions acting more or less
 * globally on code for the immpractical library.
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

#ifndef IMMPRACTICAL_H
#define IMMPRACTICAL_H

/* Auxilliary information ****************************************************/

#define EMAILADDRESS "mespar1@gmail.com"
#define WEBSITE      "http://brendelgroup.org/mespar1/"
#define DATACLEANPL  "clean_data.pl"

/* Parametric String Literals ************************************************/

#define ALFSIZE                       4 /* A,C,G,[T|U] ... no ambiguity     *
                                         * permitted since we're building a *
                                         * model, afterall!                 */
#define BUCKRATIO                   1.2 /* Ratio of adjacent buckets (b)    */
#define CERTPROB                    1.0 /* Pr(event) == 1.0                 */
#define CONFCUTOFF                 0.50 /* chi**2 confidence cutoff value   */
#define DEFRFREQ                    0.0 /* Default relative frequency (ML)  */
#define EMPTYCOUNT                    0 /* Count of a particular oligo that *
                                         * would trigger pseudocounting     */
#define EQUIPROB   (CERTPROB / ALFSIZE) /* ML of null prefix (depth -1)     */  
#ifndef FALSE
#define FALSE                         0 /* False value                      */
#endif
#define FLANK                        47 /* Hexamer flank length for quant-  *
                                         * ile specific functionality       */
#define FOPSEUDOCT                    5 /* Pseudocounts for fixed-order     */
#define HBOUND    (CERTPROB - POSCONST) /* The function to maximize has a   *
                                         * unique root in [0,1]             */
#define HYPONUM                       2 /* Rows in chi**2 contingency table */
#define INITFREQ                      0 /* Initial values for count arrays  */
#define LBOUND                 POSCONST /* The function to maximize has a   *
                                         * unique root in [0,1]             */
#define MAXDIFF                   10E-6 /* Max deviation for float compare  */
#define MAXFILENAME                1024 /* Max length for output filename   */
#define MAXLOCSTRING               1024 /* Max length for any local strings */
#define MAXITER                     100 /* Maximum number of bisections to  *
                                         * conduct. (recursion depth limit) */
#define MAXORDER                      5 /* [0..5] (pretext/history lengths) */
#define MINCOUNTUNIT                  1 /* Fundamental unit used for        *
                                         * pseudocount propagation.         */
#define MINRELIABLECHISQCT          400 /* Minimum number of occurrences of *
                                         * a history for the chi**2 model   *
                                         * to assign full probability to    */
#define NONINTERESTWEIGHT           0.0 /* default lambdahat for histories  *
                                         * not in the devel data set.       */
#define NULLPROB                    0.0 /* null probability of event        */
#ifdef THREEPERIODIC
#ifdef SHADOWSTRAND
#define NUMMODELS                     7 /* Model 0 (C g t),  (forward)      *
                                         *       1 (c G t),  (forward)      *
                                         *       2 (c g T),  (forward)      *
                                         *       3 (A c g),  (reverse)      *
                                         *       4 (a C g),  (reverse)      *
                                         *       5 (a c G),  (reverse)      *
                                         *       6 noncoding (intron)       */
#define NONCODING                     6 /* Noncoding model index            */
#else
#define NUMMODELS                     4 /* Model 0 (C g t),  (forward)      *
                                         *       1 (c G t),  (forward)      *
                                         *       2 (c g T),  (forward)      *
                                         *       3 noncoding (intron)       */
#define NONCODING                     3 /* Noncoding model index            */
#endif
#else
#define NUMMODELS                     2 /* Coding (0) and Noncoding (1)     */
#define NONCODING                     1 /* Noncoding model index            */
#endif
#define PfixOMEGA                  5E-5 /* Lower bound on adjustment for    *
                                         * NULLPROBs in PMFs (stat basis)   */
#define PI       3.14159265358979323846 /* pi--the one, the only            */
#define POSCONST                  10E-6 /* Small, pos constant (epsilon)    */
#ifdef QUANTSPEC
#define QTCT                          4 /* Number of quantiles              */
#endif
#define RELIABLECT                   50 /* Necessary counts to occur in     *
                                         * held-out set to make a D bucket  *
                                         * valid.                           */
#define SHIFTLEN                      1 /* Length to advance window by with *
                                         * quantile-specific operations     */
#ifndef TRUE
#define TRUE                          1 /* Truth value                      */
#endif
#define XI                     POSCONST /* for bottom-up DI                 */

/* Algorithm specifications **************************************************/

#define NUMALGOS   6 /* Number of training algorithms  *
                      * the library (ideally) supports */

/* Which dimension corresponds to which algorithm? */
enum {
  FIXORDX=0,
  TDDIX,
  BUDIX,
  CHI2X,
  DMMMX,
  PMMMX
};

/* The following char strings describe the training/testing *
 * algorithms currently implented in the library.           */
static const char algodesc[NUMALGOS][MAXLOCSTRING] = {
  "fixed-order (pseudocount smoothed)",
  "top-down deleted interpolation",
  "bottom-up deleted interpolation",
  "chi-square",
  "dynamically modulating",
  "partially modulating mixed model"
};

/* The following tags are appended to the binary parameter file, *
 * indicating from which statistical model the paramters were    *
 * actually generated.  Please ensure indices correspond to      *
 * those of algodesc!                                            */
static const char algotags[NUMALGOS][MAXLOCSTRING] = {
  "FIXORD",
  "TDDI",
  "BUDI",
  "CHI2",
  "DMMM",
  "PMMM"
};

/* Indices for DMMM methods */
enum {
  VARIANCEBASED=0,
  RANGEBASED
};

/* Structure definitions *****************************************************/

typedef struct {
  /* These arrays will store oligo frequencies, mononucleotides *
   * up through hexanucleotides.  These will be used to derive  *
   * max likelihood probabilities.  The NUMMODELS dimension may *
   * at first look unnecessary, but it IS required for counting *
   * three-periodically.                                        */
  int
    count1[NUMMODELS][ALFSIZE],
    count2[NUMMODELS][ALFSIZE][ALFSIZE],
    count3[NUMMODELS][ALFSIZE][ALFSIZE][ALFSIZE],
    count4[NUMMODELS][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE],
    count5[NUMMODELS][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE],
    count6[NUMMODELS][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE];
} countarrayT;

typedef struct {
  /* Arrays to store the relative frequencies (unsmoothed *
   * maximum likelihood estimates) of oligomers.          */
  double
    MLhis0[ALFSIZE],
    MLhis1[ALFSIZE][ALFSIZE],
    MLhis2[ALFSIZE][ALFSIZE][ALFSIZE],
    MLhis3[ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE],
    MLhis4[ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE],
    MLhis5[ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE];
} relfreqsT;

typedef struct {
  /* Arrays to store incrementally smoothed ML estimates in *
   * the bottom-up training procedure. The first dimension  *
   * keys the "innerorder" (i=k..0), and other dimensions   *
   * key the oligos being trained for.  Each array          *
   * corresponds to a given "order" (k).                    */ 
  double
    Psuphis0[1][ALFSIZE],
    Psuphis1[2][ALFSIZE][ALFSIZE],
    Psuphis2[3][ALFSIZE][ALFSIZE][ALFSIZE],
    Psuphis3[4][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE],
    Psuphis4[5][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE],
    Psuphis5[6][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE];
} PsupiT;

typedef struct {
  /* Arrays store final, possibly smoothed probabilities of *
   * oligomers of length 1 to MAXORDER + 1 for fast lookup  *
   * by the IMMprob function.                               */
  double
    prob1[NUMMODELS][ALFSIZE],
    prob2[NUMMODELS][ALFSIZE][ALFSIZE],
    prob3[NUMMODELS][ALFSIZE][ALFSIZE][ALFSIZE],
    prob4[NUMMODELS][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE],
    prob5[NUMMODELS][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE],
    prob6[NUMMODELS][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE];
} immprobT;

/* Structure stores fixed-order Markov model- *
 * derived probabilities for oligomers.       */
typedef struct {
  /* MAXORDER can vary from 0..5, inclusive. I use this same      *
   * data structure for any of these orders: for MAXORDER == 5    *
   * I key into the data structure using all 6 nucleotides, while *
   * for shorter oligos I "blank out" higher dimensional indices  *
   * with 0's. It's the programmer's responsibility to do this    *
   * consistently when referencing these data.                    */  
  double
    FOprob[NUMMODELS][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE][ALFSIZE];
} foprobT;

/* Structure stores arithmetic mean, standard deviation, min, and   *
 * max of hexamer probability estimates derived from partitioning   *
 * the training data according to some scheme, e.g., random splits, *
 * GC-content-based splits, etc.                                    */ 
typedef struct {
  /* Element names are self-explanatory */
  double
    mean[NUMMODELS][ALFSIZE][ALFSIZE][ALFSIZE]
                   [ALFSIZE][ALFSIZE][ALFSIZE],
    stddev[NUMMODELS][ALFSIZE][ALFSIZE][ALFSIZE]
                     [ALFSIZE][ALFSIZE][ALFSIZE],
    min[NUMMODELS][ALFSIZE][ALFSIZE][ALFSIZE]
                  [ALFSIZE][ALFSIZE][ALFSIZE],
    max[NUMMODELS][ALFSIZE][ALFSIZE][ALFSIZE]
                  [ALFSIZE][ALFSIZE][ALFSIZE],
    muless1sigma[NUMMODELS][ALFSIZE][ALFSIZE][ALFSIZE]
                           [ALFSIZE][ALFSIZE][ALFSIZE],
    muplus1sigma[NUMMODELS][ALFSIZE][ALFSIZE][ALFSIZE]
                           [ALFSIZE][ALFSIZE][ALFSIZE];
} dmprobT;

#ifdef QUANTSPEC
/* Structure stores boundaries of quantiles, e.g., *
 * min, max, median, etc., based on G+C content.   */
typedef struct {
  double codingbound[QTCT+1],
         noncodbound[QTCT+1];
} qtboundsT;
#endif

/* Public function prototypes ************************************************/

/* This function coordinates invocation of the appropriate *
 * training routine depending on which algorithm the user  *
 * wishes to employ for his/her data.                      */
int maintrain(int argc,char *argv[]);

#ifdef QUANTSPEC
/* Function to import a trained foprobT object array and a qtboundsT *
 * object (stored in a binary file) into main memory.                */
void import_fo_probsQT(char *filename,foprobT *probs,qtboundsT *bounds);

/* Function to import a trained immprobT object array and a     *
 * qtboundsT object (stored in a binary file) into main memory. */
void import_imm_probsQT(char *filename,immprobT *probs,qtboundsT *bounds);

/* Function to import a trained dmprobT object array and a qtboundsT *
 * object (stored in a binary file) into main memory.                */
void import_dm_probsQT(char *filename,dmprobT *probs,qtboundsT *bounds);

/* Function to return log-likelihood of a test input sequence, *
 * which must be null-terminated, under the model specified    *
 * using the model argument. It returns the probability of the *
 * input sequence given its first base occurs in the phase     *
 * specified by model, which parameters selected on a quantile *
 * specific basis. If testwindowgc is TRUE, this function will *
 * modulate between quantile-specific parameters as a function *
 * of window-based G+C composition.                            */
long double FOprobQT(char *input,foprobT *probs,qtboundsT *bounds,
  int model,int testwindowgc);

/* Function to return log-likelihood of a test input sequence, *
 * which must be null-terminated, under the model specified    *
 * using the model argument. It returns the probability of the *
 * input sequence given its first base occurs in the phase     *
 * specified by model, which parameters selected on a quantile *
 * specific basis. If testwindowgc is TRUE, this function will *
 * modulate between quantile-specific parameters as a function *
 * of window-based G+C composition.                            */
long double IMMprobQT(char *input,immprobT *probs,qtboundsT *bounds,
  int model,int testwindowgc);

/* Function to return log-likelihood of a test input sequence, *
 * which must be null-terminated, under the model specified    *
 * using the model argument. It returns the probability of the *
 * input sequence given its first base occurs in the phase     *
 * specified by model, which parameters selected on a quantile *
 * specific basis. If testwindowgc is TRUE, this function will *
 * modulate between quantile-specific parameters as a function *
 * of window-based G+C composition.                            */
long double DMMMprobQT(char *input,dmprobT *probs,qtboundsT *bounds,
  int cloudtype,int model,int testwindowgc);

#else

/* Function to import a trained foprobT     *
 * object (a binary file) into main memory. */
void import_fo_probs(char *filename,foprobT *probs);

/* Function to import a trained immprobT  *
 * object (binary file) into main memory. */
void import_imm_probs(char *filename,immprobT *probs);

/* Function to import a trained dmprobT     *
 * object (a binary file) into main memory. */
void import_dm_probs(char *filename,dmprobT *probs);

/* Function to return log-likelihood of a test input sequence, *
 * which must be null-terminated, under the model specified    *
 * using the model argument.  If threeper is FALSE, the        *
 * function returns the non-periodic log-probability of the    *
 * test sequence.  If non-FALSE, this function handles those   *
 * cases where model is an element of {0,1,2} if shadow is     *
 * FALSE, or {0,1,2,3,4,5} if shadow is non-FALSE; it returns  *
 * the probability of the input sequence given its first base  *
 * occurs in the phase specified by the model argument.        *
 * (see immpractical.h for phase nomenclature).                */
long double FOprob(char *input,foprobT *probs,int model,
  int threeper,int shadow);

/* Function to return log-likelihood of a test input sequence, *
 * which must be null-terminated, under the model specified    *
 * using the model argument.  If threeper is FALSE, the        *
 * function returns the non-periodic log-probability of the    *
 * test sequence.  If non-FALSE, this function handles those   *
 * cases where model is an element of {0,1,2} if shadow is     *
 * FALSE, or {0,1,2,3,4,5} if shadow is non-FALSE; it returns  *
 * the probability of the input sequence given its first base  *
 * occurs in the phase specified by the model argument.        *
 * (see immpractical.h for phase nomenclature).                */
long double IMMprob(char *input,immprobT *probs,int model,
  int threeper,int shadow);

/* Function to return log-likelihood of a test input sequence, *
 * which must be null-terminated, under the model specified    *
 * using the model argument.  If threeper is FALSE, the        *
 * function returns the non-periodic log-probability of the    *
 * test sequence.  If non-FALSE, this function handles those   *
 * cases where model is an element of {0,1,2} if shadow is     *
 * FALSE, or {0,1,2,3,4,5} if shadow is non-FALSE; it returns  *
 * the probability of the input sequence given its first base  *
 * occurs in the phase specified by the model argument.        *
 * (see immpractical.h for phase nomenclature).                */
long double DMMMprob(char *input,dmprobT *probs,int cloudtype,
  int model,int threeper,int shadow);
#endif

/* Prototypes for functions shared by (fo|dimm|chisquare)_utils.c ************/

/* Initialize the data structure for recording counts */
void init_counts(countarrayT *tal,int model);

#ifdef QUANTSPEC
/* Function to record counts occuring in input sequences,  *
 * nonperiodically. Please notice that each call requires  *
 * that the sequence stored in seq be at least 2*FLANK+    *
 * MAXORDER+1 residues long! This function, which develops *
 * the NONCODING model in the case of quantile-specific    *
 * library functionality, deposits counts of hexamers in   *
 * bins based on G+C composition of their local contexts   */
void get_countsQT(char *seq,countarrayT *tal,qtboundsT bounds);

/* This function will develop models [0-NONCODING) in the *
 * case of quantile-specific functionality. Note that all *
 * training data is assumed to be in frame 0!             */
void get_counts3percodQT(char *seq,countarrayT *tal,qtboundsT bounds);
#else
/* Function to record counts occuring in input sequences, *
 * nonperiodically. Please notice that each call requires *
 * that the sequence stored in seq be at least MAXORDER+1 *
 * residues long!                                         */
void get_counts(char *seq,countarrayT *tal,int model);

#ifdef THREEPERIODIC
/* This function will develop models [0-NONCODING) if the  *
 * library is intended to support three-periodicity.  Note *
 * that all training data is assumed to be in frame 0!     *
 * If you wish to tally counts for the reverse strand      *
 * (model indices {3,4,5}, when SHADOWSTRAND is defined),  *
 * in addition to models {0,1,2}, set the shadow argument  *
 * to TRUE.                                                */
void get_counts3percod(char *seq,countarrayT *tal,int shadow);
#endif
#endif

#ifdef EMPTYEQUIV
/* According to Jelinek 1998, the only histories of interest  *
 * are those s.t. Cd(h_k) > 0.  Given this application, we    *
 * want to assign a probability to ALL possible oligomers     *
 * that could be tested.  This subroutine will identify those *
 * oligos in the development data set that exhibit less than  *
 * EMPTYCOUNT counts, and appropriately propagate multiples   *
 * of MINCOUNTUNIT to ancestors and descendents of the        *
 * oligomer to keep counts consistent across orders. (Note:   *
 * only call this with the D tal array; develsize is not      *
 * used anymore, but its functionality may become desirable   *
 * again in the future.)                                      *
 *                                                            *
 * Note that using this approach is a slight departure from   *
 * the model presented in Jelinek's and Salzberg's reports!   */
void induce_pseudocounts(countarrayT *tal,int model,int *develsize);
#endif

#ifdef QUANTSPEC
/* Determines quantile boundaries based on overlapping windows */
void calc_quantiles(char *argv[],int algo,qtboundsT *bounds);
#endif

/* Subroutine to calculate unsmoothed maximum likelihood *
 * probabilities for oligomers.  (relative frequency)    *
 * Note: seqct and develsize are not used anymore, but   *
 * their functionality may become desirable again in     *
 * the future.                                           */
void get_rfreqs(relfreqsT *rfreqs,int model,
  countarrayT *tal,int seqct,int develsize);

/* This function oversees the building of a trained immprobT file     *
 * using the deleted interpolated Markov model approach of Jelinek or *
 * the chi**2 interpolated Markov model approach of Salzberg.         *
 * If THREEPERIODIC is set, the function oversees the building of a   *
 * trained immprobT file while taking three-periodicity into          *
 * account. (Note: Training the coding models assumes that all input  *
 * data is in phase 0!) Else, the function trains a model in terms of *
 * a binary classifier.                                               */
int imm_probmaker(int argc,char *argv[],int algo);

#endif
