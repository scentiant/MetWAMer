/* train_perceptron.c
 * Michael Sparks (mespar1@gmail.com)
 * Last modified : 20 July 2013
 *
 * This code implements train_perceptron, which trains
 * classification machinery, given a set of instances,
 * in an offline learning setting.
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

/* Include Statements */
#include <float.h>
#include <immpractical.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "classifiers.h"
#include "errors.h"

#define MAXRECORDLEN 2048 /* Max (string) length of feature vector inputs */
#define VERYLARGE (LDBL_MAX/1.075)

#define MINARGCT 6
#define BADUSE {\
  fprintf(stderr,"\a\n\
Conceptual Usage:\n\n\
  %s \\\n\
      classifier_algocode \\\n\
      learnrate \\\n\
      maxepochs \\\n\
      featvecs.dat \\\n\
      outfile_handle \\\n\
      [MCalgocode]\n\n\
  Where classifier_algocode is an integer selected from\n",argv[0]);\
  for(i=0;i<AVAILCLASSIFIERCT;++i)\
    fprintf(stderr,"    (%i) %s\n",i,classdesc[i]);\
  fprintf(stderr,"\n\
  MCalgocode can be optionally specified, allowing the output file\n\
    to be branded using the appropriate Markov chain method.\n\
    This code should be selected from\n");\
  for(i=0;i<NUMALGOS;++i)\
    if(i==PMMMX)\
      continue;\
    else\
      fprintf(stderr,"      (%i) %s\n",i,algodesc[i]);\
  fprintf(stderr,"\n");\
}

/* Main Application */
int main(int argc, char *argv[]) {
  char record[MAXRECORDLEN], /* Buffer for training instances      */
       symbol,               /* Parses token labeling data class   */
       outputname[2048];     /* Name of output file                */
  double target=1.0,         /* Desired output of machine          */
         output,             /* Actual output of machine           */
         learnrate;          /* Learning rate                      */
  featvecT fv;               /* Stores feature vector of instances */
  FILE *infile=NULL,         /* Connects to training instances     */
       *outfile=NULL;        /* Connects to parameter output file  */
  int classification,        /* TRUE -> (+), FALSE -> (-)          */
      classifiercode,        /* Designates classification machine  */
      maxepochs,             /* Max number of epochs in which to   *
                              * train neural weights               */
      mcalgocode,            /* Markov chain method identifier     */
      i;                     /* iterator variable                  */
  long double err_ltuWcurr;  /* Error given current weight vector  */
  LTUparmT ltuWcurr,         /* Weight vector from current epoch   */
           ltuWdelta;        /* Aux weight vector for linear units *
                              * trained in batch mode              */

  /* Verify that specified command line arguments are reasonable */
  if(argc<MINARGCT) {
    BADUSE;
    exit(EXIT_FAILURE);
  }

  /* Determine which classification algorithm to train for *
   * and initialize training objects as needed.            */
  classifiercode=atoi(argv[1]);
  switch(classifiercode) {
    case LINEARUNITX :
    case SIGMOIDUNITX :
      ltuWcurr.wzero=
      ltuWcurr.wmetwam=
      ltuWcurr.wmcratio=0.025; /* Small, random value :) */
      break;
    default :
      FATALERROR("\nErr (main): \
Improper classifier method requested!\n\n")
  }

  /* Parse initial learning rate.  For now, this is kept static *
   * during the training process, but a better situation would  *
   * be to implement the momentum modification trick.           */
  learnrate=atof(argv[2]);
  if(learnrate>1.0 || learnrate<0.0)
    FATALERROR("\nErr (main): \
Learning rate outside unit interval.\n\n")

  /* Parse maximum number of training epochs */
  maxepochs=atoi(argv[3]);

  /* Connect stream to training data */
  if((infile=fopen(argv[4],"rt"))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Can't seem to open %s!\n",argv[4]);
    FATALERROR(errbuff)
  }

  i=0;
  do {
    /* Compute sum of squares error as a function of this epoch's weights */
    err_ltuWcurr=0.0;

    if(classifiercode==LINEARUNITX) /* Flatten ltuWdelta */
      ltuWdelta.wzero=
      ltuWdelta.wmetwam=
      ltuWdelta.wmcratio=0.0;

    /* Iterate over training corpus */
    while(fgets(record,MAXRECORDLEN,infile)!=NULL) {
      sscanf(record,"(%c;<%le,%le>)",&symbol,&fv.metwam,&fv.mcratio);

      /* Parse class of instance and set target output appropriately */
      if(symbol=='+')
        target=+1.0;
      else if(symbol=='-')
        if(classifiercode==SIGMOIDUNITX)
          target=0.0;
        else
          target=-1.0;
      else {
        (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
\"%c\" denotes an ambiguous class demarcation!\n",symbol);
        FATALERROR(errbuff)
      }

      /* Calculate output given current state of classifier */
      output=eval_activation_unit(fv,ltuWcurr,&classification,
                                  &output,classifiercode);

      /* Update error-of-weights variable and check numeric stability. */
      err_ltuWcurr+=pow(target-output,2)/2.0;
      if(err_ltuWcurr>VERYLARGE||isinf(err_ltuWcurr)||isnan(err_ltuWcurr)) {
        #ifdef TRACKERROR
        NONFATALERROR("Warning: \
Exceptionally LARGE error!  Learning rate too large?\n")
        #endif
        err_ltuWcurr=VERYLARGE;
      }

      /* Record necessary updates */
      if(classifiercode==LINEARUNITX)
        update_weight_vector(fv,&ltuWdelta,target,output,
                             classifiercode,learnrate);
      else if(classifiercode==SIGMOIDUNITX)
        update_weight_vector(fv,&ltuWcurr,target,output,
                             classifiercode,learnrate);
      else
        FATALERROR("Err (main): Unexpected error.\n")
    } /* end instance while */

    ++i;
    #ifdef TRACKERROR
    (void)snprintf(errbuff,MAXERRLEN,"Error: \
%Le at epoch %i\n",err_ltuWcurr,i);
    NONFATALERROR(errbuff)
    #endif

    /* Step T4.2 of Mitchell 1997 */
    if(classifiercode==LINEARUNITX)
      ltuWcurr=add_weight_vectors(ltuWcurr,ltuWdelta);

    rewind(infile);

  } while(err_ltuWcurr>MAXERROFWEIGHTS && i<maxepochs);

  (void)fclose(infile);

  #ifdef TRACKERROR
  (void)snprintf(errbuff,MAXERRLEN,"Error converged to an approximate \
minimum (%Lf) after %i iterations over the dataset.\n\n",err_ltuWcurr,i);
  NONFATALERROR(errbuff)
  #endif

  /* Output learned classifier.       *
   * First, set name for output file: */
  (void)strcpy(outputname,argv[5]);
  if(argc>MINARGCT) {
    mcalgocode=atoi(argv[6]);
    if(mcalgocode>=0&&mcalgocode<(NUMALGOS)&&mcalgocode!=PMMMX) {
      (void)strcat(outputname,".");
      (void)strcat(outputname,algotags[mcalgocode]);
    }
    else
      NONFATALERROR("Warn (main): \
Inappropriate Markov chain specified (ignoring).\n")
  }
  (void)strcat(outputname,".");
  (void)strcat(outputname,classifiertags[classifiercode]);

  if((outfile=fopen(outputname,"wb"))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Can't seem to open %s!\n",outputname);
    FATALERROR(errbuff)
  }
  switch(classifiercode) {
    case LINEARUNITX :
    case SIGMOIDUNITX :
      if(fwrite(&ltuWcurr,sizeof(LTUparmT),1,outfile)!=1)
        FATALERROR("Err (main): \
Error writing model to binary output file!\n")
      break;
    default :
      FATALERROR("\nErr (main): \
Improper classifier method requested!\n\n")
  } /* end switch */
  (void)fclose(outfile);

  return(EXIT_SUCCESS);
} /* end main */
