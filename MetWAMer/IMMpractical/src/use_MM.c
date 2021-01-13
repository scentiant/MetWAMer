/* use_MM.c
 * Michael E Sparks (mespar1@gmail.com)
 *
 * This is a driver program for testing and demonstrating the use
 * of functions available in the immpractical library for
 * assessing likelihoods of input test sequences. It implements
 * methods to determine hypotheses with maximum likelihood,
 * MAP hypotheses, and implements the Genmark sequence
 * classification algorithm of Borodovsky and McIninch.
 *
 * Copyright (C) 2005,2006,2013 Michael E Sparks
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

/* Include statements */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "immpractical.h"

/* These macros will alert the user of the proper usage of    *
 * the driver. Note: algodesc array defined in immpractical.h */
#define USAGE "\a\n\
Conceptual Usage:\n\n\
  %s algorithm_type parmfile < FASTA_TESTDATA\n\n\
    Where algorithm_type is an integer selected from\n\n"

#define DMUSAGE "\a\n\
DMMM USAGE: %s %i cloud_type parms.%s < FASTA_TESTDATA\n\n\
  Where cloud_type = %i -> variance-based cloud\n\
                   = %i -> range-based cloud\n\n\
  and each sequence in FASTA_TESTDATA is contained on one line.\n\n"

#define PMUSAGE "\a\n\
PMMM USAGE: %s %i cloud_type parms.%s parms.%s < FASTA_TESTDATA\n\n\
  Where cloud_type = %i -> variance-based cloud\n\
                   = %i -> range-based cloud\n\n\
  and each sequence in FASTA_TESTDATA is contained on one line.\n\n"

#define BADUSE {                                   \
  fprintf(stderr,USAGE,argv[0]);                   \
  for(i=0;i<NUMALGOS;++i)                          \
    fprintf(stderr,"    (%i) %s\n",i,algodesc[i]); \
  fprintf(stderr,"\n    and each sequence in \
FASTA_TESTDATA is contained on one line.\n\n");    \
}

#define MAXLENGTH 1048576 /* Max linelen to accept from stdin */
#define MAX(a,b) ((a) > (b) ? (a) : (b)) 
#define VERYLARGE (LDBL_MAX / 1.075)
#define VERYSMALL (LDBL_MIN * 1.075)

/* main application */
int main(int argc,char *argv[])
{
  /* Note: NUMMODELS and the immprobT/foprobT struct def's are *
   * in immpractical.h.  Note that some of these variables may *
   * or may not be used, depending on what training algorithm  *
   * generated the parameterization to be used.                */
  char testseq[MAXLENGTH],      /* Stores test input sequences           */
       desc[MAXLENGTH],         /* Stores name of sequence               */
       *test=NULL;              /* Verifies file extension of parm file  */
  long double
     loglikerats[NUMMODELS],    /* Record log-likelihoods of each of the *
                                 * NUMMODELS models as a ratio relative  *
                                 * to the NONCODING model.               *
                                 * ( log(model X) - log(NONCODING) )     */
     loglikes[NUMMODELS],       /* Record log-likelihoods of each seq    */
     priors[NUMMODELS],         /* Stores prior probabilities            */
     scaleupfac,                /* Store average of log-odd ratios       */
     currmax=0.0,               /* Stores maximum probability from some  *
                                 * array of floats                       */
     posteriors[NUMMODELS],     /* Records posterior prob of each model  */
     constant_denom;            /* For Bayes rule; store denominator     */
  int length,                   /* Stores length of string               */
      i,                        /* Iterator variable                     */
      totalct=0,                /* Total number of test sequences        */
      hiloglikerats[NUMMODELS], /* Stores count of max log-likelihoods   */
      hiloglikeratsct=0,        /* Instances where this was feasible     */
      hiposteriors[NUMMODELS],  /* Stores count of max posteriors        */
      hiposteriorsct=0,         /* Instances where this was feasible     */
      nocando=0,                /* Bool from isnan/isinf tests           */
      algotype,                 /* Stores training algo used for parms   */
      cloudtype=-1;             /* Stores cloud definition method        */
  #if defined(THREEPERIODIC) && defined(SHADOWSTRAND)
  int genmarkcodct=0,           /* Tabulate # coding calls by genmark    */
      genmarkcandoct=0;         /* # sequences genmark could classify    */
  const double genmark_cut=0.5; /* For making classification decision    *
                                 * under the genmark algorithm           */
  long double genmark_accum;    /* Store coding hypotheses' posteriors   */
  #endif
  #ifdef QUANTSPEC
  qtboundsT
    qt_bounds,                  /* Read qtboundsT object into it         */
    qt_bounds4PMMMnoncod;       /* For PMMMs, store the qtboundsT object *
                                 * needed for the NONCODING model        */
  dmprobT *dmprobs=NULL;        /* Read in the dmprobT object array      */
  immprobT *immprobs=NULL;      /* Read in the immprobT object array     */
  foprobT *fixordprobs=NULL;    /* Read in the foprobT object array      */
  #else
  dmprobT dmprobs;              /* Read the dmprobT object into it       */
  immprobT immprobs;            /* Read the immprobT object into it      */
  foprobT fixordprobs;          /* Read the foprobT object into it       */
  #endif

  /* Verify command line parameters */
  if(argc < 2) {
    BADUSE;
    exit(EXIT_FAILURE); 
  }
  else
    algotype=atoi(argv[1]);

  /* Read in the parameter object (algotags defined in         *
   * immpractical.h, import_imm_probs defined in dimm_utils.c) */
  switch(algotype) {
    case FIXORDX :
    case TDDIX :
    case BUDIX :
    case CHI2X :
      if(argc < 3) {
        BADUSE;
        exit(EXIT_FAILURE); 
      }
      else if((test=strstr(argv[2],algotags[algotype]))==NULL) {
        fprintf(stderr,"Err (main): \
Binary parm file has improper extension!\n");
        exit(EXIT_FAILURE);
      }
      else {
        if(algotype == FIXORDX) {
          #ifdef QUANTSPEC
          if((fixordprobs=(foprobT*)malloc(sizeof(foprobT)*(int)QTCT))!=NULL)
            import_fo_probsQT(argv[2],fixordprobs,&qt_bounds);
          else {
            fprintf(stderr,"Err (main): Out of memory!\n");
            exit(EXIT_FAILURE);
          }
          #else
          import_fo_probs(argv[2],&fixordprobs);
          #endif
        }
        else { /* an IMM variant */
          #ifdef QUANTSPEC
          if((immprobs=(immprobT*)malloc(sizeof(immprobT)*(int)QTCT))!=NULL)
            import_imm_probsQT(argv[2],immprobs,&qt_bounds);
          else {
            fprintf(stderr,"Err (main): Out of memory!\n");
            exit(EXIT_FAILURE);
          }
          #else
          import_imm_probs(argv[2],&immprobs);
          #endif
        }
      }
      break;
    case DMMMX :
      /* Quick check: if user is deploying the DMMM, determine if they  *
       * want to use cloud bounded by range- or variance-based methods. */
      if(argc < 4) {
        fprintf(stderr,DMUSAGE,argv[0],DMMMX,algotags[DMMMX],
          VARIANCEBASED,RANGEBASED);
        exit(EXIT_FAILURE); 
      }
      else
        cloudtype=atoi(argv[2]);

      if((test=strstr(argv[3],algotags[DMMMX]))==NULL) {
        fprintf(stderr,"Err (main): \
Binary parm file has improper extension!\n");
        fprintf(stderr,DMUSAGE,argv[0],DMMMX,algotags[DMMMX],
          VARIANCEBASED,RANGEBASED);
        exit(EXIT_FAILURE);
      }
      else {
        #ifdef QUANTSPEC
        if((dmprobs=(dmprobT*)malloc(sizeof(dmprobT)*(int)QTCT))!=NULL)
          import_dm_probsQT(argv[3],dmprobs,&qt_bounds);
        else {
          fprintf(stderr,"Err (main): Out of memory!\n");
          exit(EXIT_FAILURE);
        }
        #else
        import_dm_probs(argv[3],&dmprobs);
        #endif
      }
      break;
    case PMMMX :
      if(argc < 5) {
        fprintf(stderr,PMUSAGE,argv[0],PMMMX,algotags[DMMMX],
          algotags[FIXORDX],VARIANCEBASED,RANGEBASED);
        exit(EXIT_FAILURE); 
      }
      else
        cloudtype=atoi(argv[2]);

      if((test=strstr(argv[3],algotags[DMMMX]))==NULL ||
         (test=strstr(argv[4],algotags[FIXORDX]))==NULL) {
        fprintf(stderr,"Err (main): \
Binary parm file has improper extension!\n");
        fprintf(stderr,PMUSAGE,argv[0],PMMMX,algotags[DMMMX],
          algotags[FIXORDX],VARIANCEBASED,RANGEBASED);
        exit(EXIT_FAILURE); 
      }
      else {
        #ifdef QUANTSPEC
        if(((dmprobs=(dmprobT*)malloc(sizeof(dmprobT)*(int)QTCT))!=NULL) &&
           ((fixordprobs=(foprobT*)malloc(sizeof(foprobT)*(int)QTCT))!=NULL)) {
          import_dm_probsQT(argv[3],dmprobs,&qt_bounds);
          import_fo_probsQT(argv[4],fixordprobs,&qt_bounds4PMMMnoncod);
        }
        else {
          fprintf(stderr,"Err (main): Out of memory!\n");
          exit(EXIT_FAILURE);
        }
        #else
        import_dm_probs(argv[3],&dmprobs);
        import_fo_probs(argv[4],&fixordprobs);
        #endif
      }
      break;
    default :
      fprintf(stderr,"\nErr (main): \
Improper algorithm requested!\n\n");
      exit(EXIT_FAILURE);
  } /* end switch */

  /* Echo invocation to output. */
  fprintf(stdout,"Invoked as:");
  for(i=0;i<argc;++i)
    fprintf(stdout," %s",argv[i]);
  fprintf(stdout,"\n\n");

  /* see immpractical.h */
  fprintf(stdout,"Note: Model nomenclature for this report\n");
  #ifdef THREEPERIODIC
  #ifndef SHADOWSTRAND
  fprintf(stdout,"Model 0 (C g t), 1st codon position, forward strand\n");
  fprintf(stdout,"      1 (c G t), 2nd codon position, forward strand\n");
  fprintf(stdout,"      2 (c g T), 3rd codon position, forward strand\n");
  fprintf(stdout,"      3 noncoding (intron)\n\n");
  #else
  fprintf(stdout,"Model 0 (C g t), 1st codon position, forward strand\n");
  fprintf(stdout,"      1 (c G t), 2nd codon position, forward strand\n");
  fprintf(stdout,"      2 (c g T), 3rd codon position, forward strand\n");
  fprintf(stdout,"      3 (A c g), 1st codon position, reverse strand\n");
  fprintf(stdout,"      4 (a C g), 2nd codon position, reverse strand\n");
  fprintf(stdout,"      5 (a c G), 3rd codon position, reverse strand\n");
  fprintf(stdout,"      6 noncoding (intron)\n\n");
  #endif
  #else
  fprintf(stdout,"Model 0 coding\n");
  fprintf(stdout,"      1 noncoding\n\n");
  #endif

  /* intialize array of priors: under the assumption that the prior    *
   * probability of noncoding model is 1/2, and whatever number of     *
   * alternative coding hypotheses must equally divvy up the remaining *
   * 1/2 probability.                                                  */
  fprintf(stdout,"Distribution of Priors:\n");
  for(i=0;i<(NUMMODELS-1);++i) {
    priors[i]=1.0 / ((NUMMODELS-1)*2);
    fprintf(stdout,"\t(model %i): %Lf\n",i,priors[i]);
  }
  priors[NONCODING]=0.5;
  fprintf(stdout,"\t(model %i): %Lf\n\n",NONCODING,priors[NONCODING]);

  /* initialize counting arrays */
  for(i=0;i<NUMMODELS;++i) {
    hiloglikerats[i]=0;
    hiposteriors[i]=0;
  }

  /* Appraise coding potential of test sequences */
  while(fgets(desc,MAXLENGTH,stdin)!=NULL) {
    /* Get sequence data, remove newline, and append null terminator *
     * Caveat emptor: the entire test sequence is on 1 line!         */
    if(fgets(testseq,MAXLENGTH,stdin)==NULL) {
      fprintf(stderr,"Err (main): Encountered a description, \
but no sequence!\n\t(%s)\n",desc);
      return(EXIT_FAILURE);
    }
    length=(int)strlen(testseq);
    testseq[--length]='\0';

    /* Compute likelihoods (Pr{data|model}) */
    switch(algotype) {
      case FIXORDX :
        /* Calculate the overall probability of          *
         * the sequence under all NUMMODELS models.      *
         * FOprob and FOprobQT are defined in fo_utils.c */
        #ifdef QUANTSPEC
        #ifdef TESTWINDOWGC
        for(i=0;i<NUMMODELS;++i)
          loglikes[i]=FOprobQT(testseq,fixordprobs,&qt_bounds,i,TRUE);
        #else
        for(i=0;i<NUMMODELS;++i)
          loglikes[i]=FOprobQT(testseq,fixordprobs,&qt_bounds,i,FALSE);
        #endif
        #else
        #ifdef THREEPERIODIC 
        /* See immpractical.h for what NONCODING *
         * means when THREEPERIODIC is defined.  */
        for(i=0;i<3;++i)
          loglikes[i]=FOprob(testseq,&fixordprobs,i,TRUE,FALSE);
        #ifdef SHADOWSTRAND
        for(i=3;i<NONCODING;++i)
          loglikes[i]=FOprob(testseq,&fixordprobs,i,TRUE,TRUE);
        #endif
        loglikes[NONCODING]=
          FOprob(testseq,&fixordprobs,NONCODING,FALSE,FALSE);
        #else
        for(i=0;i<NUMMODELS;++i)
          loglikes[i]=FOprob(testseq,&fixordprobs,i,FALSE,FALSE);
        #endif
        #endif
        break;
      case TDDIX :
      case BUDIX :
      case CHI2X :
        /* Calculate the overall probability of     *
         * the sequence under all NUMMODELS models. *
         * IMMprob is defined in immpractical.c     */
        #ifdef QUANTSPEC
        #ifdef TESTWINDOWGC
        for(i=0;i<NUMMODELS;++i)
          loglikes[i]=IMMprobQT(testseq,immprobs,&qt_bounds,i,TRUE);
        #else
        for(i=0;i<NUMMODELS;++i)
          loglikes[i]=IMMprobQT(testseq,immprobs,&qt_bounds,i,FALSE);
        #endif
        #else
        #ifdef THREEPERIODIC 
        /* See immpractical.h for what NONCODING *
         * means when THREEPERIODIC is defined.  */
        for(i=0;i<3;++i)
          loglikes[i]=IMMprob(testseq,&immprobs,i,TRUE,FALSE);
        #ifdef SHADOWSTRAND
        for(i=3;i<NONCODING;++i)
          loglikes[i]=IMMprob(testseq,&immprobs,i,TRUE,TRUE);
        #endif
        loglikes[NONCODING]=
          IMMprob(testseq,&immprobs,NONCODING,FALSE,FALSE);
        #else
        for(i=0;i<NUMMODELS;++i)
          loglikes[i]=IMMprob(testseq,&immprobs,i,FALSE,FALSE);
        #endif
        #endif
        break;
      case DMMMX :
        /* Calculate the overall probability of     *
         * the sequence under all NUMMODELS models. *
         * DMMMprob is defined in dm_utils.c        */
        #ifdef QUANTSPEC
        #ifdef TESTWINDOWGC
        for(i=0;i<NUMMODELS;++i)
          loglikes[i]=
            DMMMprobQT(testseq,dmprobs,&qt_bounds,cloudtype,i,TRUE);
        #else
        for(i=0;i<NUMMODELS;++i)
          loglikes[i]=
            DMMMprobQT(testseq,dmprobs,&qt_bounds,cloudtype,i,FALSE);
        #endif
        #else
        #ifdef THREEPERIODIC 
        /* See immpractical.h for what NONCODING *
         * means when THREEPERIODIC is defined.  */
        for(i=0;i<3;++i)
          loglikes[i]=DMMMprob(testseq,&dmprobs,cloudtype,i,TRUE,FALSE);
        #ifdef SHADOWSTRAND
        for(i=3;i<NONCODING;++i)
          loglikes[i]=DMMMprob(testseq,&dmprobs,cloudtype,i,TRUE,TRUE);
        #endif
        loglikes[NONCODING]=DMMMprob(testseq,&dmprobs,cloudtype,
                                 NONCODING,FALSE,FALSE);
        #else
        for(i=0;i<NUMMODELS;++i)
          loglikes[i]=DMMMprob(testseq,&dmprobs,cloudtype,i,FALSE,FALSE);
        #endif
        #endif
        break;
      case PMMMX :
        /* For partially modulating Markov models: Here, the coding   *
         * model(s) operates under a DMMM-type framework (fixed-order *
         * based), while the NONCODING model uses a point estimate    *
         * (fixed-order Markov model).                                */
        #ifdef QUANTSPEC
        #ifdef TESTWINDOWGC
        for(i=0;i<NONCODING;++i)
          loglikes[i]=
            DMMMprobQT(testseq,dmprobs,&qt_bounds,cloudtype,i,TRUE);
        loglikes[NONCODING]=
          FOprobQT(testseq,fixordprobs,&qt_bounds,NONCODING,TRUE);
        #else
        for(i=0;i<NONCODING;++i)
          loglikes[i]=
            DMMMprobQT(testseq,dmprobs,&qt_bounds,cloudtype,i,FALSE);
        loglikes[NONCODING]=
          FOprobQT(testseq,fixordprobs,&qt_bounds,NONCODING,FALSE);
        #endif
        #else
        #ifdef THREEPERIODIC 
        /* See immpractical.h for what NONCODING *
         * means when THREEPERIODIC is defined.  */
        for(i=0;i<3;++i)
          loglikes[i]=DMMMprob(testseq,&dmprobs,cloudtype,i,TRUE,FALSE);
        #ifdef SHADOWSTRAND
        for(i=3;i<NONCODING;++i)
          loglikes[i]=DMMMprob(testseq,&dmprobs,cloudtype,i,TRUE,TRUE);
        #endif
        #else
        loglikes[0]=DMMMprob(testseq,&dmprobs,cloudtype,0,FALSE,FALSE);
        #endif
        loglikes[NONCODING]=
          FOprob(testseq,&fixordprobs,NONCODING,FALSE,FALSE);
        #endif
        break;
      default :
        fprintf(stderr,"Err (main): Improper algorithm requested!\n");
        return(EXIT_FAILURE);
    } /* end switch */

    fprintf(stdout,"Test sequence (length %i bp): %s",length,desc);

    /* Report these log-likelihoods */
    fprintf(stdout,"Reporting log-likelihoods:\n");
    for(i=0;i<NUMMODELS;++i)
      fprintf(stdout,"\tlog Pr(data|model %d) = %Lf\n",i,loglikes[i]);

    /* Now, update loglikerats as ratios relative to NONCODING model *
     * Note: log(A) - log(B) = log(A/B)                              */
    for(i=0;i<NUMMODELS;++i)
      loglikerats[i] = loglikes[i] - loglikes[NONCODING];

    /* Report these log-likelihood ratios */
    fprintf(stdout,"Reporting log-likelihood ratios:\n");
    for(i=0;i<NUMMODELS;++i) {
      fprintf(stdout,
        "\tlog( Pr(data|model %d) / Pr(data|NONCODING) ) = %Lf\n",
        i,loglikerats[i]);
      fprintf(stdout,"\t\t[ exp(%Lf) = %E ]\n",
        loglikerats[i],exp(loglikerats[i]));
    }

    /* Record best result in hiloglikerats array, if possible. */
    if(isnan(loglikerats[0]) || isinf(loglikerats[0]))
      nocando=1;
    else {
      nocando=0;
      for(i=1,currmax=loglikerats[0];i<NUMMODELS;++i) {
        if(isnan(loglikerats[i]) || isinf(loglikerats[i])) {
          nocando=1;
          break;
        }
        else
          currmax=MAX(currmax,loglikerats[i]);
      }
    }
    if(!nocando) {
      /* Now, find which index (model) currmax *
       * (almost certainly!) came from         */
      for(i=0;i<NUMMODELS;++i) { 
        if(fabs(currmax-loglikerats[i]) < DBL_EPSILON) {
          ++hiloglikerats[i];
          break;
        }
      }
      ++hiloglikeratsct;
    }

    /* Scaling up the Bayes formula by a factor of constant scaleupfac *
     * should *help* to prevent buffer under/overflow errors.          */
    for(i=0,scaleupfac=0.0;i<NONCODING;++i)
      scaleupfac += loglikerats[i] / NUMMODELS;

    for(i=0,constant_denom=0.0;i<NUMMODELS;++i) {
      constant_denom += priors[i] * exp(loglikerats[i] - scaleupfac);

      if(constant_denom > VERYLARGE)
        fprintf(stderr,"Warning: \
constant_denom approaching long double precision upper limit!\n");
      /* A very small denominator is not expected, but... */
      if(constant_denom < VERYSMALL)
        fprintf(stderr,"Warning: \
constant_denom approaching long double precision lower limit!\n");
    }

    fprintf(stdout,"Reporting posterior probabilities:\n");
    for(i=0;i<NUMMODELS;++i) {
      posteriors[i] = (priors[i] * exp(loglikerats[i] - scaleupfac)) /
        constant_denom;

      if(posteriors[i] < VERYSMALL)
        fprintf(stderr,"Warning: \
Posterior prob approaching long double precision lower limit!\n");
      /* A very large posterior prob is not expected, but... */
      if(posteriors[i] > VERYLARGE)
        fprintf(stderr,"Warning: \
Posterior prob approaching long double precision upper limit!\n");

      fprintf(stdout,"\tPr(model %i|data) = %Lf\n",i,posteriors[i]);
    }

    /* Find greatest posterior probability (MAP hypothesis), if possible */
    if(isnan(posteriors[0]) || isinf(posteriors[0]))
      nocando=1;
    else {
      nocando=0;
      for(i=1,currmax=posteriors[0];i<NUMMODELS;++i) {
        if(isnan(posteriors[i]) || isinf(posteriors[i])) {
          nocando=1;
          break;
        }
        else
          currmax=MAX(currmax,posteriors[i]);
      }
    }
    if(!nocando) {
      /* Now, find which index (model) currmax *
       * (almost certainly!) came from         */
      for(i=0;i<NUMMODELS;++i) { 
        if(fabs(currmax-(double)posteriors[i]) < DBL_EPSILON) {
          ++hiposteriors[i];
          break;
        }
      }
      ++hiposteriorsct;
    }

    #if defined(THREEPERIODIC) && defined(SHADOWSTRAND)
    /* Genmark classification */
    genmark_accum=0.0;
    nocando=0;
    for(i=0;i<NONCODING;++i) {
      if(isnan(posteriors[i]) || isinf(posteriors[i])) {
        nocando=1;
        break;
      }
      else
        genmark_accum += posteriors[i];
    }
    if(!nocando && (isnan(posteriors[NONCODING]) || isinf(posteriors[i])))
      nocando=1;
          
    fprintf(stdout,"Genmark algorithm classification: ");
    if(nocando)
      fprintf(stdout,"Intractable\n");
    else {
      if (genmark_accum > genmark_cut) {
        fprintf(stdout,"Coding\n");
        ++genmarkcodct;
      }
      else
        fprintf(stdout,"Noncoding\n");

      ++genmarkcandoct;
    }
    #endif
    fprintf(stdout,"\n");

    ++totalct;
  } /* end while */

  #ifdef QUANTSPEC
  /* Free up whatever probability array was allocated */
  switch(algotype) {
    case FIXORDX :
      free(fixordprobs);
      break;
    case TDDIX :
    case BUDIX :
    case CHI2X :
      free(immprobs);
      break;
    case DMMMX :
    case PMMMX :
      free(dmprobs);
      if(algotype == PMMMX)
        free(fixordprobs);
      break;
    default :
      fprintf(stderr,"\nErr (main): \
Improper algorithm requested!\n\n");
      exit(EXIT_FAILURE);
  } /* end switch */
  #endif

  /* Report overall performance statistics */
  fprintf(stdout,"--Overall performance statistics:\n\n");
  fprintf(stdout,"%i sequences were tested in total.\n\n",totalct);

  for(i=0;i<NUMMODELS;++i)
    fprintf(stdout,"Model %i had greatest \
log-likelihood ratio %i times. (%.2f%%)\n",
      i,hiloglikerats[i],(hiloglikerats[i]/(double)hiloglikeratsct)*100);
  fprintf(stdout,"Could assess %i sequences \
of %i total (%.2f%%) using this method.\n",
    hiloglikeratsct,totalct,(hiloglikeratsct/(double)totalct)*100);
  fprintf(stdout,"\n");

  for(i=0;i<NUMMODELS;++i)
    fprintf(stdout,"Model %i had greatest \
posterior probability %i times. (%.2f%%)\n",
      i,hiposteriors[i],(hiposteriors[i]/(double)hiposteriorsct)*100);
  fprintf(stdout,"Could assess %i sequences \
of %i total (%.2f%%) using this method.\n",
    hiposteriorsct,totalct,(hiposteriorsct/(double)totalct)*100);
  fprintf(stdout,"\n");

  #if defined(THREEPERIODIC) && defined(SHADOWSTRAND)
  fprintf(stdout,"Genmark Coding classifications: %.2f%%\n",
    (genmarkcodct/(double)genmarkcandoct)*100);
  fprintf(stdout,"Genmark Noncoding classifications: %.2f%%\n",
    100.0-((genmarkcodct/(double)genmarkcandoct)*100));
  fprintf(stdout,"Could assess %i sequences \
of %i total (%.2f%%) using this method.\n",
    genmarkcandoct,totalct,(genmarkcandoct/(double)totalct)*100);
  #endif

  return(EXIT_SUCCESS);
} /* end main */
