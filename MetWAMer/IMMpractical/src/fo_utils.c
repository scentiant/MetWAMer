/* fo_utils.c
 * Michael E Sparks (mespar1@gmail.com)
 *
 * This file contains code for developing/using fixed-order
 * Markov models.
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
 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fo_utils.h"
#include "immpractical.h"
#include "sequence_parse.h"

/* Function definitions ******************************************************/

/* This function oversees the building of a trained *
 * foprobT file based on fixed-order Markov models. */
int fo_probmaker(int argc,char *argv[])
{
  FILE *infile=NULL,         /* For processing Fasta files           */
       *outfile=NULL;        /* Will store binary fo_probT object    */
  char *sequence=NULL,       /* For storing current sequence         */
    outputname[MAXFILENAME]; /* Stores final name of output file     */
  int seqlength,             /* Records length of current sequence   */
      model,                 /* Model iterator variable              */
      currentfile;           /* Iterator based on input file at hand */
  #ifdef QUANTSPEC
  int quantind;               /* iterator variable                   */
  countarrayT tallies[QTCT];  /* Record oligomer frequencies         */
  foprobT fo_parameters[QTCT];/* Stores final, pseudocount smoothed, *
                               * fixed-order probabilities           */
  qtboundsT qt_bounds;        /* Store quantile boundaries           */
  #else
  countarrayT tallies;       /* Record oligomer frequencies          */
  foprobT fo_parameters;     /* Stores final, pseudocount smoothed,  *
                              * fixed-order probabilities            */
  #endif

  #ifdef QUANTSPEC
  calc_quantiles(argv,FIXORDX,&qt_bounds);
  for(quantind=0;quantind<QTCT;++quantind) {
    /* Initialize countarray */
    for(model=0;model<NUMMODELS;++model)
      init_counts(&tallies[quantind],model);
  }
  #else
  /* Initialize countarray */
  for(model=0;model<NUMMODELS;++model)
    init_counts(&tallies,model);
  #endif

  /* Do the counting.                         *
   * currentfile 0 = coding                   *
   *             1 = noncoding                *
   * Note: Training the coding models assumes *
   * that all input data is in phase 0!       */
  for(currentfile=0;currentfile<2;++currentfile) {
    /* Derive counts from data */
    if((infile=fopen(argv[currentfile+2],"rt"))==NULL) {
      fprintf(stderr,"Err (fo_probmaker): Can't seem to open %s!\n",
        argv[currentfile+2]);
      exit(EXIT_FAILURE);
    }
    while ((sequence=get_fasta(infile,sequence,&seqlength))!=NULL) {
      if(seqlength >= (MAXORDER+1)) {
        switch (currentfile) {
          case 0 :
            #ifdef QUANTSPEC
            /* Quantile-specific behavior demands three-periodic modeling *
             * of both strands---locked in a NUMMODELS == 7 scheme.       */
            get_counts3percodQT(sequence,tallies,qt_bounds);
            #else
            #ifdef THREEPERIODIC
            #ifdef SHADOWSTRAND
            get_counts3percod(sequence,&tallies,TRUE);
            #else
            get_counts3percod(sequence,&tallies,FALSE);
            #endif
            #else
            get_counts(sequence,&tallies,currentfile);
            #endif
            #endif
            break;
          case 1 :
            #ifdef QUANTSPEC
            /* As NUMMODELS == 7, no need to specify NONCODING index */
            get_countsQT(sequence,tallies,qt_bounds);
            #else
            #ifdef THREEPERIODIC
            get_counts(sequence,&tallies,NONCODING);
            #else
            /* currentfile == NONCODING */
            get_counts(sequence,&tallies,currentfile);
            #endif
            #endif
            break;
          default :
            fprintf(stderr,"Err (fo_probmaker): You seem to have \
passed an invalid training file?\n");
            (void)fclose(infile);
            return(FALSE);
        }
        free(sequence);
        sequence=NULL;
      }
      else {
        fprintf(stderr,"Err (fo_probmaker): First run %s on your input!\n",
          DATACLEANPL);
        free(sequence);
        (void)fclose(infile);
        return(FALSE);
      }
    } /* end while */
    (void)fclose(infile);
  } /* end currentfile for */

  for(model=0;model<NUMMODELS;++model) {
    #ifdef QUANTSPEC
    for(quantind=0;quantind<QTCT;++quantind) {
      /* Smooth empty counts with pseudocounts. *
       * Since we're only interested in order   *
       * 5, we'll only add pseudocounts for it. */
      fo_pseudo(&tallies[quantind],model,TRUE);

      /* Record the "smoothed" probability */
      fo_rfreqs(&tallies[quantind],model,&fo_parameters[quantind]);
    }
    #else
    /* Smooth empty counts with pseudocounts. *
     * Since we're only interested in order   *
     * 5, we'll only add pseudocounts for it. */
    fo_pseudo(&tallies,model,TRUE);

    /* Record the "smoothed" probability */
    fo_rfreqs(&tallies,model,&fo_parameters);
    #endif
  }

  /* Set name for output file. */
  if((int)strlen(argv[4]) > MAXFILENAME) {
    fprintf(stderr,
      "Err (fo_probmaker): Filename exceeds %i characters!\n",MAXFILENAME);
    return(FALSE);
  }
  else {
    (void)strcpy(outputname,argv[4]);
    (void)strcat(outputname,".");
    (void)strcat(outputname,algotags[atoi(argv[1])]);
  }

  /* Write our newly trained model to the output file */
  if((outfile=fopen(outputname,"wb"))==NULL) {
    fprintf(stderr,"Err (fo_probmaker): \
Can't seem to open %s!\n",outputname);
    exit(EXIT_FAILURE);
  }
  #ifdef QUANTSPEC
  else if(fwrite(fo_parameters,sizeof(foprobT),QTCT,outfile) != QTCT) {
    fprintf(stderr,
      "Err (fo_probmaker): \
Error writing model to binary output file!\n");
    exit(EXIT_FAILURE);
  }
  else if(fwrite(&qt_bounds,sizeof(qtboundsT),1,outfile) != 1) {
    fprintf(stderr,
      "Err (fo_probmaker): \
Error writing quantile boundaries to binary output file!\n");
    exit(EXIT_FAILURE);
  }
  #else
  else if(fwrite(&fo_parameters,sizeof(foprobT),1,outfile) != 1) {
    fprintf(stderr,
      "Err (fo_probmaker): \
Error writing model to binary output file!\n");
    exit(EXIT_FAILURE);
  }
  #endif
  else
    (void)fclose(outfile);

  return(TRUE);
} /* end fo_probmaker */

/* This function adds pseudocounts for oligos that *
 * don't occur in the training data at a great     *
 * enough frequency.  It's a very primitive form   *
 * of parameter smoothing. If scaleall is set to   *
 * TRUE, the complete count matrix will be scaled  *
 * up FOPSEUDOCT units; else, the pseudocount      *
 * induction method will be more surgical.         */
void fo_pseudo(countarrayT *tal,int model,int scaleall)
{
  int i,j,k,l,m,n; /* Iterator variables */

  /* The switch statement facilitates using the  *
   * library for MAXORDER in the interval [0..5] */
  switch(MAXORDER) {
    case(0) :
      for(i=0;i<ALFSIZE;++i) {
        if(scaleall == TRUE)
          tal->count1[model][i]+=FOPSEUDOCT; 
        else{
          if (tal->count1[model][i] <= EMPTYCOUNT) {
            fprintf(stdout,"Oligo %u only \
occurs %u times in data! (model %i)\n",
            i,tal->count1[model][i],model);
            tal->count1[model][i]=FOPSEUDOCT; 
          }
        }
      }
      break;

    case(1) :
      for(i=0;i<ALFSIZE;++i) {
      for(j=0;j<ALFSIZE;++j) {
        if(scaleall == TRUE)
          tal->count2[model][i][j]+=FOPSEUDOCT; 
        else{
          if (tal->count2[model][i][j] <= EMPTYCOUNT) {
            fprintf(stdout,"Oligo %u%u only \
occurs %u times in data! (model %i)\n",
            i,j,tal->count2[model][i][j],model);
            tal->count2[model][i][j]=FOPSEUDOCT; 
          }
        }
      }}
      break;

    case(2) :
      for(i=0;i<ALFSIZE;++i) {
      for(j=0;j<ALFSIZE;++j) {
      for(k=0;k<ALFSIZE;++k) {
        if(scaleall == TRUE)
          tal->count3[model][i][j][k]+=FOPSEUDOCT; 
        else{
          if (tal->count3[model][i][j][k] <= EMPTYCOUNT) {
            fprintf(stdout,"Oligo %u%u%u only \
occurs %u times in data! (model %i)\n",
            i,j,k,tal->count3[model][i][j][k],model);
            tal->count3[model][i][j][k]=FOPSEUDOCT; 
          }
        }
      }}}
      break;

    case(3) :
      for(i=0;i<ALFSIZE;++i) {
      for(j=0;j<ALFSIZE;++j) {
      for(k=0;k<ALFSIZE;++k) {
      for(l=0;l<ALFSIZE;++l) {
        if(scaleall == TRUE)
          tal->count4[model][i][j][k][l]+=FOPSEUDOCT; 
        else{
          if (tal->count4[model][i][j][k][l] <= EMPTYCOUNT) {
            fprintf(stdout,"Oligo %u%u%u%u only \
occurs %u times in data! (model %i)\n",
            i,j,k,l,tal->count4[model][i][j][k][l],model);
            tal->count4[model][i][j][k][l]=FOPSEUDOCT; 
          }
        }
      }}}}
      break;

    case(4) :
      for(i=0;i<ALFSIZE;++i) {
      for(j=0;j<ALFSIZE;++j) {
      for(k=0;k<ALFSIZE;++k) {
      for(l=0;l<ALFSIZE;++l) {
      for(m=0;m<ALFSIZE;++m) {
        if(scaleall == TRUE)
          tal->count5[model][i][j][k][l][m]+=FOPSEUDOCT; 
        else{
          if (tal->count5[model][i][j][k][l][m] <= EMPTYCOUNT) {
            fprintf(stdout,"Oligo %u%u%u%u%u only \
occurs %u times in data! (model %i)\n",
            i,j,k,l,m,tal->count5[model][i][j][k][l][m],model);
            tal->count5[model][i][j][k][l][m]=FOPSEUDOCT; 
          }
        }
      }}}}}
      break;

    case(5) :
      for(i=0;i<ALFSIZE;++i) {
      for(j=0;j<ALFSIZE;++j) {
      for(k=0;k<ALFSIZE;++k) {
      for(l=0;l<ALFSIZE;++l) {
      for(m=0;m<ALFSIZE;++m) {
      for(n=0;n<ALFSIZE;++n) {
        if(scaleall == TRUE)
          tal->count6[model][i][j][k][l][m][n]+=FOPSEUDOCT; 
        else{
          if (tal->count6[model][i][j][k][l][m][n] <= EMPTYCOUNT) {
            fprintf(stdout,"Oligo %u%u%u%u%u%u only \
occurs %u times in data! (model %i)\n",
            i,j,k,l,m,n,tal->count6[model][i][j][k][l][m][n],model);
            tal->count6[model][i][j][k][l][m][n]=FOPSEUDOCT; 
          }
        }
      }}}}}}
      break;

    default :
      fprintf(stderr,"Err (fo_pseudo): \
Invalid MAXORDER encountered!\n");
      exit(EXIT_FAILURE);
  }

  return;
} /* end fo_pseudo */

/* This function will record the likelihood of all   *
 * possible oligomers for fixed-order Markov models. */
void fo_rfreqs(countarrayT *tal,int model,foprobT *fops)
{
  int i,j,k,l,m,n,z, /* Iterator variables */
      denom;         /* Stores total num   *
                      * oligos for order   */
  #ifdef PROBREPORT
  double pmf_test; /* Used to ascertain that a valid PMF was formed */
  #endif

  /* The switch statement facilitates using the  *
   * library for MAXORDER in the interval [0..5] */
  switch(MAXORDER) {
    case(0) :
      for(z=denom=0;z<ALFSIZE;++z)
        denom += tal->count1[model][z];
      for(n=0;n<ALFSIZE;++n) {
        if(denom == 0) {
          /* by fo_pseudo, this shouldn't happen, under the proviso *
           * that EMPTYCOUNT and FOPSEUDOCT are set to suitable     *
           * values.                                                */
          fops->FOprob[model][n][0][0][0][0][0]=DEFRFREQ;
          fprintf(stderr,"Warning (fo_rfreqs): \
fo_pseudo didn't do its job *correctly*!\n");
        }
        else
          fops->FOprob[model][n][0][0][0][0][0]=
            ((double)tal->count1[model][n])/denom;
      }
      break;

    case(1) :
      for(i=0;i<ALFSIZE;++i) {
        for(z=denom=0;z<ALFSIZE;++z)
          denom += tal->count2[model][i][z];
        for(n=0;n<ALFSIZE;++n) {
          if(denom == 0) {
            /* by fo_pseudo, this shouldn't happen, under the proviso *
             * that EMPTYCOUNT and FOPSEUDOCT are set to suitable     *
             * values.                                                */
            fops->FOprob[model][i][n][0][0][0][0]=DEFRFREQ;
            fprintf(stderr,"Warning (fo_rfreqs): \
fo_pseudo didn't do its job *correctly*!\n");
          }
          else
            fops->FOprob[model][i][n][0][0][0][0]=
              ((double)tal->count2[model][i][n])/denom;
        }
      }
      break;

    case(2) :
      for(i=0;i<ALFSIZE;++i) {
      for(j=0;j<ALFSIZE;++j) {
        for(z=denom=0;z<ALFSIZE;++z)
          denom += tal->count3[model][i][j][z];
        for(n=0;n<ALFSIZE;++n) {
          if(denom == 0) {
            /* by fo_pseudo, this shouldn't happen, under the proviso *
             * that EMPTYCOUNT and FOPSEUDOCT are set to suitable     *
             * values.                                                */
            fops->FOprob[model][i][j][n][0][0][0]=DEFRFREQ;
            fprintf(stderr,"Warning (fo_rfreqs): \
fo_pseudo didn't do its job *correctly*!\n");
          }
          else
            fops->FOprob[model][i][j][n][0][0][0]=
              ((double)tal->count3[model][i][j][n])/denom;
        }
      }}
      break;

    case(3) :
      for(i=0;i<ALFSIZE;++i) {
      for(j=0;j<ALFSIZE;++j) {
      for(k=0;k<ALFSIZE;++k) {
        for(z=denom=0;z<ALFSIZE;++z)
          denom += tal->count4[model][i][j][k][z];
        for(n=0;n<ALFSIZE;++n) {
          if(denom == 0) {
            /* by fo_pseudo, this shouldn't happen, under the proviso *
             * that EMPTYCOUNT and FOPSEUDOCT are set to suitable     *
             * values.                                                */
            fops->FOprob[model][i][j][k][n][0][0]=DEFRFREQ;
            fprintf(stderr,"Warning (fo_rfreqs): \
fo_pseudo didn't do its job *correctly*!\n");
          }
          else
            fops->FOprob[model][i][j][k][n][0][0]=
              ((double)tal->count4[model][i][j][k][n])/denom;
        }
      }}}
      break;

    case(4) :
      for(i=0;i<ALFSIZE;++i) {
      for(j=0;j<ALFSIZE;++j) {
      for(k=0;k<ALFSIZE;++k) {
      for(l=0;l<ALFSIZE;++l) {
        for(z=denom=0;z<ALFSIZE;++z)
          denom += tal->count5[model][i][j][k][l][z];
        for(n=0;n<ALFSIZE;++n) {
          if(denom == 0) {
            /* by fo_pseudo, this shouldn't happen, under the proviso *
             * that EMPTYCOUNT and FOPSEUDOCT are set to suitable     *
             * values.                                                */
            fops->FOprob[model][i][j][k][l][n][0]=DEFRFREQ;
            fprintf(stderr,"Warning (fo_rfreqs): \
fo_pseudo didn't do its job *correctly*!\n");
          }
          else
            fops->FOprob[model][i][j][k][l][n][0]=
              ((double)tal->count5[model][i][j][k][l][n])/denom;
        }
      }}}}
      break;

    case(5) :
      for(i=0;i<ALFSIZE;++i) {
      for(j=0;j<ALFSIZE;++j) {
      for(k=0;k<ALFSIZE;++k) {
      for(l=0;l<ALFSIZE;++l) {
      for(m=0;m<ALFSIZE;++m) {
        for(z=denom=0;z<ALFSIZE;++z)
          denom += tal->count6[model][i][j][k][l][m][z];
        #ifdef PROBREPORT
        pmf_test=0.0;
        #endif
        for(n=0;n<ALFSIZE;++n) {
          if(denom == 0) {
            /* by fo_pseudo, this shouldn't happen, under the proviso *
             * that EMPTYCOUNT and FOPSEUDOCT are set to suitable     *
             * values.                                                */
            fops->FOprob[model][i][j][k][l][m][n]=DEFRFREQ;
            fprintf(stderr,
              "Warning (fo_rfreqs): \
fo_pseudo didn't do its job *correctly*!\n");
          }
          else
            fops->FOprob[model][i][j][k][l][m][n]=
              ((double)tal->count6[model][i][j][k][l][m][n])/denom;
          #ifdef PROBREPORT
          pmf_test+=fops->FOprob[model][i][j][k][l][m][n];
          fprintf(stdout,"%u=%.4f,",n,fops->FOprob[model][i][j][k][l][m][n]);
          #endif
        }
        #ifdef PROBREPORT
        fprintf(stdout," Sum of prob's for \
history %i%i%i%i%i (model %i) is %.2f ",
          i,j,k,l,m,model,pmf_test);
        if( fabs(pmf_test - CERTPROB) < MAXDIFF )
          fprintf(stdout,"valid\n");
        else
          fprintf(stdout,"problematic!\n");
        #endif
      }}}}}
      break;

    default:
      fprintf(stderr,"Err (fo_rfreqs): Invalid MAXORDER encountered!\n");
      exit(EXIT_FAILURE);
  }

  return;
} /* end fo_rfreqs */

#ifdef QUANTSPEC
/* Function to return log-likelihood of a test input sequence, *
 * which must be null-terminated, under the model specified    *
 * using the model argument. It returns the probability of the *
 * input sequence given its first base occurs in the phase     *
 * specified by model, which parameters selected on a quantile *
 * specific basis. If testwindowgc is TRUE, this function will *
 * modulate between quantile-specific parameters as a function *
 * of window-based G+C composition.                            */
long double FOprobQT(char *input,foprobT *probs,qtboundsT *bounds,
  int model,int testwindowgc)
{
  int length,          /* Length of input string         */
      *nt=NULL,        /* int translation of input       */
      i,j,             /* Iterator variable              */
      GCcomp=0;        /* G+C content of input           */
  long double logprob; /* end result to return           */
  foprobT *probsptr;   /* Which element of probs to use? */

  /* Verify that model is given a reasonable value */
  if (model < 0 || model > NONCODING) {
    fprintf(stderr,"Err (FOprobQT): Invalid model value (%i)\n",model);
    exit(EXIT_FAILURE);
  }

  /* Verify that we really have a sequence to test. If input is set to NULL, *
   * this conditional will short-circuit and strlen won't get called on it.  */
  if((input==NULL) || ((length=(int)strlen(input))==0)) {
    if (input==NULL)
      fprintf(stderr,"Err (FOprobQT): No real sequence passed!\n");
    else
      fprintf(stderr,"Err (FOprobQT): Passed a zero-length test sequence!\n");
    exit(EXIT_FAILURE);
  }
  else if(length < (MAXORDER+1)) {
    fprintf(stderr,"Err (FOprobQT): Can't assess sequences < %i bases!\n",
      MAXORDER+1);
    exit(EXIT_FAILURE);
  }
  else if((nt=(int*)malloc(sizeof(int)*length))==NULL) {
    fprintf(stderr,"Err (FOprobQT): Out of memory!\n");
    exit(EXIT_FAILURE);
  }
  else /* Copy input translation into nt */
    for(i=0;i<length;++i)
      nt[i]=trans(input[i]);

  /* If sequence is shorter than a given window size, modulating *
   * parameters based on windog G+C composition is not possible. */
  if(length < (2*FLANK+MAXORDER+1))
    testwindowgc=FALSE;

  if (testwindowgc == FALSE) {
    /* Calc overall [GC] and cast to an integer *
     * precision to 100th's only.               */
    GCcomp=0;
    for(i=0;i<length;++i)
      if(input[i] == 'c' || input[i] == 'C' ||
         input[i] == 'g' || input[i] == 'G')
        ++GCcomp;
    GCcomp=(int)floor(100*(GCcomp/(float)length));
  }
  else {
    /* Calculate G+C composition of initial window. *
     * Note that, if using a window based approach  *
     * to transition probability estimate selection *
     * we use an integer to keep track of G+C comp. *
     * to minimize floating point operations. This  *
     * presumes that 2*FLANK+MAXORDER+1 == 100!     */
    GCcomp=0;
    for(i=0;i<(2*FLANK+MAXORDER+1);++i)
      if(input[i] == 'c' || input[i] == 'C' ||
         input[i] == 'g' || input[i] == 'G')
        ++GCcomp;
  }

  if(model == NONCODING) { /* homogeneous Markov model */
    if (testwindowgc == FALSE) {
      /* Select the best set of quantile-specific parameters */
      for(i=1;i<QTCT;++i)
        if(GCcomp < bounds->noncodbound[i])
          break;
      probsptr = &probs[i-1];

      /* compute likelihood */
      for(i=MAXORDER,logprob=0.0;i<length;++i)
        logprob += log(probsptr->FOprob[model]
          [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
    }
    else { /* window-based parameter modulation */
      /* Select the best set of quantile-specific   *
       * parameters for the first window. Note that *
       * GCcompint is already primed on the first   *
       * window at this point.                      */
      for(i=1;i<QTCT;++i)
        if(GCcomp < bounds->noncodbound[i])
          break;
      probsptr = &probs[i-1];

      for(i=MAXORDER,logprob=0.0;i<FLANK+MAXORDER+1;++i)
        logprob += log(probsptr->FOprob[model]
          [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);

      /* Modulate parameters until last window.  Here, i indexes *
       * the 3'-most position in the sequence no longer in the   *
       * current window.                                         */
      for(i=0;i+(2*FLANK+MAXORDER+1)<length;++i) {
        if(input[i] == 'c' || input[i] == 'C' ||
           input[i] == 'g' || input[i] == 'G')
          --GCcomp;
        if(input[i+(2*FLANK+MAXORDER+1)] == 'c' ||
           input[i+(2*FLANK+MAXORDER+1)] == 'C' ||
           input[i+(2*FLANK+MAXORDER+1)] == 'g' ||
           input[i+(2*FLANK+MAXORDER+1)] == 'G')
          ++GCcomp;

        /* Now, we're handling the i+1 through *
         * i+(2*FLANK+MAXORDER+1) window.      */
        for(j=1;j<QTCT;++j)
          if(GCcomp < bounds->noncodbound[j])
            break;
        probsptr = &probs[j-1];

        /* Assay the hexamer centered in the window above */
        logprob += log(probsptr->FOprob[model]
          [nt[i+FLANK+MAXORDER-4]]
          [nt[i+FLANK+MAXORDER-3]]
          [nt[i+FLANK+MAXORDER-2]]
          [nt[i+FLANK+MAXORDER-1]]
          [nt[i+FLANK+MAXORDER]]
          [nt[i+FLANK+MAXORDER+1]]);
      }
      /* assay remaining 3' hexamers */
      for(i=length-FLANK;i<length;++i)
        logprob += log(probsptr->FOprob[model]
          [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
    }
  }
  else { /* inhomogeneous Markov model */
    if (testwindowgc == FALSE) {
      /* Select the best set of quantile-specific parameters */
      for(i=1;i<QTCT;++i)
        if(GCcomp < bounds->codingbound[i])
          break;
      probsptr = &probs[i-1];

      /* compute likelihood. Note that MAXORDER == 5 */
      for(i=MAXORDER,logprob=0.0;i<length;++i) {
        if (model < 3) /* Forward strand */
          logprob += log(probsptr->FOprob[(i-MAXORDER+model)%3]
            [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
        else /* Reverse strand */
          logprob += log(probsptr->FOprob[((i-MAXORDER+model)%3)+3]
            [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
      }
    }
    else { /* window-based parameter modulation */
      /* Select the best set of quantile-specific   *
       * parameters for the first window. Note that *
       * GCcompint is already primed on the first   *
       * window at this point.                      */
      for(i=1;i<QTCT;++i)
        if(GCcomp < bounds->codingbound[i])
          break;
      probsptr = &probs[i-1];

      for(i=MAXORDER,logprob=0.0;i<FLANK+MAXORDER+1;++i) {
        if (model < 3) /* Forward strand */
          logprob += log(probsptr->FOprob[(i-MAXORDER+model)%3]
            [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
        else /* Reverse strand */
          logprob += log(probsptr->FOprob[((i-MAXORDER+model)%3)+3]
            [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
      }

      /* modulate parameters until last window */
      for(i=0;i+(2*FLANK+MAXORDER+1)<length;++i) {
        if(input[i] == 'c' || input[i] == 'C' ||
           input[i] == 'g' || input[i] == 'G')
          --GCcomp;
        if(input[i+(2*FLANK+MAXORDER+1)] == 'c' ||
           input[i+(2*FLANK+MAXORDER+1)] == 'C' ||
           input[i+(2*FLANK+MAXORDER+1)] == 'g' ||
           input[i+(2*FLANK+MAXORDER+1)] == 'G')
          ++GCcomp;

        for(j=1;j<QTCT;++j)
          if(GCcomp < bounds->codingbound[j])
            break;
        probsptr = &probs[j-1];

        if (model < 3) /* Forward strand */
          logprob += log(probsptr->FOprob[(i+FLANK+1+model)%3]
            [nt[i+FLANK+MAXORDER-4]]
            [nt[i+FLANK+MAXORDER-3]]
            [nt[i+FLANK+MAXORDER-2]]
            [nt[i+FLANK+MAXORDER-1]]
            [nt[i+FLANK+MAXORDER]]
            [nt[i+FLANK+MAXORDER+1]]);
        else /* Reverse strand */
          logprob += log(probsptr->FOprob[((i+FLANK+1+model)%3)+3]
            [nt[i+FLANK+MAXORDER-4]]
            [nt[i+FLANK+MAXORDER-3]]
            [nt[i+FLANK+MAXORDER-2]]
            [nt[i+FLANK+MAXORDER-1]]
            [nt[i+FLANK+MAXORDER]]
            [nt[i+FLANK+MAXORDER+1]]);
      }
      /* assay remaining 3' hexamers */
      for(i=length-FLANK;i<length;++i) {
        if (model < 3) /* Forward strand */
          logprob += log(probsptr->FOprob[(i-MAXORDER+model)%3]
            [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
        else /* Reverse strand */
          logprob += log(probsptr->FOprob[((i-MAXORDER+model)%3)+3]
            [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
      }
    }
  }

  /* back to the heap */
  free(nt);

  return(logprob);
} /* end FOprobQT */
#else
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
  int threeper,int shadow)
{
  int length,          /* Length of input string   */
      *nt=NULL,        /* int translation of input */
      i;               /* Iterator variable        */
  long double logprob; /* end result to return     */

  /* Verify that model is given a reasonable value */
  if (shadow == FALSE) {
    if (threeper == FALSE) {
      if (model < 0 || model > NONCODING) {
        fprintf(stderr,"Err (FOprob): Invalid model value (%i)\n",model);
        exit(EXIT_FAILURE);
      }
    }
    else {
      if ((model < 0 || model > 2) && model != NONCODING) {
        fprintf(stderr,"Err (FOprob): Invalid model value (%i)\n",model);
        exit(EXIT_FAILURE);
      }
    }
  }
  else {
    if (model < 3 || model > 5) {
      fprintf(stderr,"Err (FOprob): Invalid model value (%i)\n",model);
      exit(EXIT_FAILURE);
    }
  }

  /* Verify that we really have a sequence to test. If input is set to NULL, *
   * this conditional will short-circuit and strlen won't get called on it.  */
  if((input==NULL) || ((length=(int)strlen(input))==0)) {
    if (input==NULL)
      fprintf(stderr,"Err (FOprob): No real sequence passed!\n");
    else
      fprintf(stderr,"Err (FOprob): Passed a zero-length test sequence!\n");
    exit(EXIT_FAILURE);
  }
  else if(length < (MAXORDER+1)) {
    fprintf(stderr,"Err (FOprob): Can't assess sequences < %i bases!\n",
      MAXORDER+1);
    exit(EXIT_FAILURE);
  }
  else if((nt=(int*)malloc(sizeof(int)*length))==NULL) {
    fprintf(stderr,"Err (FOprob): Out of memory!\n");
    exit(EXIT_FAILURE);
  }
  else /* Copy input translation into nt */
    for(i=0;i<length;++i)
      nt[i]=trans(input[i]);

  /* compute log-likelihood of the input string */
  if (threeper != FALSE) { /* Compute 3-periodic result */
    for(i=MAXORDER,logprob=0.0;i<length;++i) {
      switch(MAXORDER) {
        case(0) :
          if (shadow == FALSE) /* Forward strand */
            logprob += log(probs->FOprob[(i-MAXORDER+model)%3]
              [nt[i]][0][0][0][0][0]);
          else /* Reverse strand */
            logprob += log(probs->FOprob[((i-MAXORDER+model)%3)+3]
              [nt[i]][0][0][0][0][0]);
          break;

        case(1) :
          if (shadow == FALSE) /* Forward strand */
            logprob += log(probs->FOprob[(i-MAXORDER+model)%3]
              [nt[i-1]][nt[i]][0][0][0][0]);
          else /* Reverse strand */
            logprob += log(probs->FOprob[((i-MAXORDER+model)%3)+3]
              [nt[i-1]][nt[i]][0][0][0][0]);
          break;

        case(2) :
          if (shadow == FALSE) /* Forward strand */
            logprob += log(probs->FOprob[(i-MAXORDER+model)%3]
              [nt[i-2]][nt[i-1]][nt[i]][0][0][0]);
          else /* Reverse strand */
            logprob += log(probs->FOprob[((i-MAXORDER+model)%3)+3]
              [nt[i-2]][nt[i-1]][nt[i]][0][0][0]);
          break;

        case(3) :
          if (shadow == FALSE) /* Forward strand */
            logprob += log(probs->FOprob[(i-MAXORDER+model)%3]
              [nt[i-3]][nt[i-2]][nt[i-1]][nt[i]][0][0]);
          else /* Reverse strand */
            logprob += log(probs->FOprob[((i-MAXORDER+model)%3)+3]
              [nt[i-3]][nt[i-2]][nt[i-1]][nt[i]][0][0]);
          break;

        case(4) :
          if (shadow == FALSE) /* Forward strand */
            logprob += log(probs->FOprob[(i-MAXORDER+model)%3]
              [nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]][0]);
          else /* Reverse strand */
            logprob += log(probs->FOprob[((i-MAXORDER+model)%3)+3]
              [nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]][0]);
          break;

        case(5) :
          if (shadow == FALSE) /* Forward strand */
            logprob += log(probs->FOprob[(i-MAXORDER+model)%3]
              [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
          else /* Reverse strand */
            logprob += log(probs->FOprob[((i-MAXORDER+model)%3)+3]
              [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
          break;

        default:
          fprintf(stderr,"Err (FOprob): Invalid MAXORDER encountered!\n");
          exit(EXIT_FAILURE);
      }
    }
  }
  else { /* Compute non-periodic result */
    for(i=MAXORDER,logprob=0.0;i<length;++i) {
      switch(MAXORDER) {
        case(0) :
          logprob += log(probs->FOprob[model][nt[i]][0][0][0][0][0]);
          break;
        case(1) :
          logprob += log(probs->FOprob[model][nt[i-1]]
                                             [nt[i]][0][0][0][0]);
          break;
        case(2) :
          logprob += log(probs->FOprob[model][nt[i-2]]
                                             [nt[i-1]]
                                             [nt[i]][0][0][0]);
          break;
        case(3) :
          logprob += log(probs->FOprob[model][nt[i-3]]
                                             [nt[i-2]]
                                             [nt[i-1]]
                                             [nt[i]][0][0]);
          break;
        case(4) :
          logprob += log(probs->FOprob[model][nt[i-4]]
                                             [nt[i-3]]
                                             [nt[i-2]]
                                             [nt[i-1]]
                                             [nt[i]][0]);
          break;
        case(5) :
          logprob += log(probs->FOprob[model][nt[i-5]]
                                             [nt[i-4]]
                                             [nt[i-3]]
                                             [nt[i-2]]
                                             [nt[i-1]]
                                             [nt[i]]);
          break;
        default:
          fprintf(stderr,"Err (FOprob): Invalid MAXORDER encountered!\n");
          exit(EXIT_FAILURE);
      }
    }
  }

  /* back to the heap */
  free(nt);

  return(logprob);
} /* end FOprob */
#endif

#ifdef QUANTSPEC
/* Function to import a trained foprobT object array and a qtboundsT *
 * object (stored in a binary file) into main memory.                */
void import_fo_probsQT(char *filename,foprobT *probs,qtboundsT *bounds)
{
  FILE *fptr=NULL; /* For connecting to binary input file */

  /* Open parameter file */
  if((fptr=fopen(filename,"rb"))==NULL) {
    fprintf(stderr,"Err (import_fo_probsQT): \
Can't open %s for binary reading!\n",filename);
    exit(EXIT_FAILURE);
  }
  /* Read in the foprobT "objects" */
  else if(fread(probs,sizeof(foprobT),QTCT,fptr) != QTCT) {
    fprintf(stderr,"Err (import_fo_probsQT): \
Error reading from %s!\n",filename);
    exit(EXIT_FAILURE);
  }
  /* Read in the qtboundsT "object" */
  else if(fread(bounds,sizeof(qtboundsT),1,fptr) != 1) {
    fprintf(stderr,"Err (import_fo_probsQT): \
Error reading from %s!\n",filename);
    exit(EXIT_FAILURE);
  }
  else
    (void)fclose(fptr);

  return;
} /* end import_fo_probsQT */
#else
/* Function to import a trained foprobT     *
 * object (a binary file) into main memory. */
void import_fo_probs(char *filename,foprobT *probs)
{
  FILE *fptr=NULL; /* For connecting to binary input file */

  /* Open parameter file */
  if((fptr=fopen(filename,"rb"))==NULL) {
    fprintf(stderr,"Err (import_fo_probs): \
Can't open %s for binary reading!\n",filename);
    exit(EXIT_FAILURE);
  }
  /* Read in the foprobT "object" */
  else if(fread(probs,sizeof(foprobT),1,fptr) != 1) {
    fprintf(stderr,"Err (import_fo_probs): \
Error reading from %s!\n",filename);
    exit(EXIT_FAILURE);
  }
  else
    (void)fclose(fptr);

  return;
} /* end import_fo_probs */
#endif
