/* dm_utils.c
 * Michael E Sparks (mespar1@gmail.com)
 *
 * This file contains code for developing/using dynamically modulating
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
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dm_utils.h"
#include "fo_utils.h"
#include "immpractical.h"
#include "sequence_parse.h"
#include "simplex.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define LPm    10          /* number of constraints in LP problem        */
#define LPn    4           /* number of nonbasic variables in LP problem */
#define LPsize (LPm+LPn+1) /* Basic LP data structure size unit          */

/* Function definitions ******************************************************/

/* This function oversees the building of a trained dmprobT *
 * file based on dynamically modulating Markov models.      */
int dm_probmaker(int argc,char *argv[])
{
  FILE *infile=NULL,          /* For processing Fasta files           */
       *outfile=NULL;         /* Will store binary dm_probT object    */
  char *sequence=NULL,        /* For storing current sequence         */
    outputname[MAXFILENAME];  /* Stores final name of output file     */
  double mean,                /* mean of dimension's probabilities    */
         stddev,              /* stddev of dimension's probabilities  */
         min,                 /* min of dimension's probabilities     */
         max;                 /* max of dimension's probabilities     */
  int seqlength,              /* Records length of current sequence   */
      model,                  /* Model iterator variable              */
      currentfile,            /* Iterator based on input file at hand */
      partition_counter,      /* Iterator for each of the data splits */
      num_parts,              /* Store number of data splits          */
      i,j,k,l,m,n;            /* Iterate through dimensions           */
  #ifdef QUANTSPEC
  int quantind;               /* iterator variable                    */
  countarrayT *tallies=NULL;  /* Record oligomer frequencies          */
  foprobT **fo_parmarray=NULL;/* Stores final, pseudocount smoothed,  *
                               * fixed-order probabilities for each   *
                               * of the training data splits.         */
  dmprobT *dm_parameters=NULL;/* Final result to write to outfile     */
  qtboundsT *qt_bounds=NULL,  /* Store quantiles per partition.       */
            qt_boundsFINAL;   /* Store quantile boundaries *averaged* *
                               * over each partition. Here, we assume *
  * that each of the num_parts partitions are independent, random     *
  * samples obtained from a source distribution.                      */
  #else
  countarrayT tallies;        /* Record oligomer frequencies          */
  foprobT *fo_parmarray=NULL; /* Stores final, pseudocount smoothed,  *
                               * fixed-order probabilities for each   *
                               * of the training data splits.         */
  dmprobT dm_parameters;      /* Final result to write to outfile     */
  #endif

  /* First, we train pseudocount-smoothed, fixed-order probability *
   * models for each of the training data partitions, a quantity   *
   * specified by argv[2]. Note that it is assumed the maintrain   *
   * function checked the validity of the command line arguments.  */
  num_parts=atoi(argv[2]);
  if(num_parts<1){
    fprintf(stderr,"Err (dm_probmaker): Invalid partition count!\n");
    exit(EXIT_FAILURE);
  }

  #ifdef QUANTSPEC
  /* Allocate storage for array of countarrayT elements */
  if((tallies=(countarrayT*)malloc(sizeof(countarrayT)*(int)QTCT))==NULL) {
    fprintf(stderr,"Err (dm_probmaker): Out of memory!\n");
    exit(EXIT_FAILURE);
  }
  /* Allocate storage for arrays of fo_probT elements */
  if((fo_parmarray=(foprobT**)malloc(sizeof(foprobT*)*num_parts))==NULL) {
    fprintf(stderr,"Err (dm_probmaker): Out of memory!\n");
    exit(EXIT_FAILURE);
  }
  for(i=0;i<num_parts;++i) {
    if((fo_parmarray[i]=(foprobT*)malloc(sizeof(foprobT)*QTCT))==NULL) {
      fprintf(stderr,"Err (dm_probmaker): Out of memory!\n");
      exit(EXIT_FAILURE);
    }
  }
  /* Allocate storage for array of dmprobT elements */
  if((dm_parameters=(dmprobT*)malloc(sizeof(dmprobT)*(int)QTCT))==NULL) {
    fprintf(stderr,"Err (dm_probmaker): Out of memory!\n");
    exit(EXIT_FAILURE);
  }
  /* Allocate storage for the array of qtboundsT elements */
  if((qt_bounds=(qtboundsT*)malloc(sizeof(qtboundsT)*num_parts))==NULL) {
    fprintf(stderr,"Err (dm_probmaker): Out of memory!\n");
    exit(EXIT_FAILURE);
  }
  #else
  /* Allocate storage for the array of fo_probT elements */
  if((fo_parmarray=(foprobT*)malloc(sizeof(foprobT)*num_parts))==NULL) {
    fprintf(stderr,"Err (dm_probmaker): Out of memory!\n");
    exit(EXIT_FAILURE);
  }
  #endif

  for(partition_counter=0;
      partition_counter<num_parts;
      ++partition_counter) {

    #ifdef QUANTSPEC
    calc_quantiles(argv,DMMMX,&qt_bounds[partition_counter]);
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
      /* open the necessary files */
      if(currentfile==0){
        if((infile=fopen(argv[3+partition_counter],"rt"))==NULL) {
          fprintf(stderr,"Err (dm_probmaker): Can't seem to open %s!\n",
            argv[3+partition_counter]);
          exit(EXIT_FAILURE);
        }
      }
      else{ /* (currentfile==1) */
        if((infile=fopen(argv[3+num_parts+partition_counter],"rt"))==NULL) {
          fprintf(stderr,"Err (dm_probmaker): Can't seem to open %s!\n",
            argv[3+num_parts+partition_counter]);
          exit(EXIT_FAILURE);
        }
      }

      /* Derive counts from data */
      while ((sequence=get_fasta(infile,sequence,&seqlength))!=NULL) {
        if(seqlength >= (MAXORDER+1)) {
          switch (currentfile) {
            case 0 :
              #ifdef QUANTSPEC
              /* Quantile-specific behavior demands three-periodic modeling *
               * of both strands---locked in a NUMMODELS == 7 scheme.       */
              get_counts3percodQT(sequence,tallies,
                qt_bounds[partition_counter]);
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
              get_countsQT(sequence,tallies,qt_bounds[partition_counter]);
              #else
              #ifdef THREEPERIODIC
              get_counts(sequence,&tallies,NONCODING);
              #else
              get_counts(sequence,&tallies,currentfile);
              #endif
              #endif
              break;
            default :
              fprintf(stderr,"Err (dm_probmaker): \
You seem to have passed an invalid training file?\n");
              (void)fclose(infile);
              return(FALSE);
          }
          free(sequence);
          sequence=NULL;
        }
        else {
          fprintf(stderr,"Err (dm_probmaker): First run %s on your input!\n",
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
        fo_rfreqs(&tallies[quantind],model,
          &fo_parmarray[partition_counter][quantind]);
      }
      #else
      /* Smooth empty counts with pseudocounts. *
       * Since we're only interested in order   *
       * 5, we'll only add pseudocounts for it. */
      fo_pseudo(&tallies,model,TRUE);

      /* Record the "smoothed" probability */
      fo_rfreqs(&tallies,model,fo_parmarray+partition_counter);
      #endif
    }
  } /* end partition_counter for */

  /* For each dimension, compute and store the mean, *
   * standard deviation, minimum, and maximum values */
  for(model=0;model<NUMMODELS;++model) {
    for(i=0;i<ALFSIZE;++i) {
    for(j=0;j<ALFSIZE;++j) {
    for(k=0;k<ALFSIZE;++k) {
    for(l=0;l<ALFSIZE;++l) {
    for(m=0;m<ALFSIZE;++m) {
    for(n=0;n<ALFSIZE;++n) {
      #ifdef QUANTSPEC
      for(quantind=0;quantind<QTCT;++quantind) {
        max=-0.1; min=1.1; /* Impossible values */
        mean=stddev=0.0;
        for(partition_counter=0;
            partition_counter<num_parts;
            ++partition_counter) {
          if( fabs(fo_parmarray[partition_counter][quantind].
            FOprob[model][i][j][k][l][m][n]) < DBL_EPSILON ||
              fabs(fo_parmarray[partition_counter][quantind].
            FOprob[model][i][j][k][l][m][n]-CERTPROB) < DBL_EPSILON) {
            fprintf(stderr,"Err (dm_probmaker): \
Expected input distributions to be pre-smoothed!\n");
            exit(EXIT_FAILURE);
          }
          mean += fo_parmarray[partition_counter][quantind].FOprob[model]
            [i][j][k][l][m][n];
          min=MIN(min,fo_parmarray[partition_counter][quantind].FOprob[model]
            [i][j][k][l][m][n]);
          max=MAX(max,fo_parmarray[partition_counter][quantind].FOprob[model]
            [i][j][k][l][m][n]);
        }
        mean /= num_parts; /* take arithmetic mean */

        /* Now for the (population) standard deviation */
        for(partition_counter=0;
            partition_counter<num_parts;
            ++partition_counter) {
          stddev +=
            (fo_parmarray[partition_counter][quantind].
              FOprob[model][i][j][k][l][m][n]-mean) *
            (fo_parmarray[partition_counter][quantind].
              FOprob[model][i][j][k][l][m][n]-mean);
        }
        stddev /= num_parts;
        stddev = sqrt(stddev);

        /* record these values */
        dm_parameters[quantind].mean[model][i][j][k][l][m][n]=mean;
        dm_parameters[quantind].stddev[model][i][j][k][l][m][n]=stddev;
        /* Since parameters input to this method must themselves  *
         * be pre-smoothed (nonsingular), min and max are         *
         * guaranteed to be probabilities in the interval (0,1)   */
        dm_parameters[quantind].min[model][i][j][k][l][m][n]=min;
        dm_parameters[quantind].max[model][i][j][k][l][m][n]=max;
        /* For mean +/- stddev (~68% of a std dist), there are *
         * no guarantees that these values will be legitimate  *
         * probabilities in (0,1)                              */
        dm_parameters[quantind].muless1sigma[model][i][j][k][l][m][n]=
          MAX(mean-stddev,POSCONST);
        dm_parameters[quantind].muplus1sigma[model][i][j][k][l][m][n]=
          MIN(mean+stddev,(CERTPROB-POSCONST));
        #ifdef DMPARMPRINT
        printf("For oligo %i%i%i%i%i%i , quantile %i: \
mean %.4f,   stddev %.4f,\n",i,j,k,l,m,n,quantind,
          dm_parameters[quantind].mean[model][i][j][k][l][m][n],
          dm_parameters[quantind].stddev[model][i][j][k][l][m][n]);
        printf("                        min %.4f,    max %.4f,\n",
          dm_parameters[quantind].min[model][i][j][k][l][m][n],
          dm_parameters[quantind].max[model][i][j][k][l][m][n]);
        printf("                        mu-sig %.4f, mu+sig %.4f\n",
          dm_parameters[quantind].muless1sigma[model][i][j][k][l][m][n],
          dm_parameters[quantind].muplus1sigma[model][i][j][k][l][m][n]);
        #endif
      } /* end quantind for */
      #else
      max=-0.1; min=1.1; /* Impossible values */
      mean=stddev=0.0;
      for(partition_counter=0;
          partition_counter<num_parts;
          ++partition_counter) {
        if( fabs(fo_parmarray[partition_counter].
          FOprob[model][i][j][k][l][m][n]) < DBL_EPSILON ||
            fabs(fo_parmarray[partition_counter].
          FOprob[model][i][j][k][l][m][n]-CERTPROB) < DBL_EPSILON) {
          fprintf(stderr,"Err (dm_probmaker): \
Expected input distributions to be pre-smoothed!\n");
          exit(EXIT_FAILURE);
        }
        mean += fo_parmarray[partition_counter].FOprob[model][i][j][k]
                                                             [l][m][n];
        min=MIN(min,fo_parmarray[partition_counter].FOprob[model][i][j][k]
                                                                 [l][m][n]);
        max=MAX(max,fo_parmarray[partition_counter].FOprob[model][i][j][k]
                                                                 [l][m][n]);
      }
      mean /= num_parts; /* take arithmetic mean */

      /* Now for the (population) standard deviation */
      for(partition_counter=0;
          partition_counter<num_parts;
          ++partition_counter) {
        stddev +=
          (fo_parmarray[partition_counter].
            FOprob[model][i][j][k][l][m][n]-mean) *
          (fo_parmarray[partition_counter].
            FOprob[model][i][j][k][l][m][n]-mean);
      }
      stddev /= num_parts;
      stddev = sqrt(stddev);

      /* record these values */
      dm_parameters.mean[model][i][j][k][l][m][n]=mean;
      dm_parameters.stddev[model][i][j][k][l][m][n]=stddev;
      /* Since parameters input to this method must themselves  *
       * be pre-smoothed (nonsingular), min and max are         *
       * guaranteed to be probabilities in the interval (0,1)   */
      dm_parameters.min[model][i][j][k][l][m][n]=min;
      dm_parameters.max[model][i][j][k][l][m][n]=max;
      /* For mean +/- stddev (~68% of a std dist), there are *
       * no guarantees that these values will be legitimate  *
       * probabilities in (0,1)                              */
      dm_parameters.muless1sigma[model][i][j][k][l][m][n]=
        MAX(mean-stddev,POSCONST);
      dm_parameters.muplus1sigma[model][i][j][k][l][m][n]=
        MIN(mean+stddev,(CERTPROB-POSCONST));
      #ifdef DMPARMPRINT
      printf("For oligo %i%i%i%i%i%i: mean %.4f,   stddev %.4f,\n",
        i,j,k,l,m,n,
        dm_parameters.mean[model][i][j][k][l][m][n],
        dm_parameters.stddev[model][i][j][k][l][m][n]);
      printf("                        min %.4f,    max %.4f,\n",
        dm_parameters.min[model][i][j][k][l][m][n],
        dm_parameters.max[model][i][j][k][l][m][n]);
      printf("                        mu-sig %.4f, mu+sig %.4f\n",
        dm_parameters.muless1sigma[model][i][j][k][l][m][n],
        dm_parameters.muplus1sigma[model][i][j][k][l][m][n]);
      #endif
      #endif
    }}}}}}
  } /*end second partition_counter for*/

  #ifdef QUANTSPEC
  /* For this application, we are expecting that each of the data *
   * partitions are independent, random samples obtained from a   *
   * distribution for a given functional class of sequences.  So, *
   * we expect the quantile boundaries to be similar in each of   *
   * them.  We take averages of each of these boundaries to store *
   * in the final qtboundsT object, which will be used in testing *
   * situations.  If for some reason the boundaries are out of    *
   * strictly increasing linear order, an error message will be   *
   * emitted, and the program will abort. One could alternatively *
   * implement some compromising fix.                             */
  for(quantind=1;quantind<QTCT;++quantind) {
    mean=0;
    for(partition_counter=0;
        partition_counter<num_parts;
        ++partition_counter)
      mean+=qt_bounds[partition_counter].codingbound[quantind];
    qt_boundsFINAL.codingbound[quantind] = mean / num_parts;

    mean=0;
    for(partition_counter=0;
        partition_counter<num_parts;
        ++partition_counter)
      mean+=qt_bounds[partition_counter].noncodbound[quantind];
    qt_boundsFINAL.noncodbound[quantind] = mean / num_parts;
  }

  /* Set 0th and QTCTth elements to min and max G+C compositions *
   * observed over all partitions, respectively.                 */
  max=-0.1; min=1.1; /* Impossible values */
  for(partition_counter=0;
      partition_counter<num_parts;
      ++partition_counter) {
    min=MIN(min,qt_bounds[partition_counter].codingbound[0]);
    max=MAX(max,qt_bounds[partition_counter].codingbound[QTCT]);
  }
  qt_boundsFINAL.codingbound[0] = min;
  qt_boundsFINAL.codingbound[QTCT] = max;

  max=-0.1; min=1.1; /* Impossible values */
  for(partition_counter=0;
      partition_counter<num_parts;
      ++partition_counter) {
    min=MIN(min,qt_bounds[partition_counter].noncodbound[0]);
    max=MAX(max,qt_bounds[partition_counter].noncodbound[QTCT]);
  }
  qt_boundsFINAL.noncodbound[0] = min;
  qt_boundsFINAL.noncodbound[QTCT] = max;

  /* confirm overall ordering */
  for(quantind=1;quantind<=QTCT;++quantind) {
    if(qt_boundsFINAL.codingbound[quantind-1] >
       qt_boundsFINAL.codingbound[quantind] ||
       fabs(qt_boundsFINAL.codingbound[quantind-1] -
            qt_boundsFINAL.codingbound[quantind]) < DBL_EPSILON) {
      fprintf(stderr,"Err (dm_probmaker): \
Error in coding data distributions!\n\
qt_boundsFINAL.codingbound[%i] (%f) >= qt_boundsFINAL.codingbound[%i] (%f)\n",
        quantind-1,qt_boundsFINAL.codingbound[quantind-1],
        quantind,qt_boundsFINAL.codingbound[quantind]);
      exit(EXIT_FAILURE);
    }

    if(qt_boundsFINAL.noncodbound[quantind-1] >
       qt_boundsFINAL.noncodbound[quantind] ||
       fabs(qt_boundsFINAL.noncodbound[quantind-1] -
            qt_boundsFINAL.noncodbound[quantind]) < DBL_EPSILON) {
      fprintf(stderr,"Err (dm_probmaker): \
Error in noncod data distributions!\n\
qt_boundsFINAL.noncodbound[%i] (%f) >= qt_boundsFINAL.noncodbound[%i] (%f)\n",
        quantind-1,qt_boundsFINAL.noncodbound[quantind-1],
        quantind,qt_boundsFINAL.noncodbound[quantind]);
      exit(EXIT_FAILURE);
    }
  }
  #endif

  /* Set name for output file. */
  if ((int)strlen(argv[(2*num_parts)+3]) > MAXFILENAME) {
    fprintf(stderr,
      "Err (dm_probmaker): Filename exceeds %i characters!\n",MAXFILENAME);
    return(FALSE);
  }
  else {
    (void)strcpy(outputname,argv[(2*num_parts)+3]);
    (void)strcat(outputname,".");
    (void)strcat(outputname,algotags[atoi(argv[1])]);
  }

  /* Write our newly trained model to the output file */
  if((outfile=fopen(outputname,"wb"))==NULL) {
    fprintf(stderr,"Err (dm_probmaker): Can't seem to open %s!\n",outputname);
    exit(EXIT_FAILURE);
  }
  #ifdef QUANTSPEC
  else if(fwrite(dm_parameters,sizeof(dmprobT),QTCT,outfile) != QTCT) {
    fprintf(stderr,
      "Err (dm_probmaker): Error writing model to binary output file!\n");
    exit(EXIT_FAILURE);
  }
  else if(fwrite(&qt_boundsFINAL,sizeof(qtboundsT),1,outfile) != 1) {
    fprintf(stderr,
      "Err (dm_probmaker): \
Error writing quantile boundaries to binary output file!\n");
    exit(EXIT_FAILURE);
  }
  #else
  else if(fwrite(&dm_parameters,sizeof(dmprobT),1,outfile) != 1) {
    fprintf(stderr,
      "Err (dm_probmaker): Error writing model to binary output file!\n");
    exit(EXIT_FAILURE);
  }
  #endif
  else
    (void)fclose(outfile);

  #ifdef QUANTSPEC
  free(tallies);
  for(i=0;i<num_parts;++i)
    free(fo_parmarray[i]);
  free(fo_parmarray);
  free(dm_parameters);
  free(qt_bounds);
  #else
  free(fo_parmarray);
  #endif

  return(TRUE);
} /* end dm_probmaker */

#ifdef QUANTSPEC
/* Function to return log-likelihood of a test input sequence, *
 * which must be null-terminated, under the model specified    *
 * using the model argument. It returns the probability of the *
 * input sequence given its first base occurs in the phase     *
 * specified by model, which parameters selected on a quantile *
 * specific basis. If testwindowgc is TRUE, this function will *
 * modulate between quantile-specific parameters as a function *
 * of window-based G+C composition.                            */
long double DMMMprobQT(char *input,dmprobT *probs,qtboundsT *bounds,
  int cloudtype,int model,int testwindowgc)
{
  int length, /* Length of test sequence instance  */
      /* stores counts of oligos in input sequence */
      testcts[NUMMODELS][ALFSIZE][ALFSIZE][ALFSIZE]
                        [ALFSIZE][ALFSIZE][ALFSIZE],
      *nt=NULL,         /* int translation of input               */
      N[LPsize],        /* Store nonbasic status of each variable */
      Ntmp[LPsize]={0,1,1,1,1,0,0,0,0,0,0,0,0,0,0}, /*  copy to N */
      status,           /* record status from simplex             */
      relevantp,        /* is some history pertinent to test seq? */
      localmodel,       /* helps index three-periodic cases       */
      centhexmod,       /* model of central hexamers              */
      GCcomp=0,         /* G+C composition                        */
      locseq[MAXORDER], /* Copy hexamers for local work           */
      i,j,k,l,m,n,
      x,y,z,aa,bb;      /* iterator variables                     */
  double **A=NULL,      /* Matrix of nonbasic var coefficients    */
         Atmp[LPsize][LPsize] =
           {  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
           },             /* This matrix will be invariant (A=Atmp)  */
         b[LPsize],       /* constraints concerning basic vars       */
         c[LPsize],       /* Nonbasic variable coefficients that     *
                           * occur in objective function.            */
         v,               /* Obj fct's current value (not used)      */
         *solution=NULL;  /* Stores information on optimal solution  */
  long double logprob=0.0;/* Stores overall log-probability of input */
  dmprobT *probsptr;      /* Which element of probs to use?          */
  #ifdef DMMCREPORT
  double pmfsum; /* Check whether DMMM-derived PMF is valid */
  #endif

  /* Verify that model is given a reasonable value */
  if (model < 0 || model > NONCODING) {
    fprintf(stderr,"Err (DMMMprobQT): Invalid model value (%i)\n",model);
    exit(EXIT_FAILURE);
  }

  /* Verify cloudtype is set reasonably */
  if (cloudtype != VARIANCEBASED && cloudtype != RANGEBASED) {
    fprintf(stderr,"Err (DMMMprobQT): \
Invalid cloudtype specified (%i)\n",cloudtype);
    exit(EXIT_FAILURE);
  }

  /* Verify that we really have a sequence to test. If input is set to NULL, *
   * this conditional will short-circuit and strlen won't get called on it.  */
  if((input==NULL) || ((length=(int)strlen(input))==0)) {
    if (input==NULL)
      fprintf(stderr,"Err (DMMMprobQT): No real sequence passed!\n");
    else
      fprintf(stderr,"Err (DMMMprobQT): \
Passed a zero-length test sequence!\n");
    exit(EXIT_FAILURE);
  }
  else if((nt=(int*)malloc(sizeof(int)*length))==NULL) {
    fprintf(stderr,"Err (DMMMprobQT): Out of memory!\n");
    exit(EXIT_FAILURE);
  }
  else /* Copy input translation into nt */
    for(i=0;i<length;++i)
      nt[i]=trans(input[i]);

  /* A is basically a 2d array--I have to develop it dynamically in *
   * order to keep my simplex implentation as generic as possible,  *
   * which is overall a very minor inconvenience, IMHO.             */
  if((A=(double**)malloc(sizeof(double*)*LPsize))!=NULL){
    for(i=0;i<LPsize;++i){
      if((A[i]=(double*)malloc(sizeof(double)*LPsize))==NULL){
        fprintf(stderr,"Err (DMMMprobQT): Out of memory!\n");
        exit(EXIT_FAILURE);
      }
    }
  }
  else {
    fprintf(stderr,"Err (DMMMprobQT): Out of memory!\n");
    exit(EXIT_FAILURE);
  }

  /* If sequence is shorter than a given window size, modulating *
   * parameters based on window G+C composition is not possible. */
  if(length < (2*FLANK+MAXORDER+1))
    testwindowgc=FALSE;

  if (testwindowgc == FALSE) {
    /* Calculate overall G+C composition of test sequence */
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

  if (testwindowgc == FALSE) {
    /* Select the best set of quantile-specific parameters */
    for(i=1;i<QTCT;++i) {
      if(model==NONCODING) {
        if(GCcomp < bounds->noncodbound[i])
          break;
      }
      else
        if(GCcomp < bounds->codingbound[i])
          break;
    }
    probsptr = &probs[i-1];

    /* The idea is to return the probability of the sequence given    *
     * the first base occurs in the frame denoted by model. z will    *
     * denote the frame the input is read from (control will exit     *
     * the loop after one iteration if model == NONCODING), and       *
     * localmodel the set of transition probability estimates to use. */
    localmodel=model;
    for(z=0;z<3;++z){
      /* Initialize count (parameter) array */
      for(i=0;i<ALFSIZE;++i)
      for(j=0;j<ALFSIZE;++j)
      for(k=0;k<ALFSIZE;++k)
      for(l=0;l<ALFSIZE;++l)
      for(m=0;m<ALFSIZE;++m)
      for(n=0;n<ALFSIZE;++n)
        testcts[0][i][j][k][l][m][n]=0;

      /* Set parameters for this test instance */
      if(model==NONCODING)
        for(i=0;i<length-5;++i)
          ++testcts[0][nt[i]][nt[i+1]][nt[i+2]][nt[i+3]][nt[i+4]][nt[i+5]];
      else
        for(i=z;i<length-5;i+=3)
          ++testcts[0][nt[i]][nt[i+1]][nt[i+2]][nt[i+3]][nt[i+4]][nt[i+5]];

      /* Find those oligomers that are relevant to this test instance */
      for(i=0;i<ALFSIZE;++i) {
      for(j=0;j<ALFSIZE;++j) {
      for(k=0;k<ALFSIZE;++k) {
      for(l=0;l<ALFSIZE;++l) {
      for(m=0;m<ALFSIZE;++m) {
        relevantp=0;
        for(n=0;n<ALFSIZE;++n) {
          if(testcts[0][i][j][k][l][m][n]>0) {
            relevantp=1;
            break;
          }
        }

        /* Executes if we've found a relevant oligomer */
        if(relevantp){
          /* Set everything up for simplex.  Note that variables x_1, *
           * x_2, x_3, and x_4 correspond to bases A,C,G,T (0,1,2,3), *
           * respectively.                                            */
          v=0;
          for(aa=0;aa<LPsize;++aa) {
            N[aa]=Ntmp[aa];
            for(bb=0;bb<LPsize;++bb)
              A[aa][bb]=Atmp[aa][bb];
            /* pad b and c with zeros prior to their "real" inits */
            b[aa]=0;
            c[aa]=0;
          }
          for(aa=1;aa<=4;++aa)
            c[aa]=testcts[0][i][j][k][l][m][aa-1];
          for(aa=5;aa<=11;aa+=2) {
            if(cloudtype == RANGEBASED) {
              b[aa]=
                probsptr->max[localmodel][i][j][k][l][m][(aa-5)/2];
              b[aa+1]=
                -probsptr->min[localmodel][i][j][k][l][m][(aa-5)/2];
            }
            else { /* (cloudtype == VARIANCEBASED) */
              b[aa]=
                probsptr->muplus1sigma[localmodel][i][j][k][l][m][(aa-5)/2];
              b[aa+1]=
                -probsptr->muless1sigma[localmodel][i][j][k][l][m][(aa-5)/2];
            }
          }
          b[13]=1;
          b[14]=-1;

          /* Solve using my simplex implementation */
          solution=NULL;
          solution=simplex(LPsize,1,N,A,b,c,&v,LPm,LPn,&status,solution);
          if(solution==NULL){
            fprintf(stderr,"Err (DMMMprobQT): \
Simplex unexpectedly failed to compute a solution!\n");
            exit(EXIT_FAILURE);
          }

          /* Increment logprob appropriately */
          #ifdef DMMCREPORT
          pmfsum=0.0;
          #endif
          for(aa=0;aa<ALFSIZE;++aa) {
            if(testcts[0][i][j][k][l][m][aa]>0)
              logprob+=testcts[0][i][j][k][l][m][aa]*log(solution[aa+1]);
            #ifdef DMMCREPORT
            printf("%i=%.4f,",aa,solution[aa+1]);
            pmfsum+=solution[aa+1];
            #endif
          }
          #ifdef DMMCREPORT
          printf("Sum of prob's for history %i%i%i%i%i (model %i) is %.2f ",
            i,j,k,l,m,model,pmfsum);
          if(fabs(pmfsum-CERTPROB) > MAXDIFF)
            printf("problematic!\n");
          else
            printf("valid\n");
          #endif

          /* return to heap */
          free(solution);
        } /* end if relevantp */
      }}}}}

      /* update localmodel */
      if(model==NONCODING)
        break;
      else if(model < 3)
        localmodel=((localmodel+1)%3);
      else /* shadow strand */
        localmodel=((localmodel+1)%3)+3;
    } /* end localmodel for */
  } /* end if */
  else { /* window-based parameter modulation */
    /* In this case, parameters for simplex are based on the *
     * oligomers encountered on a window-specific basis.     */

    /* Select the best set of quantile-specific   *
     * parameters for the first window. Note that *
     * GCcomp is already primed on the first      *
     * window at this point.                      */
    for(i=1;i<QTCT;++i) {
      if(model==NONCODING) {
        if(GCcomp < bounds->noncodbound[i])
          break;
      }
      else
        if(GCcomp < bounds->codingbound[i])
          break;
    }
    probsptr = &probs[i-1];

    /* One init is all that's needed in this section */
    for(x=0;x<NUMMODELS;++x)
      for(i=0;i<ALFSIZE;++i)
      for(j=0;j<ALFSIZE;++j)
      for(k=0;k<ALFSIZE;++k)
      for(l=0;l<ALFSIZE;++l)
      for(m=0;m<ALFSIZE;++m)
      for(n=0;n<ALFSIZE;++n)
        testcts[x][i][j][k][l][m][n]=0;

    /* Compute likelihood of 5' end of sequence */
    localmodel=model; /* gets updated at end of next for-loop */
    for(z=0;z<3;++z){
      /* Get oligomer counts from across this *
       * window, in the appropriate frame.    */
      if(model==NONCODING)
        for(i=0;i<(2*FLANK+1);++i)
          ++testcts[localmodel]
              [nt[i]][nt[i+1]][nt[i+2]][nt[i+3]][nt[i+4]][nt[i+5]];
      else
        for(i=z;i<(2*FLANK+1);i+=3)
          ++testcts[localmodel]
              [nt[i]][nt[i+1]][nt[i+2]][nt[i+3]][nt[i+4]][nt[i+5]];

      /* This loop will individually address relevant *
       * hexamers in the 5' end of the first window,  *
       * up through the first central hexamer.        */
      if(model==NONCODING)
        i=0;
      else
        i=z;
      while(i<=FLANK) { /* i updated below. Note that *
         * i is updated in such a way as to prevent   *
         * modulating localmodel in this loop.        */

        /* copy working hexamer */
        for(x=0;x<MAXORDER+1;++x)
          locseq[x]=nt[i+x];

        /* Set everything up for simplex.  Note that variables x_1, *
         * x_2, x_3, and x_4 correspond to bases A,C,G,T (0,1,2,3), *
         * respectively.                                            */
        v=0;
        for(aa=0;aa<LPsize;++aa) {
          N[aa]=Ntmp[aa];
          for(bb=0;bb<LPsize;++bb)
            A[aa][bb]=Atmp[aa][bb];
          /* pad b and c with zeros prior to their "real" inits */
          b[aa]=0;
          c[aa]=0;
        }
        for(aa=1;aa<=4;++aa)
          c[aa]=testcts[localmodel][locseq[0]]
                                   [locseq[1]]
                                   [locseq[2]]
                                   [locseq[3]]
                                   [locseq[4]][aa-1];
        for(aa=5;aa<=11;aa+=2) {
          if(cloudtype == RANGEBASED) {
            b[aa]=probsptr->max[localmodel][locseq[0]]
                                           [locseq[1]]
                                           [locseq[2]]
                                           [locseq[3]]
                                           [locseq[4]][(aa-5)/2];
            b[aa+1]=-probsptr->min[localmodel][locseq[0]]
                                              [locseq[1]]
                                              [locseq[2]]
                                              [locseq[3]]
                                              [locseq[4]][(aa-5)/2];
          }
          else { /* (cloudtype == VARIANCEBASED) */
            b[aa]=probsptr->muplus1sigma[localmodel][locseq[0]]
                                                    [locseq[1]]
                                                    [locseq[2]]
                                                    [locseq[3]]
                                                    [locseq[4]][(aa-5)/2];
            b[aa+1]=-probsptr->muless1sigma[localmodel][locseq[0]]
                                                       [locseq[1]]
                                                       [locseq[2]]
                                                       [locseq[3]]
                                                       [locseq[4]][(aa-5)/2];
          }
        }
        b[13]=1;
        b[14]=-1;

        /* Solve using my simplex implementation */
        solution=NULL;
        solution=simplex(LPsize,1,N,A,b,c,&v,LPm,LPn,&status,solution);
        if(solution==NULL){
          fprintf(stderr,"Err (DMMMprobQT): \
Simplex unexpectedly failed to compute a solution!\n");
          exit(EXIT_FAILURE);
        }

        /* Increment logprob appropriately */
        logprob+=log(solution[locseq[5]+1]);

        /* free memory */
        free(solution);

        /* Update i index */
        if(model==NONCODING)
          ++i;
        else
          i+=3; /* Hence, no need to update localmodel! */
      } /* end while */

      /* Now, update localmodel */
      if(model==NONCODING)
        break;
      else if(model < 3)
        localmodel=((localmodel+1)%3);
      else /* shadow strand */
        localmodel=((localmodel+1)%3)+3;
    } /* end localmodel (z) for */

    /* modulate parameters until last window *
     * Note that this loop begins processing *
     * hexamers central to the window that   *
     * begins at position 1, NOT 0! In other *
     * words, y indexes that position in the *
     * sequence just prior to the current    *
     * window.                               */
    localmodel=model; /* easiest just to reset this */
    for(y=0;y+(2*FLANK+MAXORDER+1)<length;++y) {
      if(input[y] == 'c' || input[y] == 'C' ||
         input[y] == 'g' || input[y] == 'G')
        --GCcomp;
      if(input[y+(2*FLANK+MAXORDER+1)] == 'c' ||
         input[y+(2*FLANK+MAXORDER+1)] == 'C' ||
         input[y+(2*FLANK+MAXORDER+1)] == 'g' ||
         input[y+(2*FLANK+MAXORDER+1)] == 'G')
        ++GCcomp;

      for(i=1;i<QTCT;++i) {
        if(model==NONCODING) {
          if(GCcomp < bounds->noncodbound[i])
            break;
        }
        else
          if(GCcomp < bounds->codingbound[i])
            break;
      }
      probsptr = &probs[i-1];

      /* Update hexamer counts across this window.  *
       * Note that in the code above that processed *
       * the 5' end of the sequence, testcts were   *
       * already accumulated in each frame.  Here,  *
       * we need only tweak those counts based on   *
       * the old/new ends of the shifting window.   *
       * Remember that localmodel indexes the model *
       * of the base just prior to the current      *
       * window.                                    */
      --testcts[localmodel]
          [nt[y]][nt[y+1]][nt[y+2]][nt[y+3]][nt[y+4]][nt[y+5]];
      if(model==NONCODING)
        ++testcts[localmodel]
            [nt[y+(2*FLANK+1)]][nt[y+(2*FLANK+2)]][nt[y+(2*FLANK+3)]]
            [nt[y+(2*FLANK+4)]][nt[y+(2*FLANK+5)]][nt[y+(2*FLANK+6)]];
      else if(model < 3)
        ++testcts[(localmodel+(2*FLANK+1))%3]
            [nt[y+(2*FLANK+1)]][nt[y+(2*FLANK+2)]][nt[y+(2*FLANK+3)]]
            [nt[y+(2*FLANK+4)]][nt[y+(2*FLANK+5)]][nt[y+(2*FLANK+6)]];
      else /* shadow strand */
        ++testcts[((localmodel+(2*FLANK+1))%3)+3]
            [nt[y+(2*FLANK+1)]][nt[y+(2*FLANK+2)]][nt[y+(2*FLANK+3)]]
            [nt[y+(2*FLANK+4)]][nt[y+(2*FLANK+5)]][nt[y+(2*FLANK+6)]];

      /* What model is the central hexamer in? */
      if(model==NONCODING)
        centhexmod = model;
      else if(model < 3)
        centhexmod = (localmodel+(FLANK+1))%3;
      else /* shadow strand */
        centhexmod = ((localmodel+(FLANK+1))%3)+3;

      /* copy working hexamer */
      for(x=0;x<MAXORDER+1;++x)
        locseq[x]=nt[y+(FLANK+1)+x];

      /* rerun simplex to estimate transition probabilities */
      v=0;
      for(aa=0;aa<LPsize;++aa) {
        N[aa]=Ntmp[aa];
        for(bb=0;bb<LPsize;++bb)
          A[aa][bb]=Atmp[aa][bb];
        /* pad b and c with zeros prior to their "real" inits */
        b[aa]=0;
        c[aa]=0;
      }
      for(aa=1;aa<=4;++aa)
        c[aa]=testcts[centhexmod][locseq[0]]
                                 [locseq[1]]
                                 [locseq[2]]
                                 [locseq[3]]
                                 [locseq[4]][aa-1];
      for(aa=5;aa<=11;aa+=2) {
        if(cloudtype == RANGEBASED) {
          b[aa]=probsptr->max[centhexmod][locseq[0]]
                                         [locseq[1]]
                                         [locseq[2]]
                                         [locseq[3]]
                                         [locseq[4]][(aa-5)/2];
          b[aa+1]=-probsptr->min[centhexmod][locseq[0]]
                                            [locseq[1]]
                                            [locseq[2]]
                                            [locseq[3]]
                                            [locseq[4]][(aa-5)/2];
        }
        else { /* (cloudtype == VARIANCEBASED) */
          b[aa]=probsptr->muplus1sigma[centhexmod][locseq[0]]
                                                  [locseq[1]]
                                                  [locseq[2]]
                                                  [locseq[3]]
                                                  [locseq[4]][(aa-5)/2];
          b[aa+1]=-probsptr->muless1sigma[centhexmod][locseq[0]]
                                                     [locseq[1]]
                                                     [locseq[2]]
                                                     [locseq[3]]
                                                     [locseq[4]][(aa-5)/2];
        }
      }
      b[13]=1;
      b[14]=-1;

      /* Solve using my simplex implementation */
      solution=NULL;
      solution=simplex(LPsize,1,N,A,b,c,&v,LPm,LPn,&status,solution);
      if(solution==NULL){
        fprintf(stderr,"Err (DMMMprobQT): \
Simplex unexpectedly failed to compute a solution!\n");
        exit(EXIT_FAILURE);
      }

      /* Increment logprob appropriately */
      logprob+=log(solution[locseq[5]+1]);
      free(solution);

      /* update localmodel */
      if(model==NONCODING)
        ;
      else if(model < 3)
        localmodel=((localmodel+1)%3);
      else /* shadow strand */
        localmodel=((localmodel+1)%3)+3;
    } /* end window stepping for */

    /* Score remaining hexamers in 3' end of final window. *
     * Note that the central hexamer in the last window    *
     * was assayed in the loop above, that appropriate     *
     * quantile-specific parameters were already selected, *
     * that the hexamer counts across the window are       *
     * recorded, and now localmodel indexes the model of   *
     * the first base of the last window.                  */
    for(y+=FLANK+1;y<length-MAXORDER;++y) {
      /* copy working hexamer */
      for(x=0;x<MAXORDER+1;++x)
        locseq[x]=nt[y+x];

      /* Calculate appropriate model for current hexamer */
      if(model==NONCODING)
        ;
      else if(model < 3)
        localmodel=((model+y)%3);
      else /* shadow strand */
        localmodel=((model+y)%3)+3;

      /* Set everything up for simplex.  Note that variables x_1, *
       * x_2, x_3, and x_4 correspond to bases A,C,G,T (0,1,2,3), *
       * respectively.                                            */
      v=0;
      for(aa=0;aa<LPsize;++aa) {
        N[aa]=Ntmp[aa];
        for(bb=0;bb<LPsize;++bb)
          A[aa][bb]=Atmp[aa][bb];
        /* pad b and c with zeros prior to their "real" inits */
        b[aa]=0;
        c[aa]=0;
      }
      for(aa=1;aa<=4;++aa)
        c[aa]=testcts[localmodel][locseq[0]]
                                 [locseq[1]]
                                 [locseq[2]]
                                 [locseq[3]]
                                 [locseq[4]][aa-1];
      for(aa=5;aa<=11;aa+=2) {
        if(cloudtype == RANGEBASED) {
          b[aa]=probsptr->max[localmodel][locseq[0]]
                                         [locseq[1]]
                                         [locseq[2]]
                                         [locseq[3]]
                                         [locseq[4]][(aa-5)/2];
          b[aa+1]=-probsptr->min[localmodel][locseq[0]]
                                            [locseq[1]]
                                            [locseq[2]]
                                            [locseq[3]]
                                            [locseq[4]][(aa-5)/2];
        }
        else { /* (cloudtype == VARIANCEBASED) */
          b[aa]=probsptr->muplus1sigma[localmodel][locseq[0]]
                                                  [locseq[1]]
                                                  [locseq[2]]
                                                  [locseq[3]]
                                                  [locseq[4]][(aa-5)/2];
          b[aa+1]=-probsptr->muless1sigma[localmodel][locseq[0]]
                                                     [locseq[1]]
                                                     [locseq[2]]
                                                     [locseq[3]]
                                                     [locseq[4]][(aa-5)/2];
        }
      }
      b[13]=1;
      b[14]=-1;

      /* Solve using my simplex implementation */
      solution=NULL;
      solution=simplex(LPsize,1,N,A,b,c,&v,LPm,LPn,&status,solution);
      if(solution==NULL){
        fprintf(stderr,"Err (DMMMprobQT): \
Simplex unexpectedly failed to compute a solution!\n");
        exit(EXIT_FAILURE);
      }

      /* Increment logprob appropriately */
      logprob+=log(solution[locseq[5]+1]);
      free(solution);
    
    } /* end 3' for */
  } /* end else */

  /* back to the heap */
  free(nt);
  for(i=0;i<LPsize;++i)
    free(A[i]);
  free(A);

  if(isinf(logprob))
    fprintf(stderr,"Warning (DMMMprobQT): \
Likelihood score grew to infinity!\n");

  return(logprob);
} /* end DMMMprobQT */
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
long double DMMMprob(char *input,dmprobT *probs,int cloudtype,
  int model,int threeper,int shadow)
{
  int length, /* Length of test sequence instance  */
      /* stores counts of oligos in input sequence */
      testcts[ALFSIZE][ALFSIZE][ALFSIZE]
             [ALFSIZE][ALFSIZE][ALFSIZE],
      *nt=NULL,          /* int translation of input               */
      N[LPsize],         /* Store nonbasic status of each variable */
      Ntmp[LPsize]={0,1,1,1,1,0,0,0,0,0,0,0,0,0,0}, /*  copy to N  */
      status,            /* record status from simplex             */
      relevantp,         /* is some history pertinent to test seq? */
      localmodel,        /* helps index three-periodic cases       */
      i,j,k,l,m,n,
      z,aa,bb;           /* iterator variables                     */
  double **A=NULL,       /* Matrix of nonbasic var coefficients    */
         Atmp[LPsize][LPsize] =
           {  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              { 0,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
           },             /* This matrix will be invariant (A=Atmp)  */
         b[LPsize],       /* constraints concerning basic vars       */
         c[LPsize],       /* Nonbasic variable coefficients that     *
                           * occur in objective function.            */
         v,               /* Obj fct's current value (not used)      */
         *solution=NULL;  /* Stores information on optimal solution  */
  long double logprob=0.0;/* Stores overall log-probability of input */
  #ifdef DMMCREPORT
  double pmfsum; /* Check whether DMMM-derived PMF is valid */
  #endif

  /* Verify that model is given a reasonable value */
  if (shadow == FALSE) {
    if (threeper == FALSE) {
      if (model < 0 || model > NONCODING) {
        fprintf(stderr,"Err (DMMMprob): Invalid model value (%i)\n",model);
        exit(EXIT_FAILURE);
      }
    }
    else {
      if ((model < 0 || model > 2) && model != NONCODING) {
        fprintf(stderr,"Err (DMMMprob): Invalid model value (%i)\n",model);
        exit(EXIT_FAILURE);
      }
    }
  }
  else {
    if (model < 3 || model > 5) {
      fprintf(stderr,"Err (DMMMprob): Invalid model value (%i)\n",model);
      exit(EXIT_FAILURE);
    }
  }

  /* Verify cloudtype is set reasonably */
  if (cloudtype != VARIANCEBASED && cloudtype != RANGEBASED) {
    fprintf(stderr,"Err (DMMMprob): \
Invalid cloudtype specified (%i)\n",cloudtype);
    exit(EXIT_FAILURE);
  }

  /* Verify that we really have a sequence to test. If input is set to NULL, *
   * this conditional will short-circuit and strlen won't get called on it.  */
  if((input==NULL) || ((length=(int)strlen(input))==0)) {
    if (input==NULL)
      fprintf(stderr,"Err (DMMMprob): No real sequence passed!\n");
    else
      fprintf(stderr,"Err (DMMMprob): Passed a 0nt test sequence!\n");
    exit(EXIT_FAILURE);
  }
  else if((nt=(int*)malloc(sizeof(int)*length))==NULL) {
    fprintf(stderr,"Err (DMMMprob): Out of memory!\n");
    exit(EXIT_FAILURE);
  }
  else /* Copy input translation into nt */
    for(i=0;i<length;++i)
      nt[i]=trans(input[i]);

  /* A is basically a 2d array--I have to develop it dynamically in *
   * order to keep my simplex implentation as generic as possible,  *
   * which is overall a very minor inconvenience, IMHO.             */
  if((A=(double**)malloc(sizeof(double*)*LPsize))!=NULL){
    for(i=0;i<LPsize;++i){
      if((A[i]=(double*)malloc(sizeof(double)*LPsize))==NULL){
        fprintf(stderr,"Err (DMMMprob): Out of memory!\n");
        exit(EXIT_FAILURE);
      }
    }
  }
  else {
    fprintf(stderr,"Err (DMMMprob): Out of memory!\n");
    exit(EXIT_FAILURE);
  }

  /* The idea is to return the probability of the sequence given *
   * the first base occurs in the frame denoted by model. z will *
   * denote the frame the input is read from, and localmodel the *
   * set of transition probability estimates to use.             */
  localmodel=model;
  for(z=0;z<3;++z){
    /* Initialize count (parameter) array */
    for(i=0;i<ALFSIZE;++i)
    for(j=0;j<ALFSIZE;++j)
    for(k=0;k<ALFSIZE;++k)
    for(l=0;l<ALFSIZE;++l)
    for(m=0;m<ALFSIZE;++m)
    for(n=0;n<ALFSIZE;++n)
      testcts[i][j][k][l][m][n]=0;

    /* Set parameters for this test instance */
    if(model==NONCODING || threeper==FALSE)
      for(i=0;i<length-5;++i)
        ++testcts[nt[i]][nt[i+1]][nt[i+2]][nt[i+3]][nt[i+4]][nt[i+5]];
    else
      for(i=z;i<length-5;i+=3)
        ++testcts[nt[i]][nt[i+1]][nt[i+2]][nt[i+3]][nt[i+4]][nt[i+5]];

    /* Find those oligomers that are relevant to this test instance */
    for(i=0;i<ALFSIZE;++i) {
    for(j=0;j<ALFSIZE;++j) {
    for(k=0;k<ALFSIZE;++k) {
    for(l=0;l<ALFSIZE;++l) {
    for(m=0;m<ALFSIZE;++m) {
      relevantp=0;
      for(n=0;n<ALFSIZE;++n) {
        if(testcts[i][j][k][l][m][n]>0) {
          relevantp=1;
          break;
        }
      }

      /* Executes if we've found a relevant oligomer */
      if(relevantp){
        /* Set everything up for simplex.  Note that variables x_1, *
         * x_2, x_3, and x_4 correspond to bases A,C,G,T (0,1,2,3), *
         * respectively.                                            */
        v=0;
        for(aa=0;aa<LPsize;++aa) {
          N[aa]=Ntmp[aa];
          for(bb=0;bb<LPsize;++bb)
            A[aa][bb]=Atmp[aa][bb];
          /* pad b and c with zeros prior to their "real" inits */
          b[aa]=0;
          c[aa]=0;
        }
        for(aa=1;aa<=4;++aa)
          c[aa]=testcts[i][j][k][l][m][aa-1];
        for(aa=5;aa<=11;aa+=2) {
          if(cloudtype == RANGEBASED) {
            b[aa]=probs->max[localmodel][i][j][k][l][m][(aa-5)/2];
            b[aa+1]=-probs->min[localmodel][i][j][k][l][m][(aa-5)/2];
          }
          else { /* (cloudtype == VARIANCEBASED) */
            b[aa]=probs->muplus1sigma[localmodel][i][j][k][l][m][(aa-5)/2];
            b[aa+1]=-probs->muless1sigma[localmodel][i][j][k][l][m][(aa-5)/2];
          }
        }
        b[13]=1;
        b[14]=-1;

        /* Solve using my simplex implementation */
        solution=NULL;
        solution=simplex(LPsize,1,N,A,b,c,&v,LPm,LPn,&status,solution);
        if(solution==NULL){
          fprintf(stderr,"Err (DMMMprob): \
Simplex unexpectedly failed to compute a solution!\n");
          exit(EXIT_FAILURE);
        }

        /* Increment logprob appropriately */
        #ifdef DMMCREPORT
        pmfsum=0.0;
        #endif
        for(aa=0;aa<ALFSIZE;++aa) {
          if(testcts[i][j][k][l][m][aa]>0)
            logprob+=testcts[i][j][k][l][m][aa]*log(solution[aa+1]);
          #ifdef DMMCREPORT
          printf("%i=%.4f,",aa,solution[aa+1]);
          pmfsum+=solution[aa+1];
          #endif
        }
        #ifdef DMMCREPORT
        printf("Sum of prob's for history %i%i%i%i%i (model %i) is %.2f ",
          i,j,k,l,m,model,pmfsum);
        if(fabs(pmfsum-CERTPROB) > MAXDIFF)
          printf("problematic!\n");
        else
          printf("valid\n");
        #endif

        /* return to heap */
        free(solution);
      } /* end if relevantp */
    }}}}}

    /* update localmodel */
    if(model==NONCODING || threeper==FALSE)
      break;
    else if(shadow==FALSE)
      localmodel=((localmodel+1)%3);
    else /* (shadow==TRUE) */
      localmodel=((localmodel+1)%3)+3;
  } /* end localmodel for */

  /* back to the heap */
  free(nt);
  for(i=0;i<LPsize;++i)
    free(A[i]);
  free(A);

  if(isinf(logprob))
    fprintf(stderr,"Warning (DMMMprob): \
Likelihood score grew to infinity!\n");

  return(logprob);
} /* end DMMMprob */
#endif

#ifdef QUANTSPEC
/* Function to import a trained dmprobT object array and a qtboundsT *
 * object (stored in a binary file) into main memory.                */
void import_dm_probsQT(char *filename,dmprobT *probs,qtboundsT *bounds)
{
  FILE *fptr=NULL; /* For connecting to binary input file */

  /* Open parameter file */
  if((fptr=fopen(filename,"rb"))==NULL) {
    fprintf(stderr,"Err (import_dm_probsQT): Can't open %s for \
binary reading!\n",filename);
    exit(EXIT_FAILURE);
  }
  /* Read in the dmprobT "objects" */
  else if(fread(probs,sizeof(dmprobT),QTCT,fptr) != QTCT) {
    fprintf(stderr,"Err (import_dm_probsQT): Error reading from %s!\n",
      filename);
    exit(EXIT_FAILURE);
  }
  /* Read in the qtboundsT "object" */
  else if(fread(bounds,sizeof(qtboundsT),1,fptr) != 1) {
    fprintf(stderr,"Err (import_dm_probsQT): Error reading from %s!\n",
      filename);
    exit(EXIT_FAILURE);
  }
  else
    (void)fclose(fptr);

  return;
} /* end import_dm_probsQT */
#else
/* Function to import a trained dmprobT     *
 * object (a binary file) into main memory. */
void import_dm_probs(char *filename,dmprobT *probs)
{
  FILE *fptr=NULL; /* For connecting to binary input file */

  /* Open parameter file */
  if((fptr=fopen(filename,"rb"))==NULL) {
    fprintf(stderr,"Err (import_dm_probs): Can't open %s for \
binary reading!\n",filename);
    exit(EXIT_FAILURE);
  }
  /* Read in the dmprobT "object" */
  else if(fread(probs,sizeof(dmprobT),1,fptr) != 1) {
    fprintf(stderr,"Err (import_dm_probs): Error reading from %s!\n",
      filename);
    exit(EXIT_FAILURE);
  }
  else
    (void)fclose(fptr);

  return;
} /* end import_dm_probs */
#endif
