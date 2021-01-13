/* immpractical.c
 * Michael E Sparks (mespar1@gmail.com)
 *
 * This file contains the core code of the IMMpractical library.
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
#include "chisquare_utils.h"
#include "dimm_utils.h"
#include "fo_utils.h"
#include "dm_utils.h"
#include "immpractical.h"
#include "sequence_parse.h"
#include "sorting.h"

/* Macro definitions *******************************************************/

#define imm_probmakerUSAGE \
  "\a\nUsage: %s %i inD-coding inD-noncod inH-coding inH-noncod out_file\n\n"

/* These two macros will alert the user of the proper usage of *
 * the driver. Note: algodesc array defined in immpractical.h  */
#define GENERALUSAGE \
  "\a\nConceptual Usage:\n\n\
  %s algorithm_type in_file-coding in_file-noncod out_file\n\n\
    Where algorithm_type is an integer selected from\n"

/* We do not list PMMM's in the training algorithm repertoire *
 * as it is not a training method. (It uses FIXORD and DMMM)  */
#define BADUSE {                                      \
  fprintf(stderr,GENERALUSAGE,argv[0]);               \
  for(i=0;i<NUMALGOS;++i)                             \
    if(i!=PMMMX)                                      \
      fprintf(stderr,"    (%i) %s\n",i,algodesc[i]);  \
  fprintf(stderr,"\n    (Note that all sequence files \
are to be in Fasta format,\n      with each sequence \
contained entirely on one line.)\n\n");               \
}

/* Local function/subroutine prototypes ************************************/

#ifdef QUANTSPEC
/* To account for local as well as global G+C compositional heterogeneity *
 * in DNA sequences, this function will develop lists of the G+C content  *
 * of 2*FLANK+6 base long windows (100nt) shifted by SHIFTLEN bases.      *
 * GCcounts is STRICTLY a 2x10x10 matrix, intented to index counts of     *
 * test windows exhibiting a G+C compositional percentage [0..100]        */
static void ctGCpercs(FILE *infile,int GCcounts[][10][10]);

/* This function computes the actual quantile boundaries *
 * based on GCcounts, and records these in bounds.       */
static void find_quantiles(int GCcounts[][10][10],
  int codingp,qtboundsT *bounds);
#endif

/* Function implementations ************************************************/

/* This function coordinates invocation of the appropriate *
 * training routine depending on which algorithm the user  *
 * wishes to employ for his/her data.                      */
int maintrain(int argc,char *argv[])
{
  int algtype, /* Which of the NUMALGOS training algorithms *
                * does the user wish to employ?             */
      status,  /* Return status of invoked function         */
      i;       /* Iterator variable                         */

  /* First, determine what training algorithm the user requests */
  if (argc < 2) {
    BADUSE;
    return(FALSE);
  }
  else
    algtype=atoi(argv[1]);

  /* If the library is intended for quantile-specific training and testing, *
   * verify that all the requisite pre-conditions are met.                  */
  #ifdef QUANTSPEC
  if(MAXORDER != 5) {
    fprintf(stderr,"Err (maintrain): Sorry, quantile-specific \
functionality requires MAXORDER == 5! (see immpractical.h)\n");
    exit(EXIT_FAILURE);
  }
  if(FLANK != 47) {
    fprintf(stderr,"Err (maintrain): \
Only FLANK == 47 is supported!\n");
    exit(EXIT_FAILURE);
  }
  #ifndef THREEPERIODIC
  fprintf(stderr,"Err (maintrain): Sorry, quantile-specific \
functionality requires THREEPERIODIC and SHADOWSTRAND be defined!\n");
  exit(EXIT_FAILURE);
  #endif
  #ifndef SHADOWSTRAND
  fprintf(stderr,"Err (maintrain): Sorry, quantile-specific \
functionality requires THREEPERIODIC and SHADOWSTRAND be defined!\n");
  exit(EXIT_FAILURE);
  #endif
  #endif

  /* Verify parameters and invoke the proper training function */
  switch(algtype) {
    case FIXORDX :
      if(argc != 5) {
        fprintf(stderr,fo_probmakerUSAGE,argv[0],algtype);
        status=FALSE;
      }
      else
        status=fo_probmaker(argc,argv);
      break;
    case TDDIX :
    case BUDIX :
    case CHI2X :
      if(argc != 7) {
        fprintf(stderr,imm_probmakerUSAGE,argv[0],algtype);
        status=FALSE;
      }
      /* pass algtype to imm_probmaker, which will *
       * distinguish between TDDI, BUDI, and CHI2  */
      else
        status=imm_probmaker(argc,argv,algtype);
      break;
    case DMMMX :
      if((int)MAXORDER != 5) {
        fprintf(stderr,"\nSorry, DMMMs (and PMMMs) require \
MAXORDER == 5!\n\n");
        status=FALSE;
      }
      else if(argc < 6 || argc != ((2*atoi(argv[2]))+4) ) {
        fprintf(stderr,dm_probmakerUSAGE,argv[0],algtype);
        status=FALSE;
      }
      else
        status=dm_probmaker(argc,argv);
      break;
    default :
      fprintf(stderr,"\nErr (maintrain): \
Improper algorithm requested!\n\n");
      status=FALSE;
  } /* end switch */

  return(status);
} /* end maintrain */

/* Initialize the data structure for recording counts */
void init_counts(countarrayT *tal,int model)
{
  int i,j,k,l,m,n; /* Iterator variables */

  for(i=0;i<ALFSIZE;++i) {
    tal->count1[model][i]=INITFREQ;
    for(j=0;j<ALFSIZE;++j) {
      tal->count2[model][i][j]=INITFREQ;
      for(k=0;k<ALFSIZE;++k) {
        tal->count3[model][i][j][k]=INITFREQ;
        for(l=0;l<ALFSIZE;++l) {
          tal->count4[model][i][j][k][l]=INITFREQ;
          for(m=0;m<ALFSIZE;++m) {
            tal->count5[model][i][j][k][l][m]=INITFREQ;
            for(n=0;n<ALFSIZE;++n) {
              tal->count6[model][i][j][k][l][m][n]=INITFREQ;
  }}}}}}

  return;
} /* end init_counts */

#ifdef QUANTSPEC
/* Function to record counts occuring in input sequences,  *
 * nonperiodically. Please notice that each call requires  *
 * that the sequence stored in seq be at least 2*FLANK+    *
 * MAXORDER+1 residues long! This function, which develops *
 * the NONCODING model in the case of quantile-specific    *
 * library functionality, deposits counts of hexamers in   *
 * bins based on G+C composition of their local contexts   */
void get_countsQT(char *seq,countarrayT *tal,qtboundsT bounds)
{
  int i,j,             /* Iterator variables                  */
      GCcomp,          /* Store total of G+C in seq           */
      calcoverall=0,   /* Should we consider overall [GC]?    */
      size,            /* record length of seq                */
      *seqi=NULL;      /* stores int translation of seq       */
  countarrayT *talptr; /* point to appropriate element of tal */

  /* Verify that we really have a sequence to test. If seq is set to NULL,  *
   * this conditional will short-circuit and strlen won't get called on it. */
  if((seq==NULL) || ((size=(int)strlen(seq))<(2*FLANK+MAXORDER+1))) {
    if (seq==NULL) {
      fprintf(stderr,"Err (get_countsQT): No real sequence passed!\n");
      exit(EXIT_FAILURE);
    }
    else {
      if(size<MAXORDER+1) { /* intractable */
        fprintf(stderr,"Err (get_counts): First run %s on your input!\n",
          DATACLEANPL);
        exit(EXIT_FAILURE);
      }
      /* If the sequence is a tad too short, we'll just *
       * calculate its overall [GC], then record counts *
       * based on that.                                 */
      calcoverall=1;
    }
  }

  if((seqi=(int*)malloc(sizeof(int)*size))==NULL) {
    fprintf(stderr,"Err (get_countsQT): Out of memory!\n");
    exit(EXIT_FAILURE);
  }
  else /* Copy over sequence */
    for(i=0;i<size;++i)
      seqi[i]=trans(seq[i]);

  /* A lower bound of 100nt on intron length is unrealistic *
   * for some species, particularly plants.  Rather than    *
   * discard these valuable data, compute the overall [GC]  *
   * and record the results to that quantile.               */
  if(calcoverall) {
    /* Calc overall [GC] and cast to an integer *
     * precision to 100th's only.               */
    for(j=0,GCcomp=0;j<size;++j)
      if(seq[j] == 'c' || seq[j] == 'C' ||
         seq[j] == 'g' || seq[j] == 'G')
        ++GCcomp;
    GCcomp=(int)floor(100*(GCcomp/(float)size));

    /* To which quantile does this [G+C] belong? Recall that, as *
     * quantile-specific functionality requires NUMMODELS == 7,  *
     * this function only deals with noncoding sequences.        */
    for(j=1;j<QTCT;++j)
      if(GCcomp < bounds.noncodbound[j])
        break;
    talptr = &tal[j-1];

    /* Do the counting */
    for(i=0;i<size;++i)
      talptr->count1[NONCODING][seqi[i]]++;

    for(i=0;i<size-1;++i)
      talptr->count2[NONCODING][seqi[i]]
                               [seqi[i+1]]++;

    for(i=0;i<size-2;++i)
      talptr->count3[NONCODING][seqi[i]]
                               [seqi[i+1]]
                               [seqi[i+2]]++;

    for(i=0;i<size-3;++i)
      talptr->count4[NONCODING][seqi[i]]
                               [seqi[i+1]]
                               [seqi[i+2]]
                               [seqi[i+3]]++;

    for(i=0;i<size-4;++i)
      talptr->count5[NONCODING][seqi[i]]
                               [seqi[i+1]]
                               [seqi[i+2]]
                               [seqi[i+3]]
                               [seqi[i+4]]++;

    for(i=0;i<size-5;++i)
      talptr->count6[NONCODING][seqi[i]]
                               [seqi[i+1]]
                               [seqi[i+2]]
                               [seqi[i+3]]
                               [seqi[i+4]]
                               [seqi[i+5]]++;

  }
  else { /* windowing is possible */
    /* As I expect 2*FLANK+MAXORDER+1 == 100, GCcomp will *
     * automatically represent a percentage.              */
    for(j=0,GCcomp=0;j<(2*FLANK+MAXORDER+1);++j)
      if(seq[j] == 'c' || seq[j] == 'C' ||
         seq[j] == 'g' || seq[j] == 'G')
        ++GCcomp;

    /* Compute G+C composition of the hexamer's local context *
     * to determine which of the QTCT bins it will be         *
     * recorded in.                                           */
    for(i=0;i<=(size-(2*FLANK+MAXORDER+1));++i) {

      /* To which quantile does this [G+C] belong? Recall that, as *
       * quantile-specific functionality requires NUMMODELS == 7,  *
       * this function only deals with noncoding sequences.        */
      for(j=1;j<QTCT;++j)
        if(GCcomp < bounds.noncodbound[j])
          break;
      talptr = &tal[j-1];

      /* Now, go ahead and store the counts */
      j=i+FLANK;
      talptr->count1[NONCODING][seqi[j]]++;

      talptr->count2[NONCODING][seqi[j]]
                               [seqi[j+1]]++;

      talptr->count3[NONCODING][seqi[j]]
                               [seqi[j+1]]
                               [seqi[j+2]]++;

      talptr->count4[NONCODING][seqi[j]]
                               [seqi[j+1]]
                               [seqi[j+2]]
                               [seqi[j+3]]++;

      talptr->count5[NONCODING][seqi[j]]
                               [seqi[j+1]]
                               [seqi[j+2]]
                               [seqi[j+3]]
                               [seqi[j+4]]++;

      talptr->count6[NONCODING][seqi[j]]
                               [seqi[j+1]]
                               [seqi[j+2]]
                               [seqi[j+3]]
                               [seqi[j+4]]
                               [seqi[j+5]]++;

      /* update GCcomp, if necessary */
      if(i+(2*FLANK+MAXORDER+1)<size) {
        if(seq[i] == 'c' || seq[i] == 'C' ||
           seq[i] == 'g' || seq[i] == 'G')
          --GCcomp;
        if(seq[i+(2*FLANK+MAXORDER+1)] == 'c' ||
           seq[i+(2*FLANK+MAXORDER+1)] == 'C' ||
           seq[i+(2*FLANK+MAXORDER+1)] == 'g' ||
           seq[i+(2*FLANK+MAXORDER+1)] == 'G')
          ++GCcomp;
      }

    } /* end for */
  } /* end else */

  free(seqi);
  return;
} /* end get_countsQT */
#else
/* Function to record counts occuring in input sequences, *
 * nonperiodically. Please notice that each call requires *
 * that the sequence stored in seq be at least MAXORDER+1 *
 * residues long!                                         */
void get_counts(char *seq,countarrayT *tal,int model)
{
  int i,          /* Iterator variable             */
      size,       /* record length of seq          */
      *seqi=NULL; /* stores int translation of seq */

  /* Verify that we really have a sequence to test. If seq is set to NULL,  *
   * this conditional will short-circuit and strlen won't get called on it. */
  if((seq==NULL) || ((size=(int)strlen(seq))<(MAXORDER+1))) {
    if (seq==NULL)
      fprintf(stderr,"Err (get_counts): No real sequence passed!\n");
    else
      fprintf(stderr,"Err (get_counts): First run %s on your input!\n",
        DATACLEANPL);
    exit(EXIT_FAILURE);
  }
  else if((seqi=(int*)malloc(sizeof(int)*size))==NULL) {
    fprintf(stderr,"Err (get_counts): Out of memory!\n");
    exit(EXIT_FAILURE);
  }
  else /* Copy over sequence */
    for(i=0;i<size;++i)
      seqi[i]=trans(seq[i]);

  /* NOTE: I don't train on every charcter in every training    *
   * sequence (only go to size-MAXORDER for all history depths. *
   * This way counts are consistent across differing orders.    */

  for(i=0;i<(size-MAXORDER);++i) {
    /* Process historydepth=0 frequencies *
     * This statement must always exectue */
    tal->count1[model][seqi[i]]++;

    if (MAXORDER > 0)
      /* Process historydepth=1 frequencies */
      tal->count2[model][seqi[i]]
                        [seqi[i+1]]++;
    if (MAXORDER > 1)
      /* Process historydepth=2 frequencies */
      tal->count3[model][seqi[i]]
                        [seqi[i+1]]
                        [seqi[i+2]]++;
    if (MAXORDER > 2)
      /* Process historydepth=3 frequencies */
      tal->count4[model][seqi[i]]
                        [seqi[i+1]]
                        [seqi[i+2]]
                        [seqi[i+3]]++;
    if (MAXORDER > 3)
      /* Process historydepth=4 frequencies */
      tal->count5[model][seqi[i]]
                        [seqi[i+1]]
                        [seqi[i+2]]
                        [seqi[i+3]]
                        [seqi[i+4]]++;
    if (MAXORDER > 4)
      /* Process historydepth=5 frequencies */
      tal->count6[model][seqi[i]]
                        [seqi[i+1]]
                        [seqi[i+2]]
                        [seqi[i+3]]
                        [seqi[i+4]]
                        [seqi[i+5]]++;
  } /* end for */

  free(seqi);
  return;
} /* end get_counts */
#endif

#ifdef QUANTSPEC
/* This function will develop models [0-NONCODING) if the  *
 * library is intended to support three-periodicity.  Note *
 * that all training data is assumed to be in frame 0!     */
void get_counts3percodQT(char *seq,countarrayT *tal,qtboundsT bounds)
{
  int GCcomp,          /* Store total of G+C in seq            */
      size,            /* record length of seq                 */
      i,j,             /* Iterator variables                   */
      shadmod,         /* stores model index for shadow        */
      *seqi=NULL,      /* stores int translation of seq        */
      *seqicomp=NULL;  /* stores complement of int translation *
                        * of seq                               */
  countarrayT *talptr; /* point to appropriate element of tal  */

  /* Verify that we really have a sequence to test. If seq is set to NULL,  *
   * this conditional will short-circuit and strlen won't get called on it. */
  if((seq==NULL) || ((size=(int)strlen(seq))<(2*FLANK+MAXORDER+1))) {
    if (seq==NULL) {
      fprintf(stderr,"Err (get_counts3percodQT): No real sequence passed!\n");
      exit(EXIT_FAILURE);
    }
    else {
      /* Sorry, but a coding sequence < 100nts *
       * seems better off ignored.  Perhaps    *
       * a different policy would be in order  *
       * if I were training on exons, but in   *
       * practice, I'm not.                    */
      fprintf(stderr,"Warning (get_counts3percodQT): \
Sequence < %i bases! Ignoring it.\n",2*FLANK+MAXORDER+1);
      return;
    }
  }
  else if((seqi=(int*)malloc(sizeof(int)*size))==NULL) {
    fprintf(stderr,"Err (get_counts3percodQT): Out of memory!\n");
    exit(EXIT_FAILURE);
  }
  else if((seqicomp=(int*)malloc(sizeof(int)*size))==NULL) {
    fprintf(stderr,"Err (get_counts3percodQT): Out of memory!\n");
    exit(EXIT_FAILURE);
  }
  else {
    for(i=0;i<size;++i) {
      seqi[i]=trans(seq[i]);
      seqicomp[i]=basecomp(trans(seq[i]));
    }
  }

  /* As I expect 2*FLANK+MAXORDER+1 == 100, GCcomp will *
   * automatically represent a percentage.              */
  for(i=0,GCcomp=0;i<(2*FLANK+MAXORDER+1);++i)
    if(seq[i] == 'c' || seq[i] == 'C' ||
       seq[i] == 'g' || seq[i] == 'G')
      ++GCcomp;

  /* Progress through windows */
  for(i=0;i<=size-(2*FLANK+MAXORDER+1);++i) {

    /* To which quantile does this [G+C] belong? */
    for(j=1;j<QTCT;++j)
      if(GCcomp < bounds.codingbound[j])
        break;
    talptr = &tal[j-1];

    /* Note that all training data is assumed to be in frame 0!!! */
    j=i+FLANK;
    talptr->count1[j%3][seqi[j]]++;

    talptr->count2[j%3][seqi[j]]
                       [seqi[j+1]]++;

    talptr->count3[j%3][seqi[j]]
                       [seqi[j+1]]
                       [seqi[j+2]]++;

    talptr->count4[j%3][seqi[j]]
                       [seqi[j+1]]
                       [seqi[j+2]]
                       [seqi[j+3]]++;

    talptr->count5[j%3][seqi[j]]
                       [seqi[j+1]]
                       [seqi[j+2]]
                       [seqi[j+3]]
                       [seqi[j+4]]++;

    talptr->count6[j%3][seqi[j]]
                       [seqi[j+1]]
                       [seqi[j+2]]
                       [seqi[j+3]]
                       [seqi[j+4]]
                       [seqi[j+5]]++;

    /* Now record counts on shadow strand */
    if(j%3==0)
      shadmod=3;
    else if(j%3==1)
      shadmod=5;
    else /*(j%3==2)*/
      shadmod=4;

    talptr->count1[shadmod][seqicomp[j+5]]++;

    talptr->count2[shadmod][seqicomp[j+5]]
                           [seqicomp[j+4]]++;

    talptr->count3[shadmod][seqicomp[j+5]]
                           [seqicomp[j+4]]
                           [seqicomp[j+3]]++;

    talptr->count4[shadmod][seqicomp[j+5]]
                           [seqicomp[j+4]]
                           [seqicomp[j+3]]
                           [seqicomp[j+2]]++;

    talptr->count5[shadmod][seqicomp[j+5]]
                           [seqicomp[j+4]]
                           [seqicomp[j+3]]
                           [seqicomp[j+2]]
                           [seqicomp[j+1]]++;

    talptr->count6[shadmod][seqicomp[j+5]]
                           [seqicomp[j+4]]
                           [seqicomp[j+3]]
                           [seqicomp[j+2]]
                           [seqicomp[j+1]]
                           [seqicomp[j]]++;

    /* update GCcomp, if necessary */
    if(i+(2*FLANK+MAXORDER+1)<size) {
      if(seq[i] == 'c' || seq[i] == 'C' ||
         seq[i] == 'g' || seq[i] == 'G')
        --GCcomp;
      if(seq[i+(2*FLANK+MAXORDER+1)] == 'c' ||
         seq[i+(2*FLANK+MAXORDER+1)] == 'C' ||
         seq[i+(2*FLANK+MAXORDER+1)] == 'g' ||
         seq[i+(2*FLANK+MAXORDER+1)] == 'G')
        ++GCcomp;
    }

  } /* end i for */

  free(seqi);
  free(seqicomp);
  return;
} /* end get_counts3percodQT */
#else
#ifdef THREEPERIODIC
/* This function will develop models [0-NONCODING) if the  *
 * library is intended to support three-periodicity.  Note *
 * that all training data is assumed to be in frame 0!     *
 * If you wish to tally counts for the reverse strand      *
 * (model indices {3,4,5}, when SHADOWSTRAND is defined),  *
 * in addition to models {0,1,2}, set the shadow argument  *
 * to TRUE.                                                */
void get_counts3percod(char *seq,countarrayT *tal,int shadow)
{
  int model,      /* Iterate through models           */
      oligolen,   /* Oligo lengths of 1..(MAXORDER+1) */
      size,       /* record length of seq             */
      i,j,        /* Iterator variables               */
      *seqi=NULL; /* stores int translation of seq    */

  /* Verify that we really have a sequence to test. If seq is set to NULL,  *
   * this conditional will short-circuit and strlen won't get called on it. */
  if((seq==NULL) || ((size=(int)strlen(seq))<(MAXORDER+1))) {
    if (seq==NULL)
      fprintf(stderr,"Err (get_counts3percod): No real sequence passed!\n");
    else
      fprintf(stderr,"Err (get_counts3percod): First run %s on your input!\n",
        DATACLEANPL);
    exit(EXIT_FAILURE);
  }
  if((seqi=(int*)malloc(sizeof(int)*size))==NULL) {
    fprintf(stderr,"Err (get_counts3percod): Out of memory!\n");
    exit(EXIT_FAILURE);
  }

  /* First, process coding strand */
  for(i=0;i<size;++i)
    seqi[i]=trans(seq[i]);

  /* Note that all training data is assumed to be in frame 0!!! */
  for(model=0;model<3;++model) { /* model in {0,1,2} */
    for(i=model;i<(size-MAXORDER);i+=3) {
      /* At first glance, this next loop may seem superfluous.  *
       * However, it is needed for MAXORDER to range from 0..5. */
      for(oligolen=1;oligolen<=(MAXORDER+1);++oligolen) {
        switch(oligolen) {
          case 1 :
            tal->count1[model][seqi[i]]++;
            break;
          case 2 :
            tal->count2[model][seqi[i]]
                              [seqi[i+1]]++;
            break;
          case 3 :
            tal->count3[model][seqi[i]]
                              [seqi[i+1]]
                              [seqi[i+2]]++;
            break;
          case 4 :
            tal->count4[model][seqi[i]]
                              [seqi[i+1]]
                              [seqi[i+2]]
                              [seqi[i+3]]++;
            break;
          case 5 :
            tal->count5[model][seqi[i]]
                              [seqi[i+1]]
                              [seqi[i+2]]
                              [seqi[i+3]]
                              [seqi[i+4]]++;
            break;
          case 6 :
            tal->count6[model][seqi[i]]
                              [seqi[i+1]]
                              [seqi[i+2]]
                              [seqi[i+3]]
                              [seqi[i+4]]
                              [seqi[i+5]]++;
            break;
          default :
            fprintf(stderr,"Err (get_counts3percod): Invalid oligolen!\n");
            exit(EXIT_FAILURE);
        } /* end switch */
      } /* end oligolen for */
    } /* end i for */
  } /* end model for */

  /* Optionally, process reverse strand */
  if (shadow == TRUE) {
    #ifndef SHADOWSTRAND
    fprintf(stderr,
      "Err (get_counts3percod): \
Did not define SHADOWSTRAND in Makefile!??\n");
    exit(EXIT_FAILURE);
    #endif

    /* take reverse complement of seq */
    for(i=size-1,j=0;i>=0;--i,++j)
      seqi[j]=basecomp(trans(seq[i]));

    /* Here, we may not know a priori which frame (model) the first *
     * base of the reverse complemented sequence will be in.        */
    if ((size%3)==0)
      /* Actually, since all data I use to train the models *
       * is based on putative peptide sequences, it will    *
       * always have length of a multiple of three.  I only *
       * expect this branch to execute, therefore.          */
      model=3;
    else if ((size%3)==2) {
      fprintf(stderr,"Warning (get_counts3percod): \
Unexpected reverse complement length!\n");
      model=4;
    }
    else { /* ((size%3)==1) */
      fprintf(stderr,"Warning (get_counts3percod): \
Unexpected reverse complement length!\n");
      model=5;
    }
    for(i=0;i<(size-MAXORDER);model=((model+1)%3)+3,++i) {
      for(oligolen=1;oligolen<=(MAXORDER+1);++oligolen) {
        switch(oligolen) {
          case 1 :
            tal->count1[model][seqi[i]]++;
            break;
          case 2 :
            tal->count2[model][seqi[i]]
                              [seqi[i+1]]++;
            break;
          case 3 :
            tal->count3[model][seqi[i]]
                              [seqi[i+1]]
                              [seqi[i+2]]++;
            break;
          case 4 :
            tal->count4[model][seqi[i]]
                              [seqi[i+1]]
                              [seqi[i+2]]
                              [seqi[i+3]]++;
            break;
          case 5 :
            tal->count5[model][seqi[i]]
                              [seqi[i+1]]
                              [seqi[i+2]]
                              [seqi[i+3]]
                              [seqi[i+4]]++;
            break;
          case 6 :
            tal->count6[model][seqi[i]]
                              [seqi[i+1]]
                              [seqi[i+2]]
                              [seqi[i+3]]
                              [seqi[i+4]]
                              [seqi[i+5]]++;
            break;
          default :
            fprintf(stderr,"Err (get_counts3percod): Invalid oligolen!\n");
            exit(EXIT_FAILURE);
        } /* end switch */
      } /* end oligolen for */
    } /* end model for */
  } /* end reverse strand */

  free(seqi);
  return;
} /* end get_counts3percod */
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
void induce_pseudocounts(countarrayT *tal,int model,int *develsize)
{
  int i,j,k,l,m,n, /* Iterator variables                 */
      currlen,     /* Tracks oligo length, 1..MAXORDER+1 */
      z;           /* Iterator for updating develsize    */

  for(currlen=1;currlen<=(MAXORDER+1);++currlen) {
    switch(currlen) {
      case 1 :
        for(i=0;i<ALFSIZE;++i) {
          if(tal->count1[model][i] <= EMPTYCOUNT) {
            fprintf(stdout,
              "History %u only occurs %u times in devel data! (model %i)\n",
              i,tal->count1[model][i],model);
            /* Update size of development data */
            for(z=currlen;z<=(MAXORDER+1);++z)
              *develsize+=(int)pow(ALFSIZE,(MAXORDER+1)-z);
            /* set pseudocounts for current oligo */
            tal->count1[model][i]=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            /* propagate pseudocounts forward */
            for(j=0;j<ALFSIZE;++j) {
              tal->count2[model][i][j]=MINCOUNTUNIT *
                (int)pow(ALFSIZE,(MAXORDER+1)-(currlen+1));
              for(k=0;k<ALFSIZE;++k) {
                tal->count3[model][i][j][k]=MINCOUNTUNIT *
                  (int)pow(ALFSIZE,(MAXORDER+1)-(currlen+2));
                for(l=0;l<ALFSIZE;++l) {
                  tal->count4[model][i][j][k][l]=MINCOUNTUNIT *
                    (int)pow(ALFSIZE,(MAXORDER+1)-(currlen+3));
                  for(m=0;m<ALFSIZE;++m) {
                    tal->count5[model][i][j][k][l][m]=MINCOUNTUNIT *
                      (int)pow(ALFSIZE,(MAXORDER+1)-(currlen+4));
                    for(n=0;n<ALFSIZE;++n) {
                      tal->count6[model][i][j][k][l][m][n]=MINCOUNTUNIT *
                        (int)pow(ALFSIZE,(MAXORDER+1)-(currlen+5));
        }}}}}}}
        break;
      case 2 :
        for(i=0;i<ALFSIZE;++i) {
        for(j=0;j<ALFSIZE;++j) {
          if(tal->count2[model][i][j] <= EMPTYCOUNT) {
            fprintf(stdout,
              "History %u%u only occurs %u times in devel data! (model %i)\n",
              i,j,tal->count2[model][i][j],model);
            /* Update size of development data */
            for(z=1;z<currlen;++z)
              *develsize+=(int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            for(z=currlen;z<=(MAXORDER+1);++z)
              *develsize+=(int)pow(ALFSIZE,(MAXORDER+1)-z);
            /* propagate pseudocounts backward */
            tal->count1[model][i]+=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            /* set pseudocounts for current oligo */
            tal->count2[model][i][j]=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            /* propagate pseudocounts forward */
            for(k=0;k<ALFSIZE;++k) {
              tal->count3[model][i][j][k]=MINCOUNTUNIT *
                (int)pow(ALFSIZE,(MAXORDER+1)-(currlen+1));
              for(l=0;l<ALFSIZE;++l) {
                tal->count4[model][i][j][k][l]=MINCOUNTUNIT *
                  (int)pow(ALFSIZE,(MAXORDER+1)-(currlen+2));
                for(m=0;m<ALFSIZE;++m) {
                  tal->count5[model][i][j][k][l][m]=MINCOUNTUNIT *
                    (int)pow(ALFSIZE,(MAXORDER+1)-(currlen+3));
                  for(n=0;n<ALFSIZE;++n) {
                    tal->count6[model][i][j][k][l][m][n]=MINCOUNTUNIT *
                      (int)pow(ALFSIZE,(MAXORDER+1)-(currlen+4));
        }}}}}}}
        break;
      case 3 :
        for(i=0;i<ALFSIZE;++i) {
        for(j=0;j<ALFSIZE;++j) {
        for(k=0;k<ALFSIZE;++k) {
          if(tal->count3[model][i][j][k] <= EMPTYCOUNT) {
            fprintf(stdout,"History %u%u%u only \
occurs %u times in devel data! (model %i)\n",
              i,j,k,tal->count3[model][i][j][k],model);
            /* Update size of development data */
            for(z=1;z<currlen;++z)
              *develsize+=(int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            for(z=currlen;z<=(MAXORDER+1);++z)
              *develsize+=(int)pow(ALFSIZE,(MAXORDER+1)-z);
            /* propagate pseudocounts backward */
            tal->count1[model][i]+=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            tal->count2[model][i][j]+=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            /* set pseudocounts for current oligo */
            tal->count3[model][i][j][k]=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            /* propagate pseudocounts forward */
            for(l=0;l<ALFSIZE;++l) {
              tal->count4[model][i][j][k][l]=MINCOUNTUNIT *
                (int)pow(ALFSIZE,(MAXORDER+1)-(currlen+1));
              for(m=0;m<ALFSIZE;++m) {
                tal->count5[model][i][j][k][l][m]=MINCOUNTUNIT *
                  (int)pow(ALFSIZE,(MAXORDER+1)-(currlen+2));
                for(n=0;n<ALFSIZE;++n) {
                  tal->count6[model][i][j][k][l][m][n]=MINCOUNTUNIT *
                    (int)pow(ALFSIZE,(MAXORDER+1)-(currlen+3));
        }}}}}}}
        break;
      case 4 :
        for(i=0;i<ALFSIZE;++i) {
        for(j=0;j<ALFSIZE;++j) {
        for(k=0;k<ALFSIZE;++k) {
        for(l=0;l<ALFSIZE;++l) {
          if(tal->count4[model][i][j][k][l] <= EMPTYCOUNT) {
            fprintf(stdout,"History %u%u%u%u only \
occurs %u times in devel data! (model %i)\n",
              i,j,k,l,tal->count4[model][i][j][k][l],model);
            /* Update size of development data */
            for(z=1;z<currlen;++z)
              *develsize+=(int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            for(z=currlen;z<=(MAXORDER+1);++z)
              *develsize+=(int)pow(ALFSIZE,(MAXORDER+1)-z);
            /* propagate pseudocounts backward */
            tal->count1[model][i]+=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            tal->count2[model][i][j]+=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            tal->count3[model][i][j][k]+=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            /* set pseudocounts for current oligo */
            tal->count4[model][i][j][k][l]=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            /* propagate pseudocounts forward */
            for(m=0;m<ALFSIZE;++m) {
              tal->count5[model][i][j][k][l][m]=MINCOUNTUNIT *
                (int)pow(ALFSIZE,(MAXORDER+1)-(currlen+1));
              for(n=0;n<ALFSIZE;++n) {
                tal->count6[model][i][j][k][l][m][n]=MINCOUNTUNIT *
                  (int)pow(ALFSIZE,(MAXORDER+1)-(currlen+2));
        }}}}}}}
        break;
      case 5 :
        for(i=0;i<ALFSIZE;++i) {
        for(j=0;j<ALFSIZE;++j) {
        for(k=0;k<ALFSIZE;++k) {
        for(l=0;l<ALFSIZE;++l) {
        for(m=0;m<ALFSIZE;++m) {
          if(tal->count5[model][i][j][k][l][m] <= EMPTYCOUNT) {
            fprintf(stdout,"History %u%u%u%u%u only \
occurs %u times in devel data! (model %i)\n",
              i,j,k,l,m,tal->count5[model][i][j][k][l][m],model);
            /* Update size of development data */
            for(z=1;z<currlen;++z)
              *develsize+=(int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            for(z=currlen;z<=(MAXORDER+1);++z)
              *develsize+=(int)pow(ALFSIZE,(MAXORDER+1)-z);
            /* propagate pseudocounts backward */
            tal->count1[model][i]+=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            tal->count2[model][i][j]+=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            tal->count3[model][i][j][k]+=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            tal->count4[model][i][j][k][l]+=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            /* set pseudocounts for current oligo */
            tal->count5[model][i][j][k][l][m]=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            /* propagate pseudocounts forward */
            for(n=0;n<ALFSIZE;++n) {
              tal->count6[model][i][j][k][l][m][n]=MINCOUNTUNIT *
                (int)pow(ALFSIZE,(MAXORDER+1)-(currlen+1));
        }}}}}}}
        break;
      case 6 :
        /* Note that for this particular currlen, with MAXORDER == 5, *
         * pow(ALFSIZE,(MAXORDER+1)-currlen) is simply 1.  I decided  *
         * to retain the code as is should a user wish to modify      *
         * MAXORDER.                                                  */
        for(i=0;i<ALFSIZE;++i) {
        for(j=0;j<ALFSIZE;++j) {
        for(k=0;k<ALFSIZE;++k) {
        for(l=0;l<ALFSIZE;++l) {
        for(m=0;m<ALFSIZE;++m) {
        for(n=0;n<ALFSIZE;++n) {
          if(tal->count6[model][i][j][k][l][m][n] <= EMPTYCOUNT) {
            fprintf(stdout,"History %u%u%u%u%u%u only \
occurs %u times in devel data! (model %i)\n",
              i,j,k,l,m,n,tal->count6[model][i][j][k][l][m][n],model);
            /* Update size of development data */
            for(z=1;z<currlen;++z)
              *develsize+=(int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            for(z=currlen;z<=(MAXORDER+1);++z)
              *develsize+=(int)pow(ALFSIZE,(MAXORDER+1)-z);
            /* propagate pseudocounts backward */
            tal->count1[model][i]+=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            tal->count2[model][i][j]+=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            tal->count3[model][i][j][k]+=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            tal->count4[model][i][j][k][l]+=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            tal->count5[model][i][j][k][l][m]+=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
            /* set pseudocounts for current oligo */
            tal->count6[model][i][j][k][l][m][n]=MINCOUNTUNIT *
              (int)pow(ALFSIZE,(MAXORDER+1)-currlen);
        }}}}}}}
        break;
      default :
        fprintf(stderr,"Err (induce_pseudocounts): Invalid oligo length\n");
        exit(EXIT_FAILURE);
  }}

  return;
} /* end induce_pseudocounts */
#endif

/* Subroutine to calculate unsmoothed maximum likelihood *
 * probabilities for oligomers.  (relative frequency)    *
 * Note: seqct and develsize are not used anymore, but   *
 * their functionality may become desirable again in     *
 * the future.                                           */
void get_rfreqs(relfreqsT *rfreqs,int model,countarrayT *tal,
  int seqct,int develsize)
{
  int i,j,k,l,m,n,z, /* Iterator variables */
      denom;         /* Stores total num   *
                      * oligos per order   */

  /* mononucleotide frequencies           *
   * These statements must always execute */
  for(z=denom=0;z<ALFSIZE;++z)
    denom += tal->count1[model][z];
  for(i=0;i<ALFSIZE;++i) {
    if(denom == 0)
      rfreqs->MLhis0[i]=DEFRFREQ;
    else
      rfreqs->MLhis0[i]=
        ((double)tal->count1[model][i])/denom;
  }

  if (MAXORDER > 0) {
    /* dinucleotide frequencies */
    for(i=0;i<ALFSIZE;++i) {
      for(z=denom=0;z<ALFSIZE;++z)
        denom += tal->count2[model][i][z];
      for(j=0;j<ALFSIZE;++j) {
        if(denom == 0)
          rfreqs->MLhis1[i][j]=DEFRFREQ;
        else
          rfreqs->MLhis1[i][j]=
            ((double)tal->count2[model][i][j])/denom;
      }
    }
  }

  if (MAXORDER > 1) {
    /* trinucleotide frequencies */
    for(i=0;i<ALFSIZE;++i) {
    for(j=0;j<ALFSIZE;++j) {
      for(z=denom=0;z<ALFSIZE;++z)
        denom += tal->count3[model][i][j][z];
      for(k=0;k<ALFSIZE;++k) {
        if(denom == 0)
          rfreqs->MLhis2[i][j][k]=DEFRFREQ;
        else
          rfreqs->MLhis2[i][j][k] =
            ((double)tal->count3[model][i][j][k])/denom;
      }
    }}
  }

  if (MAXORDER > 2) {
    /* tetranucleotide frequencies */
    for(i=0;i<ALFSIZE;++i) {
    for(j=0;j<ALFSIZE;++j) {
    for(k=0;k<ALFSIZE;++k) {
      for(z=denom=0;z<ALFSIZE;++z)
        denom += tal->count4[model][i][j][k][z];
      for(l=0;l<ALFSIZE;++l) {
        if(denom == 0)
          rfreqs->MLhis3[i][j][k][l]=DEFRFREQ;
        else
          rfreqs->MLhis3[i][j][k][l]=
            ((double)tal->count4[model][i][j][k][l])/denom;
      }
    }}}
  }

  if (MAXORDER > 3) {
    /* pentanucleotide frequencies */
    for(i=0;i<ALFSIZE;++i) {
    for(j=0;j<ALFSIZE;++j) {
    for(k=0;k<ALFSIZE;++k) {
    for(l=0;l<ALFSIZE;++l) {
      for(z=denom=0;z<ALFSIZE;++z)
        denom += tal->count5[model][i][j][k][l][z];
      for(m=0;m<ALFSIZE;++m) {
        if(denom == 0)
          rfreqs->MLhis4[i][j][k][l][m]=DEFRFREQ;
        else
          rfreqs->MLhis4[i][j][k][l][m]=
            ((double)tal->count5[model][i][j][k][l][m])/denom;
      }
    }}}}
  }

  if (MAXORDER > 4) {
    /* hexanucleotide frequencies */
    for(i=0;i<ALFSIZE;++i) {
    for(j=0;j<ALFSIZE;++j) {
    for(k=0;k<ALFSIZE;++k) {
    for(l=0;l<ALFSIZE;++l) {
    for(m=0;m<ALFSIZE;++m) {
      for(z=denom=0;z<ALFSIZE;++z)
        denom += tal->count6[model][i][j][k][l][m][z];
      for(n=0;n<ALFSIZE;++n) {
        if(denom == 0)
          rfreqs->MLhis5[i][j][k][l][m][n]=DEFRFREQ;
        else
          rfreqs->MLhis5[i][j][k][l][m][n]=
            ((double)tal->count6[model][i][j][k][l][m][n])/denom;
      }
    }}}}}
  }

  return;
} /* end get_rfreqs */

#ifdef QUANTSPEC
/* To account for local as well as global G+C compositional heterogeneity *
 * in DNA sequences, this function will develop lists of the G+C content  *
 * of 2*FLANK+6 base long windows (100nt) shifted by SHIFTLEN bases.      *
 * GCcounts is STRICTLY a 2x10x10 matrix, intented to index counts of     *
 * test windows exhibiting a G+C compositional percentage [0..100]        */
static void ctGCpercs(FILE *infile,int GCcounts[][10][10])
{
  char *sequence=NULL; /* store sequence to process */
  int seqlength,       /* length of sequence        */
      GCcomp,          /* Percent G+C in window     */
      winstosee,       /* count of windows sequence *
                        * incurs.                   */
      i,j;             /* iterator variables        */

  while((sequence=get_fasta(infile,sequence,&seqlength))!=NULL) {
    /* If the sequence is under 100nts, just compute *
     * overall [GC] and make the best of it.         */
    if(seqlength < (2*FLANK+MAXORDER+1)) {
      /* Calc overall [GC] and cast to an integer *
       * precision to 100th's only.               */
      for(i=0,GCcomp=0;i<seqlength;++i)
        if(sequence[i] == 'c' || sequence[i] == 'C' ||
           sequence[i] == 'g' || sequence[i] == 'G')
          ++GCcomp;
      GCcomp=(int)floor(100*(GCcomp/(float)seqlength));

      if(GCcomp==100)
        ++GCcounts[1][0][0];
      else
        ++GCcounts[0][GCcomp/10][GCcomp%10];
    }
    else {
      /* Determine G+C composition of first 100 bases.                      *
       * This is an IMPORTANT POINT: As I expect 2*FLANK+MAXORDER+1 == 100, *
       * GCcomp is anticipated to be a percentage, as is.                   */
      GCcomp=0;
      for(i=0;i<2*FLANK+MAXORDER+1;++i)
        if(sequence[i] == 'c' || sequence[i] == 'C' ||
           sequence[i] == 'g' || sequence[i] == 'G')
          ++GCcomp;
      if(GCcomp==100)
        ++GCcounts[1][0][0];
      else
        ++GCcounts[0][GCcomp/10][GCcomp%10];

      /* Now, for the rest */
      winstosee = ((seqlength - (2*FLANK+MAXORDER+1)) / SHIFTLEN) + 1;
      for(j=1;j<winstosee;++j) {
        /* out with the old... */
        for(i=(j-1)*SHIFTLEN;i<j*SHIFTLEN;++i)
          if(sequence[i] == 'c' || sequence[i] == 'C' ||
             sequence[i] == 'g' || sequence[i] == 'G')
            --GCcomp;
        /* ...and in with the new */
        for(i=((j-1)*SHIFTLEN)+(2*FLANK+MAXORDER+1);
            i<(j*SHIFTLEN)+(2*FLANK+MAXORDER+1);
            ++i)
          if(sequence[i] == 'c' || sequence[i] == 'C' ||
             sequence[i] == 'g' || sequence[i] == 'G')
            ++GCcomp;

        if(GCcomp==100)
          ++GCcounts[1][0][0];
        else
          ++GCcounts[0][GCcomp/10][GCcomp%10];
      }
    } /* end else */

    free(sequence);
    sequence=NULL;
  } /* end while */

  return;
} /* end ctGCpercs */

/* This function computes the actual quantile boundaries *
 * based on GCcounts, and records these in bounds.       */
static void find_quantiles(int GCcounts[][10][10],
  int codingp,qtboundsT *bounds)
{
  int GClistsize,    /* How many windows assayed */
      unitsperquant, /* # elements per quantile  */
      accum,         /* accumulator variable     */
      flag,          /* bool for triggering      */
      i,j,k;         /* iterator variable        */

  /* How many data points in total? */
  GClistsize=0;
  for(i=0;i<10;++i)
    for(j=0;j<10;++j)
      GClistsize += GCcounts[0][i][j];
  GClistsize += GCcounts[1][0][0];

  /* How large is a quantile? */
  unitsperquant=GClistsize / QTCT;

  /* Determine minimum bound. */
  flag=0;
  for(i=0;i<10 && !flag;++i)
    for(j=0;j<10 && !flag;++j)
      if(GCcounts[0][i][j]>0) {
        flag=1;
        break;
      }
  if(flag) {
    if(codingp == TRUE)
      bounds->codingbound[0]=i*10+j*1;
    else
      bounds->noncodbound[0]=i*10+j*1;
  }
  else {
    fprintf(stderr,"Err (find_quantiles): \
Min bound on G+C percentage >= 100?!!\n");
    exit(EXIT_FAILURE);
  }

  /* Determine the quantiles */
  for(flag=accum=0,k=1;i<10 && !flag;++i) {
    if(j==10) /* reset one's place index */
      j=0;
    for( ;j<10 && !flag;++j) {
      if(accum >= unitsperquant) {
        if(codingp == TRUE)
          bounds->codingbound[k]=i*10+j*1;
        else
          bounds->noncodbound[k]=i*10+j*1;
        if(++k==QTCT)
          flag=1;
        /* Surplus counts rightfully belong *
         * in the next quantile             */
        accum -= unitsperquant;
      }
      else
        accum += GCcounts[0][i][j];
    }
  }
  /* prevent possible infinite loop */
  while(k<QTCT) {
    fprintf(stderr,"Warning (find_quantiles): \
Quantile boundary (and not MAX, either!) occurs on 100%% [G+C]!??\n");
    if(codingp == TRUE)
      bounds->codingbound[k++]=100;
    else
      bounds->noncodbound[k++]=100;
  }

  /* Determine maximum bound */
  if(GCcounts[1][0][0] > 0) {
    if(codingp == TRUE)
      bounds->codingbound[QTCT]=100;
    else
      bounds->noncodbound[QTCT]=100;
  }
  else {
    flag=0;
    for(i=9;i>=0 && !flag;--i)
      for(j=9;j>=0 && !flag;--j)
        if(GCcounts[0][i][j]>0)
          flag=1;
    if(codingp == TRUE)
      bounds->codingbound[QTCT]=i*10+j*1;
    else
      bounds->noncodbound[QTCT]=i*10+j*1;
  }

  /* Some final sanity checks */
  for(k=1;k<=QTCT;++k) {
    if(codingp == TRUE) {
      if(bounds->codingbound[k] < bounds->codingbound[k-1]) {
        fprintf(stderr,"Err (find_quantiles): \
Coding quantile boundaries unordered!\n");
        exit(EXIT_FAILURE);
      }
      else if(bounds->codingbound[k] == bounds->codingbound[k-1])
        /* We consider this scenario tolerable */
        fprintf(stderr,"Warning (find_quantiles): \
Coding quantile boundaries equivalent!\n");
      else
        ;
    }
    else {
      if(bounds->noncodbound[k] < bounds->noncodbound[k-1]) {
        fprintf(stderr,"Err (find_quantiles): \
Noncoding quantile boundaries unordered!\n");
        exit(EXIT_FAILURE);
      }
      else if(bounds->noncodbound[k] == bounds->noncodbound[k-1])
        /* We consider this scenario tolerable */
        fprintf(stderr,"Warning (find_quantiles): \
Noncoding quantile boundaries equivalent!\n");
      else
        ;
    }
  }

  return;
} /* end find_quantiles */

/* Determines quantile boundaries based on overlapping windows */
void calc_quantiles(char *argv[],int algo,qtboundsT *bounds)
{
  FILE *infile=NULL;
  int i,j,k,               /* iterator variables                            */
      num_parts,           /* partitions used to define clouds (DMMMs)      */
      partition_counter,   /* iterate over num_parts                        */
      GCcounts[2][10][10]; /* Stores counts of 2*FLANK+MAXORDER+1 (== 100)  *
                            * windows exhibiting a given G+C compositional  *
                            * percentage.  The first dimension indexes 100s *
                            * the 2nd 10s, and the 3rd 1s.                  */

  if(algo == FIXORDX){
    /* Derive boundaries from coding data */
    for(i=0;i<2;++i)
      for(j=0;j<10;++j)
        for(k=0;k<10;++k)
          GCcounts[i][j][k]=0;
    if((infile=fopen(argv[2],"rt"))==NULL) {
      fprintf(stderr,"Err (calc_quantiles): Can't seem to open %s!\n",argv[2]);
      exit(EXIT_FAILURE);
    }
    ctGCpercs(infile,GCcounts);
    find_quantiles(GCcounts,TRUE,bounds);
    (void)fclose(infile);

    /* Derive boundaries from noncod data */
    for(i=0;i<2;++i)
      for(j=0;j<10;++j)
        for(k=0;k<10;++k)
          GCcounts[i][j][k]=0;
    if((infile=fopen(argv[3],"rt"))==NULL) {
      fprintf(stderr,"Err (calc_quantiles): Can't seem to open %s!\n",argv[3]);
      exit(EXIT_FAILURE);
    }
    ctGCpercs(infile,GCcounts);
    find_quantiles(GCcounts,FALSE,bounds);
    (void)fclose(infile);
  }
  else if(algo == TDDIX || algo == BUDIX || algo == CHI2X) {
    /* Derive boundaries from coding data */
    for(i=0;i<2;++i)
      for(j=0;j<10;++j)
        for(k=0;k<10;++k)
          GCcounts[i][j][k]=0;
    for(i=2;i<5;i+=2) {
      if((infile=fopen(argv[i],"rt"))==NULL) {
        fprintf(stderr,"Err (calc_quantiles): \
Can't seem to open %s!\n",argv[i]);
        exit(EXIT_FAILURE);
      }
      ctGCpercs(infile,GCcounts);
      (void)fclose(infile);
    }
    find_quantiles(GCcounts,TRUE,bounds);

    /* Derive boundaries from noncod data */
    for(i=0;i<2;++i)
      for(j=0;j<10;++j)
        for(k=0;k<10;++k)
          GCcounts[i][j][k]=0;
    for(i=3;i<6;i+=2) {
      if((infile=fopen(argv[i],"rt"))==NULL) {
        fprintf(stderr,"Err (calc_quantiles): \
Can't seem to open %s!\n",argv[i]);
        exit(EXIT_FAILURE);
      }
      ctGCpercs(infile,GCcounts);
      (void)fclose(infile);
    }
    find_quantiles(GCcounts,FALSE,bounds);
  }
  else if(algo == DMMMX) {
    num_parts=atoi(argv[2]);
    if(num_parts<1){
      fprintf(stderr,"Err (calc_quantiles): Invalid partition count!\n");
      exit(EXIT_FAILURE);
    }

    /* Derive boundaries from coding data */
    for(i=0;i<2;++i)
      for(j=0;j<10;++j)
        for(k=0;k<10;++k)
          GCcounts[i][j][k]=0;
    for(partition_counter=0;partition_counter<num_parts;++partition_counter) {
      if((infile=fopen(argv[3+partition_counter],"rt"))==NULL) {
        fprintf(stderr,"Err (calc_quantiles): Can't seem to open %s!\n",
          argv[3+partition_counter]);
        exit(EXIT_FAILURE);
      }
      ctGCpercs(infile,GCcounts);
      (void)fclose(infile);
    }
    find_quantiles(GCcounts,TRUE,bounds);

    /* Derive boundaries from noncod data */
    for(i=0;i<2;++i)
      for(j=0;j<10;++j)
        for(k=0;k<10;++k)
          GCcounts[i][j][k]=0;
    for(partition_counter=0;partition_counter<num_parts;++partition_counter) {
      if((infile=fopen(argv[3+num_parts+partition_counter],"rt"))==NULL) {
        fprintf(stderr,"Err (calc_quantiles): Can't seem to open %s!\n",
          argv[3+num_parts+partition_counter]);
        exit(EXIT_FAILURE);
      }
      ctGCpercs(infile,GCcounts);
      (void)fclose(infile);
    }
    find_quantiles(GCcounts,FALSE,bounds);
  }
  else {
    fprintf(stderr,"Err (calc_quantiles): Improper algorithm requested!\n");
    exit(EXIT_FAILURE);
  }

  return;
} /* end calc_quantiles */
#endif

/* This function oversees the building of a trained immprobT file     *
 * using the deleted interpolated Markov model approach of Jelinek or *
 * the chi**2 interpolated Markov model approach of Salzberg.         *
 * If THREEPERIODIC is set, the function oversees the building of a   *
 * trained immprobT file while taking three-periodicity into          *
 * account. (Note: Training the coding models assumes that all input  *
 * data is in phase 0!) Else, the function trains a model in terms of *
 * a binary classifier.                                               */
int imm_probmaker(int argc,char *argv[],int algo)
{
  FILE *outfile=NULL;        /* Will store binary immprobT object    */
  char *sequence=NULL,       /* For storing current sequence         */
    outputname[MAXFILENAME]; /* Stores final name of output file     */
  int seqlength,             /* Records length of current sequence   */
      model;                 /* Model iterator variable              */
  #ifdef QUANTSPEC
  int quantind;              /* iterator variable                    */
  countarrayT talliesD[QTCT],/* Record oligomer frequencies, D data  */
              talliesC[QTCT];/* Record oligomer frequencies, C data  */
  relfreqsT relfreqs[QTCT];  /* Store unsmooted ML oligo probs       */
  immprobT smoothprob[QTCT]; /* We'll build and write to outfile     */
  qtboundsT qt_bounds;       /* Store quantile boundaries            */
  #else
  countarrayT talliesD,      /* Record oligomer frequencies, D data  */
              talliesC;      /* Record oligomer frequencies, C data  */
  relfreqsT relfreqs;        /* Store unsmooted ML oligo probs       */
  immprobT smoothprob;       /* We'll build and write to outfile     */
  #endif
  #ifdef THREEPERIODIC
  int currentfile;           /* Iterator based on input file at hand */
  FILE *infile=NULL;         /* For processing all Fasta input files */
  #ifdef EMPTYEQUIV
  int stub;                  /* Pass as mandatory (unused) parameter *
                              * to the induce_pseudocounts fxn.      */
  #endif
  #else
  FILE *infileD=NULL,        /* Fasta file with development data (D) */
       *infileC=NULL;        /* Fasta file with held-out data    (C) */
  int devellen,              /* Records total length of devel data   */
      numseqs;               /* Records total number input sequences */
  #endif

  /* Verify the proper training algorithm was specified */
  if(algo != TDDIX && algo != BUDIX && algo != CHI2X) {
    fprintf(stderr,"Err (imm_probmaker): Improper algorithm requested!\n");
    exit(EXIT_FAILURE);
  }

  #ifdef QUANTSPEC
  calc_quantiles(argv,algo,&qt_bounds);
  for(quantind=0;quantind<QTCT;++quantind) {
    /* Initialize countarrays */
    for(model=0;model<NUMMODELS;++model) {
      init_counts(&talliesD[quantind],model);
      init_counts(&talliesC[quantind],model);
    }
  }
  #else
  /* Initialize countarrays */
  for(model=0;model<NUMMODELS;++model) {
    init_counts(&talliesD,model);
    init_counts(&talliesC,model);
  }
  #endif

  #ifdef THREEPERIODIC
  /* Do the counting.                                  *
   * currentfile 1 = devel coding                      *
   *             2 = devel noncoding                   *
   *             3 = heldout coding                    *
   *             4 = heldout noncoding                 *
   * Note: Training the coding models assumes that all *
   * input data is in phase 0! (see immpractical.h)    */
  for(currentfile=1;currentfile<=4;++currentfile) {
    /* Derive counts from data */
    if((infile=fopen(argv[currentfile+1],"rt"))==NULL) {
      fprintf(stderr,"Err (imm_probmaker): Can't seem to open %s!\n",
        argv[currentfile+1]);
      exit(EXIT_FAILURE);
    }
    /* get_fasta defined in sequence_parse.c */
    while ((sequence=get_fasta(infile,sequence,&seqlength))!=NULL) {
      if(seqlength >= (MAXORDER+1)) {
        switch (currentfile) {
          case 1 :
            #ifdef QUANTSPEC
            /* Quantile-specific behavior demands three-periodic modeling *
             * of both strands---locked in a NUMMODELS == 7 scheme.       */
            get_counts3percodQT(sequence,talliesD,qt_bounds);
            #else
            #ifdef SHADOWSTRAND
            get_counts3percod(sequence,&talliesD,TRUE);
            #else
            get_counts3percod(sequence,&talliesD,FALSE);
            #endif
            #endif
            break;
          case 2 :
            #ifdef QUANTSPEC
            /* As NUMMODELS == 7, no need to specify NONCODING index */
            get_countsQT(sequence,talliesD,qt_bounds);
            #else
            get_counts(sequence,&talliesD,NONCODING);
            #endif
            break;
          case 3 :
            if(algo == TDDIX || algo == BUDIX) {
              #ifdef QUANTSPEC
              get_counts3percodQT(sequence,talliesC,qt_bounds);
              #else
              #ifdef SHADOWSTRAND
              get_counts3percod(sequence,&talliesC,TRUE);
              #else
              get_counts3percod(sequence,&talliesC,FALSE);
              #endif
              #endif
            }
            else { /* (algo == CHI2X) */
              /* There is no concept of held-out data for chi**2 training, *
               * so we pool the held-out counts with the development ones  */
              #ifdef QUANTSPEC
              get_counts3percodQT(sequence,talliesD,qt_bounds);
              #else
              #ifdef SHADOWSTRAND
              get_counts3percod(sequence,&talliesD,TRUE);
              #else
              get_counts3percod(sequence,&talliesD,FALSE);
              #endif
              #endif
            }
            break;
          case 4 :
            if(algo == TDDIX || algo == BUDIX)
              #ifdef QUANTSPEC
              get_countsQT(sequence,talliesC,qt_bounds);
              #else
              get_counts(sequence,&talliesC,NONCODING);
              #endif
            else /* (algo == CHI2X) */
              #ifdef QUANTSPEC
              get_countsQT(sequence,talliesD,qt_bounds);
              #else
              get_counts(sequence,&talliesD,NONCODING);
              #endif
            break;
          default :
            fprintf(stderr,"Err (imm_probmaker): \
You seem to have passed an invalid training file?\n");
            (void)fclose(infile);
            return(FALSE);
        }
        free(sequence);
        sequence=NULL;
      }
      else {
        fprintf(stderr,"Err (imm_probmaker): First run %s on your input!\n",
          DATACLEANPL);
        free(sequence);
        (void)fclose(infile);
        return(FALSE);
      }
    } /* end while */
    (void)fclose(infile);
  } /* end currentfile for */

  /* Smooth each model in succession */
  for(model=0;model<NUMMODELS;++model) {
    #ifdef QUANTSPEC
    for(quantind=0;quantind<QTCT;++quantind) {
      #ifdef EMPTYEQUIV
      /* Make sure all possible histories are histories of interest! */
      induce_pseudocounts(&talliesD[quantind],model,&stub);
      #endif

      /* Record oligomer relative frequencies (ML estimates) from D data set *
       * The last two parameters for this function are irrelevant in a three *
       * periodic context, hence the zeros.                                  */
      get_rfreqs(&relfreqs[quantind],model,&talliesD[quantind],0,0);

      /* Parameterize the model by optimizing interpolation *
       * parameters in light of held-out data set.          */
      if(algo == TDDIX || algo == BUDIX)
        DIMM_final_probs(&smoothprob[quantind],&relfreqs[quantind],
          &talliesD[quantind],&talliesC[quantind],model,algo);
      else /* (algo == CHI2X) */
        CHISQUARE_final_probs(&smoothprob[quantind],&relfreqs[quantind],
          &talliesD[quantind],model);
    }

    #else

    #ifdef EMPTYEQUIV
    /* Make sure all possible histories are histories of interest! */
    induce_pseudocounts(&talliesD,model,&stub);
    #endif

    /* Record oligomer relative frequencies (ML estimates) from D data set *
     * The last two parameters for this function are irrelevant in a three *
     * periodic context, hence the zeros.                                  */
    get_rfreqs(&relfreqs,model,&talliesD,0,0);

    /* Parameterize the model by optimizing interpolation *
     * parameters in light of held-out data set.          */
    if(algo == TDDIX || algo == BUDIX)
      DIMM_final_probs(&smoothprob,&relfreqs,&talliesD,&talliesC,model,algo);
    else /* (algo == CHI2X) */
      CHISQUARE_final_probs(&smoothprob,&relfreqs,&talliesD,model);

    #endif
  } /* end for */

  #else

  /* Develop each model in succession */
  for(model=0;model<NUMMODELS;++model) {
    /* Populate the count arrays using the development data set. */
    if((infileD=fopen(argv[model+2],"rt"))==NULL) {
      fprintf(stderr,
        "Err (imm_probmaker): Can't seem to open %s!\n",argv[model+2]);
      exit(EXIT_FAILURE);
    }
    devellen=numseqs=0;
    while ((sequence=get_fasta(infileD,sequence,&seqlength))!=NULL) {
      if(seqlength >= (MAXORDER+1)) {
        get_counts(sequence,&talliesD,model);
        free(sequence);
        sequence=NULL;
        devellen += seqlength;
        ++numseqs;
      }
      else {
        fprintf(stderr,"Err (imm_probmaker): First run %s on your input!\n",
          DATACLEANPL);
        free(sequence);
        (void)fclose(infileD);
        return(FALSE);
      }
    } /* end while */
    (void)fclose(infileD);

    #ifdef EMPTYEQUIV
    /* Make sure all possible histories are histories of interest! */
    induce_pseudocounts(&talliesD,model,&devellen);
    #endif

    /* Record oligomer relative frequencies (ML estimates) from D data set */
    get_rfreqs(&relfreqs,model,&talliesD,numseqs,devellen);

    /* Populate the count array using the held-out data set. */
    if((infileC=fopen(argv[model+4],"rt"))==NULL) {
      fprintf(stderr,"Err (imm_probmaker): Can't seem to open %s!\n",
        argv[model+NUMMODELS+2]);
      exit(EXIT_FAILURE);
    }
    while ((sequence=get_fasta(infileC,sequence,&seqlength))!=NULL) {
      if(seqlength >= (MAXORDER+1)) {
        get_counts(sequence,&talliesC,model);
        free(sequence);
        sequence=NULL;
      }
      else {
        fprintf(stderr,"Err (imm_probmaker): First run %s on your input!\n",
          DATACLEANPL);
        free(sequence);
        (void)fclose(infileC);
        return(FALSE);
      }
    } /* end while */
    (void)fclose(infileC);

    /* Parameterize the model by optimizing interpolation *
     * parameters in light of held-out data set.          */
    if(algo == TDDIX || algo == BUDIX)
      DIMM_final_probs(&smoothprob,&relfreqs,&talliesD,&talliesC,model,algo);
    else /* (algo == CHI2X) */
      CHISQUARE_final_probs(&smoothprob,&relfreqs,&talliesD,model);
  } /* end for */
  #endif

  /* Set name for output file. */
  if ((int)strlen(argv[6]) > MAXFILENAME) {
    fprintf(stderr,
      "Err (imm_probmaker): Filename exceeds %i characters!\n",MAXFILENAME);
    return(FALSE);
  }
  else {
    (void)strcpy(outputname,argv[6]);
    (void)strcat(outputname,".");
    (void)strcat(outputname,algotags[atoi(argv[1])]);
  }

  /* Write our newly trained model to the output file */
  if((outfile=fopen(outputname,"wb"))==NULL) {
    fprintf(stderr,"Err (imm_probmaker): Can't seem to open %s!\n",outputname);
    exit(EXIT_FAILURE);
  }
  #ifdef QUANTSPEC
  else if(fwrite(smoothprob,sizeof(immprobT),QTCT,outfile) != QTCT) {
    fprintf(stderr,
      "Err (imm_probmaker): Error writing model to binary output file!\n");
    exit(EXIT_FAILURE);
  }
  else if(fwrite(&qt_bounds,sizeof(qtboundsT),1,outfile) != 1) {
    fprintf(stderr,
      "Err (fo_probmaker): \
Error writing quantile boundaries to binary output file!\n");
    exit(EXIT_FAILURE);
  }
  #else
  else if(fwrite(&smoothprob,sizeof(immprobT),1,outfile) != 1) {
    fprintf(stderr,
      "Err (imm_probmaker): Error writing model to binary output file!\n");
    exit(EXIT_FAILURE);
  }
  #endif
  else
    (void)fclose(outfile);

  return(TRUE);
} /* end imm_probmaker */

#ifdef QUANTSPEC
/* Function to return log-likelihood of a test input sequence, *
 * which must be null-terminated, under the model specified    *
 * using the model argument. It returns the probability of the *
 * input sequence given its first base occurs in the phase     *
 * specified by model, which parameters selected on a quantile *
 * specific basis. If testwindowgc is TRUE, this function will *
 * modulate between quantile-specific parameters as a function *
 * of window-based G+C composition.                            */
long double IMMprobQT(char *input,immprobT *probs,qtboundsT *bounds,
  int model,int testwindowgc)
{
  int length,              /* Length of input string         */
      *nt=NULL,            /* int translation of input       */
      i,j,                 /* Iterator variable              */
      GCcomp=0;            /* G+C content                    */
  long double logprob;     /* end result to return           */
  immprobT *probsptr;      /* Which element of probs to use? */

  /* Verify that model is given a reasonable value */
  if (model < 0 || model > NONCODING) {
    fprintf(stderr,"Err (IMMprobQT): Invalid model value (%i)\n",model);
    exit(EXIT_FAILURE);
  }

  /* Verify that we really have a sequence to test. If input is set to NULL, *
   * this conditional will short-circuit and strlen won't get called on it.  */
  if((input==NULL) || ((length=(int)strlen(input))==0)) {
    if (input==NULL)
      fprintf(stderr,"Err (IMMprobQT): No real sequence passed!\n");
    else
      fprintf(stderr,"Err (IMMprobQT): Passed a zero-length test sequence!\n");
    exit(EXIT_FAILURE);
  }
  else if(length < (MAXORDER+1)) {
    fprintf(stderr,"Err (IMMprobQT): Can't assess sequences < %i bases!\n",
      MAXORDER+1);
    exit(EXIT_FAILURE);
  }
  else if((nt=(int*)malloc(sizeof(int)*length))==NULL) {
    fprintf(stderr,"Err (IMMprobQT): Out of memory!\n");
    exit(EXIT_FAILURE);
  }
  else /* Copy input translation into nt */
    for(i=0;i<length;++i)
      nt[i]=trans(input[i]);

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

    /* Select the best set of quantile-specific parameters */
    for(i=1;i<QTCT;++i) {
      if(model == NONCODING) {
        if(GCcomp < bounds->noncodbound[i])
          break;
      }
      else
        if(GCcomp < bounds->codingbound[i])
          break;
    }
    probsptr = &probs[i-1];
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

    /* Select the best set of quantile-specific   *
     * parameters for the first window. Note that *
     * GCcompint is already primed on the first   *
     * window at this point.                      */
    for(i=1;i<QTCT;++i) {
      if(model == NONCODING) {
        if(GCcomp < bounds->noncodbound[i])
          break;
      }
      else
        if(GCcomp < bounds->codingbound[i])
          break;
    }
    probsptr = &probs[i-1];
  }

  /* Compute likelihood of the first few 5'-terminal nts  *
   * This applies for either homogeneous or inhomogeneous *
   * Markov models.                                       */
  for(i=0,logprob=0.0;i<MAXORDER;++i) {
    switch(i) {
      case 0 : /* mononucleotide */
        logprob += log(probsptr->prob1[model][nt[0]]);
        break;
      case 1 : /* dinucleotide */
        logprob += log(probsptr->prob2[model][nt[0]]
                                             [nt[1]]);
        break;
      case 2 : /* trinucleotide */
        logprob += log(probsptr->prob3[model][nt[0]]
                                             [nt[1]]
                                             [nt[2]]);
        break;
      case 3 : /* tetranucleotide */
        logprob += log(probsptr->prob4[model][nt[0]]
                                             [nt[1]]
                                             [nt[2]]
                                             [nt[3]]);
        break;
      case 4 : /* pentanucleotide */
        logprob += log(probsptr->prob5[model][nt[0]]
                                             [nt[1]]
                                             [nt[2]]
                                             [nt[3]]
                                             [nt[4]]);
        break;
      default :
        fprintf(stderr,"Err (IMMprobQT): Bad short case passed!\n");
        exit(EXIT_FAILURE);
    } /* end switch */
  }

  if(model == NONCODING) { /* homogeneous Markov model */
    if (testwindowgc == FALSE) {
      /* compute likelihood */
      for(i=MAXORDER;i<length;++i)
        logprob += log(probsptr->prob6[model]
          [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
    }
    else { /* window-based parameter modulation */
      for(i=MAXORDER;i<FLANK+MAXORDER+1;++i)
        logprob += log(probsptr->prob6[model]
          [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);

      /* modulate parameters until last window.  If the index seems  *
       * funny, it's because I'm counting "on the back of the beat". */
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
        logprob += log(probsptr->prob6[model]
          [nt[i+FLANK+MAXORDER-4]]
          [nt[i+FLANK+MAXORDER-3]]
          [nt[i+FLANK+MAXORDER-2]]
          [nt[i+FLANK+MAXORDER-1]]
          [nt[i+FLANK+MAXORDER]]
          [nt[i+FLANK+MAXORDER+1]]);
      }
      /* assay remaining 3' hexamers */
      for(i=length-FLANK;i<length;++i)
        logprob += log(probsptr->prob6[model]
          [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
    }
  }
  else { /* inhomogeneous Markov model */
    if (testwindowgc == FALSE) {
      /* compute likelihood */
      for(i=MAXORDER;i<length;++i) {
        if (model < 3) /* Forward strand */
          logprob += log(probsptr->prob6[(i-MAXORDER+model)%3]
            [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
        else /* Reverse strand */
          logprob += log(probsptr->prob6[((i-MAXORDER+model)%3)+3]
            [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
      }
    }
    else { /* window-based parameter modulation */
      for(i=MAXORDER;i<FLANK+MAXORDER+1;++i) {
        if (model < 3) /* Forward strand */
          logprob += log(probsptr->prob6[(i-MAXORDER+model)%3]
            [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
        else /* Reverse strand */
          logprob += log(probsptr->prob6[((i-MAXORDER+model)%3)+3]
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
          logprob += log(probsptr->prob6[(i+FLANK+1+model)%3]
            [nt[i+FLANK+MAXORDER-4]]
            [nt[i+FLANK+MAXORDER-3]]
            [nt[i+FLANK+MAXORDER-2]]
            [nt[i+FLANK+MAXORDER-1]]
            [nt[i+FLANK+MAXORDER]]
            [nt[i+FLANK+MAXORDER+1]]);
        else /* Reverse strand */
          logprob += log(probsptr->prob6[((i+FLANK+1+model)%3)+3]
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
          logprob += log(probsptr->prob6[(i-MAXORDER+model)%3]
            [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
        else /* Reverse strand */
          logprob += log(probsptr->prob6[((i-MAXORDER+model)%3)+3]
            [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
      }
    }
  }

  /* Back to the heap */
  free(nt);

  return(logprob);
} /* end IMMprobQT */
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
long double IMMprob(char *input,immprobT *probs,int model,
  int threeper,int shadow)
{
  int length,     /* Length of input string    */
      *nt=NULL,   /* int translation of input  */
      i;          /* Iterator variable         */
  long double logprob; /* end result to return */

  /* Verify that model is given a reasonable value */
  if (shadow == FALSE) {
    if (threeper == FALSE) {
      if (model < 0 || model > NONCODING) {
        fprintf(stderr,"Err (IMMprob): Invalid model value (%i)\n",model);
        exit(EXIT_FAILURE);
      }
    }
    else {
      if ((model < 0 || model > 2) && model != NONCODING) {
        fprintf(stderr,"Err (IMMprob): Invalid model value (%i)\n",model);
        exit(EXIT_FAILURE);
      }
    }
  }
  else {
    if (model < 3 || model > 5) {
      fprintf(stderr,"Err (IMMprob): Invalid model value (%i)\n",model);
      exit(EXIT_FAILURE);
    }
  }

  /* Verify that we really have a sequence to test. If input is set to NULL, *
   * this conditional will short-circuit and strlen won't get called on it.  */
  if((input==NULL) || ((length=(int)strlen(input))==0)) {
    if (input==NULL)
      fprintf(stderr,"Err (IMMprob): No real sequence passed!\n");
    else
      fprintf(stderr,"Err (IMMprob): Passed a 0nt test sequence!\n");
    exit(EXIT_FAILURE);
  }
  else if((nt=(int*)malloc(sizeof(int)*length))==NULL) {
    fprintf(stderr,"Err (IMMprob): Out of memory!\n");
    exit(EXIT_FAILURE);
  }
  else /* Copy input translation into nt */
    for(i=0;i<length;++i)
      nt[i]=trans(input[i]);

  /* compute log-likelihood of the input string */
  for(i=0,logprob=0.0;i<MAXORDER;++i) {
    switch(i) {
      case 0 : /* mononucleotide */
        logprob += log(probs->prob1[model][nt[0]]);
        break;
      case 1 : /* dinucleotide */
        logprob += log(probs->prob2[model][nt[0]]
                                          [nt[1]]);
        break;
      case 2 : /* trinucleotide */
        logprob += log(probs->prob3[model][nt[0]]
                                          [nt[1]]
                                          [nt[2]]);
        break;
      case 3 : /* tetranucleotide */
        logprob += log(probs->prob4[model][nt[0]]
                                          [nt[1]]
                                          [nt[2]]
                                          [nt[3]]);
        break;
      case 4 : /* pentanucleotide */
        logprob += log(probs->prob5[model][nt[0]]
                                          [nt[1]]
                                          [nt[2]]
                                          [nt[3]]
                                          [nt[4]]);
        break;
      default :
        fprintf(stderr,"Err (IMMprob): Bad short case passed!\n");
        exit(EXIT_FAILURE);
    } /* end switch */
  } /* end for */

  /* This loop processes positions such that full *
   * MAXORDER length histories are available.     */
  if (threeper != FALSE) { /* Compute 3-periodic result */
    for(i=MAXORDER;i<length;++i) {
      switch(MAXORDER) {
        case(0) :
          if (shadow == FALSE) /* Forward strand */
            logprob += log(probs->prob1[(i-MAXORDER+model)%3]
              [nt[i]]);
          else /* Reverse strand */
            logprob += log(probs->prob1[((i-MAXORDER+model)%3)+3]
              [nt[i]]);
          break;
        case(1) :
          if (shadow == FALSE) /* Forward strand */
            logprob += log(probs->prob2[(i-MAXORDER+model)%3]
              [nt[i-1]][nt[i]]);
          else /* Reverse strand */
            logprob += log(probs->prob2[((i-MAXORDER+model)%3)+3]
              [nt[i-1]][nt[i]]);
          break;
        case(2) :
          if (shadow == FALSE) /* Forward strand */
            logprob += log(probs->prob3[(i-MAXORDER+model)%3]
              [nt[i-2]][nt[i-1]][nt[i]]);
          else /* Reverse strand */
            logprob += log(probs->prob3[((i-MAXORDER+model)%3)+3]
              [nt[i-2]][nt[i-1]][nt[i]]);
          break;
        case(3) :
          if (shadow == FALSE) /* Forward strand */
            logprob += log(probs->prob4[(i-MAXORDER+model)%3]
              [nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
          else /* Reverse strand */
            logprob += log(probs->prob4[((i-MAXORDER+model)%3)+3]
              [nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
          break;
        case(4) :
          if (shadow == FALSE) /* Forward strand */
            logprob += log(probs->prob5[(i-MAXORDER+model)%3]
              [nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
          else /* Reverse strand */
            logprob += log(probs->prob5[((i-MAXORDER+model)%3)+3]
              [nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
          break;
        case(5) :
          if (shadow == FALSE) /* Forward strand */
            logprob += log(probs->prob6[(i-MAXORDER+model)%3]
              [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
          else /* Reverse strand */
            logprob += log(probs->prob6[((i-MAXORDER+model)%3)+3]
              [nt[i-5]][nt[i-4]][nt[i-3]][nt[i-2]][nt[i-1]][nt[i]]);
          break;
        default:
          fprintf(stderr,"Err (IMMprob): Invalid MAXORDER encountered!\n");
          exit(EXIT_FAILURE);
      }
    }
  }
  else { /* Compute non-periodic result */
    for(i=MAXORDER;i<length;++i) {
      switch(MAXORDER) {
        case(0) :
          logprob += log(probs->prob1[model][nt[i]]);
          break;
        case(1) :
          logprob += log(probs->prob2[model][nt[i-1]]
                                            [nt[i]]);
          break;
        case(2) :
          logprob += log(probs->prob3[model][nt[i-2]]
                                            [nt[i-1]]
                                            [nt[i]]);
          break;
        case(3) :
          logprob += log(probs->prob4[model][nt[i-3]]
                                            [nt[i-2]]
                                            [nt[i-1]]
                                            [nt[i]]);
          break;
        case(4) :
          logprob += log(probs->prob5[model][nt[i-4]]
                                            [nt[i-3]]
                                            [nt[i-2]]
                                            [nt[i-1]]
                                            [nt[i]]);
          break;
        case(5) :
          logprob += log(probs->prob6[model][nt[i-5]]
                                            [nt[i-4]]
                                            [nt[i-3]]
                                            [nt[i-2]]
                                            [nt[i-1]]
                                            [nt[i]]);
          break;
        default:
          fprintf(stderr,"Err (IMMprob): Invalid MAXORDER encountered!\n");
          exit(EXIT_FAILURE);
      }
    }
  }

  /* back to the heap */
  free(nt);

  return(logprob);
} /* end IMMprob */
#endif

#ifdef QUANTSPEC
/* Function to import a trained immprobT object array and a     *
 * qtboundsT object (stored in a binary file) into main memory. */
void import_imm_probsQT(char *filename,immprobT *probs,qtboundsT *bounds)
{
  FILE *fptr=NULL; /* For connecting to binary input file */

  /* Open parameter file */
  if((fptr=fopen(filename,"rb"))==NULL) {
    fprintf(stderr,"Err (import_imm_probsQT): Can't open %s for \
binary reading!\n",filename);
    exit(EXIT_FAILURE);
  }
  /* Read in the immprobT "objects" */
  else if(fread(probs,sizeof(immprobT),QTCT,fptr) != QTCT) {
    fprintf(stderr,"Err (import_imm_probsQT): Error reading from %s!\n",
      filename);
    exit(EXIT_FAILURE);
  }
  /* Read in the qtboundsT "object" */
  else if(fread(bounds,sizeof(qtboundsT),1,fptr) != 1) {
    fprintf(stderr,"Err (import_imm_probsQT): Error reading from %s!\n",
      filename);
    exit(EXIT_FAILURE);
  }
  else
    (void)fclose(fptr);

  return;
} /* end import_imm_probsQT */
#else
/* Function to import a trained immprobT    *
 * object (a binary file) into main memory. */
void import_imm_probs(char *filename,immprobT *probs)
{
  FILE *fptr=NULL; /* For connecting to binary input file */

  /* Open parameter file */
  if((fptr=fopen(filename,"rb"))==NULL) {
    fprintf(stderr,"Err (import_imm_probs): Can't open %s for \
binary reading!\n",filename);
    exit(EXIT_FAILURE);
  }
  /* Read in the immprobT "object" */
  else if(fread(probs,sizeof(immprobT),1,fptr) != 1) {
    fprintf(stderr,"Err (import_imm_probs): Error reading from %s!\n",
      filename);
    exit(EXIT_FAILURE);
  }
  else
    (void)fclose(fptr);

  return;
} /* end import_imm_probs */
#endif
