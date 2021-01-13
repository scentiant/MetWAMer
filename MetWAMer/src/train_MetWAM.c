/* train_MetWAM.c
 * Michael Sparks (mespar1@gmail.com)
 * Last modified : 20 July 2013
 *
 * Main implementation file for the train_MetWAM utility.
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
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "errors.h"
#include "MetWAM_utils.h"
#include "probdef.h"
#include "sequence_parse.h"

#define USAGE "\a\nUsage: %s \\\n\
          [-o first_codon_position] \\\n\
          [-e last_codon_offset (distance from sequence stop)] \\\n\
          [-l Met_site_lower_bound] \\\n\
          -h hypothesis (0 or 1) \\\n\
          -f codseqs.fas \\\n\
          -r outfile\n\
\n\
  Where -o arg is the (optional, one-based) position of the first\n\
             base in the start codon in each training sequence\n\
             from codseqs.fas [default=1].  More specifically,\n\
             scanning for Mets will not occur at positions .LT.\n\
             this value.\n\
        -e arg is the (optional, non-negative valued) offset of\n\
             the last base of the last codon in each training\n\
             sequence from codseqs.fas, relative to the last\n\
             base of the sequence [default=0].  More specifically,\n\
             scanning for Mets will not occur at positions .GT. the\n\
             downstream limit this value implies in each training\n\
             instance.\n\
        -l arg is the (optional, one-based) lower bound used for\n\
             testing whether a translation initiation signal in the\n\
             sequence is valid to score [default=1].\n\
        -h arg is the hypothesis, either 0 (for false sites) or\n\
             1 (for true sites) that the resultant WAM is expected\n\
             to characterize.\n\
        -f arg is a set of training sequences in Fasta format.\n\
        -r arg is the desired output file name.\n\
\n"

/* Function prototypes */
static void makeMetWAM(ATGparm *locATGparms,int **seq_matrix,int sitect);

/* Main Application */
int main(int argc,char *argv[]) {
  ATGparm locATGparms;       /* Stores ATGparm object         */
  char
    *charcodseq=NULL,        /* coding seq to be processed    */   
    codseqid[MAXDESCLINE+1], /* description of codseq, +1 for *
                              * terminating '\0'.             */
    *outhandle=NULL,         /* handle for outputname string  */
    outputname[MAXFILENAME]; /* Stores output file name       */
  FILE *seqfile=NULL,        /* fasta input file              */
       *outfile=NULL;        /* write locATGparms to it       */
  int option,                /* option parsed by getopt       */
      hypothesis=-1,         /* 0 or 1, hypo of TISs for WAM  */
      *codseq=NULL,          /* int translation of sequence   */
      seqlength=0,           /* stores length of codseq array */
      **seq_matrix=NULL,     /* Stores sequence data          */
      sitect=0,              /* counts number of sites        */
      firstcdsbase=1,        /* first cds base position       */
      lastcdsoff=1,          /* last cds base offset position */
      metsiglowbound=1,      /* lower bound on valid Mets     */
      i,j,k;                 /* Iterator variables            */

  /* Parse command line */
  while((option=getopt(argc,argv,":o:e:l:h:f:r:"))!=-1) {
    switch(option) {
      /* Note that, for options o and l, it is not required that l.LT.o */
      case 'o' :
        if((firstcdsbase=atoi(optarg))<1) {
          NONFATALERROR("Warn (main): \
First codon position being reset to 1\n")
          firstcdsbase=1;
        }
        break;
      case 'e' :
        if((lastcdsoff=atoi(optarg))<0) {
          NONFATALERROR("Warn (main): \
Last codon offset being reset to 0\n")
          lastcdsoff=0;
        }
        break;
      case 'l' :
        if((metsiglowbound=atoi(optarg))<1) {
          NONFATALERROR("Warn (main): \
metsiglowbound being reset to 1\n")
          metsiglowbound=1;
        }
        break;
      case 'h' :
        if((hypothesis=atoi(optarg))<0||hypothesis>1) {
          FATALERROR("Warn (main): \
Non-sensible value for hypothesis specified!\n")
        }
        break;
      case 'f' :
        /* Connect to training sequence input stream */
        if((seqfile=fopen(optarg,"rt"))==NULL) {
          (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Can't seem to open %s!\n",optarg);
          FATALERROR(errbuff)
        }
        break;
      case 'r' :
        /* We will create and write to this file later */
        outhandle=optarg;
        (void)strcpy(outputname,outhandle);
        (void)strcat(outputname,".MetWAM");
        break;
      default :
        (void)snprintf(errbuff,MAXERRLEN,USAGE,argv[0]);
        FATALERROR(errbuff)
    }
  }
  if(seqfile==NULL||outhandle==NULL||hypothesis==-1) {
    (void)snprintf(errbuff,MAXERRLEN,USAGE,argv[0]);
    FATALERROR(errbuff)
  }

  /* Process the input sequences */
  while((charcodseq=get_fasta(seqfile,charcodseq,&seqlength,codseqid))!=NULL) {
    if(seqlength==0)
      if(charcodseq==NULL)
        FATALERROR("Err (main): No real sequence passed!\n")
      else
        FATALERROR("Err (main): Passed a 0nt genomic codseq!\n")
    else if((codseq=(int*)malloc(sizeof(int)*seqlength))==NULL)
      FATALERROR("Err (main): Out of memory!\n")
    else { /* Copy charcodseq translation into codseq */
      for(i=0;i<seqlength;++i)
        codseq[i]=trans(charcodseq[i]);
      free(charcodseq);
    }

    /* sanity checks */
    if(seqlength<metsiglowbound) {
      (void)snprintf(errbuff,MAXERRLEN,"Warn (main): \
ATG signal's lower bound (%i) .GT. sequence length (%i)!??\n\
  (Skipping sequence \"%s\")\n",metsiglowbound,seqlength,codseqid);
      NONFATALERROR(errbuff)
      free(codseq);
      continue;
    }
    else if(seqlength<firstcdsbase+2+(DOWNSTREXTENT)) {
      (void)snprintf(errbuff,MAXERRLEN,"Warn (main): \
First codon can't be scored in sequence \"%s\"!?? \
(Skipping it)\n",codseqid);
      NONFATALERROR(errbuff)
      free(codseq);
      continue;
    }
    else if(seqlength-2-lastcdsoff<firstcdsbase) {
      (void)snprintf(errbuff,MAXERRLEN,"Warn (main): \
Stop codon position (%i) .LT. first codon position (%i)!??\n\
  (Skipping sequence \"%s\")\n",
        seqlength-2-lastcdsoff,firstcdsbase,codseqid);
      NONFATALERROR(errbuff)
      free(codseq);
      continue;
    }
    else if((seqlength-lastcdsoff-firstcdsbase+1)%3!=0) {
      (void)snprintf(errbuff,MAXERRLEN,"Warn (main): \
CDS length (%i) is not a multiple of three!??\n\
  (Skipping sequence \"%s\")\n",
        seqlength-lastcdsoff-firstcdsbase+1,codseqid);
      NONFATALERROR(errbuff)
      free(codseq);
      continue;
    }

    /* Find first \emph{valid} Met in sequence, which we presume to be true */
    for(i=firstcdsbase-1;i<=seqlength-lastcdsoff-3;i+=3)
      if(IS_MET_P(codseq,i)==TRUE&&
         MET_VALID_TO_SCORE_P(i,metsiglowbound-1,seqlength-1)==TRUE) {
        if(hypothesis==0) /* Building WAM for false TISs */
          break;
        if((seq_matrix=
              (int**)realloc(seq_matrix,++sitect*sizeof(int*)))==NULL||
           (seq_matrix[sitect-1]=
              (int*)malloc((STRINGSIZE)*sizeof(int)))==NULL)
          FATALERROR("Err (main): Insufficient Memory!\n")
        else
          for(j=i-(UPSTREXTENT),k=0;j<=i+2+(DOWNSTREXTENT);++j,++k)
            seq_matrix[sitect-1][k]=codseq[j];
        break;
      }

    /* Find downstream (which we presume to be false) Met's in sequence */
    if(hypothesis==0) {
      if((MetTFGAPnt)%3!=0)
        FATALERROR("Err (main): MetTFGAPnt mod 3 != 0\n")
      else
        for(i+=(MetTFGAPnt);i<=seqlength-lastcdsoff-3;i+=3)
          if(IS_MET_P(codseq,i)==TRUE&&
             MET_VALID_TO_SCORE_P(i,metsiglowbound-1,seqlength-1)==TRUE) {
            if((seq_matrix=
                  (int**)realloc(seq_matrix,++sitect*sizeof(int*)))==NULL||
               (seq_matrix[sitect-1]=
                  (int*)malloc((STRINGSIZE)*sizeof(int)))==NULL)
              FATALERROR("Err (main): Insufficient Memory!\n")
            else
              for(j=i-(UPSTREXTENT),k=0;j<=i+2+(DOWNSTREXTENT);++j,++k)
                seq_matrix[sitect-1][k]=codseq[j];
          }
    }

    free(codseq);
  } /* end while */

  /* Develop the WAMs */
  makeMetWAM(&locATGparms,seq_matrix,sitect);

  /* Close input stream and free memory */
  fclose(seqfile);
  for(i=0;i<sitect;++i)
    free(seq_matrix[i]);
  free(seq_matrix);

  /* Prepare a file for writing the trained ATGparm object to */
  if((outfile=fopen(outputname,"wb"))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,
      "Can't open %s for binary writing!\n",outputname);
    FATALERROR(errbuff)
  }
  else if(fwrite(&locATGparms,sizeof(ATGparm),1,outfile)!=1) {
    (void)snprintf(errbuff,MAXERRLEN,
      "Error exporting ATGparm data to %s!\n",outputname);
    FATALERROR(errbuff)
  }
  else
    fclose(outfile);

  return(EXIT_SUCCESS);
} /* end main */

static void makeMetWAM(
  ATGparm *locATGparms,
  int **seq_matrix,
  int sitect
) {
  int mono_ct[STRINGSIZE-1][NTALFSIZE], /* Mononuc freq */
      di_ct[STRINGSIZE-1][NTALFSIZE][NTALFSIZE], /* Dinuc freq */
      i,j,k; /* iterator variables */

  /* Inits of local variables */
  for(i=0;i<STRINGSIZE;++i)
    for(j=0;j<NTALFSIZE;++j)
      for(k=0;k<NTALFSIZE;++k)
        locATGparms->WAMtable[i][j][k]=INITVAL_FLT;
  for(i=0;i<STRINGSIZE-1;++i)
    for(j=0;j<NTALFSIZE;++j) {
      mono_ct[i][j]=INITVAL_INT;
      for(k=0;k<NTALFSIZE;++k)
        di_ct[i][j][k]=INITVAL_INT;
    }

  /* tabulate mono/dinucleotide frequencies */
  for(i=0;i<STRINGSIZE-1;++i)
    for(j=0;j<sitect;++j) {
      /* Note that the first dimension of mono_ct indexes the  *
       * position (only the first STRINGSIZE-1 are considered) *
       * that nucleotide seq_matrix[j][i] occurs in.           */
      ++mono_ct[i][seq_matrix[j][i]];        

      /* Note that the first dimension of di_ct indexes the    *
       * position that the first base of dinucleotide          *
       * [seq_matrix[j][i]][seq_matrix[j][i+1]] occurs in.     *
       * Only the first STRINGSIZE-1 positions are considered. */
      ++di_ct[i][seq_matrix[j][i]][seq_matrix[j][i+1]];
    }

  /* Record equilibrium frequencies (1st ``slot" in transition freqs). *
   * Note that only base i is ``important", but we populate all bases  *
   * indexed by j, just for consistency.                               */
  if(sitect>0)
    for(i=0;i<NTALFSIZE;++i)
      locATGparms->WAMtable[0][i][Ade]=
        ((double)mono_ct[0][i])/sitect;
  else /* We support a basal, fully uninformed parameterization. */
    for(i=0;i<NTALFSIZE;++i)
      locATGparms->WAMtable[0][i][Ade]=EQUIPROB;
  if(locATGparms->WAMtable[0][Ade][Ade]<MAXFAULT&&
     locATGparms->WAMtable[0][Cyt][Ade]<MAXFAULT&&
     locATGparms->WAMtable[0][Gua][Ade]<MAXFAULT&&
     locATGparms->WAMtable[0][Thy][Ade]<MAXFAULT)
    for(i=0;i<NTALFSIZE;++i)
      locATGparms->WAMtable[0][i][Ade]=EQUIPROB;
  else if(locATGparms->WAMtable[0][Ade][Ade]<MAXFAULT||
          locATGparms->WAMtable[0][Cyt][Ade]<MAXFAULT||
          locATGparms->WAMtable[0][Gua][Ade]<MAXFAULT||
          locATGparms->WAMtable[0][Thy][Ade]<MAXFAULT) {
    for(i=0;i<NTALFSIZE;++i)
      if(locATGparms->WAMtable[0][i][Ade]<MAXFAULT)
        locATGparms->WAMtable[0][i][Ade]=PROBMIN;
      else {
        locATGparms->WAMtable[0][i][Ade]*=1-4*PROBMIN;
        locATGparms->WAMtable[0][i][Ade]+=PROBMIN;
      }
  }
  for(i=0;i<NTALFSIZE;++i)
    for(j=0;j<NTALFSIZE;++j)
      locATGparms->WAMtable[0][i][j]=
        locATGparms->WAMtable[0][i][Ade];

  /* Populate the remaining transition frequencies.                    *
   * Note that locATGparms->WAMtable[k][i][j]                          *
   * records the frequency that dinucleotide [i][j] occurs with base j *
   * occurring in position k, relative to dinucs [i][z], z != j, also  *
   * occurring in that position in the training data.                  */
  for(k=1;k<STRINGSIZE;++k) {
    for(i=0;i<NTALFSIZE;++i) /* rows: give probability of transition to a *
                              * different base j predicated on the first  *
                              * base being i.                             */
      for(j=0;j<NTALFSIZE;++j) /* columns: base being transitioned into */
        if(mono_ct[k-1][i]!=0)
          locATGparms->WAMtable[k][i][j]=
            ((double)di_ct[k-1][i][j])/mono_ct[k-1][i];
        else
          locATGparms->WAMtable[k][i][j]=NULLPROB;

    /* Smooth transition probabilities as appropriate: */
    if(k<UPSTREXTENT||k>UPSTREXTENT+SITELEN)
      for(i=0;i<NTALFSIZE;++i) {
        /* check if all four columns are < MAXFAULT */
        if(locATGparms->WAMtable[k][i][Ade]<MAXFAULT&&
           locATGparms->WAMtable[k][i][Cyt]<MAXFAULT&&
           locATGparms->WAMtable[k][i][Gua]<MAXFAULT&&
           locATGparms->WAMtable[k][i][Thy]<MAXFAULT)
          for(j=0;j<NTALFSIZE;++j) /* Adjust all probabilities */
            locATGparms->WAMtable[k][i][j]=EQUIPROB;
        /* check if any of the four columns are < MAXFAULT */
        else if(locATGparms->WAMtable[k][i][Ade]<MAXFAULT||
                locATGparms->WAMtable[k][i][Cyt]<MAXFAULT||
                locATGparms->WAMtable[k][i][Gua]<MAXFAULT||
                locATGparms->WAMtable[k][i][Thy]<MAXFAULT) {
          for(j=0;j<NTALFSIZE;++j)
            if(locATGparms->WAMtable[k][i][j]<MAXFAULT)
              locATGparms->WAMtable[k][i][j]=PROBMIN;
            else {
              locATGparms->WAMtable[k][i][j]*=1-4*PROBMIN;
              locATGparms->WAMtable[k][i][j]+=PROBMIN;
            }
        }
      }
    else if(k==UPSTREXTENT) /* (xAtg) one column of 1's, all else 0's */
      for(i=0;i<NTALFSIZE;++i)
        for(j=0;j<NTALFSIZE;++j)
          if(j==Ade)
            locATGparms->WAMtable[k][i][j]=MAXPROB;
          else
            locATGparms->WAMtable[k][i][j]=NULLPROB;
    else if(k==UPSTREXTENT+1) /* (aTg) */
      for(i=0;i<NTALFSIZE;++i)
        for(j=0;j<NTALFSIZE;++j)
          if(i==Ade&&j==Thy)
            locATGparms->WAMtable[k][i][j]=MAXPROB;
          else
            locATGparms->WAMtable[k][i][j]=NULLPROB;
    else if(k==UPSTREXTENT+2) /* (atG) */
      for(i=0;i<NTALFSIZE;++i)
        for(j=0;j<NTALFSIZE;++j)
          if(i==Thy&&j==Gua)
            locATGparms->WAMtable[k][i][j]=MAXPROB;
          else
            locATGparms->WAMtable[k][i][j]=NULLPROB;
    else if(k==UPSTREXTENT+SITELEN) /* (atgX) one row of != 0's, all else 0 */
      for(i=0;i<NTALFSIZE;++i)
        if(i==Gua) { 
          /* check if all four columns are < MAXFAULT */
          if(locATGparms->WAMtable[k][i][Ade]<MAXFAULT&&
             locATGparms->WAMtable[k][i][Cyt]<MAXFAULT&&
             locATGparms->WAMtable[k][i][Gua]<MAXFAULT&&
             locATGparms->WAMtable[k][i][Thy]<MAXFAULT)
            for(j=0;j<NTALFSIZE;++j) /* Adjust all probabilities */
              locATGparms->WAMtable[k][i][j]=EQUIPROB;
          /* check if any of the four columns are < MAXFAULT */
          else if(locATGparms->WAMtable[k][i][Ade]<MAXFAULT||
                  locATGparms->WAMtable[k][i][Cyt]<MAXFAULT||
                  locATGparms->WAMtable[k][i][Gua]<MAXFAULT||
                  locATGparms->WAMtable[k][i][Thy]<MAXFAULT) {
            for(j=0;j<NTALFSIZE;++j)
              if(locATGparms->WAMtable[k][i][j]<MAXFAULT)
                locATGparms->WAMtable[k][i][j]=PROBMIN;
              else {
                locATGparms->WAMtable[k][i][j]*=1-4*PROBMIN;
                locATGparms->WAMtable[k][i][j]+=PROBMIN;
              }
          }
        }
        else
          for(j=0;j<NTALFSIZE;++j)
            locATGparms->WAMtable[k][i][j]=NULLPROB;
    else
      FATALERROR("Err (makeMetWAM): Unexpected WAM position!\n")
  } /* end k for */

  return;
} /* end makeMetWAM */
