/* train_MetPWM.c
 * Michael Sparks (mespar1@gmail.com)
 * Last modified : 20 July 2013
 *
 * Main implementation file for the train_MetPWM utility.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "errors.h"
#include "MetWAM_utils.h"
#include "probdef.h"
#include "sequence_parse.h"

#define ARGCT 3
#define USAGE "\a\nUsage: %s sites.fas outfile\n\
\n\
  Where sites.fas contains a set of translation initiation sites, and\n\
        outfile is the desired name for the PWM output file.\n\
\n"

/* Function prototypes */
void makeMetPWM(ATGparmPWM *locATGparmsPWM,int **seq_matrix,int sitect);

/* Main Application */
int main(int argc,char *argv[]) {
  ATGparmPWM locATGparmsPWM; /* Stores ATGparmPWM object  */
  char *sequence=NULL,       /* TIS to be processed       */
    seqid[MAXDESCLINE+1],    /* description of sequence,  *
                              * +1 for terminating '\0'.  */
    outputname[MAXFILENAME]; /* Stores output file name   */
  FILE *seqfile=NULL,        /* fasta input file          */
       *outfile=NULL;        /* write locATGparms to it   */
  int nelts,                 /* total number of sites     */
      **sites=NULL,          /* stores TISs               */
      seqlen=0,              /* stores length of sequence */
      i,j;                   /* Iterator variables        */

  /* Verify command line */
  if(argc!=ARGCT) {
    (void)snprintf(errbuff,MAXERRLEN,USAGE,argv[0]);
    FATALERROR(errbuff)
  }

  /* Learn the number of sequences in the input file... */
  if((seqfile=fopen(argv[1],"rt"))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (main): Can't open %s.\n",argv[1]);
    FATALERROR(errbuff)
  }
  nelts=0;
  while((sequence=get_fasta(seqfile,sequence,&seqlen,seqid))!=NULL) {
    if(seqlen==0)
      if(sequence==NULL)
        FATALERROR("Err (main): No real sequence passed!\n")
      else
        FATALERROR("Err (main): Passed a 0nt sequence!\n")
    else if(seqlen!=(int)(STRINGSIZE)) {
      (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Unexpected sequence length of %i (expected %i)\n",seqlen,(int)(STRINGSIZE));
      FATALERROR(errbuff)
    }
    else
      ++nelts;
    free(sequence);
  }
  /* ...and slurp these into the sites matrix */
  rewind(seqfile);
  if((sites=(int**)malloc(sizeof(int*)*nelts))==NULL)
    FATALERROR("Err (main): Out of memory.\n")
  else
    for(i=0;i<nelts;++i)
      if((sites[i]=(int*)malloc(sizeof(int)*(int)(STRINGSIZE)))==NULL)
        FATALERROR("Err (main): Out of memory.\n")
      else {
        sequence=get_fasta(seqfile,sequence,&seqlen,seqid);
        for(j=0;j<(int)(STRINGSIZE);++j)
          sites[i][j]=trans(sequence[j]);
        free(sequence);
      }

  /* Develop the PWM */
  makeMetPWM(&locATGparmsPWM,sites,nelts);

  /* Close input stream and free memory */
  fclose(seqfile);
  for(i=0;i<nelts;++i)
    free(sites[i]);
  free(sites);

  /* Prepare a file for writing the trained ATGparmPWM object to */
  (void)strcpy(outputname,argv[2]);
  (void)strcat(outputname,".MetPWM");
  if((outfile=fopen(outputname,"wb"))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,
      "Can't open %s for binary writing!\n",outputname);
    FATALERROR(errbuff)
  }
  else if(fwrite(&locATGparmsPWM,sizeof(ATGparmPWM),1,outfile)!=1) {
    (void)snprintf(errbuff,MAXERRLEN,
      "Error exporting ATGparmPWM data to %s!\n",outputname);
    FATALERROR(errbuff)
  }
  else
    fclose(outfile);

  #ifdef DBGPWM
  for(i=0;i<STRINGSIZE;++i) {
    for(j=0;j<NTALFSIZE;++j) {
      (void)snprintf(errbuff,MAXERRLEN,
        "%f  ",locATGparmsPWM.PWMtable[i][j]);
      NONFATALERROR(errbuff)
    }
    NONFATALERROR("\n")
  }
  #endif

  return(EXIT_SUCCESS);
} /* end main */

void makeMetPWM(
  ATGparmPWM *locATGparmsPWM,
  int **seq_matrix,
  int sitect
) {
  int mono_ct[STRINGSIZE][NTALFSIZE], /* Mononucleotide freq */
      i,j;                            /* iterator variables  */

  /* Inits of local variables */
  for(i=0;i<STRINGSIZE;++i)
    for(j=0;j<NTALFSIZE;++j) {
      mono_ct[i][j]=INITVAL_INT;
      locATGparmsPWM->PWMtable[i][j]=INITVAL_FLT;
    }

  /* Tabulate mononucleotide frequencies */
  for(i=0;i<STRINGSIZE;++i)
    for(j=0;j<sitect;++j)
      ++mono_ct[i][seq_matrix[j][i]];        

  /* Record position-specific nucleotide frequencies */
  for(i=0;i<STRINGSIZE;++i)
    if(sitect>0)
      for(j=0;j<NTALFSIZE;++j)
        locATGparmsPWM->PWMtable[i][j]=((double)mono_ct[i][j])/sitect;
    else /* We support a basal, fully uninformed parameterization. */
      for(j=0;j<NTALFSIZE;++j)
        locATGparmsPWM->PWMtable[i][j]=EQUIPROB;

  /* Smooth sites, if necessary */
  for(i=0;i<STRINGSIZE;++i) {
    if(i>=(UPSTREXTENT)&&i<(UPSTREXTENT)+(SITELEN))
      continue;
    else if(locATGparmsPWM->PWMtable[i][Ade]<MAXFAULT&&
            locATGparmsPWM->PWMtable[i][Cyt]<MAXFAULT&&
            locATGparmsPWM->PWMtable[i][Gua]<MAXFAULT&&
            locATGparmsPWM->PWMtable[i][Thy]<MAXFAULT)
      for(j=0;j<NTALFSIZE;++j)
        locATGparmsPWM->PWMtable[i][j]=EQUIPROB;
    else if(locATGparmsPWM->PWMtable[i][Ade]<MAXFAULT||
            locATGparmsPWM->PWMtable[i][Cyt]<MAXFAULT||
            locATGparmsPWM->PWMtable[i][Gua]<MAXFAULT||
            locATGparmsPWM->PWMtable[i][Thy]<MAXFAULT) {
      for(j=0;j<NTALFSIZE;++j)
        if(locATGparmsPWM->PWMtable[i][j]<MAXFAULT)
          locATGparmsPWM->PWMtable[i][j]=PROBMIN;
        else {
          locATGparmsPWM->PWMtable[i][j]*=1-4*PROBMIN;
          locATGparmsPWM->PWMtable[i][j]+=PROBMIN;
        }
    }
  }

  return;
} /* end makeMetWAM */
