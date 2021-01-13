/* eval_MetWAM.c
 * Michael Sparks (mespar1@gmail.com)
 * Last modified : 20 July 2013
 *
 * Computes log-likelihood ratios of Met's in test sequences.
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
#include "errors.h"
#include "MetWAM_utils.h"
#include "sequence_parse.h"

#define ARGCT 4
#define USAGE "\a\nUsage: %s x.true.MetWAM \
y.false.MetWAM seqs_to_analyze.fas\n\n"

/* Main Application */
int main(int argc, char *argv[]) {
  ATGparm locATGparms_t,     /* Stores ATGparm object (true)  */
          locATGparms_f;     /* Stores ATGparm object (false) */
  char *charcodseq=NULL,     /* coding seq to be processed    */   
    codseqid[MAXDESCLINE+1]; /* description of codseq, +1 for *
                              * terminating '\0'.             */
  FILE *parmfile=NULL,       /* Connect to parameter object   */
       *seqfile=NULL;        /* Connect to sequence stream    */
  int *codseq=NULL,          /* int translation of sequence   */
      *codseqrev=NULL,       /* reverse complement of codseq  */
      seqlength=0,           /* stores length of codseq array */
      h,i,j;                 /* iterator variables            */
  long double loglikerat;    /* stores log-likelihood ratio   */

  /* Verify command line */
  if(argc!=ARGCT) {
    fprintf(stderr,USAGE,argv[0]);
    exit(EXIT_FAILURE);
  }

  /* Read in the true site parameter object */
  if((parmfile=fopen(argv[1],"rb"))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,
      "Can't open %s for binary reading!\n",argv[1]);
    FATALERROR(errbuff)
  }
  else if(fread(&locATGparms_t,sizeof(ATGparm),1,parmfile)!=1) {
    (void)snprintf(errbuff,MAXERRLEN,
      "Error importing ATGparm data from %s!\n",argv[1]);
    FATALERROR(errbuff)
  }
  else
    fclose(parmfile);

  /* Read in the false site parameter object */
  if((parmfile=fopen(argv[2],"rb"))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,
      "Can't open %s for binary reading!\n",argv[2]);
    FATALERROR(errbuff)
  }
  else if(fread(&locATGparms_f,sizeof(ATGparm),1,parmfile)!=1) {
    (void)snprintf(errbuff,MAXERRLEN,
      "Error importing ATGparm data from %s!\n",argv[2]);
    FATALERROR(errbuff)
  }
  else
    fclose(parmfile);

  /* Process sequences piped from seqfile */
  if((seqfile=fopen(argv[3],"rt"))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,"Can't open %s!\n",argv[3]);
    FATALERROR(errbuff)
  }
  while((charcodseq=get_fasta(seqfile,charcodseq,
                              &seqlength,codseqid))!=NULL) {
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

    fprintf(stdout,"\nProcessing: \"%s\"\n",codseqid);

    /* Compute log-likelihood ratios of Met's *
     * encountered in each reading frame.     */
    for(h=0;h<3;++h) {
      fprintf(stdout,"  Frame %i:\n",h+1);
      for(i=h;i<seqlength-2;i+=3)
        if(IS_MET_P(codseq,i)==TRUE&&
           MET_VALID_TO_SCORE_P(i,0,seqlength-1)==TRUE) {
          loglikerat=calc_llkr(locATGparms_t,locATGparms_f,
                               &codseq[i-(UPSTREXTENT)]);
          fprintf(stdout,"    %-8i: %+.4Le\n",i+1,loglikerat);
        }
    }

    /* Derive the reverse complement... */
    if((codseqrev=(int*)malloc(sizeof(int)*seqlength))==NULL)
      FATALERROR("Err (main): Out of memory!\n")
    else
      for(i=0,j=seqlength-1;j>=0;++i,--j)
        codseqrev[i]=basecomp(codseq[j]);
    free(codseq);
    if(i!=seqlength)
      FATALERROR("Err (main): codseq and codseqrev lengths differ?\n")

    /* ...and score its Met's, too. */
    for(;h<6;++h) {
      fprintf(stdout,"  Frame %i:\n",h+1);
      for(i=h-3;i<seqlength-2;i+=3)
        if(IS_MET_P(codseqrev,i)==TRUE&&
           MET_VALID_TO_SCORE_P(i,0,seqlength-1)==TRUE) {
          loglikerat=calc_llkr(locATGparms_t,locATGparms_f,
                               &codseqrev[i-(UPSTREXTENT)]);
          fprintf(stdout,"    %-8i: %+.4Le\n",(seqlength-1-i)+1,loglikerat);
        }
    }
    free(codseqrev);

    fprintf(stdout,"\n--\n");
  } /* end while */
  fclose(seqfile);

  return(EXIT_SUCCESS);
} /* end main */
