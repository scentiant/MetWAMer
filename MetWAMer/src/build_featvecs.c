/* build_featvecs.c
 * Michael Sparks (mespar1@gmail.com)
 * Last modified : 20 July 2013
 *
 * Main implementation file for the build_featvecs utility.
 *
 * Copyright (c) 2007,2008,2013 Michael E Sparks
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
#include <immpractical.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "errors.h"
#include "MetWAM_utils.h"
#include "probdef.h"
#include "sequence_parse.h"

#define USAGE {\
  int w;\
  (void)snprintf(errbuff,MAXERRLEN,\
    "\a\nUsage: %s \\\n\
          [-o first_codon_position] \\\n\
          [-e last_codon_offset (distance from sequence stop)] \\\n\
          [-l Met_site_lower_bound] \\\n\
          -m foo.true.MetWAM \\\n\
          -x foo.false.MetWAM \\\n\
          -f codseqs.fas \\\n\
          -a algocode\n\
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
        -m arg denotes a trained MetWAM object characterizing true TISs.\n\
        -x arg denotes a trained MetWAM object characterizing false TISs.\n\
        -f arg is a set of training sequences in Fasta format.\n\
        -a arg is the Markov chain algorithm code, selected from\n",argv[0]);\
  NONFATALERROR(errbuff)\
  for(w=0;w<NUMALGOS;++w) {\
    if(w==PMMMX) /* PMMMs are not applicable here */\
      continue;\
    (void)snprintf(errbuff,MAXERRLEN,\
"             (%i) %s\n",w,algodesc[w]);\
    NONFATALERROR(errbuff)\
  }\
  NONFATALERROR("\
             The Markov chain parameter object file should be the final\n\
             argument of the command line.  (If using dynamically\n\
             modulating Markov models, the final two arguments\n\
             must be the cloud type (0 for variance-based, 1 for range-\n\
             based) and finally the DMMM parameter object file.)\n\n")\
}

/* Main Application */
int main(int argc,char *argv[]) {
  ATGparm locATGparms[2];       /* Stores t and nil ATGparm objects */
  char
    *charcodseq=NULL,           /* coding seq to be processed       */
    codseqid[MAXDESCLINE+1],    /* description of codseq, +1 for    *
                                 * terminating '\0'.                */
    *metwam_t_obj=NULL,         /* pointer to MetWAM filename (t)   */
    *metwam_f_obj=NULL,         /* pointer to MetWAM filename (nil) */
    upseq[CONTENTSWATHLEN+1],   /* Seq to run MC over               */
    downseq[CONTENTSWATHLEN+1], /* Seq to run MC over               */
    *sufftest=NULL;             /* Verifies MC file extension       */
  dmprobT dmprobs;              /* Read the dmprobT object into it  */
  FILE *seqfile=NULL,           /* fasta input file                 */
       *parmfile=NULL;          /* connects to MetWAM object        */
  foprobT fixordprobs;          /* Read the foprobT object into it  */
  immprobT immprobs;            /* Read the immprobT object into it */
  int option,                   /* option parsed by getopt          */
      *codseq=NULL,             /* int translation of sequence      */
      seqlength=0,              /* stores length of codseq array    */
      cdslength=0,              /* stores length of the cds proper  */
      firstcdsbase=1,           /* first cds base position          */
      lastcdsoff=1,             /* last cds base offset position    */
      metsiglowbound=1,         /* lower bound on valid Mets        */
      mcalgocode=-1,            /* code for Markov chain method     */
      cloudtype=-1,             /* Stores cloud definition method   */
      i,j;                      /* Iterator variables               */
  long double
    Metscore,                   /* weighted log-likelihood ratio    */
    upscore,                    /* posterior prob of upstream seq   */
    downscore,                  /* posterior prob of downstream seq */
    updenom,                    /* Pr(Data) for upstream seq        */
    downdenom;                  /* Pr(Data) for downstream seq      */

  /* Parse command line */
  while((option=getopt(argc,argv,":o:e:l:m:x:f:a:"))!=-1) {
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
      case 'm' :
        metwam_t_obj=optarg;
        break;
      case 'x' :
        metwam_f_obj=optarg;
        break;
      case 'f' :
        /* Connect to training sequence input stream */
        if((seqfile=fopen(optarg,"rt"))==NULL) {
          (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Can't seem to open %s!\n",optarg);
          FATALERROR(errbuff)
        }
        break;
      case 'a' :
        mcalgocode=atoi(optarg);
        break;
      default :
        USAGE
        exit(EXIT_FAILURE);
    }
  }
  if(metwam_t_obj==NULL||metwam_f_obj==NULL||seqfile==NULL||
     mcalgocode<0||mcalgocode>=(NUMALGOS)||mcalgocode==PMMMX) {
    USAGE
    exit(EXIT_FAILURE);
  }

  /* Read in the Markov chain parameter file(s) */
  switch(mcalgocode) {
    case FIXORDX :
    case TDDIX :
    case BUDIX :
    case CHI2X :
      if(optind==argc) { /* no parameter file(s) specified */
        USAGE
        exit(EXIT_FAILURE); 
      }
      else if((sufftest=strstr(argv[optind],algotags[mcalgocode]))==NULL) {
        (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Binary parm file %s has improper extension!\n",argv[optind]);
        FATALERROR(errbuff)
      }
      else
        if(mcalgocode==FIXORDX)
          import_fo_probs(argv[optind],&fixordprobs);
        else /* an IMM variant */
          import_imm_probs(argv[optind],&immprobs);
      break;
    case DMMMX :
      /* Quick check: if user is deploying the DMMM, determine if they  *
       * want to use cloud bounded by range- or variance-based methods. */
      if(argc-optind<2) {
        USAGE
        exit(EXIT_FAILURE); 
      }
      else
        cloudtype=atoi(argv[optind++]);

      if((sufftest=strstr(argv[optind],algotags[DMMMX]))==NULL) {
        (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Binary parm file %s has improper extension!\n",argv[optind]);
        FATALERROR(errbuff)
      }
      else
        import_dm_probs(argv[optind],&dmprobs);
      break;
    case PMMMX :
    default :
      FATALERROR("\nErr (main): \
Improper algorithm requested!\n\n")
  } /* end switch */

  /* Read in the true MetWAM parameter object... */
  if((parmfile=fopen(metwam_t_obj,"rb"))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,
      "Can't open %s for binary reading!\n",metwam_t_obj);
    FATALERROR(errbuff)
  }
  else if(fread(&locATGparms[0],sizeof(ATGparm),1,parmfile)!=1) {
    (void)snprintf(errbuff,MAXERRLEN,
      "Error importing ATGparm data from %s!\n",metwam_t_obj);
    FATALERROR(errbuff)
  }
  else
    fclose(parmfile);
  /* ...and the false. */
  if((parmfile=fopen(metwam_f_obj,"rb"))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,
      "Can't open %s for binary reading!\n",metwam_f_obj);
    FATALERROR(errbuff)
  }
  else if(fread(&locATGparms[1],sizeof(ATGparm),1,parmfile)!=1) {
    (void)snprintf(errbuff,MAXERRLEN,
      "Error importing ATGparm data from %s!\n",metwam_f_obj);
    FATALERROR(errbuff)
  }
  else
    fclose(parmfile);

  /* Process the input sequences */
  while((charcodseq=get_fasta(seqfile,charcodseq,&seqlength,codseqid))!=NULL) {
    if(seqlength==0)
      if(charcodseq==NULL)
        FATALERROR("Err (main): No real sequence passed!\n")
      else
        FATALERROR("Err (main): Passed a 0nt genomic codseq!\n")
    else if((codseq=(int*)malloc(sizeof(int)*seqlength))==NULL)
      FATALERROR("Err (main): Out of memory!\n")
    else /* Copy charcodseq translation into codseq */
      for(i=0;i<seqlength;++i)
        codseq[i]=trans(charcodseq[i]);

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
    else if((cdslength=seqlength-lastcdsoff-firstcdsbase+1)%3!=0) {
      (void)snprintf(errbuff,MAXERRLEN,"Warn (main): \
CDS length (%i) is not a multiple of three!??\n\
  (Skipping sequence \"%s\")\n",cdslength,codseqid);
      NONFATALERROR(errbuff)
      free(codseq);
      continue;
    }

    /* Find first \emph{valid} Met in sequence, which we presume to be true */
    for(i=firstcdsbase-1;i<=seqlength-lastcdsoff-3;i+=3)
      if(IS_MET_P(codseq,i)==TRUE&&
         MET_VALID_TO_SCORE_P(i,metsiglowbound-1,seqlength-1)==TRUE&&
         MET_FLANKS_VALID_TO_SCORE_P(i,0,seqlength-1)==TRUE) {
        /* Score the Met (signal) */
        Metscore=calc_llkr(locATGparms[0],locATGparms[1],
                           &codseq[i-(UPSTREXTENT)]);
        Metscore+=
          log(pow(((double)(cdslength-(i-firstcdsbase+1)))/cdslength,3));

        /* upstream flank (content) */
        strncpy(upseq,&charcodseq[i-(UPSTREXTENT)-(CONTENTSWATHLEN)],
          CONTENTSWATHLEN);
        upseq[CONTENTSWATHLEN]='\0';

        /* downstream flank (content) */
        strncpy(downseq,&charcodseq[i+3+(DOWNSTREXTENT)],CONTENTSWATHLEN);
        downseq[CONTENTSWATHLEN]='\0';

        /* compute probability ratio of flanks under coding hypothesis */
        updenom=downdenom=upscore=downscore=0.0;
        switch(mcalgocode) {
          case FIXORDX :
            for(j=0;j<3;++j) {
              updenom+=CODINGPRIOR*
                exp(FOprob(upseq,&fixordprobs,j,TRUE,FALSE));
              downdenom+=CODINGPRIOR*
                exp(FOprob(downseq,&fixordprobs,j,TRUE,FALSE));
            }
            updenom+=NONCODPRIOR*
              exp(FOprob(upseq,&fixordprobs,NONCODING,TRUE,FALSE));
            downdenom+=NONCODPRIOR*
              exp(FOprob(downseq,&fixordprobs,NONCODING,TRUE,FALSE));
            for(j=0;j<3;++j) {
              upscore+=CODINGPRIOR*
                exp(FOprob(upseq,&fixordprobs,j,TRUE,FALSE))/updenom;
              downscore+=CODINGPRIOR*
                exp(FOprob(downseq,&fixordprobs,j,TRUE,FALSE))/downdenom;
            }
            break;
          case TDDIX :
          case BUDIX :
          case CHI2X :
            for(j=0;j<3;++j) {
              updenom+=CODINGPRIOR*
                exp(IMMprob(upseq,&immprobs,j,TRUE,FALSE));
              downdenom+=CODINGPRIOR*
                exp(IMMprob(downseq,&immprobs,j,TRUE,FALSE));
            }
            updenom+=NONCODPRIOR*
              exp(IMMprob(upseq,&immprobs,NONCODING,TRUE,FALSE));
            downdenom+=NONCODPRIOR*
              exp(IMMprob(downseq,&immprobs,NONCODING,TRUE,FALSE));
            for(j=0;j<3;++j) {
              upscore+=CODINGPRIOR*
                exp(IMMprob(upseq,&immprobs,j,TRUE,FALSE))/updenom;
              downscore+=CODINGPRIOR*
                exp(IMMprob(downseq,&immprobs,j,TRUE,FALSE))/downdenom;
            }
            break;
          case DMMMX :
            for(j=0;j<3;++j) {
              updenom+=CODINGPRIOR*
                exp(DMMMprob(upseq,&dmprobs,cloudtype,j,TRUE,FALSE));
              downdenom+=CODINGPRIOR*
                exp(DMMMprob(downseq,&dmprobs,cloudtype,j,TRUE,FALSE));
            }
            updenom+=NONCODPRIOR*
              exp(DMMMprob(upseq,&dmprobs,cloudtype,NONCODING,TRUE,FALSE));
            downdenom+=NONCODPRIOR*
              exp(DMMMprob(downseq,&dmprobs,cloudtype,NONCODING,TRUE,FALSE));
            for(j=0;j<3;++j) {
              upscore+=CODINGPRIOR*
                exp(DMMMprob(upseq,&dmprobs,cloudtype,j,TRUE,FALSE))/
                updenom;
              downscore+=CODINGPRIOR*
                exp(DMMMprob(downseq,&dmprobs,cloudtype,j,TRUE,FALSE))/
                downdenom;
            }
            break;
          case PMMMX :
          default :
            fprintf(stderr,"Err (main): Improper algorithm requested!\n");
            return(EXIT_FAILURE);
        } /* end switch */

        /* Report positive instances to stdout */
        fprintf(stdout,"(+;<%+.4Le,%+.4f>)\n",
          Metscore,log(downscore/upscore));
        break;
      } /* end if */

    /* Find downstream (which we presume to be false) Met's in sequence */
    if((MetTFGAPnt)%3!=0)
      FATALERROR("Err (main): MetTFGAPnt mod 3 != 0\n")
    else
      for(i+=(MetTFGAPnt);i<=seqlength-lastcdsoff-3;i+=3)
        if(IS_MET_P(codseq,i)==TRUE&&
           MET_VALID_TO_SCORE_P(i,metsiglowbound-1,seqlength-1)==TRUE&&
           MET_FLANKS_VALID_TO_SCORE_P(i,0,seqlength-1)==TRUE) {
          /* Score the Met (signal) */
          Metscore=calc_llkr(locATGparms[0],locATGparms[1],
                             &codseq[i-(UPSTREXTENT)]);
          Metscore+=
            log(pow(((double)(cdslength-(i-firstcdsbase+1)))/cdslength,3));

          /* upstream flank (content) */
          strncpy(upseq,&charcodseq[i-(UPSTREXTENT)-(CONTENTSWATHLEN)],
            CONTENTSWATHLEN);
          upseq[CONTENTSWATHLEN]='\0';

          /* downstream flank (content) */
          strncpy(downseq,&charcodseq[i+3+(DOWNSTREXTENT)],CONTENTSWATHLEN);
          downseq[CONTENTSWATHLEN]='\0';

          /* compute likelihood ratio of flanks under coding hypothesis */
          updenom=downdenom=upscore=downscore=0.0;
          switch(mcalgocode) {
            case FIXORDX :
              for(j=0;j<3;++j) {
                updenom+=CODINGPRIOR*
                  exp(FOprob(upseq,&fixordprobs,j,TRUE,FALSE));
                downdenom+=CODINGPRIOR*
                  exp(FOprob(downseq,&fixordprobs,j,TRUE,FALSE));
              }
              updenom+=NONCODPRIOR*
                exp(FOprob(upseq,&fixordprobs,NONCODING,TRUE,FALSE));
              downdenom+=NONCODPRIOR*
                exp(FOprob(downseq,&fixordprobs,NONCODING,TRUE,FALSE));
              for(j=0;j<3;++j) {
                upscore+=CODINGPRIOR*
                  exp(FOprob(upseq,&fixordprobs,j,TRUE,FALSE))/updenom;
                downscore+=CODINGPRIOR*
                  exp(FOprob(downseq,&fixordprobs,j,TRUE,FALSE))/downdenom;
              }
              break;
            case TDDIX :
            case BUDIX :
            case CHI2X :
              for(j=0;j<3;++j) {
                updenom+=CODINGPRIOR*
                  exp(IMMprob(upseq,&immprobs,j,TRUE,FALSE));
                downdenom+=CODINGPRIOR*
                  exp(IMMprob(downseq,&immprobs,j,TRUE,FALSE));
              }
              updenom+=NONCODPRIOR*
                exp(IMMprob(upseq,&immprobs,NONCODING,TRUE,FALSE));
              downdenom+=NONCODPRIOR*
                exp(IMMprob(downseq,&immprobs,NONCODING,TRUE,FALSE));
              for(j=0;j<3;++j) {
                upscore+=CODINGPRIOR*
                  exp(IMMprob(upseq,&immprobs,j,TRUE,FALSE))/updenom;
                downscore+=CODINGPRIOR*
                  exp(IMMprob(downseq,&immprobs,j,TRUE,FALSE))/downdenom;
              }
              break;
            case DMMMX :
              for(j=0;j<3;++j) {
                updenom+=CODINGPRIOR*
                  exp(DMMMprob(upseq,&dmprobs,cloudtype,j,TRUE,FALSE));
                downdenom+=CODINGPRIOR*
                  exp(DMMMprob(downseq,&dmprobs,cloudtype,j,TRUE,FALSE));
              }
              updenom+=NONCODPRIOR*
                exp(DMMMprob(upseq,&dmprobs,cloudtype,NONCODING,TRUE,FALSE));
              downdenom+=NONCODPRIOR*
                exp(DMMMprob(downseq,&dmprobs,cloudtype,NONCODING,TRUE,FALSE));
              for(j=0;j<3;++j) {
                upscore+=CODINGPRIOR*
                  exp(DMMMprob(upseq,&dmprobs,cloudtype,j,TRUE,FALSE))/
                  updenom;
                downscore+=CODINGPRIOR*
                  exp(DMMMprob(downseq,&dmprobs,cloudtype,j,TRUE,FALSE))/
                  downdenom;
              }
              break;
            case PMMMX :
            default :
              fprintf(stderr,"Err (main): Improper algorithm requested!\n");
              return(EXIT_FAILURE);
          } /* end switch */

          /* Report negative instances to stderr */
          fprintf(stderr,"(-;<%+.4Le,%+.4f>)\n",
            Metscore,log(downscore/upscore));
        } /* end if */

    free(codseq);
    free(charcodseq);
  } /* end while */

  /* Close input stream */
  fclose(seqfile);

  return(EXIT_SUCCESS);
} /* end main */
