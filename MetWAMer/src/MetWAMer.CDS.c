/* MetWAMer.CDS.c
 *
 * Michael Sparks (michael.sparks2@usda.gov)
 * Last modified : 21 December 2020
 *   (more extensive error reporting for CDSs whose
 *    lengths aren't multiples of 3. MES)
 *
 * This utility predicts TISs in CDS sequences.
 *
 * Copyright (c) 2007,2013 Michael E Sparks
 * All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */

/* Include Statements */
#include <getopt.h>
#include <immpractical.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calc_medoids.h"
#include "classifiers.h"
#include "contact.h"
#include "errors.h"
#include "index_utils.h"
#include "libxml2_addendum.h"
#include "MetWAM_utils.h"
#include "sequence_parse.h"

#define USAGE {\
  int w;\
  (void)snprintf(errbuff,MAXERRLEN,\
   "\a\nUsage: %s \\\n\
          [-b] \\\n\
          [-s] \\\n\
          [-k cluster_count] \\\n\
          [-d deployment_method] \\\n\
          [-l medoids.xml] \\\n\
          -a algorithm.index \\\n\
          [-w foo.MetPWM] \\\n\
          -m foo.true.MetWAM \\\n\
          -x foo.false.MetWAM \\\n\
          -g test_seqs.fas \\\n\
          [-n perceptron.foo] \\\n\
          [-t perceptron_type.bar] \\\n\
          [-c markovchain.index]\n\n\
  Where -b prints bug reporting information and aborts\n\
        -s enables static cluster-specific parameter deployment\n\
             (use parameters based on first in-frame ATG; default=false)\n\
        -k arg denotes the number of cluster-specific parameters\n\
             to use [default: 1]\n\
        -d arg denotes the parameter deployment strategy to use if k>1,\n\
             selected from\n",argv[0]);\
  NONFATALERROR(errbuff)\
  for(w=0;w<DEPLOYCOUNT;++w) {\
    (void)snprintf(errbuff,MAXERRLEN,\
"             (%i) %s\n",w,deploydesc[w]);\
    NONFATALERROR(errbuff)\
  }\
  NONFATALERROR("\
        -l arg denotes medoids; this option is mandatory if -k is used\n\
        -a arg denotes the index of the translation initiation\n\
             codon identification method to use, selected from\n")\
  for(w=0;w<METHODCOUNT;++w) {\
    (void)snprintf(errbuff,MAXERRLEN,\
"             (%i) %s\n",w,methoddesc[w]);\
    NONFATALERROR(errbuff)\
  }\
  NONFATALERROR("\
        -w arg is only relevant for stratified parameter usage, in which\n\
             the filename handle for the array of such cluster-specific\n\
             position weight matrix object files, available on the local\n\
             filesystem, must be specified. Each MetPWM object file must\n\
             have a suffix corresponding to the id attribute for it,\n\
             as set in the medoids data file.\n\
        -m arg denotes a trained MetWAM object characterizing true TISs\n\
             .OR. if k is used, the filename handle for the array of such\n\
             cluster-specific, true site-characterizing object files,\n\
             available on the local filesystem.  In this latter case,\n\
             each MetWAM object file must have a suffix corresponding to\n\
             the id attribute for it, as set in the medoids data file.\n\
        -x arg denotes a trained MetWAM object characterizing false TISs.\n\
             Regardless of whether cluster-specific parameter deployment\n\
             is to be used, only a single WAM characterizes false sites.\n\
        -g arg denotes the Fasta-formatted CDS sequence file (see manual)\n\
        -n arg is a trained perceptron filename .OR. if k is used, the\n\
             filename handle for the array of such cluster-specific object\n\
             files, available on the local filesystem.  In this latter case,\n\
             each LTUparmT object file must have a suffix corresponding to\n\
             the id attribute for it, as set in the medoids data file, and\n\
             all files must correspond to the same neural element type.\n\
        -t arg is the code for the type of perceptron(s), selected from\n")\
  for(w=0;w<AVAILCLASSIFIERCT;++w) {\
    (void)snprintf(errbuff,MAXERRLEN,\
"             (%i) %s\n",w,classdesc[w]);\
    NONFATALERROR(errbuff)\
  }\
  NONFATALERROR("\
        -c arg is the Markov chain algorithm code, selected from\n")\
  for(w=0;w<NUMALGOS;++w) {\
    if(w==PMMMX)\
      continue;\
    (void)snprintf(errbuff,MAXERRLEN,\
"             (%i) %s\n",w,algodesc[w]);\
    NONFATALERROR(errbuff)\
  }\
  NONFATALERROR("\
             This is required if a content-based translation initiation\n\
             method is called for. The Markov chain parameter object file\n\
             should be the final argument of the command line.  (If using\n\
             dynamically modulating Markov models, the final two arguments\n\
             must be the cloud type (0 for variance-based, 1 for range-\n\
             based) and finally the DMMM parameter object file.)\n\n")\
}

/* Main Application */
int main(int argc, char *argv[]) {
  ATGparmPWM
    *locATGparmsPWM=NULL;     /* Stores ATGparmPWM object          */
  ATGparm *locATGparms;       /* Stores ATGparm object             */
  char
    *charcodseq=NULL,         /* coding seq to be processed        */
    codseqid[MAXDESCLINE+1],  /* description of codseq             */
    *medoidxml=NULL,          /* Stores medoid information         */
    *metpwmobj=NULL,          /* pointer to MetPWM filename        */
    *metwam_t_obj=NULL,       /* pointer to MetWAM filename (true) */
    *metwam_f_obj=NULL,       /* pointer to MetWAM filename (nil)  */
    *neuronobj=NULL,          /* pointer to perceptron filename    */
    *cdspath=NULL,            /* pointer to CDS seq filename       */
    **medoids=NULL,           /* Stores the k-many medoids         */
    *sufftest=NULL,           /* Verifies MC file extension        */
    locname[MAXFILENAME],     /* for storing appended filenames    */
    locnum[MAXFILENAME];      /* string representations of ints    */
  dmprobT *dmprobs=NULL;      /* Read in the dmprobT object(s)     */
  fcwllkrauxT fcwllkrparms;   /* relevant data for fcwllkr         */
  FILE *wamparmfile=NULL,     /* connects to trained MetWAM object */
       *pwmparmfile=NULL,     /* connects to trained MetPWM object */
       *perceptfile=NULL,     /* connects to trained neuron file   */
       *locifile=NULL;        /* genomic loci processed by PASIF   */
  foprobT *fixordprobs=NULL;  /* Read in the foprobT object(s)     */
  immprobT *immprobs=NULL;    /* Read in the immprobT object(s)    */
  int *codseq=NULL,           /* int translation of sequence       */
      use_static_p=0,         /* Static or dynamic parm deployment */
      deploy_meth=-1,         /* hamming, pwm or array             */
      seqlength=0,            /* stores length of codseq array     */
      result,                 /* TIS position or NOTISPRED         */
      option,                 /* option parsed by getopt           */
      transinitmeth=-1,       /* Key into appropriate method       */
      use_stratified_p=FALSE, /* stratified parameter deployment?  */
      kint,                   /* int of k attribute                */
      stringsizeint,          /* int of stringsize attribute       */
      medoidseqlen,           /* length of parsed medoid sequence  */
      mcalgocode=-1,          /* code for Markov chain method      */
      cloudtype=-1,           /* Stores cloud definition method    */
      classifiercode=-1,      /* Designates classification machine */
      k=1,                    /* stores number of clusters         */
      medoidid,               /* identifier of medoid              */
      i;                      /* iterator variable                 */
  LTUparmT *perceptron=NULL;  /* trained perceptron                */
  void *MCparms=NULL;         /* Markov chain parameter object     */
  xmlChar *medoidseq;         /* stores the medoid sequence        */
  xmlDocPtr medxmldoc;        /* xml document object               */
  xmlNodePtr myNodePtr;       /* for traversing DOM trees          */
  xmlXPathContextPtr
    xpcontext;                /* eval context of xpath query       */
  xmlXPathObjectPtr
    xpresults;                /* results of xpath query            */

  /* Parse command line */
  while((option=getopt(argc,argv,":bsk:d:l:a:w:m:x:g:n:t:c:"))!=-1) {
    switch(option) {
      case 'b' :
        BUGGERS_OFF
        return(EXIT_SUCCESS);
        break;
      case 's' :
        use_static_p=1;
        break;
      case 'k' :
        use_stratified_p=TRUE;
        if((k=atoi(optarg))<2) {
          (void)snprintf(errbuff,MAXERRLEN,
            "Err (main): Nonsensible value of k (= %i < 2) passed!\n",k);
          FATALERROR(errbuff)
        }
        break;
      case 'd' :
        if((deploy_meth=atoi(optarg))<0||deploy_meth>=(int)(DEPLOYCOUNT)) {
          USAGE
          exit(EXIT_FAILURE);
        }
        break;
      case 'l' :
        medoidxml=optarg;
        break;
      case 'a' :
        transinitmeth=atoi(optarg);
        break;
      case 'w' :
        metpwmobj=optarg;
        break;
      case 'm' :
        metwam_t_obj=optarg;
        break;
      case 'x' :
        metwam_f_obj=optarg;
        break;
      case 'g' :
        cdspath=optarg;
        break;
      case 'n' :
        neuronobj=optarg;
        break;
      case 't' :
        classifiercode=atoi(optarg);
        break;
      case 'c' :
        mcalgocode=atoi(optarg);
        break;
      default :
        USAGE
        return(EXIT_FAILURE);
    }
  }
  if(metwam_t_obj==NULL||metwam_f_obj==NULL||cdspath==NULL||
     transinitmeth<0||transinitmeth>=METHODCOUNT) {
    USAGE
    return(EXIT_FAILURE);
  }
  else if(use_stratified_p==TRUE) {
    if(medoidxml==NULL)
      FATALERROR("Err (main): \
k was specified, but a medoid file was not!??\n")
    else if(deploy_meth==-1)
      FATALERROR("Err (main): \
k was specified, but a deployment method was not!??\n")
    else if(deploy_meth==PWMMODX&&metpwmobj==NULL)
      FATALERROR("Err (main): \
k was specified under PWM, but a PWM file handle was not!??\n")
    if(use_static_p==1)
      deploy_meth+=(DEPLOYCOUNT);
    cluster_parm_selection=deploy_meth;
  }

  /* For content-based methods, a Markov chain method must be specified. */
  if(transinitmeth==MFCWLLKRX||transinitmeth==NFCWLLKRX) {
    if(mcalgocode<0||mcalgocode>=NUMALGOS||mcalgocode==PMMMX||
       optind==argc) {
      USAGE
      exit(EXIT_FAILURE);
    }
    else
      switch(mcalgocode) {
        case FIXORDX :
        case TDDIX :
        case BUDIX :
        case CHI2X :
          if((sufftest=strstr(argv[optind],algotags[mcalgocode]))==NULL) {
            (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Binary parm file %s has improper extension!\n",argv[optind]);
            FATALERROR(errbuff)
          }
          else {
            if(mcalgocode==FIXORDX) {
              if((fixordprobs=(foprobT*)malloc(sizeof(foprobT)))==NULL)
                FATALERROR("Err (main): Insufficient memory\n")
              else {
                import_fo_probs(argv[optind],fixordprobs);
                MCparms=fixordprobs;
              }
            }
            else {/* an IMM variant */
              if((immprobs=(immprobT*)malloc(sizeof(immprobT)))==NULL)
                FATALERROR("Err (main): Insufficient memory\n")
              else {
                import_imm_probs(argv[optind],immprobs);
                MCparms=immprobs;
              }
            }
          }
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
          else {
            if((dmprobs=(dmprobT*)malloc(sizeof(dmprobT)))==NULL)
              FATALERROR("Err (main): Insufficient memory\n")
            else {
              import_dm_probs(argv[optind],dmprobs);
              MCparms=dmprobs;
            }
          }
          break;
        case PMMMX :
        default :
          FATALERROR("\nErr (main): \
Improper algorithm requested!\n\n")
      } /* end switch */
  }

  /* Check that neuronal information was passed for NFCWLLKR */
  if(transinitmeth==NFCWLLKRX) {
    if(neuronobj==NULL||
       classifiercode<0||classifiercode>=AVAILCLASSIFIERCT) {
      USAGE
      exit(EXIT_FAILURE);
    }
    else if((perceptron=(LTUparmT*)malloc(sizeof(LTUparmT)*k))==NULL)
      FATALERROR("Err (main): Insufficient memory\n")
    else
      for(i=0;i<k;++i) {
        strcpy(locname,neuronobj);
        if(use_stratified_p==TRUE) {
          snprintf(locnum,MAXFILENAME,".%i",i+1);
          strcat(locname,locnum);
        }

        if((perceptfile=fopen(locname,"rb"))==NULL) {
          (void)snprintf(errbuff,MAXERRLEN,
            "Err (main): Can't open %s for binary reading!\n",locname);
          FATALERROR(errbuff)
        }
        else if(fread(&perceptron[i],sizeof(LTUparmT),1,perceptfile)!=1) {
          (void)snprintf(errbuff,MAXERRLEN,
            "Err (main): Error importing LTUparmT data from %s!\n",locname);
          FATALERROR(errbuff)
        }
        else
          fclose(perceptfile);
      } /* end for */
  }

  /* Import medoid xml document and parse medoids, *
   * if a stratified approach is requested.        */
  if(use_stratified_p==TRUE) {
    xmlInitParser();
    LIBXML_TEST_VERSION
    if((medxmldoc=xmlParseFile(medoidxml))==NULL) {
      (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Error parsing medoid xml document %s!\n",medoidxml);
      FATALERROR(errbuff)
    }
    else if((xpcontext=xmlXPathNewContext(medxmldoc))==NULL)
      FATALERROR("Err (main): Error setting xpath eval context!\n")
    else if(reg_ns(xpcontext,BAD_CAST "k=http://www.kmedoids.org\0")!=TRUE)
      FATALERROR("Err (main): \
Error registering namespaces in xpath eval context object!\n")
    else if((xpresults=xmlXPathEvalExpression(
        BAD_CAST "/k:kmedoids",xpcontext))==NULL)
      FATALERROR("Err (main): Error executing xpath query\n")
    else if(xpresults->nodesetval&&xpresults->nodesetval->nodeNr!=1) {
      (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Medoid xml document %s is nonsensible \
(contains %i kmedoids elts)!\n",
        medoidxml,xpresults->nodesetval->nodeNr);
      FATALERROR(errbuff)
    }
    else {
      myNodePtr=xpresults->nodesetval->nodeTab[0];
      RECORD_ATTR_VAL_AS_INT(myNodePtr,stringsizeint,"stringsize")
      RECORD_ATTR_VAL_AS_INT(myNodePtr,kint,"k")
      if(kint!=k) {
        (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Medoid xml document %s differs in k (%i, expected %i)!\n",medoidxml,kint,k);
        FATALERROR(errbuff)
      }
      else if((medoids=(char**)malloc(sizeof(char*)*k))==NULL)
        FATALERROR("Err (main): Insufficient memory\n")
      else
        for(i=0;i<k;++i)
          if((medoids[i]=(char*)malloc(sizeof(char)*(stringsizeint+1)))==NULL)
            FATALERROR("Err (main): Insufficient memory\n")

      DESCEND_TO_CHILD(myNodePtr,"medoid")
      for(i=0;myNodePtr!=NULL;++i) {
        RECORD_ATTR_VAL_AS_INT(myNodePtr,medoidid,"id")
        --medoidid;

        medoidseq=xmlNodeListGetString(medxmldoc,myNodePtr->children,1);
        if((medoidseqlen=xmlStrlen(medoidseq))!=stringsizeint) {
          (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Medoid in %s has incorrect length (%i, expected %i)!\n",
            medoidxml,medoidseqlen,stringsizeint);
          FATALERROR(errbuff)
        }
        else {
          strncpy(medoids[medoidid],(char*)medoidseq,stringsizeint);
          medoids[medoidid][stringsizeint]='\0';
        }

        xmlFree(medoidseq);
        FIND_SIBLING(myNodePtr,"medoid")
      }
      if(i!=k) {
        (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Medoid xml document %s hasn't k (%i) medoids!??\n",medoidxml,k);
        FATALERROR(errbuff)
      }

      /* cleanup */
      xmlXPathFreeObject(xpresults);
      xmlXPathFreeContext(xpcontext);
      xmlFreeDoc(medxmldoc);
      xmlCleanupParser();
    }
  }

  /* Read in the MetPWM parameter objects, if appropriate */
  if(deploy_meth==PWMMODX||deploy_meth==PWMSTATX) {
    if((locATGparmsPWM=(ATGparmPWM*)malloc(sizeof(ATGparmPWM)*k))==NULL)
      FATALERROR("Err (main): Insufficient memory\n")
    for(i=0;i<k;++i) {
      strcpy(locname,metpwmobj);
      snprintf(locnum,MAXFILENAME,".%i",i+1);
      strcat(locname,locnum);

      if((pwmparmfile=fopen(locname,"rb"))==NULL) {
        (void)snprintf(errbuff,MAXERRLEN,
          "Err (main): Can't open %s for binary reading!\n",locname);
        FATALERROR(errbuff)
      }
      else if(fread(&locATGparmsPWM[i],sizeof(ATGparmPWM),1,pwmparmfile)!=1) {
        (void)snprintf(errbuff,MAXERRLEN,
          "Err (main): Error importing ATGparm data from %s!\n",locname);
        FATALERROR(errbuff)
      }
      else
        fclose(pwmparmfile);
    }
  }

  /* Read in the true MetWAM parameter object(s)... */
  if((locATGparms=(ATGparm*)malloc(sizeof(ATGparm)*(k+1)))==NULL)
    FATALERROR("Err (main): Insufficient memory\n")
  for(i=0;i<k;++i) {
    strcpy(locname,metwam_t_obj);
    if(use_stratified_p==TRUE) {
      snprintf(locnum,MAXFILENAME,".%i",i+1);
      strcat(locname,locnum);
    }

    if((wamparmfile=fopen(locname,"rb"))==NULL) {
      (void)snprintf(errbuff,MAXERRLEN,
        "Err (main): Can't open %s for binary reading!\n",locname);
      FATALERROR(errbuff)
    }
    else if(fread(&locATGparms[i],sizeof(ATGparm),1,wamparmfile)!=1) {
      (void)snprintf(errbuff,MAXERRLEN,
        "Err (main): Error importing ATGparm data from %s!\n",locname);
      FATALERROR(errbuff)
    }
    else
      fclose(wamparmfile);
  }
  /* ...and the false MetWAM object. */
  if((wamparmfile=fopen(metwam_f_obj,"rb"))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,
      "Err (main): Can't open %s for binary reading!\n",metwam_f_obj);
    FATALERROR(errbuff)
  }
  else if(fread(&locATGparms[k],sizeof(ATGparm),1,wamparmfile)!=1) {
    (void)snprintf(errbuff,MAXERRLEN,
      "Err (main): Error importing ATGparm data from %s!\n",metwam_f_obj);
    FATALERROR(errbuff)
  }
  else
    fclose(wamparmfile);

  /* Regardless of whether the ?fcwllkr methods will be *
   * used, we will populate the fcwllkrparms variable.  */
  fcwllkrparms.markovchainmeth=mcalgocode;
  fcwllkrparms.dmcloudindex=cloudtype;
  fcwllkrparms.neurontype=classifiercode;
  fcwllkrparms.localneuron=perceptron;
  fcwllkrparms.markovchainparms=MCparms;

  /* Connect to testing (CDS) sequence input stream */
  if((locifile=fopen(cdspath,"rt"))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,"Can't open %s!\n",cdspath);
    FATALERROR(errbuff)
  }
  /* Predict TISs */
  while((charcodseq=get_fasta(locifile,charcodseq,
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

    if((result=searchforbestinitincds(transinitmeth,k,medoids,
        locATGparms,locATGparmsPWM,&fcwllkrparms,codseq,seqlength))
        !=(NOTMULT3))
      printf("%s\n%i\n",codseqid,result);
    else {
      (void)snprintf(errbuff,MAXERRLEN,"\tIgnoring %s\n\n",codseqid);
      NONFATALERROR(errbuff)
    }

    free(codseq);
  } /* end while */
  fclose(locifile);

  /* cleanup memory and return to OS */
  if(transinitmeth==MFCWLLKRX||transinitmeth==NFCWLLKRX)
    free(MCparms);
  if(transinitmeth==NFCWLLKRX)
    free(perceptron);
  if(use_stratified_p==TRUE) {
    for(i=0;i<k;++i)
      free(medoids[i]);
    free(medoids);
    if(deploy_meth==PWMMODX||deploy_meth==PWMSTATX)
      free(locATGparmsPWM);
  }
  free(locATGparms);
  return(EXIT_SUCCESS);
} /* end main */
