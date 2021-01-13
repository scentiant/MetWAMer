/* MetWAM_utils.c
 * Michael Sparks (michael.sparks2@usda.gov)
 * Last modified : 21 December 2020
 *   (more extensive error reporting for CDSs whose
 *    lengths aren't multiples of 3. MES)
 *
 * Functions germane to calculations derived from a trained MetWAM.
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
#include <immpractical.h>
#include <libxml/parser.h>
#include <libxml/relaxng.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "errors.h"
#include "index_utils.h"
#include "libxml2_addendum.h"
#include "probdef.h"
#include "MetWAM_utils.h"
#include "sequence_parse.h"

/* Function prototypes *******************************************************/

static long double calc_llk(ATGparm parms,int *seq);
static int *derivemaximalorfwflanks(int transinitmeth,char *genomicseq,
  exoncoorsT *structure,int *orf_start,int *orf_stop);
static int searchforbestinit(int transinitmeth,int strata_ct,char **medoids,
  ATGparm *ATGparms,ATGparmPWM *ATGparmsPWM,fcwllkrauxT *fcwllkrparms,
  int *sequence,exoncoorsT *coords,int proc_pasif_p,long double *bestTISscore,
  int *Met_seen);
static int transinit_llkr(int strata_ct,char **medoids,ATGparm *ATGparms,
  ATGparmPWM *ATGparmsPWM,int *coords_orig,int *sequence,int seqlen,
  int start_orig,long double *bestTISscore,int *Met_seen);
static int transinit_wllkr(int strata_ct,char **medoids,ATGparm *ATGparms,
  ATGparmPWM *ATGparmsPWM,int *coords_orig,int *sequence,int seqlen,
  int start_orig,long double *bestTISscore,
  long double *bestTISscorenoweight,int *Met_seen);
static int transinit_bayes(int strata_ct,char **medoids,ATGparm *ATGparms,
  ATGparmPWM *ATGparmsPWM,int *coords_orig,int *sequence,int seqlen,
  int start_orig,long double *bestTISscore,int *best_is_true,int *Met_seen);
static int transinit_fcwllkr(int transinitmeth,int strata_ct,char **medoids,
  ATGparm *ATGparms,ATGparmPWM *ATGparmsPWM,fcwllkrauxT *fcwllkrparms,
  int *coords_orig,int *sequence,int seqlen,int start_orig,
  long double *bestTISscore,int *fc_class,int *Met_seen);
static double eval_pwm_score(ATGparmPWM PWM,char *site);
 
/* Function definitions ******************************************************/

/* Function to coordinate prediction of translation *
 * initiation sites based on <orf_entry> data in    *
 * gthXML documents.                                */
void predicttransinitsites(
  int transinitmeth,
  int strata_ct,
  char **medoids,
  ATGparm *ATGparms,
  ATGparmPWM *ATGparmsPWM,
  fcwllkrauxT *fcwllkrparms,
  char *template,
  const char *xmldoc,
  const char *relaxngdoc
) {
  char methodbuff[BUFFSIZE]; /* stores prediction method info     */
  int
    *maximalorfwflanks=NULL, /* will store non-stop reading frame *
                              * and appropriate flanksing seqs    */
    rf_start=-1,             /* overall reading frame start       */
    orf_start=-1,            /* open reading frame start          */
    orf_stop=-1,             /* open reading frame stop           */
    Met_seen=0,              /* Met seen in max reading frame?    */
    i;                       /* iterator variable                 */
  long double bestTISscore;  /* Score of best TIS, or default     *
                              * value if none exists.             */
  exoncoorsT
    *currexon=NULL, /* Pointers for traversing linked list */
    *prevexon=NULL,
    *structure=NULL;
  const xmlChar *registerednamespaces=BAD_CAST
    "g=http://www.genomethreader.org/GTH_output/ \
     gh=http://www.genomethreader.org/GTH_output/header/ \
     gppatp=http://www.genomethreader.org/GTH_output/PGL_module/\
predicted_gene_location/AGS_information/\
three_phase_translation/probable_ORFs/\0"; /* namespaces */
  xmlDocPtr gthxmldoc;          /* xml document object         */
  xmlNodePtr myNodePtr;         /* for traversing trees        */
  xmlXPathContextPtr xpcontext; /* eval context of xpath query */
  xmlXPathObjectPtr xpresults;  /* results of xpath query      */

  /* Import gthXML document */
  xmlInitParser();
  LIBXML_TEST_VERSION
  if((gthxmldoc=xmlParseFile(xmldoc))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (predicttransinitsites): \
Error parsing gthXML document %s!\n",xmldoc);
    FATALERROR(errbuff)
  }

  /* Validate gthXML document, if RELAX NG grammar was specified */
  if(relaxngdoc!=NULL&&validateagainstRNG(relaxngdoc,gthxmldoc)==FALSE) {
    xmlFreeDoc(gthxmldoc);
    xmlCleanupParser();
    free(template);
    (void)snprintf(errbuff,MAXERRLEN,"Err (predicttransinitsites): \
gthXML file %s is incorrect per %s!\n",xmldoc,relaxngdoc);
    FATALERROR(errbuff)
  }

  /* Set eval context and register relevant namespaces */
  if((xpcontext=xmlXPathNewContext(gthxmldoc))==NULL)
    FATALERROR("Err (predicttransinitsites): \
Error setting xpath eval context!\n")
  else if(reg_ns(xpcontext,registerednamespaces)!=TRUE)
    FATALERROR("Err (predicttransinitsites): \
Error registering namespaces in xpath eval context object!\n")

  /* Sanity check: Verify that at most one template was  *
   * used in producing the specified gthXML document, if *
   * the document contains a header element.             */
  if((xpresults=xmlXPathEvalExpression(
       BAD_CAST "/g:GTH_output/gh:header/\
gh:gDNA_template_files/gh:temp_name",xpcontext))==NULL)
    FATALERROR("Err (predicttransinitsites): \
Error executing xpath query\n")
  else if(xpresults->nodesetval&&xpresults->nodesetval->nodeNr>1) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (predicttransinitsites): \
One template per gthXML doc only ---\n\
  %s has %i listed!  (Please see the included Readme.pdf)\n",
      xmldoc,xpresults->nodesetval->nodeNr);
    xmlXPathFreeObject(xpresults);
    xmlXPathFreeContext(xpcontext);
    xmlFreeDoc(gthxmldoc);
    xmlCleanupParser();
    free(template);
    FATALERROR(errbuff)
  }
  else
    xmlXPathFreeObject(xpresults);

  if((xpresults=xmlXPathEvalExpression(
       BAD_CAST "//gppatp:orf_entry",xpcontext))==NULL)
    FATALERROR("Err (predicttransinitsites): \
Error executing xpath query\n")

  /* Process each relevant node in turn */
  for(i=0;i<xpresults->nodesetval->nodeNr;++i) {
    myNodePtr=xpresults->nodesetval->nodeTab[i];
    DESCEND_TO_CHILD(myNodePtr,"id_line")
    DESCEND_TO_CHILD(myNodePtr,"orf_info")
    DESCEND_TO_CHILD(myNodePtr,"exon_boundaries")
    DESCEND_TO_CHILD(myNodePtr,"exon")

    /* Build up an exoncoorsT linked list, representing the gene structure.  *
     * Note that coordinates parsed from gthXML documents will be one-based, *
     * but we manipulate them here using a zero-based system.                */
    if((currexon=(exoncoorsT*)malloc(sizeof(exoncoorsT)))==NULL)
      FATALERROR("Err (predicttransinitsites): Insufficient memory!\n")
    else {
      RECORD_ATTR_VAL_AS_INT(myNodePtr,currexon->start,"start")
      --currexon->start;
      RECORD_ATTR_VAL_AS_INT(myNodePtr,currexon->stop,"stop")
      --currexon->stop;
      structure=prevexon=currexon;
      orf_start=currexon->start;
      FIND_SIBLING(myNodePtr,"exon")
    }
    while(myNodePtr!=NULL) {
      if((currexon=(exoncoorsT*)malloc(sizeof(exoncoorsT)))==NULL)
        FATALERROR("Err (predicttransinitsites): Insufficient memory!\n")
      else {
        RECORD_ATTR_VAL_AS_INT(myNodePtr,currexon->start,"start")
        --currexon->start;
        RECORD_ATTR_VAL_AS_INT(myNodePtr,currexon->stop,"stop")
        --currexon->stop;
        prevexon->nextelt=currexon;
        prevexon=currexon;
        FIND_SIBLING(myNodePtr,"exon")
      }
    }
    currexon->nextelt=NULL;
    orf_stop=currexon->stop;

    /* Derive sequence to search over, which will include    *
     * terminal overhangs.  Note that orf_start and orf_stop *
     * will be modified to reflect the maximal orf (ignoring *
     * the added flanks).                                    */
    if((maximalorfwflanks=derivemaximalorfwflanks(transinitmeth,
         template,structure,&orf_start,&orf_stop))==NULL) {
      /* There was not enough excess 5'- and/or 3'-terminal overhang  *
       * to permit scoring Methionines, i.e., we had boundary issues. */
      freestructure(structure);
      continue;
    }
    else { /* update the linked list to reflect the maximal orf */
      structure->start=rf_start=orf_start;
      currexon->stop=orf_stop;
    }

    /* Identify the best translation start bound */
    orf_start=searchforbestinit(transinitmeth,strata_ct,medoids,
      ATGparms,ATGparmsPWM,fcwllkrparms,maximalorfwflanks,structure,
      FALSE,&bestTISscore,&Met_seen);

    /* Record specifics on the TIS prediction method. */
    (void)strncpy(methodbuff,"",BUFFSIZE);
    if(strata_ct>1) { /* cluster-specific parm deployment strategy */
      if(cluster_parm_selection==HAMMODX)
        (void)strncpy(methodbuff,
          ", hamming-distance based (modulating) ; ",
          BUFFSIZE);
      else if(cluster_parm_selection==PWMMODX)
        (void)strncpy(methodbuff,
          ", position weight matrix-distance based (modulating) ; ",
          BUFFSIZE);
      else if(cluster_parm_selection==WAMMODX)
        (void)strncpy(methodbuff,
          ", position weight matrix-distance based (modulating) ; ",
          BUFFSIZE);
      else if(cluster_parm_selection==HAMSTATX)
        (void)strncpy(methodbuff,
          ", position weight matrix-distance based (static) ; ",
          BUFFSIZE);
      else if(cluster_parm_selection==PWMSTATX)
        (void)strncpy(methodbuff,
          ", position weight matrix-distance based (static) ; ",
          BUFFSIZE);
      else if(cluster_parm_selection==WAMSTATX)
        (void)strncpy(methodbuff,
          ", position weight matrix-distance based (static) ; ",
          BUFFSIZE);
      else
        FATALERROR("Err (predicttransinitsites): \
Bad parameter deployment method specified.\n")
    }
    if(transinitmeth==LLKRX)
      (void)strncat(methodbuff,
        "log-likelihood ratios",
        BUFFSIZE);
    else if(transinitmeth==WLLKRX)
      (void)strncat(methodbuff,
        "weighted log-likelihood ratios",
        BUFFSIZE);
    else if(transinitmeth==BAYESX)
      (void)strncat(methodbuff,
        "Bayesian",
        BUFFSIZE);
    else if(transinitmeth==MFCWLLKRX)
      (void)strncat(methodbuff,
        "multiplicative flank contrasting with log-likelihood ratios",
        BUFFSIZE);
    else if(transinitmeth==NFCWLLKRX)
      (void)strncat(methodbuff,
        "neuron-based flank contrasting with log-likelihood ratios",
        BUFFSIZE);
    else
      FATALERROR("Err (predicttransinitsites): \
Bad algorithm specified.\n")

    /* Update XML document tree */
    myNodePtr=xpresults->nodesetval->nodeTab[i];
    if((myNodePtr=xmlNewChild(myNodePtr,NULL,
        BAD_CAST "MetWAMer_annot",NULL))==NULL)
      FATALERROR("Err (predicttransinitsites): \
Unable to create a new child node!\n")
    else {
      METHOD_PRINTER(myNodePtr,strata_ct,methodbuff)
      /* Report start/stop codons, on a one-based scale */
      SET_INT_AS_ATTR_VAL(myNodePtr,"orf_amino_bound",rf_start+1)
      SET_INT_AS_ATTR_VAL(myNodePtr,"start_codon",orf_start+1)
      if(orf_start<orf_stop)
        SET_INT_AS_ATTR_VAL(myNodePtr,"stop_codon",orf_stop-1)
      else
        SET_INT_AS_ATTR_VAL(myNodePtr,"stop_codon",orf_stop+3)
      SET_REAL_AS_ATTR_VAL(myNodePtr,"best_TIS_score",(float)bestTISscore)
      SET_INT_AS_ATTR_VAL(myNodePtr,"Met_seen_p",Met_seen)
    }
    /* Add whitespace formatting manually (much safer!). */
    xmlNodeAddContent(myNodePtr->prev,BAD_CAST "  ");
    ASCEND_TO_PARENT(myNodePtr)
    if((myNodePtr=xmlNewText(BAD_CAST "\n            "))==NULL)
      FATALERROR("Err (predicttransinitsites): \
Unable to create a new child node!\n")
    else if((myNodePtr=xmlAddChild(
             xpresults->nodesetval->nodeTab[i],myNodePtr))==NULL)
      FATALERROR("Err (predicttransinitsites): \
Unable to splice in a new child node!\n")

    /* Cleanup memory */
    free(maximalorfwflanks);
    freestructure(structure);
  } /* end <orf_entry> for */

  /* Output results of post-processing */
  xmlDocFormatDump(UPDATEDXMLOUTSTREAM,gthxmldoc,0);

  xmlXPathFreeObject(xpresults);
  xmlXPathFreeContext(xpcontext);
  xmlFreeDoc(gthxmldoc);
  xmlCleanupParser();

  return;
} /* end predicttransinitsites */

/* Function to coordinate prediction of translation   *
 * initiation sites based on <isoform> data, produced *
 * specifically by the PASIF system.                  */
void predicttransinitsites4PASIF(
  int transinitmeth,
  int strata_ct,
  char **medoids,
  ATGparm *ATGparms,
  ATGparmPWM *ATGparmsPWM,
  fcwllkrauxT *fcwllkrparms,
  FILE *loci,
  const char *xmldoc,
  const char *relaxngdoc
) { 
  char *charloc=NULL,          /* genomic locus to be processed     */
       locusid[MAXDESCLINE+1], /* description of locus, +1 for '\0' */
       *idtag,                 /* used to recover original locus    */
       methodbuff[BUFFSIZE];   /* stores prediction method info     */
  long double bestTISscore;    /* Score of best TIS, or default     *
                                * value if none exists.             */
  exoncoorsT
    *currexon=NULL,            /* for traversing linked lists       */
    *prevexon=NULL,
    *structure=NULL;
  int
    locuslength=0,             /* stores length of locus sequence   */
    *maximalorfwflanks=NULL,   /* will store non-stop reading frame *
                                * and appropriate flanksing seqs    */
    rf_start=-1,               /* overall reading frame start       */
    orf_start=-1,              /* open reading frame start          */
    orf_stop=-1,               /* open reading frame stop           */
    Met_seen=0,                /* Met seen in max reading frame?    */
    i;                         /* iterator variable                 */
  const xmlChar *registerednamespaces=BAD_CAST
    "p=http://www.pasif.org/PASIF_output/\0"; /* namespaces */
  xmlDocPtr pasifxmldoc;        /* xml document object         */
  xmlNodePtr myNodePtr,         /* for traversing trees        */
             myNodePtr2;
  xmlXPathContextPtr xpcontext; /* eval context of xpath query */
  xmlXPathObjectPtr xpresults;  /* results of xpath query      */

  /* Import PASIF xml document */
  xmlInitParser();
  LIBXML_TEST_VERSION
  if((pasifxmldoc=xmlParseFile(xmldoc))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (predicttransinitsites4PASIF): \
Error parsing gthXML document %s!\n",xmldoc);
    FATALERROR(errbuff)
  }

  /* Validate PASIF xml document, if RELAX NG grammar was specified */
  if(relaxngdoc!=NULL&&validateagainstRNG(relaxngdoc,pasifxmldoc)==FALSE) {
    xmlFreeDoc(pasifxmldoc);
    xmlCleanupParser();
    (void)snprintf(errbuff,MAXERRLEN,"Err (predicttransinitsites4PASIF): \
PASIF xml file %s is incorrect per %s!\n",xmldoc,relaxngdoc);
    FATALERROR(errbuff)
  }

  /* Set eval context and register relevant namespaces */
  if((xpcontext=xmlXPathNewContext(pasifxmldoc))==NULL)
    FATALERROR("Err (predicttransinitsites4PASIF): \
Error setting xpath eval context!\n")
  else if(reg_ns(xpcontext,registerednamespaces)!=TRUE)
    FATALERROR("Err (predicttransinitsites4PASIF): \
Error registering namespaces in xpath eval context object!\n")

  if((xpresults=xmlXPathEvalExpression(
       BAD_CAST "//p:locus",xpcontext))==NULL)
    FATALERROR("Err (predicttransinitsites4PASIF): \
Error executing xpath query\n")

  /* Process each relevant node in turn */
  for(i=0;i<xpresults->nodesetval->nodeNr;++i) {
    myNodePtr=xpresults->nodesetval->nodeTab[i];
    RECORD_ATTR_VAL_AS_CHAR_STRING(myNodePtr,idtag,"id")

    /* Obtain the corresponding genomic locus */
    while((charloc=get_fasta(loci,charloc,&locuslength,locusid))!=NULL)
      if(strcmp(locusid+1,idtag))
        free(charloc);
      else
        break;

    if(charloc==NULL) {
      (void)snprintf(errbuff,MAXERRLEN,"Err (predicttransinitsites4PASIF): \
Could not find \"%s\" in loci file!\n",idtag);
      xmlXPathFreeObject(xpresults);
      xmlXPathFreeContext(xpcontext);
      xmlFreeDoc(pasifxmldoc);
      xmlCleanupParser();
      free(idtag);
      free(charloc);
      fclose(loci);
      FATALERROR(errbuff)
    }
    else
      rewind(loci);

    /* It is possible that no isoforms were predicted at the given locus. */
    ATTEMPT_DESCENT_TO_CHILD(myNodePtr,"isoform")
    if(myNodePtr==NULL) {
      free(idtag);
      free(charloc);
      continue;
    }
    do { /* process each predicted isoform */
      myNodePtr2=myNodePtr;
      DESCEND_TO_CHILD(myNodePtr2,"exon")

      /* Build up an exoncoorsT linked list, representing the gene *
       * structure.  Note that coordinates parsed from PASIF xml   *
       * documents will be one-based but we manipulate them here   *
       * using a zero-based system.                                */
      if((currexon=(exoncoorsT*)malloc(sizeof(exoncoorsT)))==NULL)
        FATALERROR("Err (predicttransinitsites4PASIF): \
Insufficient memory!\n")
      else {
        RECORD_ATTR_VAL_AS_INT(myNodePtr2,currexon->start,"start")
        --currexon->start;
        RECORD_ATTR_VAL_AS_INT(myNodePtr2,currexon->stop,"stop")
        --currexon->stop;
        structure=prevexon=currexon;
        orf_start=currexon->start;
        FIND_SIBLING(myNodePtr2,"exon")
      }
      while(myNodePtr2!=NULL) {
        if((currexon=(exoncoorsT*)malloc(sizeof(exoncoorsT)))==NULL)
          FATALERROR("Err (predicttransinitsites4PASIF): \
Insufficient memory!\n")
        else {
          RECORD_ATTR_VAL_AS_INT(myNodePtr2,currexon->start,"start")
          --currexon->start;
          RECORD_ATTR_VAL_AS_INT(myNodePtr2,currexon->stop,"stop")
          --currexon->stop;
          prevexon->nextelt=currexon;
          prevexon=currexon;
          FIND_SIBLING(myNodePtr2,"exon")
        }
      }
      currexon->nextelt=NULL;
      orf_stop=currexon->stop;

      /* Derive sequence to search over, which will include    *
       * terminal overhangs.  Note that orf_start and orf_stop *
       * will be modified to reflect the maximal orf (ignoring *
       * the added flanks).                                    */
      if((maximalorfwflanks=derivemaximalorfwflanks(transinitmeth,
           charloc,structure,&orf_start,&orf_stop))==NULL) {
        /* There was not enough excess 5'- and/or 3'-terminal overhang  *
         * to permit scoring Methionines, i.e., we had boundary issues. */
        freestructure(structure);
        FIND_SIBLING(myNodePtr,"isoform")
        continue;
      }
      else { /* update the linked list to reflect the maximal orf */
        structure->start=rf_start=orf_start;
        currexon->stop=orf_stop;
      }

      /* Identify the best translation start bound */
      orf_start=searchforbestinit(transinitmeth,strata_ct,medoids,
        ATGparms,ATGparmsPWM,fcwllkrparms,maximalorfwflanks,structure,
        TRUE,&bestTISscore,&Met_seen);

      /* Record specifics on the TIS prediction method. */
      (void)strncpy(methodbuff,"",BUFFSIZE);
      if(strata_ct>1) { /* cluster-specific parm deployment strategy */
        if(cluster_parm_selection==HAMMODX)
          (void)strncpy(methodbuff,
            ", hamming-distance based (modulating) ; ",
            BUFFSIZE);
        else if(cluster_parm_selection==PWMMODX)
          (void)strncpy(methodbuff,
            ", position weight matrix-distance based (modulating) ; ",
            BUFFSIZE);
        else if(cluster_parm_selection==WAMMODX)
          (void)strncpy(methodbuff,
            ", position weight matrix-distance based (modulating) ; ",
            BUFFSIZE);
        else if(cluster_parm_selection==HAMSTATX)
          (void)strncpy(methodbuff,
            ", position weight matrix-distance based (static) ; ",
            BUFFSIZE);
        else if(cluster_parm_selection==PWMSTATX)
          (void)strncpy(methodbuff,
            ", position weight matrix-distance based (static) ; ",
            BUFFSIZE);
        else if(cluster_parm_selection==WAMSTATX)
          (void)strncpy(methodbuff,
            ", position weight matrix-distance based (static) ; ",
            BUFFSIZE);
        else
          FATALERROR("Err (predicttransinitsites4PASIF): \
Bad parameter deployment method specified.\n")
      }
      if(transinitmeth==LLKRX)
        (void)strncat(methodbuff,
          "log-likelihood ratios",
          BUFFSIZE);
      else if(transinitmeth==WLLKRX)
        (void)strncat(methodbuff,
          "weighted log-likelihood ratios",
          BUFFSIZE);
      else if(transinitmeth==BAYESX)
        (void)strncat(methodbuff,
          "Bayesian",
          BUFFSIZE);
      else if(transinitmeth==MFCWLLKRX)
        (void)strncat(methodbuff,
          "multiplicative flank contrasting with log-likelihood ratios",
          BUFFSIZE);
      else if(transinitmeth==NFCWLLKRX)
        (void)strncat(methodbuff,
          "neuron-based flank contrasting with log-likelihood ratios",
          BUFFSIZE);
      else
        FATALERROR("Err (predicttransinitsites4PASIF): \
Bad algorithm specified.\n")

      /* Update XML document tree */
      myNodePtr2=myNodePtr;
      if((myNodePtr2=xmlNewChild(myNodePtr2,NULL,
          BAD_CAST "MetWAMer_annot",NULL))==NULL)
        FATALERROR("Err (predicttransinitsites4PASIF): \
Unable to create a new child node!\n")
      else {
        METHOD_PRINTER(myNodePtr2,strata_ct,methodbuff)
        /* Report start/stop codons, on a one-based scale */
        SET_INT_AS_ATTR_VAL(myNodePtr2,"orf_amino_bound",rf_start+1)
        SET_INT_AS_ATTR_VAL(myNodePtr2,"start_codon",orf_start+1)
        if(orf_start<orf_stop)
          SET_INT_AS_ATTR_VAL(myNodePtr2,"stop_codon",orf_stop-1)
        else {
          (void)snprintf(errbuff,MAXERRLEN,
            "Warn (predicttransinitsites4PASIF): \
PASIF xml file has %i as start and %i as stop?\n\
  (Expected forward-sense only!)\n",orf_start+1,orf_stop-1);
          NONFATALERROR(errbuff)
          SET_INT_AS_ATTR_VAL(myNodePtr2,"stop_codon",orf_stop+3)
        }
        SET_REAL_AS_ATTR_VAL(myNodePtr2,"best_TIS_score",(float)bestTISscore)
        SET_INT_AS_ATTR_VAL(myNodePtr2,"Met_seen_p",Met_seen)
      }
      /* Add whitespace formatting manually (much safer!). */
      xmlNodeAddContent(myNodePtr2->prev,BAD_CAST "  ");
      ASCEND_TO_PARENT(myNodePtr2)
      if((myNodePtr2=xmlNewText(BAD_CAST "\n    "))==NULL)
        FATALERROR("Err (predicttransinitsites4PASIF): \
Unable to create a new child node!\n")
      else if((myNodePtr2=xmlAddChild(myNodePtr,myNodePtr2))==NULL)
        FATALERROR("Err (predicttransinitsites4PASIF): \
Unable to splice in a new child node!\n")

      /* cleanup isoform-specific memory */
      free(maximalorfwflanks);
      freestructure(structure);

      FIND_SIBLING(myNodePtr,"isoform")
    } while(myNodePtr!=NULL);

    /* cleanup locus-specific memory */
    free(idtag);
    free(charloc);
  } /* end <locus> for */

  /* Output results of post-processing */
  xmlDocFormatDump(UPDATEDXMLOUTSTREAM,pasifxmldoc,0);

  xmlXPathFreeObject(xpresults);
  xmlXPathFreeContext(xpcontext);
  xmlFreeDoc(pasifxmldoc);
  xmlCleanupParser();

  return;
} /* end predicttransinitsites4PASIF */

/* Calculate log-likelihood of a putative translation initiation    *
 * codon, whose first base (A) corresponds to *(seq+(UPSTREXTENT)). */
static long double calc_llk(
  ATGparm parms,
  int *seq
) {
  int i;               /* iterator variable     */
  long double loglike; /* Stores log-likelihood */

  /* Calculate likelihood under the specified hypothesis.  The Ade used *
   * to compute the first slot's probability might just as well be Cyt, *
   * Gua, or Thy!                                                       */
  loglike=log(parms.WAMtable[0][seq[0]][Ade]);
  for(i=1;i<STRINGSIZE;++i)
    loglike+=log(parms.WAMtable[i][seq[i-1]][seq[i]]);

  return(loglike);
} /* end calc_llk */

/* Calculate log-likelihood ratio of a putative translation *
 * initiation codon, whose first base (A) corresponds to    *
 * *(seq+(UPSTREXTENT)).                                    */
long double calc_llkr(
  ATGparm parms_true,
  ATGparm parms_false,
  int *seq
) {
  long double loglikes[HYPOTHESISCT]; /* Store log-likelihoods */

  /* Sanity check */
  if(*(seq+(UPSTREXTENT))!=Ade||
     *(seq+(UPSTREXTENT)+1)!=Thy||
     *(seq+(UPSTREXTENT)+2)!=Gua)
    NONFATALERROR("Warning (calc_llkr): Scoring non-ATG site!??\n")

  /* Calculate likelihood under each hypothesis. */
  loglikes[TRUEMODIND]=calc_llk(parms_true,seq);
  loglikes[FALSEMODIND]=calc_llk(parms_false,seq);

  /* Return the likelihood ratio.      *
   * Note: log(A/B) := log(A) - log(B) */
  return(loglikes[TRUEMODIND]-loglikes[FALSEMODIND]);
} /* end calc_llkr */

/* Returns TRUE if gthXML doc validates against RELAX NG schema, FALSE o.w. */
int validateagainstRNG(
  const char *relaxngdoc,
  xmlDocPtr gthxmldoc
) {
  int valstat; /* Store status of RNG schema validation */
  xmlDocPtr rngdoc;
  xmlRelaxNGParserCtxtPtr rngctxt;
  xmlRelaxNGPtr rngparsed;
  xmlRelaxNGValidCtxtPtr rngctxtptr;

  if((rngdoc=xmlParseFile(relaxngdoc))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (validateagainstRNG): \
Error parsing RELAX NG document %s!\n",relaxngdoc);
    FATALERROR(errbuff)
  }
  if((rngctxt=xmlRelaxNGNewDocParserCtxt(rngdoc))==NULL)
    FATALERROR("Err (validateagainstRNG): \
Can't create parser context for RELAX NG schema.\n")
  if((rngparsed=xmlRelaxNGParse(rngctxt))==NULL)
    FATALERROR("Err (validateagainstRNG): \
Can't create a parsed RELAX NG schema.\n")
  if((rngctxtptr=xmlRelaxNGNewValidCtxt(rngparsed))==NULL)
    FATALERROR("Err (validateagainstRNG): \
Can't create ptr to parser context for RELAX NG schema.\n")

  valstat=xmlRelaxNGValidateDoc(rngctxtptr,gthxmldoc);
  xmlRelaxNGFreeValidCtxt(rngctxtptr);
  xmlRelaxNGFree(rngparsed);
  xmlRelaxNGFreeParserCtxt(rngctxt);
  xmlFreeDoc(rngdoc);

  return(valstat==0?TRUE:FALSE);
} /* end validateagainstRNG */

/* Function to modify xpcontext such that registerednamespaces *
 * will be registered, and thus available to XPath queries.    *
 * registerednamespaces format is                              *
 * abbrev1=verbatim1 abbrev2=verbatim2 ... abbrev3=verbatim3   */
int reg_ns(
  xmlXPathContextPtr xpcontext,
  const xmlChar *registerednamespaces
) {
  xmlChar *loc_regns, /* copy of registerednamespaces */
          *help_ptr,  /* navigate loc_regns with it   */
          *abbrev,    /* alias of full namespace URL  */
          *verbatim;  /* full URL of namespace        */

  if((loc_regns=xmlStrdup(registerednamespaces))==NULL)
    FATALERROR("Err (reg_ns): \
Error copying registerednamespaces!\n")

  for(help_ptr=loc_regns;help_ptr!=NULL; ) {
    if(help_ptr[0]==' ') { /* scroll past whitespace */
      ++help_ptr;
      continue;
    }

    abbrev=help_ptr;
    if((help_ptr=(xmlChar*)xmlStrchr(help_ptr,'='))==NULL) {
      xmlFree(loc_regns);
      NONFATALERROR("Warning (reg_ns): \
Blatantly nonconformant crud in registerednamespaces!\n")
      return(FALSE);
    }
    else
      *help_ptr='\0';

    verbatim=++help_ptr;
    if((help_ptr=(xmlChar*)xmlStrchr(help_ptr,' '))!=NULL)
      *(help_ptr++)='\0';

    if(xmlXPathRegisterNs(xpcontext,abbrev,verbatim)!=0) {
      (void)snprintf(errbuff,MAXERRLEN,"Err (reg_ns): \
Couldn't register %s=%s\n",abbrev,verbatim);
      FATALERROR(errbuff)
    }
  }

  xmlFree(loc_regns);
  return(TRUE);
} /* end reg_ns */

/* Returns an int array of the maximally extended reading frame, *
 * with upstream and downstream flanking sequences; of lengths   *
 * prependlen and appendlen, respectively; if possible.          */
static int *derivemaximalorfwflanks(
  int transinitmeth,     /* Different methods require differing   *
                          * amounts of flanking sequences.        */
  char *genomicseq,      /* Genomic template, whose first elt's   *
                          * index is 0, last strlen(genomicseq)-1 */
  exoncoorsT *structure, /* stores gene structure info            */
  int *orf_start,        /* overall orf start bound, 0-based      */
  int *orf_stop          /* overall orf stop bound, 0-based       */
) {
  char strand='^';          /* Watson (+) or Crick (-)        */
  exoncoorsT
    *help_ptr=structure;    /* used for traversing structure  */
  int prependlen=-1,        /* N-terminus nt flank            */
      appendlen=-1,         /* C-terminus nt flank            */
      *orf=NULL,            /* stores int translation of orf  */
      orf_pos=-1,           /* tracks current position in orf */
      codon[3],             /* for testing codon identity     */
      genomicseqlen=
        strlen(genomicseq), /* length of genomic template     */
      i,j;                  /* iterator variables             */

  /* Learn the CDS's orientation */
  if(*orf_start<*orf_stop) /* Watson strand */
    strand='+';
  else if(*orf_start>*orf_stop) /* Crick strand */
    strand='-';
  else
    FATALERROR("Err (derivemaximalorfwflanks): \
A one-base coding sequence would seem impossible?!!\n")

  /* Based on the translation initiation detection method, we need *
   * to set the appropriate lengths for flanking sequences, which  *
   * will be added to the ends of the orf variable.                */
  prependlen=(int)(UPSTREXTENT);
  appendlen=(int)(DOWNSTREXTENT);
  if(transinitmeth==MFCWLLKRX||transinitmeth==NFCWLLKRX) {
    prependlen+=(int)(CONTENTSWATHLEN);
    appendlen+=(int)(CONTENTSWATHLEN);
  }

  /* Search upstream of CDS, in-frame, until a stop codon is encountered. *
   * Similarly, extend the CDS downstream of *orf_stop, if possible.      *
   * Note that, prependlen and appendlen are enacted here, to make it     *
   * possible to score all possible Methionine residues in the CDS.       */
  if(strand=='+') { /* Watson strand */
    for(i=*orf_start-3;i>=0;i-=3) {
      for(j=0;j<3;++j)
        codon[j]=trans(genomicseq[i+j]);
      if(IS_STOP_P(codon,0)==TRUE) {
        *orf_start=i+3;
        break;
      }
    }
    if((*orf_start-=prependlen)<0)
      return(NULL);

    for(i=*orf_stop-2;i<genomicseqlen-2;i+=3) {
      for(j=0;j<3;++j)
        codon[j]=trans(genomicseq[i+j]);
      if(IS_STOP_P(codon,0)==TRUE) {
        *orf_stop=i+2;
        break;
      }
    }
    if((*orf_stop+=appendlen)>=genomicseqlen)
      return(NULL);
  }
  else { /* Crick strand */
    for(i=*orf_start+3;i<genomicseqlen;i+=3) {
      for(j=0;j<3;++j)
        codon[j]=basecomp(trans(genomicseq[i-j]));
      if(IS_STOP_P(codon,0)==TRUE) {
        *orf_start=i-3;
        break;
      }
    }
    if((*orf_start+=prependlen)>=genomicseqlen)
      return(NULL);

    for(i=*orf_stop+2;i>=2;i-=3) {
      for(j=0;j<3;++j)
        codon[j]=basecomp(trans(genomicseq[i-j]));
      if(IS_STOP_P(codon,0)==TRUE) {
        *orf_stop=i-2;
        break;
      }
    }
    if((*orf_stop-=appendlen)<0)
      return(NULL);
  }

  /* Macro to ammend a fragment to the orfptr argument */
  #define COPYFRAG(gseq,orfptr,orfpos,start,stop) { \
    char *fragment=NULL; \
    int fraglen=abs((start)-(stop))+1, \
        q; \
    if((fragment=parseseqfrag(gseq,(start),(stop)))==NULL) \
      FATALERROR("Err (derivemaximalorfwflanks): \
Indexed substring outside range!??\n") \
    else if((orfptr=(int*)realloc( \
              orfptr,sizeof(int)*(fraglen+orfpos+1)))==NULL) \
      FATALERROR("Err (derivemaximalorfwflanks): Insufficient Memory!\n") \
    if((start)==(stop)&&strand=='-') \
      TAKE_COMPLEMENT(fragment[0],fragment[0]) \
    for(q=0;q<fraglen;++q) \
      orf[++orfpos]=trans(fragment[q]); \
    free(fragment); \
  }

  /* We assume that, from head to tail, the exoncoorsT linked *
   * list elements are arranged in the order of exons in the  *
   * actual CDS they represent.                               */
  if(help_ptr->nextelt==NULL) /* single-exon CDS */
    COPYFRAG(genomicseq,orf,orf_pos,*orf_start,*orf_stop)
  else {
    COPYFRAG(genomicseq,orf,orf_pos,*orf_start,help_ptr->stop)

    for(help_ptr=help_ptr->nextelt;
        help_ptr->nextelt!=NULL;
        help_ptr=help_ptr->nextelt)
      COPYFRAG(genomicseq,orf,orf_pos,help_ptr->start,help_ptr->stop)

    COPYFRAG(genomicseq,orf,orf_pos,help_ptr->start,*orf_stop)
  }

  /* Revert overall start/stop bounds to reflect the maximal ORF, proper */
  if(strand=='+') { /* Watson strand */
    *orf_start+=prependlen;
    *orf_stop-=appendlen;
  }
  else { /* Crick strand */
    *orf_start-=prependlen;
    *orf_stop+=appendlen;
  }

  return(orf);
} /* end derivemaximalorfwflanks */

/* Finds the most likely translation initiation site from the  *
 * amino-terminus of the reading frame, and returns this value *
 * using a zero-based coordinate system.  If no "good" site    *
 * exists, it simply returns the start of the maximal reading  *
 * frame as calculated from, e.g., derivemaximalorfwflanks.    */
static int searchforbestinit(
  int transinitmeth,
  int strata_ct,
  char **medoids,
  ATGparm *ATGparms,
  ATGparmPWM *ATGparmsPWM,
  fcwllkrauxT *fcwllkrparms,
  int *sequence,
  exoncoorsT *coords,
  int proc_pasif_p,
  long double *bestTISscore,
  int *Met_seen
) {
  exoncoorsT
    *helper=NULL;      /* for traversing linked list  */
  int
    prependlen=-1,     /* N-terminus nt flank         */
    appendlen=-1,      /* C-terminus nt flank         */
    *coords_orig=NULL, /* mapping derived coordinates *
                        * to the original ones        */
    seqlen=-1,         /* sequence length             */
    watson_p=-1,       /* Forward-sense strand?       */
    start_codon=-1,    /* Index to return to caller   */
    begin,end,         /* index orf fragment bounds   */
    i,j,               /* iterator variables          */
    nouse_int;         /* dummy variable, used only   *
                        * for interface consistency   */
  long double
    nouse_dbl;         /* dummy variable, used only   *
                        * for interface consistency   */

  /* Learn orientation of orf */
  for(helper=coords;helper->nextelt!=NULL;helper=helper->nextelt)
    ;
  if(coords->start<helper->stop)
    watson_p=TRUE;
  else if(coords->start>helper->stop)
    watson_p=FALSE;
  else
    FATALERROR("Err (searchforbestinit): \
A one-base coding sequence would seem impossible?!!\n")

  /* Based on the translation initiation detection method, we    *
   * need to set the appropriate lengths for flanking sequences. */
  prependlen=(int)(UPSTREXTENT);
  appendlen=(int)(DOWNSTREXTENT);
  if(transinitmeth==MFCWLLKRX||transinitmeth==NFCWLLKRX) {
    prependlen+=(int)(CONTENTSWATHLEN);
    appendlen+=(int)(CONTENTSWATHLEN);
  }

  /* Learn length of rf, plus flanking sequences, in nucleotides. *
   * Recall that the overall boundaries in coords do not reflect  *
   * the flanking sequences, but rather the reading frame proper, *
   * embedded in sequence (implicitly, via genomic coordinates).  */
  helper=coords;
  if(helper->nextelt==NULL) /* single-exon CDS */
    if(watson_p==TRUE)
      seqlen=(helper->stop+appendlen)-(helper->start-prependlen)+1;
    else
      seqlen=(helper->start+prependlen)-(helper->stop-appendlen)+1;
  else {
    if(watson_p==TRUE)
      seqlen=helper->stop-(helper->start-prependlen)+1;
    else
      seqlen=(helper->start+prependlen)-helper->stop+1;
    for(helper=helper->nextelt;helper->nextelt!=NULL;helper=helper->nextelt)
      seqlen+=abs(helper->start-helper->stop)+1;
    if(watson_p==TRUE)
      seqlen+=(helper->stop+appendlen)-helper->start+1;
    else
      seqlen+=helper->start-(helper->stop-appendlen)+1;
  }

  /* Verify that the maximal reading frame's length is a multiple of 3. *
   * Because the PASIF algorithm does not take predicted translation    *
   * products into account, the three-periodicity of its results cannot *
   * be guaranteed.                                                     */
  if(proc_pasif_p==FALSE&&(seqlen-prependlen-appendlen)%3!=0) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (searchforbestinit): \
Maximal reading frame's length is not a multiple of 3!?? (length = %i)\n",
      seqlen-prependlen-appendlen);
    FATALERROR(errbuff)
  }

  /* Build a parallel array storing original gDNA coordinates */
  if((coords_orig=(int*)malloc(sizeof(int)*seqlen))==NULL)
    FATALERROR("Err (searchforbestinit): Insufficient Memory!\n")
  else {
    for(helper=coords,i=0;helper!=NULL;helper=helper->nextelt) {
      /* First, establish fragment begin/end bounds */
      if(helper==coords) /* first exon */
        if(helper->nextelt==NULL) /* single-exon CDS */
          if(watson_p==TRUE) {
            begin=helper->start-prependlen;
            end=helper->stop+appendlen;
          }
          else {
            begin=helper->start+prependlen;
            end=helper->stop-appendlen;
          }
        else {
          if(watson_p==TRUE)
            begin=helper->start-prependlen;
          else
            begin=helper->start+prependlen;
          end=helper->stop;
        }
      else if(helper->nextelt==NULL) { /* last exon */
        begin=helper->start;
        if(watson_p==TRUE)
          end=helper->stop+appendlen;
        else
          end=helper->stop-appendlen;
      }
      else { /* internal exon */
        begin=helper->start;
        end=helper->stop;
      }

      /* Now, slurp these original coordinates into a local buffer */
      if(watson_p==TRUE)
        for(j=begin;j<=end;++j)
          coords_orig[i++]=j;
      else
        for(j=begin;j>=end;--j)
          coords_orig[i++]=j;
    }
    if(i!=seqlen)
      FATALERROR("Err (searchforbestinit): Inconsistent Coordinates!\n")
  }

  /* Call the appropriate calculation method */
  if(transinitmeth==LLKRX)
    start_codon=transinit_llkr(strata_ct,medoids,ATGparms,ATGparmsPWM,
      coords_orig,sequence,seqlen,coords->start,bestTISscore,Met_seen);
  else if(transinitmeth==WLLKRX)
    start_codon=transinit_wllkr(strata_ct,medoids,ATGparms,ATGparmsPWM,
      coords_orig,sequence,seqlen,coords->start,bestTISscore,&nouse_dbl,
      Met_seen);
  else if(transinitmeth==BAYESX)
    start_codon=transinit_bayes(strata_ct,medoids,ATGparms,ATGparmsPWM,
      coords_orig,sequence,seqlen,coords->start,bestTISscore,&nouse_int,
      Met_seen);
  else if(transinitmeth==MFCWLLKRX||transinitmeth==NFCWLLKRX)
    start_codon=transinit_fcwllkr(transinitmeth,strata_ct,medoids,
      ATGparms,ATGparmsPWM,fcwllkrparms,coords_orig,sequence,seqlen,
      coords->start,bestTISscore,&nouse_int,Met_seen);
  else
    FATALERROR("Err (searchforbestinit): Bad algorithm specified.\n")

  free(coords_orig);
  return(start_codon);
} /* end searchforbestinit */

/* Finds the most likely translation initiation site in a CDS   *
 * and returns this value using a one-based coordinate system.  *
 * If no "good" site exists, NOTISPRED is returned.  If the     *
 * coding sequence is not a multiple of 3, NOTMULT3 is returned */
int searchforbestinitincds(
  int transinitmeth,
  int strata_ct,
  char **medoids,
  ATGparm *ATGparms,
  ATGparmPWM *ATGparmsPWM,
  fcwllkrauxT *fcwllkrparms,
  int *sequence,
  int seqlen
) {
  int prependlen,         /* N-terminus nt flank         */
      appendlen,          /* C-terminus nt flank         */
      start_codon=0,      /* Predicted TIS               */
      *coords_orig=NULL,  /* mapping derived coordinates *
                           * to the original ones        */
      i,                  /* iterator variables          */
      Met_seen,           /* was a TIS predicted?        *
                           * prediction routines         */
      best_is_true_p,     /* BAYES, was MAP hypo true?   */
      fc_verdict;         /* ?FCWLLKR, was TIS pred?     */
  long double
    bestTISscore,         /* "Real" score of best TIS    */
    bestTISscorenoweight; /* WLLKR, non-weighted score   */ 

  /* Any CDS passed to this routine MUST have upstream *
   * and downstrextent amendments: 96 5 CDS 3 96       */
  prependlen=(int)(UPSTREXTENT)+(int)(CONTENTSWATHLEN);
  appendlen=(int)(DOWNSTREXTENT)+(int)(CONTENTSWATHLEN);

  /* Verify that the CDS length is a multiple of three */
  if((seqlen-prependlen-appendlen)%3!=0) {
    (void)snprintf(errbuff,MAXERRLEN,"Warn (searchforbestinitincds): \
CDS length (%int) is not a multiple of 3!??\n",
      seqlen-prependlen-appendlen);
    NONFATALERROR(errbuff)
    return(NOTMULT3);
  }

  /* Build an array storing original CDS coordinates */
  if((coords_orig=(int*)malloc(sizeof(int)*seqlen))==NULL)
    FATALERROR("Err (searchforbestinitincds): Insufficient Memory!\n")
  for(i=0;i<seqlen;++i)
    coords_orig[i]=i+1;

  if(transinitmeth==LLKRX) {
    start_codon=transinit_llkr(
      strata_ct,medoids,ATGparms,ATGparmsPWM,
      coords_orig+(int)(CONTENTSWATHLEN),
      sequence+(int)(CONTENTSWATHLEN),
      seqlen-(2*(int)(CONTENTSWATHLEN)),
      prependlen+1,&bestTISscore,&Met_seen);
    if(Met_seen==FALSE||bestTISscore<0.0)
      start_codon=NOTISPRED;
  }
  else if(transinitmeth==WLLKRX) {
    start_codon=transinit_wllkr(
      strata_ct,medoids,ATGparms,ATGparmsPWM,
      coords_orig+(int)(CONTENTSWATHLEN),
      sequence+(int)(CONTENTSWATHLEN),
      seqlen-(2*(int)(CONTENTSWATHLEN)),
      prependlen+1,&bestTISscore,
      &bestTISscorenoweight,&Met_seen);
    if(Met_seen==FALSE||bestTISscorenoweight<0.0)
      start_codon=NOTISPRED;
  }
  else if(transinitmeth==BAYESX) {
    start_codon=transinit_bayes(
      strata_ct,medoids,ATGparms,ATGparmsPWM,
      coords_orig+(int)(CONTENTSWATHLEN),
      sequence+(int)(CONTENTSWATHLEN),
      seqlen-(2*(int)(CONTENTSWATHLEN)),
      prependlen+1,&bestTISscore,
      &best_is_true_p,&Met_seen);
    if(Met_seen==FALSE||best_is_true_p==FALSE)
      start_codon=NOTISPRED;
  }
  else if(transinitmeth==MFCWLLKRX||transinitmeth==NFCWLLKRX) {
    start_codon=transinit_fcwllkr(transinitmeth,strata_ct,medoids,
      ATGparms,ATGparmsPWM,fcwllkrparms,coords_orig,sequence,seqlen,
      prependlen,&bestTISscore,&fc_verdict,&Met_seen);
    if(Met_seen==FALSE||fc_verdict==FALSE)
      start_codon=NOTISPRED;
  }
  else
    FATALERROR("Err (searchforbestinitincds): Bad algorithm specified.\n")

  free(coords_orig);
  return(start_codon);
} /* end searchforbestinitincds */

/* Predicts the best translation initiation site on the basis of *
 * log-likelihood ratios and hard constraints on protein length. */
static int transinit_llkr(
  int strata_ct,
  char **medoids,
  ATGparm *ATGparms,
  ATGparmPWM *ATGparmsPWM,
  int *coords_orig,
  int *sequence,
  int seqlen,
  int start_orig,
  long double *bestTISscore,
  int *Met_seen
) {
  char *currmet=NULL;   /* Stores sequence of init site       */
  double shortdist,     /* distance of closest medoid         */
         currdist;      /* used to compute shortdist          */
  int bestind=-1,       /* tracks best start site yet seen    */
      Met_seen_p=FALSE, /* Has candidate start-Met been seen? */
      start_codon,      /* Index to return to caller          */
      medoidlen,        /* length of sites in medoids         */
      closemedoidind=0, /* index of closest medoid            */
      i,j;              /* iterator variables                 */
  long double
    loglikerat,    /* log-likelihood ratios       */
    bestscore=0.0; /* score of best site yet seen */

  /* Find the best start-Methionine */
  for(i=UPSTREXTENT;i<seqlen-(DOWNSTREXTENT)-2;i+=3)
    if(IS_MET_P(sequence,i)==TRUE&&
       PEPTIDE_LENGTH_CONSTRAINTS_MET_P(seqlen-i-(DOWNSTREXTENT))==TRUE) {
      if(MET_VALID_TO_SCORE_P(i,0,seqlen-1)==FALSE) /* sanity check */
        FATALERROR("Err (transinit_llkr): Met indexed incorrectly!??\n")

      if(strata_ct==1)
        loglikerat=calc_llkr(ATGparms[0],ATGparms[strata_ct],
                             &sequence[i-(UPSTREXTENT)]);
      else { /* stratified approach */
        medoidlen=strlen(medoids[0]);
        if((currmet=(char*)malloc(sizeof(char)*(medoidlen+1)))==NULL)
          FATALERROR("Err (transinit_llkr): Insufficient memory.\n")

        /* find the closest of the medoids */
        for(j=0;j<medoidlen;++j)
          RECOVER_BASE(sequence[i-(UPSTREXTENT)+j],currmet[j])
        currmet[medoidlen]='\0';

        /* Determine the correct parameters to score the site with */
        switch(cluster_parm_selection) {
          case HAMSTATX :
            if(Met_seen_p==TRUE)
              break;
          case HAMMODX :
            shortdist=eval_dist_kernel(medoids[0],currmet);
            closemedoidind=0;
            for(j=1;j<strata_ct;++j) {
              currdist=eval_dist_kernel(medoids[j],currmet);
              if(currdist<shortdist) {
                shortdist=currdist;
                closemedoidind=j;
              }
            }
            break;
          case PWMSTATX :
            if(Met_seen_p==TRUE)
              break;
          case PWMMODX :
            shortdist=eval_pwm_score(ATGparmsPWM[0],currmet);
            closemedoidind=0;
            for(j=1;j<strata_ct;++j) {
              currdist=eval_pwm_score(ATGparmsPWM[j],currmet);
              if(currdist>shortdist) {
                shortdist=currdist;
                closemedoidind=j;
              }
            }
            break;
          case WAMSTATX :
            if(Met_seen_p==TRUE)
              break;
          case WAMMODX :
            shortdist=calc_llk(ATGparms[0],&sequence[i-(UPSTREXTENT)]);
            closemedoidind=0;
            for(j=1;j<strata_ct;++j) {
              currdist=calc_llk(ATGparms[j],&sequence[i-(UPSTREXTENT)]);
              if(currdist>shortdist) {
                shortdist=currdist;
                closemedoidind=j;
              }
            }
            break;
          default :
            fprintf(stderr,"Err (transinit_llkr): \
Improper cluster deployment strategy requested!\n");
            return(EXIT_FAILURE);
        } /* end switch */
        #ifdef SHOWCLOSEMEDOID
        fprintf(stdout,"=%i %i\n",i+(CONTENTSWATHLEN)+1,closemedoidind);
        #endif

        loglikerat=calc_llkr(ATGparms[closemedoidind],ATGparms[strata_ct],
                             &sequence[i-(UPSTREXTENT)]);
        free(currmet);
      }

      if((Met_seen_p==TRUE&&loglikerat>bestscore)||Met_seen_p==FALSE) {
        bestscore=loglikerat;
        bestind=i;
        if(Met_seen_p==FALSE)
          Met_seen_p=TRUE;
      }
    }

  /* We require the predicted start-Met's log-likelihood *
   * ratio to favor its actually being a true site!      */
  if(Met_seen_p==FALSE||bestscore<0.0)
    start_codon=start_orig;
  else
    start_codon=coords_orig[bestind];

  /* If no ATG was seen, we still report a null score. */
  *bestTISscore=bestscore;
  *Met_seen=Met_seen_p;

  return(start_codon);
} /* end transinit_llkr */

/* Predicts the best translation initiation site on the basis *
 * of weighted log-likelihood ratios and hard constraints on  *
 * protein length.                                            */
static int transinit_wllkr(
  int strata_ct,
  char **medoids,
  ATGparm *ATGparms,
  ATGparmPWM *ATGparmsPWM,
  int *coords_orig,
  int *sequence,
  int seqlen,
  int start_orig,
  long double *bestTISscore,
  long double *bestTISscorenoweight,
  int *Met_seen
) {
  char *currmet=NULL;   /* Stores sequence of init site       */
  double shortdist,     /* distance of closest medoid         */
         currdist;      /* used to compute shortdist          */
  int bestind=-1,       /* tracks best start site yet seen    */
      Met_seen_p=FALSE, /* Has candidate start-Met been seen? */
      start_codon,      /* Index to return to caller          */
      medoidlen,        /* length of sites in medoids         */
      closemedoidind=0, /* index of closest medoid            */
      cdslength=0,      /* stores length of the cds proper    */
      i,j;              /* iterator variables                 */
  long double
    loglikerat,     /* log-likelihood ratios         */
    llkrweight,     /* weight for loglikerat         */
    wloglikerat,    /* weighted log-likelihood ratio */
    bestwscore=0.0, /* best weighted score yet seen  */
    bestscore=0.0;  /* score of best site yet seen   */

  /* Find the best start-Methionine */
  cdslength=seqlen-(UPSTREXTENT)-(DOWNSTREXTENT);
  #ifdef DBGWLLKR
  fprintf(stderr,"length(maximal_reading_frame)=%i\n",cdslength);
  #endif
  for(i=UPSTREXTENT;i<seqlen-(DOWNSTREXTENT)-2;i+=3)
    if(IS_MET_P(sequence,i)==TRUE&&
       PEPTIDE_LENGTH_CONSTRAINTS_MET_P(seqlen-i-(DOWNSTREXTENT))==TRUE) {
      if(MET_VALID_TO_SCORE_P(i,0,seqlen-1)==FALSE) /* sanity check */
        FATALERROR("Err (transinit_wllkr): Met indexed incorrectly!??\n")

      if(strata_ct==1)
        loglikerat=calc_llkr(ATGparms[0],ATGparms[strata_ct],
                             &sequence[i-(UPSTREXTENT)]);
      else { /* stratified approach */
        medoidlen=strlen(medoids[0]);
        if((currmet=(char*)malloc(sizeof(char)*(medoidlen+1)))==NULL)
          FATALERROR("Err (transinit_wllkr): Insufficient memory.\n")

        /* find the closest of the medoids */
        for(j=0;j<medoidlen;++j)
          RECOVER_BASE(sequence[i-(UPSTREXTENT)+j],currmet[j])
        currmet[medoidlen]='\0';

        /* Determine the correct parameters to score the site with */
        switch(cluster_parm_selection) {
          case HAMSTATX :
            if(Met_seen_p==TRUE)
              break;
          case HAMMODX :
            shortdist=eval_dist_kernel(medoids[0],currmet);
            closemedoidind=0;
            for(j=1;j<strata_ct;++j) {
              currdist=eval_dist_kernel(medoids[j],currmet);
              if(currdist<shortdist) {
                shortdist=currdist;
                closemedoidind=j;
              }
            }
            break;
          case PWMSTATX :
            if(Met_seen_p==TRUE)
              break;
          case PWMMODX :
            shortdist=eval_pwm_score(ATGparmsPWM[0],currmet);
            closemedoidind=0;
            for(j=1;j<strata_ct;++j) {
              currdist=eval_pwm_score(ATGparmsPWM[j],currmet);
              if(currdist>shortdist) {
                shortdist=currdist;
                closemedoidind=j;
              }
            }
            break;
          case WAMSTATX :
            if(Met_seen_p==TRUE)
              break;
          case WAMMODX :
            shortdist=calc_llk(ATGparms[0],&sequence[i-(UPSTREXTENT)]);
            closemedoidind=0;
            for(j=1;j<strata_ct;++j) {
              currdist=calc_llk(ATGparms[j],&sequence[i-(UPSTREXTENT)]);
              if(currdist>shortdist) {
                shortdist=currdist;
                closemedoidind=j;
              }
            }
            break;
          default :
            fprintf(stderr,"Err (transinit_wllkr): \
Improper cluster deployment strategy requested!\n");
            return(EXIT_FAILURE);
        } /* end switch */
        #ifdef SHOWCLOSEMEDOID
        fprintf(stdout,"=%i %i\n",i+(CONTENTSWATHLEN)+1,closemedoidind);
        #endif

        loglikerat=calc_llkr(ATGparms[closemedoidind],ATGparms[strata_ct],
                             &sequence[i-(UPSTREXTENT)]);
        free(currmet);
      }

      /* Compute weights that empirically approximate our beliefs as to *
       * where in the maximal reading frame a start-Met should occur.   */
      llkrweight=pow(((double)(cdslength-(i-(UPSTREXTENT))))/cdslength,3);
      wloglikerat=log(llkrweight)+loglikerat;
      if(isinf(wloglikerat)||isnan(wloglikerat)) {
        NONFATALERROR("Warn (transinit_wllkr): Buffer underflow\n")
        continue;
      }
      #ifdef DBGWLLKR
      fprintf(stderr,"position(Met)=%-4i, relative_length(protein)=%.2f, \
Weight=%Lf, LLKR=%+Lf, WLLKR=%+Lf\n",
        i-(UPSTREXTENT)+1,((double)cdslength-(i-(UPSTREXTENT)))/cdslength,
        llkrweight,loglikerat,wloglikerat);
      #endif
      if((Met_seen_p==TRUE&&wloglikerat>bestwscore)||Met_seen_p==FALSE) {
        bestwscore=wloglikerat;
        bestscore=loglikerat;
        bestind=i;
        if(Met_seen_p==FALSE)
          Met_seen_p=TRUE;
      }
    }
  #ifdef DBGWLLKR
  if(Met_seen_p==FALSE||bestscore<0.0)
    fprintf(stderr,"position(WLLKR_Met)=DNE\n\n");
  else
    fprintf(stderr,"position(WLLKR_Met)=%i\n\n",
      bestind-(UPSTREXTENT)+1);
  #endif

  /* We require the predicted start-Met's log-likelihood *
   * ratio to favor its actually being a true site!      */
  if(Met_seen_p==FALSE||bestscore<0.0)
    start_codon=start_orig;
  else
    start_codon=coords_orig[bestind];

  /* If no ATG was seen, we still report a null score. */
  *bestTISscore=bestwscore;
  *bestTISscorenoweight=bestscore;
  *Met_seen=Met_seen_p;

  return(start_codon);
} /* end transinit_wllkr */

/* Predicts the best translation initiation site using a     *
 * Bayesian approach and hard constraints on protein length. */
static int transinit_bayes(
  int strata_ct,
  char **medoids,
  ATGparm *ATGparms,
  ATGparmPWM *ATGparmsPWM,
  int *coords_orig,
  int *sequence,
  int seqlen,
  int start_orig,
  long double *bestTISscore,
  int *best_is_true,
  int *Met_seen
) {
  char *currmet=NULL;   /* Stores sequence of init site   */
  double shortdist,     /* distance of closest medoid     */
         currdist;      /* used to compute shortdist      */
  int metct,            /* Count of in-frame Methionines  */
      *metpos=NULL,     /* Array of start-Met positions   */
      bestind,          /* for finding MAP hypothesis     */
      start_codon,      /* Index to return to caller      */
      medoidlen,        /* length of sites in medoids     */
      closemedoidind=0, /* index of closest medoid        */
      cdslength=0,      /* stores length of the cds       */
      i,j,k;            /* iterator variables             */
  long double
    *metpriors=NULL,    /* Stores prior distribution      */
    normconst,          /* Prior normalization constant   */
    probsum,            /* Verify prior dist is a PMF     */
    *metllks=NULL,      /* Stores likelihoods of Mets     */
    *numerator=NULL,    /* Stores numerator of Bayes rule */
    bestscore=0.0;      /* for finding MAP hypothesis     */

  /* First, learn the count of in-frame Mets */
  for(i=UPSTREXTENT,metct=0;i<seqlen-(DOWNSTREXTENT)-2;i+=3)
    if(IS_MET_P(sequence,i)==TRUE&&
       PEPTIDE_LENGTH_CONSTRAINTS_MET_P(seqlen-i-(DOWNSTREXTENT))==TRUE) {
      if(MET_VALID_TO_SCORE_P(i,0,seqlen-1)==FALSE) /* sanity check */
        FATALERROR("Err (transinit_bayes): Met indexed incorrectly!??\n")
      else
        ++metct;
    }
  if(metct==0) {
    #ifdef DBGBAYES
    fprintf(stderr,"position(MAP_hypo)=DNE (no viable Met's present)\n");
    #endif

    /* If no ATG was seen, we still report a null score. */
    *bestTISscore=bestscore;
    *Met_seen=FALSE;
    *best_is_true=FALSE;

    return(start_orig);
  }

  /* Allocate arrays for Met positions, priors, *
   * log-likelihoods, and numerators...         */
  if((metpos=(int*)malloc(sizeof(int)*metct))==NULL||
     (metpriors=(long double*)malloc(sizeof(long double)*(metct<<1)))==NULL||
     (metllks=(long double*)malloc(sizeof(long double)*(metct<<1)))==NULL||
     (numerator=(long double*)malloc(sizeof(long double)*(metct<<1)))==NULL)
    FATALERROR("Err (transinit_bayes): Insufficient memory.\n")

  /* ...and populate all (but the numerators) accordingly. */
  cdslength=seqlen-(UPSTREXTENT)-(DOWNSTREXTENT);
  for(i=UPSTREXTENT,j=0;i<seqlen-(DOWNSTREXTENT)-2;i+=3)
    if(IS_MET_P(sequence,i)==TRUE&&
       PEPTIDE_LENGTH_CONSTRAINTS_MET_P(seqlen-i-(DOWNSTREXTENT))==TRUE) {
      if(MET_VALID_TO_SCORE_P(i,0,seqlen-1)==FALSE) /* sanity check */
        FATALERROR("Err (transinit_bayes): Met indexed incorrectly!??\n")

      metpos[j]=i;
      metpriors[j]=
        (cdslength-(metpos[j]-(UPSTREXTENT)))/(long double)cdslength;
      metpriors[j+metct]=1.0-metpriors[j];

      if(strata_ct==1) {
        /* site under true hypothesis */
        metllks[j]=calc_llk(ATGparms[0],&sequence[i-(UPSTREXTENT)]);
        /* site under false hypothesis */
        metllks[(j++)+metct]=calc_llk(ATGparms[strata_ct],
                                    &sequence[i-(UPSTREXTENT)]);
      }
      else { /* stratified approach */
        medoidlen=strlen(medoids[0]);
        if((currmet=(char*)malloc(sizeof(char)*(medoidlen+1)))==NULL)
          FATALERROR("Err (transinit_bayes): Insufficient memory.\n")

        /* find the closest of the medoids */
        for(k=0;k<medoidlen;++k)
          RECOVER_BASE(sequence[i-(UPSTREXTENT)+k],currmet[k])
        currmet[medoidlen]='\0';

        /* Determine the correct parameters to score the site with */
        switch(cluster_parm_selection) {
          case HAMSTATX :
            if(j>0)
              break;
          case HAMMODX :
            shortdist=eval_dist_kernel(medoids[0],currmet);
            closemedoidind=0;
            for(k=1;k<strata_ct;++k) {
              currdist=eval_dist_kernel(medoids[k],currmet);
              if(currdist<shortdist) {
                shortdist=currdist;
                closemedoidind=k;
              }
            }
            break;
          case PWMSTATX :
            if(j>0)
              break;
          case PWMMODX :
            shortdist=eval_pwm_score(ATGparmsPWM[0],currmet);
            closemedoidind=0;
            for(k=1;k<strata_ct;++k) {
              currdist=eval_pwm_score(ATGparmsPWM[k],currmet);
              if(currdist>shortdist) {
                shortdist=currdist;
                closemedoidind=k;
              }
            }
            break;
          case WAMSTATX :
            if(j>0)
              break;
          case WAMMODX :
            shortdist=calc_llk(ATGparms[0],&sequence[i-(UPSTREXTENT)]);
            closemedoidind=0;
            for(k=1;k<strata_ct;++k) {
              currdist=calc_llk(ATGparms[k],&sequence[i-(UPSTREXTENT)]);
              if(currdist>shortdist) {
                shortdist=currdist;
                closemedoidind=k;
              }
            }
            break;
          default :
            fprintf(stderr,"Err (transinit_bayes): \
Improper cluster deployment strategy requested!\n");
            return(EXIT_FAILURE);
        } /* end switch */
        #ifdef SHOWCLOSEMEDOID
        fprintf(stdout,"=%i %i\n",i+(CONTENTSWATHLEN)+1,closemedoidind);
        #endif

        /* site under true hypothesis */
        metllks[j]=calc_llk(ATGparms[closemedoidind],
                            &sequence[i-(UPSTREXTENT)]);
        /* site under false hypothesis */
        metllks[(j++)+metct]=calc_llk(ATGparms[strata_ct],
                                    &sequence[i-(UPSTREXTENT)]);

        free(currmet);
      }
    }
  if(j!=metct) /* sanity check */
    FATALERROR("Err (transinit_bayes): Inconsistency detected!??\n")

  /* The prior array must be transformed into a valid PMF */
  for(i=0,normconst=0.0;i<metct<<1;++i)
    normconst+=metpriors[i];
  for(i=0,probsum=0.0;i<metct<<1;++i) {
    metpriors[i]/=normconst;
    probsum+=metpriors[i];
  }
  if(fabs(probsum-1.0)>MAXFAULT)
    FATALERROR("Err (transinit_bayes): Priors are not a PMF!\n")

  /* Because identification of the MAP hypothesis can be *
   * done by finding argmax_h of the numerator of Bayes  *
   * rule, derivation of the full posterior distribution *
   * is unnecessary.                                     */
  for(i=0;i<metct<<1;++i)
    numerator[i]=log(metpriors[i])+metllks[i];

  #ifdef DBGBAYES
  fprintf(stderr,"length(maximal_reading_frame)=%i\n",cdslength);
  fprintf(stderr,"position(Met)\trel_len(protein)\t\
Prior_t\tLLK_t\tlog(Prior_t)+LLK_t\t\
Prior_nil\tLLK_nil\tlog(Prior_nil)+LLK_nil\n");
  for(i=0;i<metct;++i)
    fprintf(stderr,"%-4i\t%.2Lf\t%.2Lf\t%+Lf\t%+Lf\t%.2Lf\t%+Lf\t%+Lf\n",
      metpos[i]-(UPSTREXTENT)+1,
      (cdslength-(metpos[i]-(UPSTREXTENT)))/(long double)cdslength,
      metpriors[i],metllks[i],numerator[i],
      metpriors[i+metct],metllks[i+metct],numerator[i+metct]);
  #endif

  /* Determine the MAP hypothesis */
  bestscore=numerator[0];
  bestind=0;
  for(i=1;i<metct<<1;++i)
    if(numerator[i]>bestscore) {
      bestscore=numerator[i];
      bestind=i;
    }
  #ifdef DBGBAYES
  if(bestind<metct)
    fprintf(stderr,"position(MAP_Met)=%i t\n\n",
      metpos[bestind]-(UPSTREXTENT)+1);
  else
    fprintf(stderr,"position(MAP_Met)=%i nil (DNE)\n\n",
      metpos[bestind-metct]-(UPSTREXTENT)+1);
  #endif

  /* Identify the best start codon, which must be favored *
   * as actually being a true start-Methionine!           */
  if(bestind<metct) {
    start_codon=coords_orig[metpos[bestind]];
    *best_is_true=TRUE;
  }
  else {
    start_codon=start_orig;
    *best_is_true=FALSE;
  }

  *bestTISscore=bestscore;
  *Met_seen=TRUE;

  free(numerator);
  free(metllks);
  free(metpriors);
  free(metpos);
  return(start_codon);
} /* end transinit_bayes */

/* Predicts the best translation initiation site on the basis *
 * of a combination of weighted log-likelihood ratios and     *
 * ratios of downstream to upstream sequence log-likelihoods  *
 * under a coding hypothesis, and hard constraints on protein *
 * length.                                                    */
static int transinit_fcwllkr(
  int transinitmeth,
  int strata_ct,
  char **medoids,
  ATGparm *ATGparms,
  ATGparmPWM *ATGparmsPWM,
  fcwllkrauxT *fcwllkrparms,
  int *coords_orig,
  int *sequence,
  int seqlen,
  int start_orig,
  long double *bestTISscore,
  int *fc_class,
  int *Met_seen
) {
  char *currmet=NULL,           /* Stores sequence of init site       */
    upseq[CONTENTSWATHLEN+1],   /* Seq to run MC over                 */
    downseq[CONTENTSWATHLEN+1]; /* Seq to run MC over                 */
  double shortdist,             /* distance of closest medoid         */
         currdist,              /* used to compute shortdist          */
         neuronscore;           /* fct(dot_prod(weights,features))    */
  featvecT features;            /* represenation of instance to score */
  int bestind=-1,               /* tracks best start site yet seen    */
      Met_seen_p=FALSE,         /* Has candidate start-Met been seen? */
      start_codon,              /* Index to return to caller          */
      verdict,                  /* records classification result      */
      bestclass=FALSE,          /* verdict of best bestnscore seen    */
      medoidlen,                /* length of sites in medoids         */
      closemedoidind=0,         /* index of closest medoid            */
      cdslength=0,              /* stores length of the cds           */
      i,j,k;                    /* iterator variables                 */
  long double
    loglikerat,                 /* log-likelihood ratios              */
    llkrweight,                 /* weight for loglikerat              */
    wloglikerat,                /* weighted log-likelihood ratio      */
    bestnscore=0.0,             /* best neuronal score yet seen       */
    multscore,                  /* wloglikerat+downscore-upscore      */
    bestmscore=0.0,             /* best multiplicative score yet seen */
    upscore,                    /* posterior prob of upstream seq     */
    downscore,                  /* posterior prob of downstream seq   */
    updenom,                    /* Pr(Data) for upstream seq          */
    downdenom;                  /* Pr(Data) for downstream seq        */

  /* Find the best start-Methionine */
  cdslength=seqlen-(UPSTREXTENT)-(DOWNSTREXTENT)-2*(CONTENTSWATHLEN);
  #ifdef DBGFCWLLKR
  fprintf(stderr,"length(maximal_reading_frame)=%i\n",cdslength);
  #endif
  for(i=(UPSTREXTENT)+(CONTENTSWATHLEN);
      i<seqlen-(DOWNSTREXTENT)-(CONTENTSWATHLEN)-2;
      i+=3)
    if(IS_MET_P(sequence,i)==TRUE&&
       PEPTIDE_LENGTH_CONSTRAINTS_MET_P(seqlen-i-
         (DOWNSTREXTENT)-(CONTENTSWATHLEN))==TRUE) {
      if(MET_VALID_TO_SCORE_P(i,(CONTENTSWATHLEN),seqlen-1)==FALSE||
         MET_FLANKS_VALID_TO_SCORE_P(i,0,seqlen-1)==FALSE) /* sanity check */
        FATALERROR("Err (transinit_fcwllkr): Met indexed incorrectly!??\n")

      if(strata_ct==1)
        loglikerat=calc_llkr(ATGparms[0],ATGparms[strata_ct],
                             &sequence[i-(UPSTREXTENT)]);
      else { /* stratified approach */
        medoidlen=strlen(medoids[0]);
        if((currmet=(char*)malloc(sizeof(char)*(medoidlen+1)))==NULL)
          FATALERROR("Err (transinit_fcwllkr): Insufficient memory.\n")

        /* find the closest of the medoids */
        for(j=0;j<medoidlen;++j)
          RECOVER_BASE(sequence[i-(UPSTREXTENT)+j],currmet[j])
        currmet[medoidlen]='\0';

        /* Determine the correct parameters to score the site with */
        switch(cluster_parm_selection) {
          case HAMSTATX :
            if(Met_seen_p==TRUE)
              break;
          case HAMMODX :
            shortdist=eval_dist_kernel(medoids[0],currmet);
            closemedoidind=0;
            for(j=1;j<strata_ct;++j) {
              currdist=eval_dist_kernel(medoids[j],currmet);
              if(currdist<shortdist) {
                shortdist=currdist;
                closemedoidind=j;
              }
            }
            break;
          case PWMSTATX :
            if(Met_seen_p==TRUE)
              break;
          case PWMMODX :
            shortdist=eval_pwm_score(ATGparmsPWM[0],currmet);
            closemedoidind=0;
            for(j=1;j<strata_ct;++j) {
              currdist=eval_pwm_score(ATGparmsPWM[j],currmet);
              if(currdist>shortdist) {
                shortdist=currdist;
                closemedoidind=j;
              }
            }
            break;
          case WAMSTATX :
            if(Met_seen_p==TRUE)
              break;
          case WAMMODX :
            shortdist=calc_llk(ATGparms[0],&sequence[i-(UPSTREXTENT)]);
            closemedoidind=0;
            for(j=1;j<strata_ct;++j) {
              currdist=calc_llk(ATGparms[j],&sequence[i-(UPSTREXTENT)]);
              if(currdist>shortdist) {
                shortdist=currdist;
                closemedoidind=j;
              }
            }
            break;
          default :
            fprintf(stderr,"Err (transinit_fcllkr): \
Improper cluster deployment strategy requested!\n");
            return(EXIT_FAILURE);
        } /* end switch */
        #ifdef SHOWCLOSEMEDOID
        fprintf(stdout,"=%i %i\n",i+1,closemedoidind);
        #endif

        loglikerat=calc_llkr(ATGparms[closemedoidind],ATGparms[strata_ct],
                             &sequence[i-(UPSTREXTENT)]);
        free(currmet);
      }

      /* Compute weights that empirically approximate our beliefs as to *
       * where in the maximal reading frame a start-Met should occur.   */
      llkrweight=pow(
        ((double)(cdslength-(i-(UPSTREXTENT)-(CONTENTSWATHLEN))))/cdslength,3);
      wloglikerat=log(llkrweight)+loglikerat;
      if(isinf(wloglikerat)||isnan(wloglikerat)) {
        NONFATALERROR("Warn (transinit_fcwllkr): Buffer underflow\n")
        continue;
      }

      /* upstream flank (content) */
      for(j=i-(UPSTREXTENT)-(CONTENTSWATHLEN),k=0;k<CONTENTSWATHLEN;)
        RECOVER_BASE(sequence[j++],upseq[k++])
      upseq[CONTENTSWATHLEN]='\0';

      /* downstream flank (content) */
      for(j=i+3+(DOWNSTREXTENT),k=0;k<CONTENTSWATHLEN;)
        RECOVER_BASE(sequence[j++],downseq[k++])
      downseq[CONTENTSWATHLEN]='\0';

      /* compute likelihood ratio of flanks under coding hypothesis */
      updenom=downdenom=upscore=downscore=0.0;
      switch(fcwllkrparms->markovchainmeth) {
        case FIXORDX :
          for(j=0;j<3;++j) {
            updenom+=CODINGPRIOR*
              exp(FOprob(upseq,fcwllkrparms->markovchainparms,
                         j,TRUE,FALSE));
            downdenom+=CODINGPRIOR*
              exp(FOprob(downseq,fcwllkrparms->markovchainparms,
                         j,TRUE,FALSE));
          }
          updenom+=NONCODPRIOR*
            exp(FOprob(upseq,fcwllkrparms->markovchainparms,
                       NONCODING,TRUE,FALSE));
          downdenom+=NONCODPRIOR*
            exp(FOprob(downseq,fcwllkrparms->markovchainparms,
                       NONCODING,TRUE,FALSE));
          for(j=0;j<3;++j) {
            upscore+=CODINGPRIOR*
              exp(FOprob(upseq,fcwllkrparms->markovchainparms,
                         j,TRUE,FALSE))/updenom;
            downscore+=CODINGPRIOR*
              exp(FOprob(downseq,fcwllkrparms->markovchainparms,
                         j,TRUE,FALSE))/downdenom;
          }
          break;
        case TDDIX :
        case BUDIX :
        case CHI2X :
          for(j=0;j<3;++j) {
            updenom+=CODINGPRIOR*
              exp(IMMprob(upseq,fcwllkrparms->markovchainparms,
                          j,TRUE,FALSE));
            downdenom+=CODINGPRIOR*
              exp(IMMprob(downseq,fcwllkrparms->markovchainparms,
                          j,TRUE,FALSE));
          }
          updenom+=NONCODPRIOR*
            exp(IMMprob(upseq,fcwllkrparms->markovchainparms,
                        NONCODING,TRUE,FALSE));
          downdenom+=NONCODPRIOR*
            exp(IMMprob(downseq,fcwllkrparms->markovchainparms,
                        NONCODING,TRUE,FALSE));
          for(j=0;j<3;++j) {
            upscore+=CODINGPRIOR*
              exp(IMMprob(upseq,fcwllkrparms->markovchainparms,
                          j,TRUE,FALSE))/updenom;
            downscore+=CODINGPRIOR*
              exp(IMMprob(downseq,fcwllkrparms->markovchainparms,
                          j,TRUE,FALSE))/downdenom;
          }
          break;
        case DMMMX :
          for(j=0;j<3;++j) {
            updenom+=CODINGPRIOR*
              exp(DMMMprob(upseq,fcwllkrparms->markovchainparms,
                           fcwllkrparms->dmcloudindex,j,TRUE,FALSE));
            downdenom+=CODINGPRIOR*
              exp(DMMMprob(downseq,fcwllkrparms->markovchainparms,
                           fcwllkrparms->dmcloudindex,j,TRUE,FALSE));
          }
          updenom+=NONCODPRIOR*
            exp(DMMMprob(upseq,fcwllkrparms->markovchainparms,
                         fcwllkrparms->dmcloudindex,NONCODING,TRUE,FALSE));
          downdenom+=NONCODPRIOR*
            exp(DMMMprob(downseq,fcwllkrparms->markovchainparms,
                         fcwllkrparms->dmcloudindex,NONCODING,TRUE,FALSE));
          for(j=0;j<3;++j) {
            upscore+=CODINGPRIOR*
              exp(DMMMprob(upseq,fcwllkrparms->markovchainparms,
                           fcwllkrparms->dmcloudindex,j,TRUE,FALSE))/
              updenom;
            downscore+=CODINGPRIOR*
              exp(DMMMprob(downseq,fcwllkrparms->markovchainparms,
                           fcwllkrparms->dmcloudindex,j,TRUE,FALSE))/
              downdenom;
          }
          break;
        case PMMMX :
        default :
          fprintf(stderr,"Err (transinit_fcwllkr): \
Improper algorithm requested!\n");
          return(EXIT_FAILURE);
      } /* end switch */

      /* update the featvecT structure */
      features.metwam=wloglikerat;
      features.mcratio=log(downscore/upscore);

      if(transinitmeth==MFCWLLKRX) {
        multscore=features.metwam+features.mcratio;
        #ifdef DBGFCWLLKR
        fprintf(stderr,"position(Met)=%-4i, relative_length(protein)=%.2f, \
Weight=%Lf, LLKR=%+Lf, WLLKR=%+Lf, MCRATIO=%+.4f, log(score)=%+Lf\n",
          i-(UPSTREXTENT)-(CONTENTSWATHLEN)+1,
          ((double)cdslength-(i-(UPSTREXTENT)-(CONTENTSWATHLEN)))/cdslength,
          llkrweight,loglikerat,wloglikerat,
          features.mcratio,multscore);
        #endif
        if((Met_seen_p==TRUE&&multscore>bestmscore)||Met_seen_p==FALSE) {
          bestmscore=multscore;
          bestind=i;
          if(Met_seen_p==FALSE)
            Met_seen_p=TRUE;
        }
      }
      else if(transinitmeth==NFCWLLKRX) {
        /* invoke the threshold logic unit */
        if(strata_ct==1)
          eval_activation_unit(features,
            *fcwllkrparms->localneuron,&verdict,
            &neuronscore,fcwllkrparms->neurontype);
        else
          eval_activation_unit(features,
            fcwllkrparms->localneuron[closemedoidind],
            &verdict,&neuronscore,fcwllkrparms->neurontype);
        #ifdef DBGFCWLLKR
        fprintf(stderr,"position(Met)=%-4i, relative_length(protein)=%.2f, \
Weight=%Lf, LLKR=%+Lf, WLLKR=%+Lf, MCRATIO=%+.4f, \
neuron_score=%.4f, status=%i\n",
          i-(UPSTREXTENT)-(CONTENTSWATHLEN)+1,
          ((double)cdslength-(i-(UPSTREXTENT)-(CONTENTSWATHLEN)))/cdslength,
          llkrweight,loglikerat,wloglikerat,
          features.mcratio,neuronscore,verdict);
        #endif
        if((Met_seen_p==TRUE&&neuronscore>bestnscore)||Met_seen_p==FALSE) {
          bestnscore=neuronscore;
          bestclass=verdict;
          bestind=i;
          if(Met_seen_p==FALSE)
            Met_seen_p=TRUE;
        }
      }
      else {
        (void)snprintf(errbuff,MAXERRLEN,"Err (transinit_fcwllkr): \
An unrecognized method was specified (indexed as %i)!\n",transinitmeth);
        FATALERROR(errbuff)
      }
    } /* end if */

  if(transinitmeth==MFCWLLKRX) {
    /* We require the predicted start-Met's score  *
     * ratio to favor both its actually being a    *
     * true site and the downstream sequence swath *
     * to have greater coding likelihood than that *
     * upstream.                                   */
    if(Met_seen_p==FALSE||bestmscore<MFCTHRESH) {
      start_codon=start_orig;
      *fc_class=FALSE;
    }
    else {
      start_codon=coords_orig[bestind];
      *fc_class=TRUE;
    }

    #ifdef DBGFCWLLKR
    if(Met_seen_p==FALSE||bestmscore<MFCTHRESH)
      fprintf(stderr,"position(MFCWLLKR_Met)=DNE\n\n");
    else
      fprintf(stderr,"position(MFCWLLKR_Met)=%i\n\n",
        bestind-(UPSTREXTENT)-(CONTENTSWATHLEN)+1);
    #endif

    /* If no ATG was seen, we still report a null score. */
    *bestTISscore=bestmscore;
    *Met_seen=Met_seen_p;
  }
  else { /* (transinitmeth==NFCWLLKRX) */
    /* We require the predicted start-Met to be *
     * classified as true by the perceptron.    */
    if(Met_seen_p==FALSE||bestclass==FALSE) {
      start_codon=start_orig;
      *fc_class=FALSE;
    }
    else {
      start_codon=coords_orig[bestind];
      *fc_class=TRUE;
    }

    #ifdef DBGFCWLLKR
    if(Met_seen_p==FALSE||bestclass==FALSE)
      fprintf(stderr,"position(NFCWLLKR_Met)=DNE\n\n");
    else
      fprintf(stderr,"position(NFCWLLKR_Met)=%i\n\n",
        bestind-(UPSTREXTENT)-(CONTENTSWATHLEN)+1);
    #endif

    /* If no ATG was seen, we still report a null score. */
    *bestTISscore=bestnscore;
    *Met_seen=Met_seen_p;
  }

  return(start_codon);
} /* end transinit_fcwllkr */

/* Free memory in the linked list representing a structure */
void freestructure(exoncoorsT *structure) {
  exoncoorsT
    *current,*next; /* Pointers for traversing linked list */

  current=structure;
  while(current->nextelt!=NULL) {
    next=current->nextelt;
    free(current);
    current=next;
  }
  free(current);
  structure=NULL;

  return;
} /* end freestructure */

/* Computes the distance between two data items, *
 * namely the number of nucleotides differing    *
 * between the pair, i.e., edit distance.        */
double eval_dist_kernel(
  char *dat1,
  char *dat2
) {
  int dist=0, /* value to return to caller */
      i;      /* iterator variable         */

  for(i=0;i<STRINGSIZE;++i)
    if(dat1[i]!=dat2[i])
      ++dist;

  return((double)dist);
} /* end eval_dist_kernel */

/* Computes the log-score of the argument, given a PWM */
static double eval_pwm_score(
  ATGparmPWM PWM,
  char *site
) {
  double score=0.0; /* log-score of site */
  int    i;         /* iterator variable */

  for(i=0;i<STRINGSIZE;++i)
    score+=log(PWM.PWMtable[i][trans(site[i])]);

  return(score);
} /* end eval_pwm_score */
