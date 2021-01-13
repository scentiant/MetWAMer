/* parse_tigrxml_cds.c
 * Michael Sparks (mespar1@gmail.com)
 * Last modified : 20 July 2013
 *
 * This utility parses a TIGR/TAIR xml document,
 * and generates training sequences that can be used to
 * parameterize a start-Met WAM, potentially with content-
 * based sensing (see CONTENTSWATHLEN).
 * It trains on the basis of annotated CDS's of curated
 * gene models, s.t. they begin with a Methionine, and
 * encode a protein at least MINAALEN residues in length.
 * Note that, if using output from this tool as input to
 * train_MetWAM, the user must specify an offset of
 * prependlen+1 (see below, near ``Verify command line'').
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
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <stdio.h>
#include <string.h>
#include "errors.h"
#include "index_utils.h"
#include "libxml2_addendum.h"
#include "MetWAM_utils.h"

#define TRAINSEQOUTSTREAM stdout /* Where to write results */
#define MINAALEN             100 /* Min len of peptides to *
                                  * output for training    */

#define MINARGCT 3
#define USAGE "\a\nUsage: %s foo.tigrxml foo.fas.ind [c]\n\
\n\
  Where foo.tigrxml is a TIGR xml document for template foo,\n\
        foo.fas.ind is the indexFasSeq'ed index of foo.fas, and\n\
        c can be optionally specified to extract adequate upstream\n\
          and downstream sequences to train for content-based methods.\n\
\n"

/* Prototypes */
exoncoorsT *sort_exon_order(exoncoorsT *genstruct,char strand);

/* Main Application */
int main(int argc, char *argv[]) {
  char
    *seqcomplete=NULL, /* Stores imported template sequence  */
    *seqfrag=NULL,     /* stores parsed sequence fragments   */
    *seqprint=NULL,    /* stores training sequence to print  */
    strand='^';        /* Watson (+) or Crick (-)?           */
  exoncoorsT 
    *exon=NULL,        /* stores exon start/stop coordinates */
    *headptr=NULL,
    *prevptr=NULL,
    *currptr=NULL;
  int
    prependlen=-1,      /* N-terminus nt flank      */
    appendlen=-1,       /* C-terminus nt flank      */
    peplen,             /* length, in AAs, of PPS   */
    seqct=1,            /* training seq count       */
    genstart,           /* gene starting coordinate */
    genstop,            /* gene stopping coordinate */
    exonstart,          /* exon starting coordinate */
    exonstop,           /* exon stopping coordinate */
    currpos=-1,         /* index for seqprint       */
    i;                  /* iterator variable        */
  xmlChar *pepseq;      /* stores length of PPS's   */
  xmlDocPtr tigrxmldoc; /* xml document object      */
  xmlNodePtr            /* for DOM tree traversal   */
    modelParent,
    nodePtrGeneric1,
    nodePtrGeneric2;
  xmlXPathContextPtr /* eval context of xpath query */
    xpcontext;
  xmlXPathObjectPtr /* results of xpath query */
    xpresults;

  /* Verify command line */
  if(argc<MINARGCT) {
    (void)snprintf(errbuff,MAXERRLEN,USAGE,argv[0]);
    FATALERROR(errbuff)
  }
  else if(argc==MINARGCT) {
    prependlen=(int)(UPSTREXTENT);
    appendlen=(int)(DOWNSTREXTENT);
  }
  else {
    prependlen=(int)((UPSTREXTENT)+(CONTENTSWATHLEN));
    appendlen=(int)((DOWNSTREXTENT)+(CONTENTSWATHLEN));
  }

  /* Import TIGR/TAIR xml document */
  xmlInitParser();
  LIBXML_TEST_VERSION
  if((tigrxmldoc=xmlParseFile(argv[1]))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Error parsing TIGR xml document %s!\n",argv[1]);
    FATALERROR(errbuff)
  }

  /* Import indexed template sequence */
  if((seqcomplete=slurpindex(argv[2],seqcomplete))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Error parsing index from %s!\n",argv[2]);
    FATALERROR(errbuff)
  }

  /* Set eval context and register relevant namespaces */
  if((xpcontext=xmlXPathNewContext(tigrxmldoc))==NULL)
    FATALERROR("Err (main): Error setting xpath eval context!\n")

  /* Obtain node set */
  if((xpresults=xmlXPathEvalExpression(
       BAD_CAST "//PROTEIN_CODING/TU/MODEL[@CURATED='1']",xpcontext))==NULL)
    FATALERROR("Err (main): Error executing xpath query\n")

  /* Assuming parent points to a COORDSET element, this macro  *
   * will record the start (END5) and stop (END3) coordinates. */
  #define SLURP_COORDSET(doc,parent,start,stop) { \
    xmlChar *coordchar; \
    xmlNodePtr helper; \
    helper=parent; \
    DESCEND_TO_CHILD(helper,"END5") \
    coordchar=xmlNodeListGetString(doc,helper->children,1); \
    start=atoi((char*)coordchar); \
    free(coordchar); \
    helper=parent; \
    DESCEND_TO_CHILD(helper,"END3") \
    coordchar=xmlNodeListGetString(doc,helper->children,1); \
    stop=atoi((char*)coordchar); \
    free(coordchar); \
  }

  /* Process each relevant node in turn */
  for(i=0;i<xpresults->nodesetval->nodeNr;++i) {
    modelParent=nodePtrGeneric1=xpresults->nodesetval->nodeTab[i];

    /* Verify protein sequence is of adequate length */
    ATTEMPT_DESCENT_TO_CHILD(nodePtrGeneric1,"PROTEIN_SEQUENCE")
    if(nodePtrGeneric1==NULL)
      continue;
    pepseq=xmlNodeListGetString(tigrxmldoc,nodePtrGeneric1->children,1);
    peplen=xmlStrlen(pepseq);
    if((char)pepseq[peplen-1]=='*')
      --peplen;
    xmlFree(pepseq);
    if(peplen<MINAALEN)
      continue;

    /* Determine orientation of gene */
    nodePtrGeneric1=modelParent;
    DESCEND_TO_CHILD(nodePtrGeneric1,"COORDSET")
    SLURP_COORDSET(tigrxmldoc,nodePtrGeneric1,genstart,genstop)
    if(genstart<genstop)
      strand='+';
    else if(genstart>genstop)
      strand='-';
    else
      NONFATALERROR("Warn (main): One-nucleotide CDS encountered?!!\n")

    /* Record CDS coordinates of individual exons. Because TIGR xml *
     * does not guarantee that exons are recorded in any particular *
     * order, we will need to sort these, prior to parsing the CDS. */
    nodePtrGeneric1=modelParent;
    DESCEND_TO_CHILD(nodePtrGeneric1,"EXON")
    while(nodePtrGeneric1!=NULL) {
      nodePtrGeneric2=nodePtrGeneric1;
      ATTEMPT_DESCENT_TO_CHILD(nodePtrGeneric2,"CDS")
      if(nodePtrGeneric2!=NULL) { /* exon with coding bases */
        DESCEND_TO_CHILD(nodePtrGeneric2,"COORDSET")
        SLURP_COORDSET(tigrxmldoc,nodePtrGeneric2,exonstart,exonstop)
        if((currptr=(exoncoorsT*)malloc(sizeof(exoncoorsT)))==NULL)
          FATALERROR("Err (main): Insufficient memory.\n")
        else {
          if(headptr==NULL)
            headptr=currptr;
          currptr->start=exonstart;
          currptr->stop=exonstop;
          currptr->nextelt=NULL;
          if(prevptr!=NULL)
            prevptr->nextelt=currptr;
          prevptr=currptr;
        }
      }
      FIND_SIBLING(nodePtrGeneric1,"EXON")
    } /* end while */

    /* Sort exoncoorsT linked list */
    exon=headptr=sort_exon_order(headptr,strand);

    /* Parse amino-terminus flanking residues.  Note that it is the *
     * responsibility of training/testing code to ensure that Met   *
     * codons occurring in the first prependlen residues of the     *
     * parsed sequence not be considered; see the macro,            *
     * MET_VALID_TO_SCORE_P in MetWAM_utils.h to facilitate this.   */
    if(strand=='+')
      exon->start-=prependlen;
    else /* (strand=='-') */
      exon->start+=prependlen;
    if(exon->nextelt==NULL) { /* single-exon CDS */
      if(strand=='+')
        exon->stop+=appendlen;
      else
        exon->stop-=appendlen;
      if((seqfrag=parseseqfrag(seqcomplete,exon->start-1,exon->stop-1))==NULL)
        FATALERROR("Err (main): Couldn't parse sequence fragment\n")
      else if(exon->start==exon->stop&&strand=='-')
        TAKE_COMPLEMENT(seqfrag[0],seqfrag[0])
      if(strlen(seqfrag)>=prependlen+3&&seqfrag[prependlen]=='A'&&
         seqfrag[prependlen+1]=='T'&&seqfrag[prependlen+2]=='G')
        fprintf(TRAINSEQOUTSTREAM,"> %i\n%s\n",seqct++,seqfrag);
      free(seqfrag);
    }
    else {
      if((seqfrag=parseseqfrag(seqcomplete,exon->start-1,exon->stop-1))==NULL)
        FATALERROR("Err (main): Couldn't parse sequence fragment\n")
      else if(exon->start==exon->stop&&strand=='-')
        TAKE_COMPLEMENT(seqfrag[0],seqfrag[0])
      if((seqprint=(char*)malloc(sizeof(char)*
            (abs(exon->start-exon->stop)+2)))==NULL)
        FATALERROR("Err (main): Insufficient memory\n")
      else {
        (void)snprintf(seqprint,abs(exon->start-exon->stop)+2,"%s",seqfrag);
        currpos=abs(exon->start-exon->stop)+1;
      }
      free(seqfrag);
      /* Parse interior exons */
      exon=exon->nextelt;
      while(exon->nextelt!=NULL) {
        /* Note that, on first iteration, exon will record the 2nd exon */
        if((seqfrag=parseseqfrag(seqcomplete,
              exon->start-1,exon->stop-1))==NULL)
          FATALERROR("Err (main): Couldn't parse sequence fragment\n")
        else if(exon->start==exon->stop&&strand=='-')
          TAKE_COMPLEMENT(seqfrag[0],seqfrag[0])
        if((seqprint=(char*)realloc(seqprint,sizeof(char)*
              (abs(exon->start-exon->stop)+2+currpos)))==NULL)
          FATALERROR("Err (main): Insufficient memory\n")
        else {
          (void)snprintf(seqprint+currpos,abs(exon->start-exon->stop)+2,
                  "%s",seqfrag);
          currpos+=abs(exon->start-exon->stop)+1;
        }
        free(seqfrag);
        exon=exon->nextelt;
      }
      /* Parse carboxyl-terminus flanking residues */
      if(strand=='+')
        exon->stop+=appendlen;
      else /* (strand=='-') */
        exon->stop-=appendlen;
      if((seqfrag=parseseqfrag(seqcomplete,exon->start-1,exon->stop-1))==NULL)
        FATALERROR("Err (main): Couldn't parse sequence fragment\n")
      else if(exon->start==exon->stop&&strand=='-')
        TAKE_COMPLEMENT(seqfrag[0],seqfrag[0])
      if((seqprint=(char*)realloc(seqprint,sizeof(char)*
            (abs(exon->start-exon->stop)+2+currpos)))==NULL)
        FATALERROR("Err (main): Insufficient memory\n")
      else
        (void)snprintf(seqprint+currpos,abs(exon->start-exon->stop)+2,
                "%s",seqfrag);
      free(seqfrag);

      if(strlen(seqprint)>=prependlen+3&&seqprint[prependlen]=='A'&&
         seqprint[prependlen+1]=='T'&&seqprint[prependlen+2]=='G')
        fprintf(TRAINSEQOUTSTREAM,"> %i\n%s\n",seqct++,seqprint);
      free(seqprint);
    }
    freestructure(headptr);
    headptr=currptr=prevptr=NULL;
  } /* end for */

  /* cleanup */
  xmlXPathFreeObject(xpresults);
  xmlXPathFreeContext(xpcontext);
  xmlFreeDoc(tigrxmldoc);
  xmlCleanupParser();
  free(seqcomplete);

  return(EXIT_SUCCESS);
} /* end main */

/* Sort exons in genstruct, returning a pointer to head node of list */
exoncoorsT *sort_exon_order(
  exoncoorsT *genstruct,
  char strand
) {
  exoncoorsT
    *ptr=genstruct,
    **ptrarray=NULL;
  int exonct=0,
      i,j;

  /* Learn the count of exons in the structure */
  do {
    ++exonct;
    ptr=ptr->nextelt;
  } while(ptr!=NULL);

  /* Build and sort an array of pointers to exons, based *
   * on increasing order of the exon start positions.    */
  if((ptrarray=(exoncoorsT**)malloc(sizeof(exoncoorsT*)*exonct))==NULL)
    FATALERROR("Err (sort_exon_order): Insufficient memory\n")
  else
    for(i=0,ptr=genstruct;i<exonct;++i,ptr=ptr->nextelt)
      ptrarray[i]=ptr;

  /* bubble sort */
  for(i=0;i<exonct;++i)
    for(j=exonct-1;j>i;--j)
      if(ptrarray[j]->start<ptrarray[j-1]->start) {
        ptr=ptrarray[j];
        ptrarray[j]=ptrarray[j-1];
        ptrarray[j-1]=ptr;
      }

  /* Reconfigure links in the list */
  if(strand=='+') {
    for(i=0;i<exonct-1;++i)
      ptrarray[i]->nextelt=ptrarray[i+1];
    ptrarray[i]->nextelt=NULL;
    genstruct=ptrarray[0];
  }
  else if(strand=='-') {
    for(i=exonct-1;i>0;--i)
      ptrarray[i]->nextelt=ptrarray[i-1];
    ptrarray[i]->nextelt=NULL;
    genstruct=ptrarray[exonct-1];
  }
  else
    FATALERROR("Err (sort_exon_order): Ambiguous strand specified\n")

  free(ptrarray);
  return(genstruct);
} /* end sort_exon_order */
