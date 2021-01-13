/* parse_pps_codseqs.c
 * Michael Sparks (mespar1@gmail.com)
 * Last modified : 20 July 2013
 *
 * This utility parses a gthXML document,
 * and generates training sequences that can be used to
 * parameterize a start-Met WAM, potentially with content-
 * based sensing (see CONTENTSWATHLEN).
 * It trains on the basis of the first <orf_entry> element
 * only, as this is typically the longest, and thus, most
 * likely correct (caveat emptor), reading frame.
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
#include <libxml/parser.h>
#include <libxml/relaxng.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <stdio.h>
#include <stdlib.h>
#include "errors.h"
#include "index_utils.h"
#include "libxml2_addendum.h"
#include "MetWAM_utils.h"

#define TRAINSEQOUTSTREAM stdout /* Where to write results */
#define MINAALEN             100 /* Min len of peptides to *
                                  * output for training    */

#define USAGE "\a\nUsage: %s \\\n\
          [-c] \\\n\
          -f bar.gthxml \\\n\
          -g bar.fas.ind \\\n\
          [-r gthxml.rng]\n\
\n\
  Where -c indicates that adequate upstream and downstream sequences\n\
             should be extracted to train for content-based methods,\n\
        -f arg is a gthXML file,\n\
        -g arg is the indexFasSeq'ed index of the genomic sequence\n\
             cognate to the specified gthXML document, and\n\
        -r arg is the RELAX NG schema for the gthXML grammar (optional).\n\
\n"

/* Main Application */
int main(int argc, char *argv[]) {
  char *seqcomplete=NULL,      /* Stores imported template sequence  */
       strand,                 /* Watson (+) or Crick (-)?           */
       *seqfrag,               /* stores parsed sequence fragments   */
       *gthxmlfile=NULL,       /* pointer to gthXML filename         */
       *gdnapath=NULL,         /* pointer to gDNA index filename     */
       *rngschema=NULL;        /* pointer to RNG Schema filename     */
  exoncoorsT exon;             /* stores exon start/stop coordinates */
  int option,                  /* option parsed by getopt            */
      seqct=1,                 /* tracks training instance count     */
      aalen,                   /* length, in AAs, of PPS             */
      prependlen=-1,           /* N-terminus nt flank                */
      appendlen=-1,            /* C-terminus nt flank                */
      content_sensors_p=FALSE, /* Extract data with expectation to   *
                                * train a content-based model?       */
      i;                       /* iterator variable                  */
  const xmlChar *registerednamespaces=BAD_CAST
    "gppatp=http://www.genomethreader.org/GTH_output/PGL_module/\
predicted_gene_location/AGS_information/\
three_phase_translation/probable_ORFs/\0";
  xmlChar *aalenchar;           /* stores length of PPS's      */
  xmlDocPtr gthxmldoc;          /* xml document object         */
  xmlNodePtr myNodePtr,         /* for traversing DOM trees    */
             aalenNodePtr;      /* to find amino acid length   */
  xmlXPathContextPtr xpcontext; /* eval context of xpath query */
  xmlXPathObjectPtr xpresults;  /* results of xpath query      */

  /* Parse command line */
  while((option=getopt(argc,argv,":cf:g:r:"))!=-1) {
    switch(option) {
      case 'c' :
        content_sensors_p=TRUE;
        break;
      case 'f' :
        gthxmlfile=optarg;
        break;
      case 'g' :
        gdnapath=optarg;
        break;
      case 'r' :
        rngschema=optarg;
        break;
      default :
        (void)snprintf(errbuff,MAXERRLEN,USAGE,argv[0]);
        FATALERROR(errbuff)
    }
  }
  if(gthxmlfile==NULL||gdnapath==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,USAGE,argv[0]);
    FATALERROR(errbuff)
  }

  /* Determine amount of flanking sequence to extract */
  prependlen=(int)(UPSTREXTENT);
  appendlen=(int)(DOWNSTREXTENT);
  if(content_sensors_p==TRUE) {
    prependlen+=(int)(CONTENTSWATHLEN);
    appendlen+=(int)(CONTENTSWATHLEN);
  }

  /* Import gthXML document */
  xmlInitParser();
  LIBXML_TEST_VERSION
  if((gthxmldoc=xmlParseFile(gthxmlfile))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Error parsing gthXML document %s!\n",gthxmlfile);
    FATALERROR(errbuff)
  }

  /* Import indexed template sequence */
  if((seqcomplete=slurpindex(gdnapath,seqcomplete))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Error parsing index from %s!\n",gdnapath);
    FATALERROR(errbuff)
  }

  /* Validate the gthXML document, if desired */
  if(rngschema!=NULL&&validateagainstRNG(rngschema,gthxmldoc)==FALSE) {
    xmlFreeDoc(gthxmldoc);
    xmlCleanupParser();
    free(seqcomplete);
    (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
gthXML file %s is incorrect per %s!\n",gthxmlfile,rngschema);
    FATALERROR(errbuff)
  }

  /* Set eval context and register relevant namespaces */
  if((xpcontext=xmlXPathNewContext(gthxmldoc))==NULL)
    FATALERROR("Err (main): Error setting xpath eval context!\n")
  else if(reg_ns(xpcontext,registerednamespaces)!=TRUE)
    FATALERROR("Err (main): \
Error registering namespaces in xpath eval context object!\n")

  /* Obtain node set */
  if((xpresults=xmlXPathEvalExpression(
       BAD_CAST "//gppatp:orf_entry[1]",xpcontext))==NULL)
    FATALERROR("Err (main): Error executing xpath query\n")
  /* Process each relevant node in turn */
  for(i=0;i<xpresults->nodesetval->nodeNr;++i) {
    myNodePtr=xpresults->nodesetval->nodeTab[i];
    DESCEND_TO_CHILD(myNodePtr,"id_line")

    DESCEND_TO_CHILD(myNodePtr,"gDNA")
    RECORD_ATTR_VAL_AS_CHAR(myNodePtr,strand,"strand")

    FIND_SIBLING(myNodePtr,"orf_info")

    aalenNodePtr=myNodePtr;
    DESCEND_TO_CHILD(aalenNodePtr,"number_encoded_amino_acids")

    /* Verify PPS is of adequate length */
    aalenchar=xmlNodeListGetString(gthxmldoc,aalenNodePtr->children,1);
    aalen=atoi((char*)aalenchar);
    xmlFree(aalenchar);
    if(aalen<MINAALEN)
      continue;

    fprintf(TRAINSEQOUTSTREAM,"> %i\n",seqct++);
    DESCEND_TO_CHILD(myNodePtr,"exon_boundaries")
    DESCEND_TO_CHILD(myNodePtr,"exon")

    /* Parse amino-terminus flanking residues */
    RECORD_ATTR_VAL_AS_INT(myNodePtr,exon.start,"start")
    RECORD_ATTR_VAL_AS_INT(myNodePtr,exon.stop,"stop")
    if(strand=='+')
      exon.start-=prependlen;
    else if(strand=='-')
      exon.start+=prependlen;
    else
      FATALERROR("Err (main): gthXML file is non-conformant\n")

    FIND_SIBLING(myNodePtr,"exon")
    if(myNodePtr==NULL) { /* single-exon CDS */
      if(strand=='+')
        exon.stop+=appendlen;
      else /* (strand=='-') */
        exon.stop-=appendlen;
      if((seqfrag=parseseqfrag(seqcomplete,exon.start-1,exon.stop-1))==NULL) {
        (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Couldn't parse sequence fragment (%i to %i)\n",exon.start-1,exon.stop-1);
        FATALERROR(errbuff)
      }
      else if(exon.start==exon.stop&&strand=='-')
        TAKE_COMPLEMENT(seqfrag[0],seqfrag[0])
      fprintf(TRAINSEQOUTSTREAM,"%s",seqfrag);
      free(seqfrag);
    }
    else {
      if((seqfrag=parseseqfrag(seqcomplete,exon.start-1,exon.stop-1))==NULL) {
        (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Couldn't parse sequence fragment (%i to %i)\n",exon.start-1,exon.stop-1);
        FATALERROR(errbuff)
      }
      else if(exon.start==exon.stop&&strand=='-')
        TAKE_COMPLEMENT(seqfrag[0],seqfrag[0])
      fprintf(TRAINSEQOUTSTREAM,"%s",seqfrag);
      free(seqfrag);
      /* Parse interior exons */
      while(myNodePtr!=NULL) {
        /* Note that, on first iteration, exon will record the 2nd exon */
        RECORD_ATTR_VAL_AS_INT(myNodePtr,exon.start,"start")
        RECORD_ATTR_VAL_AS_INT(myNodePtr,exon.stop,"stop")

        FIND_SIBLING(myNodePtr,"exon")
        if(myNodePtr!=NULL) {
          if((seqfrag=parseseqfrag(seqcomplete,
                                   exon.start-1,exon.stop-1))==NULL) {
            (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Couldn't parse sequence fragment (%i to %i)\n",exon.start-1,exon.stop-1);
            FATALERROR(errbuff)
          }
          else if(exon.start==exon.stop&&strand=='-')
            TAKE_COMPLEMENT(seqfrag[0],seqfrag[0])
          fprintf(TRAINSEQOUTSTREAM,"%s",seqfrag);
          free(seqfrag);
        }
      }
      /* Parse carboxyl-terminus flanking residues */
      if(strand=='+')
        exon.stop+=appendlen;
      else /* (strand=='-') */
        exon.stop-=appendlen;
      if((seqfrag=parseseqfrag(seqcomplete,exon.start-1,exon.stop-1))==NULL) {
        (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Couldn't parse sequence fragment (%i to %i)\n",exon.start-1,exon.stop-1);
        FATALERROR(errbuff)
      }
      else if(exon.start==exon.stop&&strand=='-')
        TAKE_COMPLEMENT(seqfrag[0],seqfrag[0])
      fprintf(TRAINSEQOUTSTREAM,"%s",seqfrag);
      free(seqfrag);
    }
    fprintf(TRAINSEQOUTSTREAM,"\n");
  }

  /* cleanup */
  xmlXPathFreeObject(xpresults);
  xmlXPathFreeContext(xpcontext);
  xmlFreeDoc(gthxmldoc);
  xmlCleanupParser();
  free(seqcomplete);

  return(EXIT_SUCCESS);
} /* end main */
