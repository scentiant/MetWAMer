/* calc_medoids.c
 * Michael Sparks (mespar1@gmail.com)
 * Last modified : 20 July 2013
 *
 * This utility computes the user-specified number of medoids
 * in a dataset of same-length nucleotide sequences by optimizing
 * over a number of kmedoids runs.  These best results encountered
 * during this operation are written to the specified output file.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calc_medoids.h"
#include "cluster.h"
#include "errors.h"
#include "MetWAM_utils.h"
#include "sequence_parse.h"

#define EMMAXRUN     100 /* Max number of passes of the EM algorithm *
                          * during k-medoids clustering.             */
#define KMEDOIDCALL   30 /* Number of attempts to run the kmedoid    *
                          * algorithm.  The best clustering found,   *
                          * i.e., that with minimal error, in this   *
                          * many attempts will be returned.          */

#define ARGCT 4
#define USAGE "\a\nUsage: %s k sites.fas outfile\n\
\n\
  Where k is the number of distinct clusters to produce,\n\
        sites.fas contains a set of translation initiation sites, and\n\
        outfile is the desired name for the output file.\n\
\n"

/* Structure to record the numclusts distinct medoids,  *
 * each of length stringsize, into a concatemer of such *
 * strings, medoid_concat.                              */
typedef struct medoidlist {
  int stringsize,
      numclusts;
  char *medoid_concat;
} medoidlistT;

/* Prototypes */
static void record_medoids(medoidlistT *medoids,int nelts,
  int nclusts,int *labels,int *labelsbest,char **sites);

/* Main Application */
int main(int argc,char *argv[]) {
  char *sequence=NULL,       /* genomic locus to be processed */   
       seqid[MAXDESCLINE+1], /* description of locus, +1 for  *
                              * terminating '\0'.             */
       outname[MAXFILENAME], /* Stores output file name       */
       **sites=NULL;         /* stores translation init sites */
  double **kernel=NULL,      /* stores the evaluated kernel   */
         currerr,            /* stores error of clustering    */
         besterr;            /* stores error of clustering    */
  FILE *seqfile=NULL,        /* connects to input stream      */
       *outfile=NULL;        /* connects to output stream     */
  int nclusts,               /* stores k, the cluster count   */
      nelts,                 /* total number of sites         */
      seqlen,                /* length of a seq from seqfile  */
      *labels=NULL,          /* array of cluster assignments  */
      *labelsbest=NULL,      /* array of cluster assignments  */
      optclustfreq,          /* optimal clustering frequency  */
      i,j,k;                 /* iterator variables            */
  medoidlistT medoidsloc;    /* stores the identified medoids */

  /* Verify command line */
  if(argc!=ARGCT) {
    (void)snprintf(errbuff,MAXERRLEN,USAGE,argv[0]);
    FATALERROR(errbuff)
  }

  /* Parse nclusts (k) */
  if((nclusts=atoi(argv[1]))<2)
    FATALERROR("Err (main): Nonsensible value for k specified.\n")

  /* Learn the number of sequences in the input file... */
  if((seqfile=fopen(argv[2],"rt"))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (main): Can't open %s.\n",argv[2]);
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
  if((sites=(char**)malloc(sizeof(char*)*nelts))==NULL)
    FATALERROR("Err (main): Out of memory.\n")
  else
    for(i=0;i<nelts;++i)
      if((sites[i]=(char*)malloc(sizeof(char)*((int)(STRINGSIZE)+1)))==NULL)
        FATALERROR("Err (main): Out of memory.\n")
      else {
        sequence=get_fasta(seqfile,sequence,&seqlen,seqid);
        for(j=0;j<(int)(STRINGSIZE);++j)
          sites[i][j]=sequence[j];
        sites[i][STRINGSIZE]='\0';
        free(sequence);
      }
  fclose(seqfile);

  /* Allocate space for the cached, evaluated kernel function... */
  if((kernel=(double**)malloc(sizeof(double*)*nelts))==NULL)
    FATALERROR("Err (main): Out of memory.\n")
  else
    for(kernel[0]=NULL,i=1;i<nelts;++i)
      if((kernel[i]=(double*)malloc(sizeof(double)*i))==NULL)
        FATALERROR("Err (main): Out of memory.\n")
  /* ...and populate its lower triangle */
  for(i=1;i<nelts;++i) /* rows */
    for(j=0;j<i;++j) /* cols */
      kernel[i][j]=eval_dist_kernel(sites[i],sites[j]);

  /* Allocate the (working) class label array */
  if((labels=(int*)malloc(sizeof(int)*nelts))==NULL||
     (labelsbest=(int*)malloc(sizeof(int)*nelts))==NULL)
    FATALERROR("Err (main): Out of memory.\n")

  /* Find an optimal kmedoids-based clustering *
   * over up to KMEDOIDCALL iterations.        */
  if((medoidsloc.medoid_concat=(char*)malloc(sizeof(char)*
      (nclusts*(int)(STRINGSIZE))))==NULL)
    FATALERROR("Err (main): Out of memory.\n")
  do {
    (void)kmedoids(nclusts,nelts,kernel,EMMAXRUN,
      labels,&besterr,&optclustfreq);
    record_medoids(&medoidsloc,nelts,nclusts,
      labels,labelsbest,sites);
  } while(optclustfreq<=0);
  for(i=1;i<(int)(KMEDOIDCALL); ) {
    (void)kmedoids(nclusts,nelts,kernel,EMMAXRUN,
      labels,&currerr,&optclustfreq);
    if(optclustfreq>0) {
      ++i;
      if(currerr<besterr) {
        record_medoids(&medoidsloc,nelts,nclusts,
          labels,labelsbest,sites);
        besterr=currerr;
      }
    }
  }

  #ifdef SHOWCLUSTERING
  for(i=0;i<nelts;++i)
    fprintf(stdout,"%i\t%i\t%s\n",i+1,labelsbest[i]+1,sites[i]);
  #endif

  /* Write learned medoids to output stream */
  (void)strcpy(outname,argv[3]);
  (void)strcat(outname,".medoids.xml");
  if((outfile=fopen(outname,"wt"))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,
      "Err (main): Can't open %s for writing!\n",outname);
    FATALERROR(errbuff)
  }
  fprintf(outfile,"<?xml \
version=\"1.0\" \
encoding=\"ISO-8859-1\"?>\n");
  fprintf(outfile,"<kmedoids \
xmlns=\"http://www.kmedoids.org\" \
stringsize=\"%i\" \
k=\"%i\" \
error=\"%.2f\">\n",medoidsloc.stringsize,medoidsloc.numclusts,besterr);
  for(i=k=0;i<nclusts;++i) {
    fprintf(outfile,"  <medoid id=\"%i\">",i+1);
    for(j=0;j<(int)(STRINGSIZE);++j)
      fprintf(outfile,"%c",medoidsloc.medoid_concat[k++]);
    fprintf(outfile,"</medoid>\n");
  }
  fprintf(outfile,"</kmedoids>\n");
  fclose(outfile);

  /* cleanup memory */
  free(medoidsloc.medoid_concat);
  free(labels);
  free(labelsbest);
  for(i=1;i<nelts;++i)
    free(kernel[i]);
  free(kernel);
  for(i=0;i<nelts;++i)
    free(sites[i]);
  free(sites);

  return(EXIT_SUCCESS);
} /* end main */

/* Routine to identify the k (nclusts) distinct medoids from   *
 * the labels array, recording these in the medoids structure. */
static void record_medoids(
  medoidlistT *medoids,
  int nelts,
  int nclusts,
  int *labels,
  int *labelsbest,
  char **sites
) {
  int *clusts=NULL, /* records k distinct centroids */
      clustssize=0, /* size of clusts array         */
      new_clust_p,  /* predicate for growing clusts */
      i,j,k;        /* iterator variables           */

  /* Copy the labels array to labelsbest */
  for(i=0;i<nelts;++i)
    labelsbest[i]=labels[i];

  /* Learn the k distinct medoids */
  for(i=0;i<nelts;++i) {
    new_clust_p=TRUE;
    for(j=0;j<clustssize;++j)
      if(clusts[j]==labels[i]) {
        new_clust_p=FALSE;
        break;
      }
    if(new_clust_p==TRUE) { /* push medoid onto list */
      if((clusts=(int*)realloc(clusts,sizeof(int)*++clustssize))==NULL)
        FATALERROR("Err (record_medoids): Out of memory.\n")
      else
        clusts[clustssize-1]=labels[i];
    }
  }
  /* Sanity check */
  if(clustssize!=nclusts)
    FATALERROR("Err (record_medoids): Error from kmedoids function!\n")

  /* Record in medoids structure */
  medoids->stringsize=(int)(STRINGSIZE);
  medoids->numclusts=nclusts;
  for(i=k=0;i<medoids->numclusts;++i)
    for(j=0;j<(int)(STRINGSIZE);++j)
      medoids->medoid_concat[k++]=sites[clusts[i]][j];

  free(clusts);
  return;
} /* end record_medoids */
