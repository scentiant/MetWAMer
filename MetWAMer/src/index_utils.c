/* index_utils.c
 * Michael E Sparks (mespar1@gmail.com)
 * Last modified: 20 July 2013
 *
 * Functions used by indexFasSeq and related utilities.
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
#include "errors.h"
#include "index_utils.h"
#include "sequence_parse.h"

/* Function to import an indexFasSeq'ed sequence from the filesystem */
char *slurpindex(
  const char *infile,
  char *seq
) {
  FILE *fptr=NULL;
  int length;

  if(seq!=NULL)
    FATALERROR("Err (slurpindex): Passed seq as non-NULL?!!\n")
  else if((fptr=fopen(infile,"rb"))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (slurpindex): \
Can't open %s for binary reading!\n",infile);
    FATALERROR(errbuff)
  }
  else if(fread(&length,sizeof(int),1,fptr)!=1) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (slurpindex): \
Error reading sequence length from %s!\n",infile);
    FATALERROR(errbuff)
  }
  else if((seq=(char*)malloc(sizeof(char)*(length+1)))==NULL)
    FATALERROR("Err (slurpindex): Out of memory!\n")
  else if(fread(seq,length*sizeof(char),1,fptr)!=1) {
    free(seq);
    (void)snprintf(errbuff,MAXERRLEN,"Err (slurpindex): \
Error reading sequence from %s!\n",infile);
    FATALERROR(errbuff)
  }
  else
    (void)fclose(fptr);

  seq[length]='\0';
  return(seq);
} /* end slurpindex */

/* Function to efficiently parse substrings *
 * of an index in a one-shot manner.  start *
 * and stop should be supplied using zero-  *
 * based coordinates.                       */
char *oneshotparse(
  const char *infile,
  char *seq,
  int start,
  int stop
) {
  char c;
  FILE *fptr=NULL;
  int indexlen,
      fraglen,
      i,j,k;

  if(seq!=NULL)
    FATALERROR("Err (oneshotparse): Passed seq as non-NULL?!!\n")
  else if((fptr=fopen(infile,"rb"))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (oneshotparse): \
Can't open %s for binary reading!\n",infile);
    FATALERROR(errbuff)
  }
  else if(fread(&indexlen,sizeof(int),1,fptr)!=1) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (oneshotparse): \
Error reading sequence length from %s!\n",infile);
    FATALERROR(errbuff)
  }

  fraglen=abs(start-stop)+1;
  if(fraglen>indexlen||start<0||stop<0||start>=indexlen||stop>=indexlen)
    FATALERROR("Err (oneshotparse): Indexed substring outside range\n")
  /* Cue up the appropriate starting point in the file */
  if(start<=stop)
    fseek(fptr,start,SEEK_CUR);
  else
    fseek(fptr,stop,SEEK_CUR);
  /* Read in the fragment */
  if((seq=(char*)malloc(sizeof(char)*(fraglen+1)))==NULL)
    FATALERROR("Err (oneshotparse): Insufficient memory.\n")
  else if(fread(seq,sizeof(char)*fraglen,1,fptr)!=1) {
    (void)snprintf(errbuff,MAXERRLEN,"Err (oneshotparse): \
Error reading sequence substring from %s!\n",infile);
    FATALERROR(errbuff)
  }
  else {
    seq[fraglen]='\0';
    fclose(fptr);
  }

  /* Ensure sequence is in upper case, and *
   * take reverse complement, if desired.  */
  for(i=0;i<fraglen;++i) {
    seq[i]=(char)toupper((int)seq[i]);
    /* Recall that, if start .EQ. stop, and the user *
     * wants the reverse complementary base, he/she  *
     * must post-process the result, externally!     */
    if(start>stop)
      TAKE_COMPLEMENT(seq[i],seq[i])
  }
  if(start>stop) /* reverse the complement */
    for(j=0,k=fraglen-1;k>=fraglen/2;++j,--k) {
      /* This swap is pointless at the center of an odd-length   *
       * string, but it's cheaper to make one wasteful operation *
       * than checking whether fraglen%2==1 every iteration.     */
      c=seq[k];
      seq[k]=seq[j];
      seq[j]=c;
    }

  return(seq);
} /* end oneshotparse */

/* Function to parse a sequence fragment from hostseq, *
 * using a ZERO-BASED coordinate system.  While for a  *
 * one-shot fragment parse, an fseek() operation would *
 * be more efficient, this code was designed with      *
 * repeated extractions in mind, for which importing   *
 * the data into memory is quicker.                    */
char *parseseqfrag(
  char *hostseq,
  int start,
  int stop
) {
  char *frag=NULL;
  int hostlen,
      fraglen,
      i,j;

  hostlen=(int)strlen(hostseq);
  fraglen=abs(start-stop)+1;

  if(fraglen>hostlen||start<0||stop<0||start>=hostlen||stop>=hostlen)
    return(NULL);
  else if((frag=(char*)malloc(sizeof(char)*(fraglen+1)))==NULL)
    FATALERROR("Err (parseseqfrag): Out of memory!\n")

  if(start<=stop)
    for(i=start,j=0;i<=stop;++i,++j)
      frag[j]=(char)toupper((int)hostseq[i]);
  else
    for(i=start,j=0;i>=stop;--i,++j)
      TAKE_COMPLEMENT(hostseq[i],frag[j])
  frag[j]='\0';
  return(frag);
} /* end parseseqfrag */

/* Returns a char string of the translated cds */
char *translate_cds(
  char *cds,
  char gencode[NTALFSIZE][NTALFSIZE][NTALFSIZE]
) {
  char *peptide=NULL;
  int cdslen,
      *cdsint=NULL,
      peplen,
      i,j;

  if((cdslen=strlen(cds))<1)
    FATALERROR("Err (translate_cds): No CDS to translate!\n")
  else if((cdsint=(int*)malloc(sizeof(int)*cdslen))==NULL)
    FATALERROR("Err (translate_cds): Insufficient memory.\n")
  else if((peptide=(char*)malloc(sizeof(char)*
           (peplen=cdslen/3+(cdslen%3?1:0)+1)))==NULL)
    FATALERROR("Err (translate_cds): Insufficient memory.\n")

  /* Translate the cds using supplied gencode */
  for(i=0;i<cdslen;++i)
    cdsint[i]=trans(cds[i]);
  for(i=0,j=0;i<cdslen-2;i+=3)
    peptide[j++]=gencode[cdsint[i]][cdsint[i+1]][cdsint[i+2]];
  if(cdslen%3)
    peptide[j++]='?';
  peptide[j]='\0';

  free(cdsint);
  return(peptide);
} /* end translate_cds */
