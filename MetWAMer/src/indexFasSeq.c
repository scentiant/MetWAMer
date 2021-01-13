/* Michael E Sparks (mespar1@gmail.com)
 * 20 July 2013
 *
 * indexFasSeq.c
 *
 * This code takes as input the first fasta formatted sequence
 * from a specified input file and either
 * 1) builds an index of it in binary .OR.
 * 2) parses substrings from the indexed sequence on Watson
 *    or Crick strands.  (If a single-base is specified, it will
 *    be returned on the Watson strand; the user must post-process
 *    this result if the complementary base is intended.)
 * See USAGE below for details
 *
 * Copyright (c) 2005,2007,2013 Michael E Sparks
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

#define USAGE "\a\n\
    `%s file2index.fas` -> build index named \"file2index.fas.ind\" \n\
  .OR. \n\
    `%s indexedfile.ind start stop` -> parse substrings\n\
        (if start .LE. stop, use Watson; else, use Crick)\n\
  .OR. \n\
    `%s indexedfile.ind start stop 1` -> parse substrings .AND. translate\n\
        (if start .LE. stop, use Watson; else, use Crick)\n\n"

int main(int argc, char *argv[]) {
  char *seqcomplete=NULL,
       *seqfrag=NULL,
       desc[MAXDESCLINE+1],
       indname[255],
       *peptide=NULL;
  FILE *fptr=NULL;
  int start,
      stop,
      length;

  /* "Operator overloading" results in context-specific behavior */
  if(argc!=2&&argc!=4&&argc!=5) {
    (void)snprintf(errbuff,MAXERRLEN,USAGE,argv[0],argv[0],argv[0]);
    FATALERROR(errbuff)
  }
  else if(argc>2) { /* Parse requested sequence fragment from index */
    /* The user requests start/stop coordinates using a one-based *
     * system, but parseseqfrag operates with a zero-based one.   *
     * For this utility, if the user wants to parse a single base *
     * from the reverse-sense strand, they will have to post-     *
     * process their results.                                     */
    start=atoi(argv[2])-1;
    stop=atoi(argv[3])-1;
    if((seqfrag=oneshotparse(argv[1],seqfrag,start,stop))==NULL)
      FATALERROR("Err (main): Indexed substring outside range!\n")
    else
      if(argc==4)
        fprintf(stdout,"%s\n",seqfrag);
      else {
        INIT_STD_GEN_CODE
        peptide=translate_cds(seqfrag,stdgencode);
        fprintf(stdout,"%s\n",peptide);
        free(peptide);
      }
    free(seqfrag);
  }
  else if(argc==2) { /* Create a sequence index of FIRST seq in argv[1] */
    if((fptr=fopen(argv[1],"rt"))==NULL||
       (seqcomplete=get_fasta(fptr,seqcomplete,&length,desc))==NULL) {
      (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Error parsing sequence from %s!\n",argv[1]);
      FATALERROR(errbuff)
    }
    else
      fclose(fptr);

    (void)strcpy(indname,argv[1]);
    (void)strcat(indname,".ind");
    if((fptr=fopen(indname,"wb"))==NULL||
       fwrite(&length,sizeof(int),1,fptr)!=1||
       fwrite(seqcomplete,sizeof(char),length,fptr)!=length) {
      (void)snprintf(errbuff,MAXERRLEN,"Err (main): \
Failure manipulating binary file %s.\n",indname);
      FATALERROR(errbuff)
    }
    else {
      free(seqcomplete);
      fclose(fptr);
    }
  }
  else
    FATALERROR("Err (main): An unexpected error was encountered.\n")

  return(EXIT_SUCCESS);
}
