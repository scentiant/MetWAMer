/* Michael E Sparks (mespar1@iastate.edu)
 * Last modified: 17 July 2007
 *
 * indexFasSeq.c
 *
 * This code takes as input the first fasta formatted sequence
 * from a specified input file and either
 * 1) builds an index of it in binary .OR.
 * 2) parses substrings from the indexed sequence on Watson
 *    or Crick strands.
 * See USAGE below for details
 *
 * Copyright (c) 2005,2007 Michael E Sparks
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

#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Note that, if the user wants to take the reverse complement of a *
 * single base, they will be required to post-process the results.  */
#define USAGE "\a\n  `%s file2index.fas` -> \
build index named \"file2index.fas.ind\" \n\
    .OR. \n\
  `%s indexedfile.ind start stop` -> parse substrings\n\
     (if start <= stop, use Watson; else, use Crick)\n\n"

char *get_fasta(FILE *file,char *seq,int *seqlength);

int main (int argc, char *argv[]) {
  char *seq=NULL,
       indname[255];
  FILE *fptr=NULL;
  int length,
      fraglen,
      start,
      stop,
      i,
      c;

  if(argc!=2&&argc!=4) {
    fprintf(stdout,USAGE,argv[0],argv[0]);
    exit(EXIT_FAILURE);
  }
  else if(argc==4) { /* Input file, start stop */
    if((fptr=fopen(argv[1],"rb"))==NULL) {
      fprintf(stderr,"Error: Can't open %s for binary reading!\n",argv[1]);
      exit(EXIT_FAILURE);
    }
    else if(fread(&length,sizeof(int),1,fptr)!=1) {
      fprintf(stderr,"Error: Error reading from %s!\n",argv[1]);
      exit(EXIT_FAILURE);
    }
    start=atoi(argv[2])-1;
    stop=atoi(argv[3])-1;
    fraglen=abs(start-stop)+1;
    if(fraglen>length||start<0||stop<0||start>=length||stop>=length) {
      free(seq);
      fprintf(stderr,"Indexed substring outside range\n");
      exit(EXIT_FAILURE);
    }
    /* Cue up the appropriate starting point in the file */
    if(start<=stop)
      fseek(fptr,start,SEEK_CUR);
    else
      fseek(fptr,stop,SEEK_CUR);
    /* Read in the fragment */
    if((seq=(char*)malloc(sizeof(char)*(fraglen+1)))==NULL) {
      fprintf(stderr,"Error: Out of memory!\n");
      exit(EXIT_FAILURE);
    }
    else if(fread(seq,sizeof(char)*fraglen,1,fptr)!=1) {
      free(seq);
      fprintf(stderr,"Error: Error reading from %s!\n",argv[1]);
      exit(EXIT_FAILURE);
    }
    else {
      seq[fraglen]='\0';
      (void)fclose(fptr);
    }
    /* Report the fragment to the stdout stream */
    if(start<=stop)
      for(i=0;i<fraglen;++i) {
        c=toupper((int)seq[i]);
        fprintf(stdout,"%c",(char)c);
      }
    else
      for(i=fraglen-1;i>=0;--i) {
        switch(seq[i]) {
          case 'A' :
          case 'a' :
            c='T';
            break;
          case 'C' :
          case 'c' :
            c='G';
            break;
          case 'G' :
          case 'g' :
            c='C';
            break;
          case 'T' :
          case 't' :
            c='A';
            break;
          default :
            c=toupper((int)seq[i]);
        }
        fprintf(stdout,"%c",(char)c);
      }
    free(seq);
  }
  else if(argc==2) {
    /* Input sequence */
    if((fptr=fopen(argv[1],"rt"))==NULL||
       (seq=get_fasta(fptr,seq,&length))==NULL) {
      fprintf(stderr,"An error was encountered.\n");
      exit(EXIT_FAILURE);
    }
    else
      fclose(fptr);
    /* Write sequence to binary output file */
    (void)strcpy(indname,argv[1]);
    (void)strcat(indname,".ind");
    if((fptr=fopen(indname,"wb"))==NULL||
       fwrite(&length,sizeof(int),1,fptr)!=1||
       fwrite(seq,sizeof(char),length,fptr)!=length) {
      free(seq);
      fprintf(stderr,"Failure handling binary file.\n");
      exit(EXIT_FAILURE);
    }
    else {
      fclose(fptr);
      free(seq);
    }
  }
  else {
    fprintf(stderr,"An error was encountered.\n");
    exit(EXIT_FAILURE);
  }

  return(EXIT_SUCCESS);
}

/* get_fasta:  A general function to parse FASTA entries  *
 * Input :                                                *
 *   A file pointer: This must be connected to a stream   *
 *                   containing Fasta-formatted sequences *
 *   A char pointer: This must be set to NULL.  The       *
 *                   function will attempt to allocate    *
 *                   sufficient memory to store the seq   *
 *   An int pointer: This is a single integer.  The       *
 *                   function will store the length of    *
 *                   the current sequence in it.          */
char *get_fasta(FILE *file,char *seq,int *seqlength) {
  int c,
      length;
  fpos_t fpos;

  /* Verify comment line */
  if ( ((c=fgetc(file)) != (int)'>') || (bool)(feof(file)) ) {
    if ((bool)feof(file))
      return(NULL);
    else { /* Improper format */
      fprintf(stderr,"Error: This requires FASTA formatted input!\n");
      return(NULL);
    }
  }
  else {
    /* Scroll past the comment line and ignore it */
    while ( ((c=fgetc(file)) != (int)'\n') && (bool)(c != EOF) )
      ;
    /* We should always expect sequence after the comment line */
    if (c == EOF) {
      fprintf(stderr,"Error: A comment, but...no sequence!??\n");
      return(NULL);
    }
  }
  /* Determine sequence length */
  (void)fgetpos(file,&fpos);
  for(length=0;((c=fgetc(file)) != (int)'>') && (bool)(c != EOF); )
    if(!isspace(c))
      ++length;

  /* Copy data into sequence pointer */
  if ((seq=(char*)malloc(sizeof(char)*length))!=NULL) {
    for((void)fsetpos(file,&fpos),length=0;
        ((c=fgetc(file)) != (int)'>') && (bool)(c != EOF); )
      if(!isspace(c))
        seq[length++]=(char)c;

    *seqlength=length;
    return(seq);
  }
  else {
    fprintf(stderr,"Error: Out of memory!\n");
    return(NULL);
  }
} /* end get_fasta */
