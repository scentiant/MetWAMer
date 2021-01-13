/* print_MetWAM.c
 * Michael Sparks (mespar1@gmail.com)
 * Last modified : 20 July 2013
 *
 * Prints values stored in the specified ATGparm object to stderr.
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

#define ARGCT 2
#define USAGE "\a\nUsage: %s x.MetWAM\n\n"

/* Main Application */
int main(int argc, char *argv[]) {
  ATGparm locATGparms; /* Stores ATGparm object */
  FILE *infile;        /* File pointer          */
  int i,j,k;           /* iterator variables    */

  /* Verify command line */
  if(argc!=ARGCT) {
    fprintf(stderr,USAGE,argv[0]);
    exit(EXIT_FAILURE);
  }

  /* Read in the parameter object */
  if((infile=fopen(argv[1],"rb"))==NULL) {
    (void)snprintf(errbuff,MAXERRLEN,
      "Can't open %s for binary reading!\n",argv[1]);
    FATALERROR(errbuff)
  }
  else if(fread(&locATGparms,sizeof(ATGparm),1,infile)!=1) {
    (void)snprintf(errbuff,MAXERRLEN,
      "Error importing ATGparm data from %s!\n",argv[1]);
    FATALERROR(errbuff)
  }
  else
    fclose(infile);

  /* Print values to stderr */
  for(i=0;i<STRINGSIZE;++i) {
    for(j=0;j<NTALFSIZE;++j)
      for(k=0;k<NTALFSIZE;++k) {
        (void)snprintf(errbuff,MAXERRLEN,"%.4f",
          locATGparms.WAMtable[i][j][k]);
        NONFATALERROR(errbuff)
        if(k<NTALFSIZE-1)
          NONFATALERROR(" ")
        else
          NONFATALERROR("\n")
      }
    NONFATALERROR("\n")
  }

  return(EXIT_SUCCESS);
} /* end main */
