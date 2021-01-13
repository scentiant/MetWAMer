/* BSSM_printer.c
 * Michael Sparks (mespar1@iastate.edu)
 * Last modified : 3 April 2007
 *
 * This is a driver program for debugging the parameterization routines.
 *
 * Copyright (c) 2003,2004,2007 Michael E Sparks
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
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include "bssm_utils.h"

#define ARGCT 2
#define USAGE "Usage: %s x.bssm\n"

/* Main Application */
int main(int argc, char *argv[]) {
  Bssmparm my_bssmparm; /* Stores model parameterization */
  FILE *f_ptr1;         /* File pointers                 */

  /* Verify command line specification of a parameter file to inspect */
  if (argc != ARGCT) {
    fprintf(stderr,USAGE,argv[0]);
    exit(EXIT_FAILURE);
  }

  /* Determine if the BSSM parameter file *
   * already exists in the filesystem     */
  if ((f_ptr1=fopen(argv[1],"rb"))!=NULL) {
    if (fread(&my_bssmparm,sizeof(Bssmparm),1,f_ptr1) != 1) {
      fprintf(stderr,"Error importing Bssmparm data from %s!\n",argv[1]);
      fclose(f_ptr1);
      exit(EXIT_FAILURE);
    }
  }
  else {
    fprintf(stderr,"Can't open %s for binary reading!\n",argv[1]);
    exit(EXIT_FAILURE);
  }

  /* Close files */
  fclose(f_ptr1);

  /* Print debugging information to stderr */
  echo_bssm(&my_bssmparm);

  return(EXIT_SUCCESS);
} /* end main */
