/* train_MC.c
 * Michael E Sparks (mespar1@gmail.com)
 *
 * This utility allows the user to interface with the IMMpractical
 * library, allowing for development of Markov chain transition
 * probability matrices.
 *
 * Copyright (c) 2007 Michael E Sparks
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

/* Included Libraries */
#include <immpractical.h>
#include <stdio.h>
#include <stdlib.h>

/* Main application */
int main(int argc,char *argv[])
{
  int test;

  /* maintrain is defined in immpractical.c */
  if((test=maintrain(argc,argv))==TRUE)
    return(EXIT_SUCCESS);
  else
    return(EXIT_FAILURE);
}
