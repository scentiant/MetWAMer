/* probdef.h
 *
 * Michael Sparks <mespar1@iastate.edu>
 * Last modified : 3 April 2007
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

#ifndef PROBDEF_H
#define PROBDEF_H

#define ALFSIZE        4  /* Cardinality of nucleotide alphabet  */
#define HYPOTHESIS2    2  /* true (T) or false (F) sites         */
#define HYPOTHESIS7    7  /* T1, T2, T0, F1, F2, F0 and Fi       */
#define STRINGSIZE   102  /* Dinuc + 50nt flanks                 */
#define TERMCOUNT      2  /* Donor, Acceptor                     */
#define MAXSPLICESIG  50  /* Maximal splice site window extent   */
#define WINSIZE      100  /* 50nt on the left and right of GT|AG */

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

typedef struct {
  int hypothesisnum, /* number of hypotheses, either *
                      * HYPOTHESIS2 or HYPOTHESIS7   */
      wsizedonorleft, 
      wsizedonorright,
      wsizeacceptorleft,
      wsizeacceptorright;
  union {
    float hypo2table[TERMCOUNT][HYPOTHESIS2][WINSIZE+2][ALFSIZE][ALFSIZE];
    float hypo7table[TERMCOUNT][HYPOTHESIS7][WINSIZE+2][ALFSIZE][ALFSIZE];
  } hypotables;
} Bssmmodel;
    
typedef struct {
  int gtmodelset, /* Set to either TRUE or FALSE */
      gcmodelset; 
  Bssmmodel gtmodel,
            gcmodel;
} Bssmparm;

#endif 
