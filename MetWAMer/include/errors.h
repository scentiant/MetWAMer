/* errors.h
 * Michael E Sparks (michael.sparks2@usda.gov)
 * Last modified: 21 December 2020
 *   (increased MAXERRLEN value in response to format-truncation
 *    warnings issued by gcc version 10.1.0. MES)
 *
 * Macros for advising about errors and (potentially)
 * correcting them.
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

#ifndef ERRORS_H
#define ERRORS_H

#define MAXERRLEN 1048576 /* Max length of error message + 1 */

char errbuff[MAXERRLEN]; /* (Global) buffer for storing error messages */

/* Print a descriptive message for a recoverable error. *
 * This macro can be called if the user simply wants    *
 * to print a message to stderr, too.                   */
#define NONFATALERROR(message) { \
  fprintf(stderr,"%s",message); \
  fflush(stderr); \
}

/* Print a descriptive message for a nonrecoverable error, *
 * then terminate the program.                             */
#define FATALERROR(message) { \
  fprintf(stderr,"%s",message); \
  exit(EXIT_FAILURE); \
}

#endif
