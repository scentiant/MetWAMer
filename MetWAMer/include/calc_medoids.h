/* calc_medoids.h
 * Michael Sparks (mespar1@gmail.com)
 * Last modified : 17 August 2007
 *
 * Functions exported from the calc_medoids implementation file.
 *
 * Copyright (c) 2007 Michael E Sparks
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

#ifndef CALC_MEDOIDS_H
#define CALC_MEDOIDS_H

/* Computes the distance between two data items, *
 * namely the number of nucleotides differing    *
 * between the pair, i.e., edit distance.        */
double eval_dist_kernel(
  char *dat1,
  char *dat2
);

#endif
