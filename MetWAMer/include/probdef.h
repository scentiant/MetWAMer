/* probdef.h
 *
 * Michael Sparks <mespar1@gmail.com>
 * Last modified : 17 August 2007
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

/* Kozak sequence:
 *
 * A  C  C  A  T  G  G
 *      -1 +1
 *
 * G or A works well in position -3, G at +4.
 * (See Fig 17.1 of Weaver's Molecular Biology, 2ed.)
 *
 * Per a seminar given by Nikolai Alexandrov, 25 July 2007,
 * in maize, this was seen to be [GA]CCATGGCG
 */

#ifndef PROBDEF_H
#define PROBDEF_H

#include <immpractical.h>

#define INITVAL_FLT    0.0
#define INITVAL_INT      0
#define PROBMIN     0.0500 /* Min prob in PMF dists of WAM */
#define MAXFAULT    0.0005
#define MAXPROB     1.0000 /* Max prob in PMF dists of WAM */

/* The following definitions are already specified *
 * as needed in the immpractical header.           */
//#define EQUIPROB    0.2500
//#define NULLPROB    0.0000

#endif 
