/*  bssm_utils.h
 *  Michael Sparks (mespar1@iastate.edu)
 *  Last modified : 3 April 2007
 *
 *  This is a collection of functions associated with
 *  manipulating Bssmparm objects for training purposes.
 *
 *  Input data files are--strictly!--named as follows, using
 *  an obvious schema:
 *    F0_don  F1_don  F2_don  Fi_don  T0_don  T1_don  T2_don
 *    F0_acc  F1_acc  F2_acc  Fi_acc  T0_acc  T1_acc  T2_acc
 *
 *  Phase is denoted as follows:
 *    1 -> C O D |
 *    2 -> C | O D
 *    0 -> C O | D
 *
 * Copyright (c) 2003,2004,2006,2007 Michael E Sparks
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

#ifndef BSSM_UTILS_H
#define BSSM_UTILS_H

#include "probdef.h"

#define INITVAL_FLT    0.0
#define INITVAL_INT      0
#define PROBMIN     0.0500
#define EQUIPROB    0.2500
#define MAXPROB     1.0000
#define NULLPROB    0.0000
#define MAXFAULT    0.0005

/* init_bssm : Function to initialize all elements of the *
 *             Bssmparm object to INITVAL                */
void init_bssm(Bssmparm *bssm);

/*  build_bssm : Function to update our BSSM parameterization *
 *               file.                                        */
void build_bssm(int **seq_matrix, int num_entries, Bssmparm *bssm,
                char *type, int dim1, int dim2);

/* parse_dimensions : Function to determine dimensions of the *
 *                    Bssmparm object that should be updated */
void parse_dimensions(char *filename, int *dim1, int *dim2);

/* echo_bssm : Function to assist with debugging by printing *
 *             the model parameterization to stderr stream   */
void echo_bssm(Bssmparm *bssm);

#endif
