/* fo_utils.h
 * Michael E Sparks (mespar1@gmail.com)
 *
 * This file contains code for developing/using fixed-order
 * Markov models.
 *
 * Copyright (C) 2005,2006 Michael E Sparks
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef FO_UTILS_H
#define FO_UTILS_H

#include "immpractical.h"

#define fo_probmakerUSAGE \
  "\a\nUsage: %s %i in-coding in-noncod out_file\n\n"

/* This function oversees the building of a trained *
 * foprobT file based on fixed-order Markov models. */
int fo_probmaker(int argc,char *argv[]);

/* This function adds pseudocounts for oligos that *
 * don't occur in the training data at a great     *
 * enough frequency.  It's a very primitive form   *
 * of parameter smoothing. If scaleall is set to   *
 * TRUE, the complete count matrix will be scaled  *
 * up FOPSEUDOCT units; else, the pseudocount      *
 * induction method will be more surgical.         */
void fo_pseudo(countarrayT *tal,int model,int scaleall);

/* This function will record the likelihood of all   *
 * possible oligomers for fixed-order Markov models. */
void fo_rfreqs(countarrayT *tal,int model,foprobT *fops);

#endif
