/* dimm_utils.h
 * Michael E Sparks (mespar1@gmail.com)
 *
 * This file contains code specific to training deleted interpolated
 * Markov models.  Top-down and bottom-up deleted interpolated
 * Markov models are implemented.
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

#ifndef DIMM_UTILS_H
#define DIMM_UTILS_H
 
#include "immpractical.h"

/* Subroutine to coordinate development of final, smoothed oligo *
 * likelihoods by the top-down or bottom-up training methods.    *
 * Set topdown to TRUE for top-down, FASLE, for bottom-up.       *
 * Note that ind_listT structures are defined in sorting.h       */
void DIMM_final_probs(immprobT *finalprob,relfreqsT *rfreqs,
  countarrayT *talD,countarrayT *talC,int model,int algo);

#endif
