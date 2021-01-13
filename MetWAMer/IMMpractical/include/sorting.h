/* sorting.h
 * Michael E Sparks (mespar1@gmail.com)
 * 
 * This file contains code implementing sorting routines needed
 * by functions in dimm_utils.c
 *
 * Copyright (C) 2005  Michael E Sparks
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

#ifndef SORTING_H
#define SORTING_H

/* This definition is the element to make arrays for storing *
 * history counts in a key-value manner for quicksorting.    */
typedef struct ind_list {
  int key[MAXORDER], /* For keying by history (0..MAXORDER -> 5'..3') */
      ct_pret;       /* Store the key's (history's) count             */
} ind_listT;

/* Function to record tally data into an ind_listT array */
void copy_data(int order,int model,
  countarrayT *talH,ind_listT *index_list);

/* Here's a straightforward implementation of *
 * quicksort for ind_listT arrays.            */
void quicksort(ind_listT *ind_list,int left,int right);

#endif
