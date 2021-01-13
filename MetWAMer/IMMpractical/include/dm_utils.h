/* dm_utils.h
 * Michael E Sparks (mespar1@gmail.com)
 *
 * This file contains code for developing/using dynamically modulating
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

#ifndef DM_UTILS_H
#define DM_UTILS_H

#include "immpractical.h"

#define dm_probmakerUSAGE \
  "\a\nUsage: %s %i num_partitions\n\
  in-coding.part(1) in-coding.part(2) ... in-coding.part(num_partitions)\n\
  in-noncod.part(1) in-noncod.part(2) ... in-noncod.part(num_partitions) \
out_file\n\n"

/* This function oversees the building of a trained dmprobT *
 * file based on dynamically modulating Markov models.      */
int dm_probmaker(int argc,char *argv[]);

#endif
