/* sequence_parse.h
 * Michael Sparks (mespar1@gmail.com)
 *
 * Contains code for generic functions meant to manipulate
 * biological sequences stored in common format files.
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

#ifndef SEQUENCE_PARSE_H
#define SEQUENCE_PARSE_H

/* get_fasta:  A general function to parse FASTA entries  *
 * Input :                                                *
 *   A file pointer: This must be connected to a stream   *
 *                   containing Fasta-formatted sequences *
 *   A char pointer: This must be set to NULL.  The       *
 *                   function will attempt to allocate    *
 *                   sufficient memory to store the seq   *
 *   An int pointer: This is a single integer.  The       *
 *                   function will store the length of    *
 *                   the current sequence in it.          */
char *get_fasta(FILE *file,char *seq,int *seqlength);

/* Function returns the integer translation of a character     *
 * in the alphabet [AaCcGgTtUu].  Depending on preprocessing   *
 * conditions, any other non-whitespace character encountered  *
 * will either throw and error or a random residue is returned *
 * via random_residue.                                         *
 *                                                             *
 * Transliteration scheme                                      *
 *     A   -> 0                                                *
 *     C   -> 1                                                *
 *     G   -> 2                                                *
 *     T/U -> 3                                                */
int trans(char c);

/* Return complement of argument (base) */
int basecomp(int x);

#endif
