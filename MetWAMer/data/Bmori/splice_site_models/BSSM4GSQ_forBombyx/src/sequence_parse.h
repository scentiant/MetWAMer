/* sequence_parse.h
 * Michael Sparks (mespar1@iastate.edu)
 * Last modified : 12 July 2004
 *
 * Code for parsing FASTA formatted sequences.
 *
 * Please Note: Transliteration scheme
 *   A   -> 0
 *   C   -> 1
 *   G   -> 2
 *   T/U -> 3
 *
 * Copyright (c) 2004 Michael E Sparks
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

#ifndef SEQUENCE_PARSE_H
#define SEQUENCE_PARSE_H

/* String literals */
#ifndef TRUE
#define TRUE  1
#endif

#ifndef FALSE
#define FALSE 0
#endif

/* get_fasta:  A general function to parse FASTA entries  *
 * Input :                                                *
 *   A file pointer: This must be connected to a stream   *
 *                   containing Fasta-formatted sequences *
 *   An int pointer: This must have sufficient space      *
 *                   allocated to store the integer       *
 *                   translation of the sequence in the   *
 *                   file.                                *
 *   An int pointer: This is a single integer.  The       *
 *                   function will store the length of    *
 *                   the current sequence in it.          *
 *   const int     : size of alphabet, e.g, 4 for nt's    *
 * Output : Returns TRUE if a sequence has been           *
 *            encountered.                                *
 *          Returns FALSE if the end of file has been     *
 *            reached.                                    */
int get_fasta(FILE *file,int *seq,int *seqlength,const int alphsize);

/* Returns a random value in the interval [0,ALFSIZE-1] */
int random_residue(const int alphsize);

#endif
