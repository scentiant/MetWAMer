/* sorting.c
 * Michael E Sparks (mespar1@gmail.com)
 *
 * This file contains code implementing sorting routines needed
 * by functions in dimm_utils.c
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "immpractical.h"
#include "sorting.h"

/* Local function prototypes *************************************************/

/* Function that quicksort depends on */
static int partition(ind_listT *ind_list,int left,int right);

/* Function definitions ******************************************************/

/* Here's a straightforward implementation of *
 * quicksort for ind_listT arrays.            */
void quicksort(ind_listT *ind_list,int left,int right)
{
  int pivot;

  if(left < right) {
    pivot=partition(ind_list,left,right);
    quicksort(ind_list,left,pivot-1);
    quicksort(ind_list,pivot+1,right);
  }
  else
    return;
} /* end quicksort */

/* Function that quicksort depends on */
static int partition(ind_listT *ind_list,int left,int right)
{
  ind_listT temp; /* For swapping elements */
  int compare2me; /* For comparing element *
                   * values (ct_pret).     */
  int i,j;        /* iterator variables    */

  compare2me=ind_list[right].ct_pret;
  i=left-1;

  for(j=left;j<=right-1;++j) {
    if(ind_list[j].ct_pret <= compare2me) {
      temp=ind_list[++i];
      ind_list[i]=ind_list[j];
      ind_list[j]=temp;
    }
  }

  temp=ind_list[i+1];
  ind_list[i+1]=ind_list[right];
  ind_list[right]=temp;

  return(i+1);
} /* end partition */

/* Function to record tally data into an ind_listT array */
void copy_data(int order,int model,countarrayT *talH,
  ind_listT *index_list)
{
  int i,j,k,l,m; /* Iterator variables. */
  const int powof2=(int)pow(ALFSIZE,2.0), /* Useful for indexing */
            powof3=(int)pow(ALFSIZE,3.0), /* array elements.     */
            powof4=(int)pow(ALFSIZE,4.0);

  if(index_list == NULL) {
    fprintf(stderr,"Err (copy_data): called with null ind_listT pointer!\n");
    exit(EXIT_FAILURE);
  }

  switch (order) {
    /* Case 0 need not be handled, as this function (called from *
     * dimm_utils.c) will never be called with order < 1         */
    case 1 :
      for(i=0;i<ALFSIZE;++i) {
        index_list[i].ct_pret=
          talH->count1[model][i];
        index_list[i].key[0]=i;
      }
      break;
    case 2 :
      for(i=0;i<ALFSIZE;++i) {
      for(j=0;j<ALFSIZE;++j) {
        index_list[(ALFSIZE*i)+j].ct_pret=
          talH->count2[model][i][j];
        index_list[(ALFSIZE*i)+j].key[0]=i;
        index_list[(ALFSIZE*i)+j].key[1]=j;
      }}
      break;
    case 3 :
      for(i=0;i<ALFSIZE;++i) {
      for(j=0;j<ALFSIZE;++j) {
      for(k=0;k<ALFSIZE;++k) {
        index_list[(powof2*i)+(ALFSIZE*j)+k].ct_pret=
          talH->count3[model][i][j][k];
        index_list[(powof2*i)+(ALFSIZE*j)+k].key[0]=i;
        index_list[(powof2*i)+(ALFSIZE*j)+k].key[1]=j;
        index_list[(powof2*i)+(ALFSIZE*j)+k].key[2]=k;
      }}}
      break;
    case 4 :
      for(i=0;i<ALFSIZE;++i) {
      for(j=0;j<ALFSIZE;++j) {
      for(k=0;k<ALFSIZE;++k) {
      for(l=0;l<ALFSIZE;++l) {
        index_list[(powof3*i)+(powof2*j)+(ALFSIZE*k)+l].ct_pret=
          talH->count4[model][i][j][k][l];
        index_list[(powof3*i)+(powof2*j)+(ALFSIZE*k)+l].key[0]=i;
        index_list[(powof3*i)+(powof2*j)+(ALFSIZE*k)+l].key[1]=j;
        index_list[(powof3*i)+(powof2*j)+(ALFSIZE*k)+l].key[2]=k;
        index_list[(powof3*i)+(powof2*j)+(ALFSIZE*k)+l].key[3]=l;
      }}}}
      break;
    case 5 :
      for(i=0;i<ALFSIZE;++i) {
      for(j=0;j<ALFSIZE;++j) {
      for(k=0;k<ALFSIZE;++k) {
      for(l=0;l<ALFSIZE;++l) {
      for(m=0;m<ALFSIZE;++m) {
        index_list[(powof4*i)+(powof3*j)+(powof2*k)+(ALFSIZE*l)+m].ct_pret=
          talH->count5[model][i][j][k][l][m];
        index_list[(powof4*i)+(powof3*j)+(powof2*k)+(ALFSIZE*l)+m].key[0]=i;
        index_list[(powof4*i)+(powof3*j)+(powof2*k)+(ALFSIZE*l)+m].key[1]=j;
        index_list[(powof4*i)+(powof3*j)+(powof2*k)+(ALFSIZE*l)+m].key[2]=k;
        index_list[(powof4*i)+(powof3*j)+(powof2*k)+(ALFSIZE*l)+m].key[3]=l;
        index_list[(powof4*i)+(powof3*j)+(powof2*k)+(ALFSIZE*l)+m].key[4]=m;
      }}}}}
      break;
    default :
      fprintf(stderr,"Err (copy_data): Invalid case encountered!\n");
      exit(EXIT_FAILURE);
  } /* end switch */

  return;
} /* end copy_data */
