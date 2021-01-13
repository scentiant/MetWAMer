/* simplex.h
 * Michael E Sparks (mespar1@gmail.com)
 *
 * Code to facilitate solution of general linear programming
 * problems using Dantzig's simplex algorithm, coupled with
 * Bland's rule to escape the (rare) cycling condition.  This
 * implementation is based on the exposition given in Chapter
 * 29 of CLRS.
 *
 * Copyright (C) 2006 Michael E Sparks
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

#ifndef SIMPLEX_H
#define SIMPLEX_H

#define SUCCESS    0 /* LP solution solved                        */
#define UNBOUNDED  1 /* LP problem's objective function unbounded */
#define INFEASIBLE 2 /* Original LP problem infeasible            */

/* This routine assumes that all arguments passed (except status,  *
 * which is modified directly and indicates whether the LP problem *
 * specified was optimized successfully or was determined to be    *
 * either infeasible or having an unbounded objective function,    *
 * and result, which stores information on the optimal solution    *
 * and must be passed as NULL) are properly initialized for some   *
 * instance of a linear programming problem. They will not be      *
 * modified in this code. Also, totvarct must equal m+n+1 (+1 b/c  *
 * the auxilliary objective function is set to -x_0; all other     *
 * variables must be indexed starting with 1. The function will    *
 * return a pointer to an array of totvarct+1 doubles, of which    *
 * the first totvarct elements correspond to the coefficients of   *
 * variables x_0 up through x_{totvarct-1}, and the totvarct'th    *
 * element will give the objective function's optimal value--if    *
 * for some reason the LP problem was not computable, it returns   *
 * NULL.                                                           */
double *simplex(int totvarct,int actvarstart,int *Ntmp,double **Atmp,
  double *btmp,double *ctmp,double *vtmp,int m,int n,
  int *status,double *result);

#endif
