/* simplex.c
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "simplex.h"

/* Porting this code to 64-bit platforms caused errors in my test     *
 * cases (Really, it *found* errors in my test cases! :o) ).  To keep *
 * things consistent, rather than use float.h's DBL_EPSILON, I        *
 * define my own constant for determining float equality. If someone  *
 * ever wanted to use this file for some high-precision application,  *
 * one could set this macro to DBL_EPSILON and include float.h.       */
#define FLOATEQL 1E-9 

/* Local function prototypes *************************************************/

/* This routine performs the pivot operation.  Since pointers  *
 * are passed as args, it adjusts the relevant data structures *
 * directly.                                                   */
static void pivot(int startind,int size,int *N,
  double **A,double *b,double *c,
  double *v,int l, int e);

#ifdef SIMPLEXDBG
/* Function to print relevant data structures to stdout *
 * for debugging purposes, should such be necessary.    */
static void print_data(int startind,int size,int *N,
  double **A,double *b,double *c,
  double *v);
#endif

/* This function determines if any of the objective function's *
 * nonbasic variables are .GT. 0, meaning another iteration of *
 * the simplex algorithm is required.                          */
static int scanforpositive(int startind,int size,int *N,double *c);

/* Why is this not in math.h?!! */
static int factorial(int arg);

/* n choose r */
static int combinations(int n,int r);

/* This routine will set status to a default value of SUCCESS, *
 * then return its remaining arguments unmodified if they give *
 * a feasible initial basic solution.  Else, it will attempt   *
 * to derive one, set status accordingly, and return.          */
static void init_simplex(int totvarct,int actvarstart,int *N,
  double **A,double *b,double *c,double *v,int m,int n,
  int *status);

/* Function definitions ******************************************************/

/* This routine performs the pivot operation.  Since pointers  *
 * are passed as args, it adjusts the relevant data structures *
 * directly.                                                   */
static void pivot(int startind,int size,int *N,
  double **A,double *b,double *c,
  double *v,int l, int e) {
  int i,j; /* Iterator variables */

  /* compute coefficients of eqn for new basic var x_e */
  if(fabs(A[l][e]) < FLOATEQL){ /* prevent division-by-zero */
    fprintf(stderr,"Err (pivot): Did not specify A[%i][%i] correctly\n",l,e);
    exit(EXIT_FAILURE);
  }
  if(fabs(b[e]=b[l]/A[l][e])<FLOATEQL)
    b[e]=0.0;
  for(j=startind;j<size;++j)
    if(N[j] && j != e)
      if(fabs(A[e][j]=A[l][j]/A[l][e])<FLOATEQL)
        A[e][j]=0.0;
  if(fabs(A[e][l]=1.0/A[l][e])<FLOATEQL)
    A[e][l]=0.0;

  /* compute coefficients of remaining constraints */
  for(i=startind;i<size;++i) {
    if(!N[i] && i != l) {
      if(fabs(b[i]=b[i]-A[i][e]*b[e])<FLOATEQL)
        b[i]=0.0;
      for(j=startind;j<size;++j)
        if(N[j] && j != e)
          if(fabs(A[i][j]=A[i][j]-A[i][e]*A[e][j])<FLOATEQL)
            A[i][j]=0.0;
      if(fabs(A[i][l]=-A[i][e]*A[e][l])<FLOATEQL)
        A[i][l]=0.0;
    }
  }

  /* compute the objective fct */
  if(fabs(*v=*v+c[e]*b[e])<FLOATEQL)
    *v=0.0;
  for(j=startind;j<size;++j)
    if(N[j] && j != e)
      if(fabs(c[j]=c[j]-c[e]*A[e][j])<FLOATEQL)
        c[j]=0.0;
  if(fabs(c[l]=-c[e]*A[e][l])<FLOATEQL)
    c[l]=0.0;

  /* update N (and, indirectly, B) */
  N[e]=0;
  N[l]=1;

  return;
} /* end pivot */

#ifdef SIMPLEXDBG
/* Function to print relevant data structures to stdout *
 * for debugging purposes, should such be necessary.    */
static void print_data(int startind,int size,int *N,
  double **A,double *b,double *c,
  double *v) {
  int i,j;

  printf("N is\n");
  for(i=startind;i<size;++i)
    printf("%i,",N[i]);
  printf("\n");

  printf("A is\n");
  for(i=startind;i<size;++i) {
    for(j=startind;j<size;++j)
      printf("%.2f,",A[i][j]);
    printf("\n");
  }

  printf("b is\n");
  for(i=startind;i<size;++i)
    printf("%.2f,",b[i]);
  printf("\n");

  printf("c is\n");
  for(i=startind;i<size;++i)
    printf("%.2f,",c[i]);
  printf("\n");

  printf("v is %.2f\n",*v);
  printf("\n");

  return;
} /* end print_data */
#endif

/* This function determines if any of the objective function's *
 * nonbasic variables are .GT. 0, meaning another iteration of *
 * the simplex algorithm is required.                          */
static int scanforpositive(int startind,int size,int *N,double *c){
  int i;

  for(i=startind;i<size;++i)
    /* if nonbasic, != 0.0, and positive */
    if (N[i] && fabs(c[i]) > FLOATEQL && c[i] > 0)
      return 1;

  return 0;
} /* end scanforpositive */

/* Why is this not in math.h?!! */
static int factorial(int arg){
  int result=1,
      i;
  if(arg==0)
    return result;
  else
    for(i=arg;i>1;--i)
      result*=i;
  return result;
} /* end factorial */

/* n choose r */
static int combinations(int n,int r){
  return( factorial(n) / (factorial(n-r)*factorial(r)) );
} /* end combinations */

/* This routine will set status to a default value of SUCCESS, *
 * then return its remaining arguments unmodified if they give *
 * a feasible initial basic solution.  Else, it will attempt   *
 * to derive one, set status accordingly, and return.          */
static void init_simplex(int totvarct,int actvarstart,int *N,
  double **A,double *b,double *c,double *v,int m,int n,
  int *status){
  int *Norig=NULL, /* Store nonbasic vars before calling pivot     */
      feasible=1,  /* boolean variable for determining feasibility */
      e,           /* index of entry variable                      */
      l,           /* index of leaving variable                    */
      iterlimit,   /* limit on iterations until cycling            */
      i,z;         /* iterator variable                            */
  double *cprime,  /* for auxilliary objective function            */
         vprime,   /* value of aux obj fct                         */
         *delta=NULL, /* how does change in a given basic          *
                       * variable impact the obj fct?              */
         lowest;   /* scratch variable for finding least value     */

  /* Things must go wrong for this to change */
  *status=SUCCESS;

  /* Determine if current solution is feasible */
  for(i=actvarstart;i<totvarct;++i)
    if(!N[i])
      if(fabs(b[i]) > FLOATEQL && b[i] < 0)
        feasible=0;

  if(feasible) /* nothing for init_simplex to do */
    return;

  /* find index l of minimal constraint in b, which *
   * will be needed for the first call to pivot.    */
  lowest=INFINITY; /* initial smallest value */
  l=-1;            /* infeasible index value */
  for(i=totvarct-1;i>=actvarstart;--i) {
    if(!N[i]){
      if(b[i]<lowest || fabs(b[i]-lowest)<FLOATEQL){
        lowest=b[i];
        l=i;
      }
    }
  }
  if(l==-1){
    fprintf(stderr,"Err (init_simplex): Detected infeasible lower limit!\n");
    exit(EXIT_FAILURE);
  }

  ++n; /* Added x_0 to the mix */
  --actvarstart; /* will now be 0 */
  if(actvarstart!=0){
    fprintf(stderr,"Err (init_simplex): actvarstart != 0 in aux routine?\n");
    exit(EXIT_FAILURE);
  }

  /* Allocate space for delta array */
  if((delta=(double*)malloc(sizeof(double)*totvarct))==NULL) {
    fprintf(stderr,"Err (simplex): Out of memory!\n");
    exit(EXIT_FAILURE);
  }

  /* Formulate auxilliary linear program */
  N[0]=1;
  c[0]=0.0;
  vprime=0.0;
  if( ((cprime=(double*)malloc(sizeof(double)*totvarct))!=NULL) &&
      ((Norig=(int*)malloc(sizeof(int)*totvarct))!=NULL)
    ) {
    cprime[0]=-1.0; /* z' = -x_0 */
    for(i=1;i<totvarct;++i)
      cprime[i]=0.0;
    /* Norig will be used later in this function */
  }
  else {
    fprintf(stderr,"Err (init_simplex): Out of memory!\n");
    exit(EXIT_FAILURE);
  }
  /* set coefficients for x_0 in the set of constraints */
  for(i=actvarstart;i<totvarct;++i)
    if(!N[i])
      A[i][0]=-1.0;

  /* Note that if .GT. n+m choose m iterations  *
   * occur, the algorithm is cycling. :( This   *
   * should not occur, as Bland's rule is       *
   * implemented below.                         */
  iterlimit=combinations(n+m,m);
  #ifdef SIMPLEXDBG
  printf("Cycling info (in init_simplex): \
n = %i, m = %i, n+m choose m = %i\n\n",
    n,m,iterlimit);
  #endif

  #ifdef SIMPLEXDBG
  /* Print values for debugging purposes */
  printf("Data prior to pivoting on x_0 x_%i element:\n",l);
  print_data(actvarstart,totvarct,N,A,b,cprime,&vprime);
  #endif
  /* Call pivot, exchanging x_l and x_0 */
  pivot(actvarstart,totvarct,N,A,b,cprime,&vprime,l,0);
  #ifdef SIMPLEXDBG
  /* Print values for debugging purposes */
  printf("Data subsequent to pivoting on x_0 x_%i element:\n",l);
  print_data(actvarstart,totvarct,N,A,b,cprime,&vprime);
  #endif

  /* lines 2-11 of the simplex algorithm, modified *
   * for the auxilliary obj fct problem.           */
  for(z=0;scanforpositive(actvarstart,totvarct,N,cprime);++z) {
    if(z>iterlimit){
      fprintf(stderr,"Err (init_simplex): \
Simplex is cycling despite Bland's rule!\n");
      exit(EXIT_FAILURE);
    }

    /* choose an index e in N with cprime_e > 0 */
    for(i=actvarstart,e=-1;i<totvarct;++i) {
      if(N[i] && fabs(cprime[i]) > FLOATEQL && cprime[i] > 0) {
        e=i; /* found column of pivot! */
        break;
      }
    }
    if(e==-1){
      fprintf(stderr,"Err (init_simplex): \
Failed to find e in N with c_e > 0\n");
      exit(EXIT_FAILURE);
    }

    /* foreach index in B */
    for(i=actvarstart;i<totvarct;++i) {
      if(!N[i]) {
        if(fabs(A[i][e]) > FLOATEQL && A[i][e] > 0) { /* Aie > 0 */
          if(fabs(delta[i]=b[i]/A[i][e])<FLOATEQL)
            delta[i]=0.0;
        }
        else
          delta[i]=INFINITY;
        #ifdef SIMPLEXDBG
        printf("delta of %i is %.4f\n",i,delta[i]);
        #endif
      }
    }
    #ifdef SIMPLEXDBG
    printf("\n");
    #endif
    /* find index l of minimal delta element */
    lowest=INFINITY; /* initial smallest value */
    l=-1;            /* infeasible index value */
    /* This variable iterates such that it applies Bland's *
     * rule for cycling prevention of simplex              */
    for(i=totvarct-1;i>=actvarstart;--i) {
      if(!N[i]){
        if(fabs(delta[i]-lowest)<FLOATEQL || delta[i]<lowest){
          lowest=delta[i];
          l=i;
        }
      }
    }
    if(l==-1){
      fprintf(stderr,"Err (init_simplex): Detected infeasible lower limit!\n");
      exit(EXIT_FAILURE);
    }
    #ifdef SIMPLEXDBG
    printf("e is %i, l is %i, delta is %f\n\n",e,l,delta[l]);
    #endif

    /* found row of pivot! */
    if(isinf(delta[l])) {
      #ifdef SIMPLEXDBG
      printf("Auxilliary LP obj fct unbounded.\n");
      #endif
      *status=UNBOUNDED;
      break;
    }
    else {
      /* Make copy of N prior to its modification in pivot */
      for(i=actvarstart;i<totvarct;++i)
        Norig[i]=N[i];

      /* pivot modifies its args as a side effect */
      pivot(actvarstart,totvarct,N,A,b,cprime,&vprime,l,e);

      /* Perform necessary substitutions on original objective function *
       * This can be done with the following pseudocode:                *
       * v = v + c_e * b_e                                              *
       * for j in Norig - {e}                                           *
       *   do c_j = c_j - c_e * a_ej                                    *
       * c_l = -c_e * a_el                                              */
      if(fabs(*v+=c[e]*b[e])<FLOATEQL)
        *v=0.0;
      for(i=actvarstart;i<totvarct;++i)
        if(Norig[i] && i != e)
          if(fabs(c[i]+=-c[e]*A[e][i])<FLOATEQL)
            c[i]=0.0;
      if(fabs(c[l]=-c[e]*A[e][l])<FLOATEQL)
        c[l]=0.0;
    }

    #ifdef SIMPLEXDBG
    /* Print values for debugging purposes */
    print_data(actvarstart,totvarct,N,A,b,cprime,&vprime);
    #endif
  } /* end for */

  if(*status!=UNBOUNDED && !N[0]){ /* implies the original LP *
                                    * problem is infeasible   */
    #ifdef SIMPLEXDBG
    printf("Original LP problem is infeasible.\n");
    #endif
    *status=INFEASIBLE;
  }

  free(delta);
  free(cprime);
  free(Norig);

  return;
} /* end init_simplex */

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
  int *status,double *result){
  int *N=NULL,     /* vector giving nonbasic variable status   *
                    * (and conversely, basic)                  */
      e,           /* index of entry variable                  */
      l,           /* index of leaving variable                */
      iterlimit,   /* limit on iterations until cycling        */
      i,j,z;       /* iterator variables                       */
  double **A=NULL, /* matrix of nonbasic variable coefficients *
                    * in the constraint equations              */
         *b=NULL,  /* vector of constraint limits              */
         *c=NULL,  /* vector of nonbasic variable coefficients *
                    * in the objective function.               */
         v,        /* current solution of the objective fct    */
         *delta=NULL, /* how does change in a given basic      *
                       * variable impact the obj fct?          */
         lowest;   /* scratch variable for finding least value */

  /* Make certain totvarct is set reasonably */
  if(totvarct != (m+n+1)){
    fprintf(stderr,"Err (simplex): totvarct != (m+n+1)\n");
    exit(EXIT_FAILURE);
  }

  /* Make certain result is set to NULL */
  if(result!=NULL) {
    fprintf(stderr,"Err (simplex): Passed result as non-NULL!\n");
    exit(EXIT_FAILURE);
  }

  /* Allocate space for delta array */
  if((delta=(double*)malloc(sizeof(double)*totvarct))==NULL) {
    fprintf(stderr,"Err (simplex): Out of memory!\n");
    exit(EXIT_FAILURE);
  }

  /* initialize local data structures */
  if (
       ((N=(int*)malloc(sizeof(int)*totvarct))!=NULL) &&
       ((A=(double**)malloc(sizeof(double*)*totvarct))!=NULL) &&
       ((b=(double*)malloc(sizeof(double)*totvarct))!=NULL) &&
       ((c=(double*)malloc(sizeof(double)*totvarct))!=NULL)
     ) {
    for(i=0;i<totvarct;++i){
      N[i]=Ntmp[i];
      if((A[i]=(double*)malloc(sizeof(double)*totvarct))!=NULL){
        for(j=0;j<totvarct;++j)
          if(fabs(A[i][j]=Atmp[i][j])<FLOATEQL)
            A[i][j]=0.0;
      }
      else{
        fprintf(stderr,"Err (simplex): Out of memory!\n");
        exit(EXIT_FAILURE);
      }

      if(fabs(b[i]=btmp[i])<FLOATEQL)
        b[i]=0.0;
      if(fabs(c[i]=ctmp[i])<FLOATEQL)
        c[i]=0.0;
    }
    if(fabs(v=*vtmp)<FLOATEQL)
      v=0.0;
  }
  else {
    fprintf(stderr,"Err (simplex): Out of memory!\n");
    exit(EXIT_FAILURE);
  }

  /* First, find a reasonable starting point */
  #ifdef SIMPLEXDBG
  printf("Entering init_simplex...\n\n");
  print_data(actvarstart,totvarct,N,A,b,c,&v);
  #endif
  init_simplex(totvarct,actvarstart,N,A,b,c,&v,m,n,status);
  #ifdef SIMPLEXDBG
  printf("Leaving init_simplex...\n\n");
  #endif

  if(*status != SUCCESS) { /* init_simplex indicated that *
                            * we cannot proceed.          */
    /* release storage */
    free(delta);
    free(N);
    for(i=0;i<totvarct;++i)
      free(A[i]);
    free(A);
    free(b);
    free(c);

    return(NULL);
  }

  /* Note that if .GT. n+m choose m iterations  *
   * occur, the algorithm is cycling. :( This   *
   * should not occur, as Bland's rule is       *
   * implemented below.                         */
  iterlimit=combinations(n+m,m);
  #ifdef SIMPLEXDBG
  printf("Cycling info: n = %i, m = %i, n+m choose m = %i\n\n",
    n,m,iterlimit);
  #endif

  /* lines 2-11 of the simplex algorithm */
  #ifdef SIMPLEXDBG
  print_data(actvarstart,totvarct,N,A,b,c,&v);
  #endif
  for(z=0;scanforpositive(actvarstart,totvarct,N,c);++z) {
    if(z>iterlimit){
      fprintf(stderr,"Err (simplex): \
Simplex is cycling despite Bland's rule!\n");
      exit(EXIT_FAILURE);
    }

    /* choose an index e in N with c_e > 0 */
    for(i=actvarstart,e=-1;i<totvarct;++i) {
      if(N[i] && fabs(c[i]) > FLOATEQL && c[i] > 0) {
        e=i; /* found column of pivot! */
        break;
      }
    }
    if(e==-1){
      fprintf(stderr,"Err (simplex): \
Failed to find e in N with c_e > 0\n");
      exit(EXIT_FAILURE);
    }

    /* foreach index in B */
    for(i=actvarstart;i<totvarct;++i) {
      if(!N[i]) {
        if(fabs(A[i][e]) > FLOATEQL && A[i][e] > 0) { /* Aie > 0 */
          if(fabs(delta[i]=b[i]/A[i][e])<FLOATEQL)
            delta[i]=0.0;
        }
        else
          delta[i]=INFINITY;
        #ifdef SIMPLEXDBG
        printf("delta of %i is %f\n",i,delta[i]);
        #endif
      }
    }
    #ifdef SIMPLEXDBG
    printf("\n");
    #endif
    /* find index l of minimal delta element */
    lowest=INFINITY; /* initial smallest value */
    l=-1;            /* infeasible index value */
    /* This variable iterates such that it applies Bland's *
     * rule for cycling prevention of simplex              */
    for(i=totvarct-1;i>=actvarstart;--i) {
      if(!N[i]){
        if(fabs(delta[i]-lowest)<FLOATEQL || delta[i]<lowest){
          lowest=delta[i];
          l=i;
        }
      }
    }
    if(l==-1){
      fprintf(stderr,"Err (simplex): Detected infeasible lower limit!\n");
      exit(EXIT_FAILURE);
    }
    #ifdef SIMPLEXDBG
    printf("e is %i, l is %i, delta is %f\n\n",e,l,delta[l]);
    #endif

    /* found row of pivot! */
    if(isinf(delta[l])) {
      #ifdef SIMPLEXDBG
      printf("Original LP obj fct unbounded.\n");
      #endif
      *status=UNBOUNDED;
      break;
    }
    else
      /* pivot modifies its args as a side effect */
      pivot(actvarstart,totvarct,N,A,b,c,&v,l,e);

    #ifdef SIMPLEXDBG
    /* Print values for debugging purposes */
    print_data(actvarstart,totvarct,N,A,b,c,&v);
    #endif
  } /* end for */

  /* Now for the readoff */
  if(*status == SUCCESS){
    if((result=(double*)malloc(sizeof(double)*(totvarct+1)))!=NULL) {
      for(i=actvarstart;i<totvarct;++i) {
        if(N[i] || fabs(b[i])<FLOATEQL)
          result[i]=0.0;
        else
          result[i]=b[i];
        #ifdef SIMPLEXDBG
        printf("x_%i is %.2f\n",i,result[i]);
        #endif
      }
      if(fabs(v)<FLOATEQL)
        result[i]=0.0;
      else
        result[i]=v;
      #ifdef SIMPLEXDBG
      printf("Objective fct is maximized at %.2f\n\n",result[i]);
      #endif
    }
    else {
      fprintf(stderr,"Err (simplex): Out of memory!\n");
      exit(EXIT_FAILURE);
    }
  }

  /* release storage */
  free(delta);
  free(N);
  for(i=0;i<totvarct;++i)
    free(A[i]);
  free(A);
  free(b);
  free(c);

  return(result);
} /* end simplex */
