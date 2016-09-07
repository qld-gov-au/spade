
/**************************************************************************
**
** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/


/* iter.h  14/09/93 */

/* 

  Structures for iterative methods

*/

#ifndef ITERHH

#define ITERHH

/* RCS id: $Id: iter.h,v 1.2 1994/03/08 05:48:27 des Exp $  */


#include	"sparse.h"


/* basic structure for iterative methods */

/* type Fun_Ax for functions to get y = A*x */
#ifdef ANSI_C
typedef MeVEC  *(*Fun_Ax)(void *,MeVEC *,MeVEC *);
#else
typedef MeVEC *(*Fun_Ax)();
#endif


/* type ITER */
typedef struct Iter_data {
   int shared_x;   /* if TRUE then x is shared and it will not be free'd */ 
   int shared_b;   /* if TRUE then b is shared and it will not be free'd */
   unsigned k;   /* no. of direction (search) vectors; =0 - none */
   int limit;    /* upper bound on the no. of iter. steps */
   int steps;    /* no. of iter. steps done */
   Real eps;     /* accuracy required */
   
   MeVEC *x;       /* input: initial guess;
		    output: approximate solution */
   MeVEC *b;       /* right hand side of the equation A*x = b */

   Fun_Ax   Ax;		 /* function computing y = A*x */
   void *A_par;         /* parameters for Ax */

   Fun_Ax  ATx;		 /* function  computing y = A^T*x;
					       T = transpose */
   void *AT_par;         /* parameters for ATx */

   Fun_Ax  Bx; /* function computing y = B*x; B - preconditioner */
   void *B_par;         /* parameters for Bx */

   Fun_Ax  BTx; /* function computing y = B^T*x; B - preconditioner */
   void *BT_par;         /* parameters for BTx */

#ifdef ANSI_C

#ifdef PROTOTYPES_IN_STRUCT
   void (*info)(struct Iter_data *, double, MeVEC *,MeVEC *);
            /* function giving some information for a user;
	       nres - a norm of a residual res */
   
   int (*stop_crit)(struct Iter_data *, double, MeVEC *,MeVEC *);
           /* stopping criterion:
	      nres - a norm of res;
	      res - residual;
	    if returned value == TRUE then stop;
	    if returned value == FALSE then continue; */
#else
   void (*info)();
   int  (*stop_crit)();
#endif /* PROTOTYPES_IN_STRUCT */

#else

   void (*info)();
            /* function giving some information for a user */
   
   int (*stop_crit)();
           /* stopping criterion:
	    if returned value == TRUE then stop;
	    if returned value == FALSE then continue; */

#endif /* ANSI_C */

   Real init_res;   /* the norm of the initial residual */

}  ITER;


#define INULL   (ITER *)NULL

/* type Fun_info */
#ifdef ANSI_C
typedef void (*Fun_info)(ITER *, double, MeVEC *,MeVEC *);
#else
typedef void (*Fun_info)();
#endif

/* type Fun_stp_crt */
#ifdef ANSI_C
typedef int (*Fun_stp_crt)(ITER *, double, MeVEC *,MeVEC *);
#else
typedef int (*Fun_stp_crt)();
#endif



/* macros */
/* default values */

#define ITER_LIMIT_DEF  1000
#define ITER_EPS_DEF    1e-6

/* other macros */

/* set ip->Ax=fun and ip->A_par=fun_par */
#define iter_Ax(ip,fun,fun_par) \
  (ip->Ax=(Fun_Ax)(fun),ip->A_par=(void *)(fun_par),0)
#define iter_ATx(ip,fun,fun_par) \
  (ip->ATx=(Fun_Ax)(fun),ip->AT_par=(void *)(fun_par),0)
#define iter_Bx(ip,fun,fun_par) \
  (ip->Bx=(Fun_Ax)(fun),ip->B_par=(void *)(fun_par),0)
#define iter_BTx(ip,fun,fun_par) \
  (ip->BTx=(Fun_Ax)(fun),ip->BT_par=(void *)(fun_par),0)

/* save free macro */
#define ITER_FREE(ip)  (iter_free(ip), (ip)=(ITER *)NULL)


/* prototypes from iter0.c */

#ifdef ANSI_C
/* standard information */
void iter_std_info(const ITER *ip,double nres,MeVEC *res,MeVEC *Bres);
/* standard stopping criterion */
int iter_std_stop_crit(const ITER *ip, double nres, MeVEC *res,MeVEC *Bres);

/* get, resize and free ITER variable */
ITER *iter_get(int lenb, int lenx);
ITER *iter_resize(ITER *ip,int lenb,int lenx);
int iter_free(ITER *ip);

void iter_dump(FILE *fp,ITER *ip);

/* copy ip1 to ip2 copying also elements of x and b */
ITER *iter_copy(const ITER *ip1, ITER *ip2);
/* copy ip1 to ip2 without copying elements of x and b */
ITER *iter_copy2(ITER *ip1,ITER *ip2);

/* functions for generating sparse matrices with random elements */
SPMeMAT	*iter_gen_sym(int n, int nrow);
SPMeMAT	*iter_gen_nonsym(int m,int n,int nrow,double diag);
SPMeMAT	*iter_gen_nonsym_posdef(int n,int nrow);

#else

void iter_std_info();
int iter_std_stop_crit();
ITER *iter_get();
int iter_free();
ITER *iter_resize();
void iter_dump();
ITER *iter_copy();
ITER *iter_copy2();
SPMeMAT	*iter_gen_sym();
SPMeMAT	*iter_gen_nonsym();
SPMeMAT	*iter_gen_nonsym_posdef();

#endif

/* prototypes from iter.c */

/* different iterative procedures */
#ifdef ANSI_C
MeVEC  *iter_cg(ITER *ip);
MeVEC  *iter_cg1(ITER *ip);
MeVEC  *iter_spcg(SPMeMAT *A,SPMeMAT *LLT,MeVEC *b,double eps,MeVEC *x,int limit,
		int *steps);
MeVEC  *iter_cgs(ITER *ip,MeVEC *r0);
MeVEC  *iter_spcgs(SPMeMAT *A,SPMeMAT *B,MeVEC *b,MeVEC *r0,double eps,MeVEC *x,
		 int limit, int *steps);
MeVEC  *iter_lsqr(ITER *ip);
MeVEC  *iter_splsqr(SPMeMAT *A,MeVEC *b,double tol,MeVEC *x,
		  int limit,int *steps);
MeVEC  *iter_gmres(ITER *ip);
MeVEC  *iter_spgmres(SPMeMAT *A,SPMeMAT *B,MeVEC *b,double tol,MeVEC *x,int k,
		   int limit, int *steps);
MeMAT  *iter_arnoldi_iref(ITER *ip,Real *h,MeMAT *Q,MeMAT *H);
MeMAT  *iter_arnoldi(ITER *ip,Real *h,MeMAT *Q,MeMAT *H);
MeMAT  *iter_sparnoldi(SPMeMAT *A,MeVEC *x0,int k,Real *h,MeMAT *Q,MeMAT *H);
MeVEC  *iter_mgcr(ITER *ip);
MeVEC  *iter_spmgcr(SPMeMAT *A,SPMeMAT *B,MeVEC *b,double tol,MeVEC *x,int k,
		  int limit, int *steps);
void	iter_lanczos(ITER *ip,MeVEC *a,MeVEC *b,Real *beta2,MeMAT *Q);
void    iter_splanczos(SPMeMAT *A,int m,MeVEC *x0,MeVEC *a,MeVEC *b,Real *beta2,
		       MeMAT *Q);
MeVEC  *iter_lanczos2(ITER *ip,MeVEC *evals,MeVEC *err_est);
MeVEC  *iter_splanczos2(SPMeMAT *A,int m,MeVEC *x0,MeVEC *evals,MeVEC *err_est);
MeVEC  *iter_cgne(ITER *ip);
MeVEC  *iter_spcgne(SPMeMAT *A,SPMeMAT *B,MeVEC *b,double eps,MeVEC *x,
		  int limit,int *steps);
#else
MeVEC  *iter_cg();
MeVEC  *iter_cg1();
MeVEC  *iter_spcg();
MeVEC  *iter_cgs();
MeVEC  *iter_spcgs();
MeVEC  *iter_lsqr();
MeVEC  *iter_splsqr();
MeVEC  *iter_gmres();
MeVEC  *iter_spgmres();
MeMAT  *iter_arnoldi_iref();
MeMAT  *iter_arnoldi();
MeMAT  *iter_sparnoldi();
MeVEC  *iter_mgcr();
MeVEC  *iter_spmgcr();
void  iter_lanczos();
void  iter_splanczos();
MeVEC  *iter_lanczos2();
MeVEC  *iter_splanczos2();
MeVEC  *iter_cgne();
MeVEC  *iter_spcgne();

#endif


#endif  /* ITERHH */
