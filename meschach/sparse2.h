
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


/* Sparse matrix factorise/solve header */
/* RCS id: $Id: sparse2.h,v 1.4 1994/01/13 05:33:46 des Exp $ */



#ifndef SPARSE2H

#define SPARSE2H

#include "sparse.h"


#ifdef ANSI_C
SPMeMAT	*spCHfactor(SPMeMAT *A), *spICHfactor(SPMeMAT *A), *spCHsymb(SPMeMAT *A);
MeVEC	*spCHsolve(SPMeMAT *CH, const MeVEC *b, MeVEC *x);

SPMeMAT	*spLUfactor(SPMeMAT *A,PERM *pivot,double threshold);
SPMeMAT	*spILUfactor(SPMeMAT *A,double theshold);
MeVEC	*spLUsolve(const SPMeMAT *LU,PERM *pivot, const MeVEC *b,MeVEC *x),
	*spLUTsolve(SPMeMAT *LU,PERM *pivot, const MeVEC *b,MeVEC *x);

SPMeMAT	*spBKPfactor(SPMeMAT *, PERM *, PERM *, double);
MeVEC	*spBKPsolve(SPMeMAT *, PERM *, PERM *, const MeVEC *, MeVEC *);

MeVEC	*pccg(MeVEC *(*A)(),void *A_par,MeVEC *(*M_inv)(),void *M_par,MeVEC *b,
						double tol,MeVEC *x);
MeVEC	*sp_pccg(SPMeMAT *,SPMeMAT *,MeVEC *,double,MeVEC *);
MeVEC	*cgs(MeVEC *(*A)(),void *A_par,MeVEC *b,MeVEC *r0,double tol,MeVEC *x);
MeVEC	*sp_cgs(SPMeMAT *,MeVEC *,MeVEC *,double,MeVEC *);
MeVEC	*lsqr(MeVEC *(*A)(),MeVEC *(*AT)(),void *A_par,MeVEC *b,double tol,MeVEC *x);
MeVEC	*sp_lsqr(SPMeMAT *,MeVEC *,double,MeVEC *);
int	cg_set_MeMemaxiter(int);

void	lanczos(MeVEC *(*A)(),void *A_par,int m,MeVEC *x0,MeVEC *a,MeVEC *b,
						Real *beta_m1,MeMAT *Q);
void	sp_lanczos(SPMeMAT *,int,MeVEC *,MeVEC *,MeVEC *,Real *,MeMAT *);
MeVEC	*lanczos2(MeVEC *(*A)(),void *A_par,int m,MeVEC *x0,MeVEC *evals,
						MeVEC *err_est);
MeVEC	*sp_lanczos2(SPMeMAT *,int,MeVEC *,MeVEC *,MeVEC *);
extern  void    scan_to(SPMeMAT *,IMeVEC *,IMeVEC *,IMeVEC *,int);
extern  row_elt  *chase_col(const SPMeMAT *,int,int *,int *,int);
extern  row_elt  *chase_past(const SPMeMAT *,int,int *,int *,int);
extern  row_elt  *bump_col(const SPMeMAT *,int,int *,int *);

#else
extern SPMeMAT	*spCHfactor(), *spICHfactor(), *spCHsymb();
extern MeVEC	*spCHsolve();

extern SPMeMAT	*spLUfactor();
extern SPMeMAT	*spILUfactor();
extern MeVEC	*spLUsolve(), *spLUTsolve();

extern SPMeMAT	*spBKPfactor();
extern MeVEC	*spBKPsolve();

extern MeVEC	*pccg(), *sp_pccg(), *cgs(), *sp_cgs(), *lsqr(), *sp_lsqr();
extern int	cg_set_MeMemaxiter();

void	lanczos(), sp_lanczos();
MeVEC	*lanczos2(), *sp_lanczos2();
extern  void    scan_to();
extern  row_elt  *chase_col();
extern  row_elt  *chase_past();
extern  row_elt  *bump_col();

#endif


#endif
