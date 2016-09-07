
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


/*
	2nd header file for Meschach's complex routines.
	This file contains declarations for complex factorisation/solve
	routines.

*/


#ifndef ZMeMATRIX2H
#define ZMeMATRIX2H

#include "zmatrix.h"

#ifdef ANSI_C
extern ZMeVEC	*zUsolve(ZMeMAT *matrix, ZMeVEC *b, ZMeVEC *out, double diag);
extern ZMeVEC	*zLsolve(ZMeMAT *matrix, ZMeVEC *b, ZMeVEC *out, double diag);
extern ZMeVEC	*zUAsolve(ZMeMAT *U, ZMeVEC *b, ZMeVEC *out, double diag);
extern ZMeVEC	*zDsolve(ZMeMAT *A, ZMeVEC *b, ZMeVEC *x);
extern ZMeVEC	*zLAsolve(ZMeMAT *L, ZMeVEC *b, ZMeVEC *out, double diag);

extern ZMeVEC	*zhhvec(ZMeVEC *,int,Real *,ZMeVEC *,complex *);
extern ZMeVEC	*zhhtrvec(ZMeVEC *,double,int,ZMeVEC *,ZMeVEC *);
extern ZMeMAT	*zhhtrrows(ZMeMAT *,int,int,ZMeVEC *,double);
extern ZMeMAT	*zhhtrcols(ZMeMAT *,int,int,ZMeVEC *,double);
extern ZMeMAT	*_zhhtrcols(ZMeMAT *,int,int,ZMeVEC *,double,ZMeVEC *);
extern ZMeMAT     *zHfactor(ZMeMAT *,ZMeVEC *);
extern ZMeMAT     *zHQunpack(ZMeMAT *,ZMeVEC *,ZMeMAT *,ZMeMAT *);

extern ZMeMAT	*zQRfactor(ZMeMAT *A, ZMeVEC *diag);
extern ZMeMAT	*zQRCPfactor(ZMeMAT *A, ZMeVEC *diag, PERM *px);
extern ZMeVEC	*_zQsolve(ZMeMAT *QR, ZMeVEC *diag, ZMeVEC *b, ZMeVEC *x, ZMeVEC *tmp);
extern ZMeMAT	*zmakeQ(ZMeMAT *QR, ZMeVEC *diag, ZMeMAT *Qout);
extern ZMeMAT	*zmakeR(ZMeMAT *QR, ZMeMAT *Rout);
extern ZMeVEC	*zQRsolve(ZMeMAT *QR, ZMeVEC *diag, ZMeVEC *b, ZMeVEC *x);
extern ZMeVEC	*zQRAsolve(ZMeMAT *QR, ZMeVEC *diag, ZMeVEC *b, ZMeVEC *x);
extern ZMeVEC	*zQRCPsolve(ZMeMAT *QR,ZMeVEC *diag,PERM *pivot,ZMeVEC *b,ZMeVEC *x);
extern ZMeVEC	*zUmlt(ZMeMAT *U, ZMeVEC *x, ZMeVEC *out);
extern ZMeVEC	*zUAmlt(ZMeMAT *U, ZMeVEC *x, ZMeVEC *out);
extern double	zQRcondest(ZMeMAT *QR);

extern ZMeVEC	*zLsolve(ZMeMAT *, ZMeVEC *, ZMeVEC *, double);
extern ZMeMAT	*zset_col(ZMeMAT *, int, ZMeVEC *);

extern ZMeMAT	*zLUfactor(ZMeMAT *A, PERM *pivot);
extern ZMeVEC	*zLUsolve(ZMeMAT *A, PERM *pivot, ZMeVEC *b, ZMeVEC *x);
extern ZMeVEC	*zLUAsolve(ZMeMAT *LU, PERM *pivot, ZMeVEC *b, ZMeVEC *x);
extern ZMeMAT	*zm_inverse(ZMeMAT *A, ZMeMAT *out);
extern double	zLUcondest(ZMeMAT *LU, PERM *pivot);

extern void	zgivens(complex, complex, Real *, complex *);
extern ZMeMAT	*zrot_rows(ZMeMAT *A, int i, int k, double c, complex s,
			   ZMeMAT *out);
extern ZMeMAT	*zrot_cols(ZMeMAT *A, int i, int k, double c, complex s,
			   ZMeMAT *out);
extern ZMeVEC	*rot_zvec(ZMeVEC *x, int i, int k, double c, complex s,
			  ZMeVEC *out);
extern ZMeMAT	*zschur(ZMeMAT *A,ZMeMAT *Q);
/* extern ZMeMAT	*schur_vecs(ZMeMAT *T,ZMeMAT *Q,X_re,X_im) */
#else
extern ZMeVEC	*zUsolve(), *zLsolve(), *zUAsolve(), *zDsolve(), *zLAsolve();

extern ZMeVEC	*zhhvec();
extern ZMeVEC	*zhhtrvec();
extern ZMeMAT	*zhhtrrows();
extern ZMeMAT     *zhhtrcols();
extern ZMeMAT     *_zhhtrcols();
extern ZMeMAT     *zHfactor();
extern ZMeMAT     *zHQunpack();


extern ZMeMAT	*zQRfactor(), *zQRCPfactor();
extern ZMeVEC	*_zQsolve();
extern ZMeMAT	*zmakeQ(), *zmakeR();
extern ZMeVEC	*zQRsolve(), *zQRAsolve(), *zQRCPsolve();
extern ZMeVEC	*zUmlt(), *zUAmlt();
extern double	zQRcondest();

extern ZMeVEC	*zLsolve();
extern ZMeMAT	*zset_col();

extern ZMeMAT	*zLUfactor();
extern ZMeVEC	*zLUsolve(), *zLUAsolve();
extern ZMeMAT	*zm_inverse();
extern double	zLUcondest();

extern void	zgivens();
extern ZMeMAT	*zrot_rows(), *zrot_cols();
extern ZMeVEC	*rot_zvec();
extern ZMeMAT	*zschur();
/* extern ZMeMAT	*schur_vecs(); */
#endif /* ANSI_C */

#endif /* ZMeMATRIX2H */

