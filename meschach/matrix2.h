
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
	Header file for ``matrix2.a'' library file
*/


#ifndef MeMATRIX2H
#define MeMATRIX2H

#include "matrix.h"

/* Unless otherwise specified, factorisation routines overwrite the
   matrix that is being factorised */

#ifndef ANSI_C

extern	MeMAT	*BKPfactor(), *CHfactor(), *LUfactor(), *QRfactor(),
		*QRCPfactor(), *LDLfactor(), *Hfactor(), *MCHfactor(),
		*m_inverse();
extern	double	LUcondest(), QRcondest();
extern	MeMAT	*makeQ(), *makeR(), *makeHQ(), *makeH();
extern	MeMAT	*LDLupdate(), *QRupdate();

extern	MeVEC	*BKPsolve(), *CHsolve(), *LUsolve(), *_Qsolve(), *QRsolve(),
		*LDLsolve(), *Usolve(), *Lsolve(), *Dsolve(), *LTsolve(),
		*UTsolve(), *LUTsolve(), *QRCPsolve();

extern  BAND    *bdLUfactor(), *bdLDLfactor();
extern  MeVEC     *bdLUsolve(), *bdLDLsolve();

extern	MeVEC	*hhvec();
extern	MeVEC	*hhtrvec();
extern	MeMAT	*hhtrrows();
extern	MeMAT	*hhtrcols(), *_hhtrcols();

extern	void	givens();
extern	MeVEC	*rot_vec();	/* in situ */
extern	MeMAT	*rot_rows();	/* in situ */
extern	MeMAT	*rot_cols();	/* in situ */


/* eigenvalue routines */
extern	MeVEC	*trieig(), *symmeig();
extern	MeMAT	*schur();
extern	void	schur_evals();
extern	MeMAT	*schur_vecs();

/* singular value decomposition */
extern	MeVEC	*bisvd(), *svd();

/* matrix powers and exponent */
MeMAT  *_m_pow();
MeMAT  *m_pow();
MeMAT  *m_exp(), *_m_exp();
MeMAT  *m_poly();

/* FFT */
void fft();
void ifft();


#else

                 /* forms Bunch-Kaufman-Parlett factorisation for
                        symmetric indefinite matrices */
extern	MeMAT	*BKPfactor(MeMAT *A,PERM *pivot,PERM *blocks),
                 /* Cholesky factorisation of A
                        (symmetric, positive definite) */
		*CHfactor(MeMAT *A),
                /* LU factorisation of A (with partial pivoting) */ 
                *LUfactor(MeMAT *A,PERM *pivot),
                /* QR factorisation of A; need dim(diag) >= # rows of A */
		*QRfactor(MeMAT *A,MeVEC *diag),
                /* QR factorisation of A with column pivoting */
		*QRCPfactor(MeMAT *A,MeVEC *diag,PERM *pivot),
                /* L.D.L^T factorisation of A */
		*LDLfactor(MeMAT *A), 
                /* Hessenberg factorisation of A -- for schur() */
                *Hfactor(MeMAT *A,MeVEC *diag1,MeVEC *diag2),
                /* modified Cholesky factorisation of A;
                        actually factors A+D, D diagonal with no
                        diagonal entry in the factor < sqrt(tol) */
                *MCHfactor(MeMAT *A,double tol),
		*m_inverse(const MeMAT *A,MeMAT *out);

                /* returns condition estimate for A after LUfactor() */
extern	double	LUcondest(const MeMAT *A, PERM *pivot),
                /* returns condition estimate for Q after QRfactor() */
                QRcondest(const MeMAT *A);

/* Note: The make..() and ..update() routines assume that the factorisation
        has already been carried out */

     /* Qout is the "Q" (orthongonal) matrix from QR factorisation */
extern	MeMAT	*makeQ(const MeMAT *QR,const MeVEC *diag,MeMAT *Qout),
                /* Rout is the "R" (upper triangular) matrix
                        from QR factorisation */
		*makeR(const MeMAT *A,MeMAT *Rout),
                /* Qout is orthogonal matrix in Hessenberg factorisation */
		*makeHQ(MeMAT *A,MeVEC *diag1,MeVEC *diag2,MeMAT *Qout),
                /* Hout is the Hessenberg matrix in Hessenberg factorisation */
		*makeH(const MeMAT *A,MeMAT *Hout);

                /* updates L.D.L^T factorisation for A <- A + alpha.u.u^T */
extern	MeMAT	*LDLupdate(MeMAT *A,MeVEC *u,double alpha),
                /* updates QR factorisation for QR <- Q.(R+u.v^T)
		   Note: we need explicit Q & R matrices,
                        from makeQ() and makeR() */
		*QRupdate(MeMAT *Q,MeMAT *R,MeVEC *u,MeVEC *v);

/* Solve routines assume that the corresponding factorisation routine
        has already been applied to the matrix along with auxiliary
        objects (such as pivot permutations)

        These solve the system A.x = b,
        except for LUTsolve and QRTsolve which solve the transposed system
                                A^T.x. = b.
        If x is NULL on entry, then it is created.
*/

extern	MeVEC	*BKPsolve(const MeMAT *A,PERM *pivot,const PERM *blocks,
			  const MeVEC *b,MeVEC *x),
		*CHsolve(const MeMAT *A,const MeVEC *b,MeVEC *x),
		*LDLsolve(const MeMAT *A,const MeVEC *b,MeVEC *x),
		*LUsolve(const MeMAT *A, PERM *pivot, const MeVEC *b,MeVEC *x),
		*_Qsolve(const MeMAT *A, const MeVEC *diag, const MeVEC *b, 
			 MeVEC *x, MeVEC *tmp),
		*QRsolve(const MeMAT *A, const MeVEC *diag, const MeVEC *b,MeVEC *x),
    		*QRTsolve(const MeMAT *A,const MeVEC *,const MeVEC *b,MeVEC *x),


     /* Triangular equations solve routines;
        U for upper triangular, L for lower traingular, D for diagonal
        if diag_val == 0.0 use that values in the matrix */

		*Usolve(const MeMAT *A,const MeVEC *b,MeVEC *x,double diag_val),
		*Lsolve(const MeMAT *A,const MeVEC *b,MeVEC *x,double diag_val),
		*Dsolve(const MeMAT *A,const MeVEC *b,MeVEC *x),
		*LTsolve(const MeMAT *A,const MeVEC *b,MeVEC *x,double diag_val),
		*UTsolve(const MeMAT *A,const MeVEC *b,MeVEC *x,double diag_val),
                *LUTsolve(const MeMAT *A,PERM *pivot,const MeVEC *b, MeVEC *x),
                *QRCPsolve(const MeMAT *QR,const MeVEC *diag,PERM *pivot,
			   const MeVEC *b,MeVEC *x);

extern  BAND    *bdLUfactor(BAND *A,PERM *pivot),
                *bdLDLfactor(BAND *A);
extern  MeVEC     *bdLUsolve(const BAND *A,PERM *pivot,const MeVEC *b,MeVEC *x),
                *bdLDLsolve(const BAND *A,const MeVEC *b,MeVEC *x);



extern	MeVEC	*hhvec(const MeVEC *,unsigned int,Real *,MeVEC *,Real *);
extern	MeVEC	*hhtrvec(const MeVEC *,double,unsigned int,const MeVEC *,MeVEC *);
extern	MeMAT	*hhtrrows(MeMAT *,unsigned int,unsigned int,const MeVEC *,double);
extern	MeMAT	*hhtrcols(MeMAT *,unsigned int,unsigned int,const MeVEC *,double);
extern	MeMAT	*_hhtrcols(MeMAT *,unsigned int,unsigned int,const MeVEC *,double,MeVEC *);

extern	void	givens(double,double,Real *,Real *);
extern	MeVEC	*rot_vec(const MeVEC *,unsigned int,unsigned int,
			 double,double,MeVEC *); /* in situ */
extern	MeMAT	*rot_rows(const MeMAT *,unsigned int,unsigned int,
			  double,double,MeMAT *); /* in situ */
extern	MeMAT	*rot_cols(const MeMAT *,unsigned int,unsigned int,
			  double,double,MeMAT *); /* in situ */


/* eigenvalue routines */

               /* compute eigenvalues of tridiagonal matrix
                  with diagonal entries a[i], super & sub diagonal entries
                  b[i]; eigenvectors stored in Q (if not NULL) */
extern	MeVEC	*trieig(MeVEC *a,MeVEC *b,MeMAT *Q),
                 /* sets out to be vector of eigenvectors; eigenvectors
                   stored in Q (if not NULL). A is unchanged */
		*symmeig(const MeMAT *A,MeMAT *Q,MeVEC *out);

               /* computes real Schur form = Q^T.A.Q */
extern	MeMAT	*schur(MeMAT *A,MeMAT *Q);
         /* computes real and imaginary parts of the eigenvalues
                        of A after schur() */
extern	void	schur_evals(MeMAT *A,MeVEC *re_part,MeVEC *im_part);
          /* computes real and imaginary parts of the eigenvectors
                        of A after schur() */
extern	MeMAT	*schur_vecs(MeMAT *T,MeMAT *Q,MeMAT *X_re,MeMAT *X_im);


/* singular value decomposition */

        /* computes singular values of bi-diagonal matrix with
                   diagonal entries a[i] and superdiagonal entries b[i];
                   singular vectors stored in U and V (if not NULL) */
MeVEC	*bisvd(MeVEC *a,MeVEC *b,MeMAT *U,MeMAT *V),
               /* sets "out" to be vector of singular values;
                   singular vectors stored in U and V */
	*svd(MeMAT *A,MeMAT *U,MeMAT *V,MeVEC *out);

/* matrix powers and exponent */
MeMAT  *_m_pow(const MeMAT *A, int p, MeMAT *tmp,MeMAT *out);
MeMAT  *m_pow(const MeMAT *A, int p, MeMAT *out);
MeMAT  *m_exp(MeMAT *,double,MeMAT *);
MeMAT  *_m_exp(MeMAT *A, double eps, MeMAT *out, int *q_out, int *j_out);
MeMAT  *m_poly(const MeMAT *,const MeVEC *,MeMAT *);

/* FFT */
void fft(MeVEC *,MeVEC *);
void ifft(MeVEC *,MeVEC *);

#endif


#endif
