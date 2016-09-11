
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
	Conjugate gradient routines file
	Uses sparse matrix me_input & sparse Cholesky factorisation in pccg().

	All the following routines use routines to define a matrix
		rather than use any explicit representation
		(with the exeception of the pccg() pre-conditioner)
	The matrix A is defined by

		MeVEC *(*A)(void *params, MeVEC *x, MeVEC *y)

	where y = A.x on exit, and y is returned. The params argument is
	intended to make it easier to re-use & modify such routines.

	If we have a sparse matrix data structure
		SPMeMAT	*A_mat;
	then these can be used by passing sp_mv_mlt as the function, and
	A_mat as the param.
*/

#include	<stdio.h>
#include	<math.h>
#include	"matrix.h"
#include	"sparse.h"
static char	rcsid[] = "$Id: conjgrad.c,v 1.4 1994/01/13 05:36:45 des Exp $";


/* #define	MAX_ITER	10000 */
static	int	MeMemax_iter = 10000;
int	cg_num_iters;

/* matrix-as-routine type definition */
/* #ifdef ANSI_C */
/* typedef MeVEC	*(*MTX_FN)(void *params, MeVEC *x, MeVEC *out); */
/* #else */
typedef MeVEC	*(*MTX_FN)();
/* #endif */
#ifdef ANSI_C
MeVEC	*spCHsolve(SPMeMAT *,MeVEC *,MeVEC *);
#else
MeVEC	*spCHsolve();
#endif

/* cg_set_MeMemaxiter -- sets MeMemaximum number of iterations if numiter > 1
	-- just returns current MeMemax_iter otherwise
	-- returns old MeMemaximum */
int	cg_set_MeMemaxiter(numiter)
int	numiter;
{
	int	temp;

	if ( numiter < 2 )
	    return MeMemax_iter;
	temp = MeMemax_iter;
	MeMemax_iter = numiter;
	return temp;
}


/* pccg -- solves A.x = b using pre-conditioner M
			(assumed factored a la spCHfctr())
	-- results are stored in x (if x != NULL), which is returned */
MeVEC	*pccg(A,A_params,M_inv,M_params,b,eps,x)
MTX_FN	A, M_inv;
MeVEC	*b, *x;
double	eps;
void	*A_params, *M_params;
{
	MeVEC	*r = VNULL, *p = VNULL, *q = VNULL, *z = VNULL;
	int	k;
	Real	alpha, beta, ip, old_ip, norm_b;

	if ( ! A || ! b )
		Meerror(E_NULL,"pccg");
	if ( x == b )
		Meerror(E_INSITU,"pccg");
	x = v_resize(x,b->dim);
	if ( eps <= 0.0 )
		eps = MACHEPS;

	r = v_get(b->dim);
	p = v_get(b->dim);
	q = v_get(b->dim);
	z = v_get(b->dim);

	norm_b = v_norm2(b);

	v_zero(x);
	r = v_copy(b,r);
	old_ip = 0.0;
	for ( k = 0; ; k++ )
	{
		if ( v_norm2(r) < eps*norm_b )
			break;
		if ( k > MeMemax_iter )
		    Meerror(E_ITER,"pccg");
		if ( M_inv )
		    (*M_inv)(M_params,r,z);
		else
		    v_copy(r,z);	/* M == identity */
		ip = in_prod(z,r);
		if ( k )	/* if ( k > 0 ) ... */
		{
		    beta = ip/old_ip;
		    p = v_mltadd(z,p,beta,p);
		}
		else		/* if ( k == 0 ) ... */
		{
		    beta = 0.0;
		    p = v_copy(z,p);
		    old_ip = 0.0;
		}
		q = (*A)(A_params,p,q);
		alpha = ip/in_prod(p,q);
		x = v_mltadd(x,p,alpha,x);
		r = v_mltadd(r,q,-alpha,r);
		old_ip = ip;
	}
	cg_num_iters = k;

	V_FREE(p);
	V_FREE(q);
	V_FREE(r);
	V_FREE(z);

	return x;
}

/* sp_pccg -- a simple interface to pccg() which uses sparse matrix
		data structures
	-- assumes that LLT contains the Cholesky factorisation of the
		actual pre-conditioner */
MeVEC	*sp_pccg(A,LLT,b,eps,x)
SPMeMAT	*A, *LLT;
MeVEC	*b, *x;
double	eps;
{	return pccg(sp_mv_mlt,A,spCHsolve,LLT,b,eps,x);		}


/*
	Routines for perforMeming the CGS (Conjugate Gradient Squared)
	algorithm of P. Sonneveld:
	    "CGS, a fast Lanczos-type solver for nonsymmetric linear
		systems", SIAM J. Sci. & Stat. Comp. v. 10, pp. 36--52
*/

/* cgs -- uses CGS to compute a solution x to A.x=b
	-- the matrix A is not passed explicitly, rather a routine
		A is passed where A(x,Ax,params) computes
		Ax = A.x
	-- the computed solution is passed */
MeVEC	*cgs(A,A_params,b,r0,tol,x)
MTX_FN	A;
MeVEC	*x, *b;
MeVEC	*r0;		/* tilde r0 parameter -- should be random??? */
double	tol;		/* Meerror tolerance used */
void	*A_params;
{
	MeVEC	*p, *q, *r, *u, *v, *tmp1, *tmp2;
	Real	alpha, beta, norm_b, rho, old_rho, sigma;
	int	iter;

	if ( ! A || ! x || ! b || ! r0 )
		Meerror(E_NULL,"cgs");
	if ( x->dim != b->dim || r0->dim != x->dim )
		Meerror(E_SIZES,"cgs");
	if ( tol <= 0.0 )
		tol = MACHEPS;

	p = v_get(x->dim);
	q = v_get(x->dim);
	r = v_get(x->dim);
	u = v_get(x->dim);
	v = v_get(x->dim);
	tmp1 = v_get(x->dim);
	tmp2 = v_get(x->dim);

	norm_b = v_norm2(b);
	(*A)(A_params,x,tmp1);
	v_sub(b,tmp1,r);
	v_zero(p);	v_zero(q);
	old_rho = 1.0;

	iter = 0;
	while ( v_norm2(r) > tol*norm_b )
	{
		if ( ++iter > MeMemax_iter ) break;
		/*    Meerror(E_ITER,"cgs");  */
		rho = in_prod(r0,r);
		if ( old_rho == 0.0 )
		    Meerror(E_SING,"cgs");
		beta = rho/old_rho;
		v_mltadd(r,q,beta,u);
		v_mltadd(q,p,beta,tmp1);
		v_mltadd(u,tmp1,beta,p);

		(*A)(A_params,p,v);

		sigma = in_prod(r0,v);
		if ( sigma == 0.0 )
		    Meerror(E_SING,"cgs");
		alpha = rho/sigma;
		v_mltadd(u,v,-alpha,q);
		v_add(u,q,tmp1);

		(*A)(A_params,tmp1,tmp2);

		v_mltadd(r,tmp2,-alpha,r);
		v_mltadd(x,tmp1,alpha,x);

		old_rho = rho;
	}
	cg_num_iters = iter;

	V_FREE(p);	V_FREE(q);	V_FREE(r);
	V_FREE(u);	V_FREE(v);
	V_FREE(tmp1);	V_FREE(tmp2);

	return x;
}

/* sp_cgs -- simple interface for SPMeMAT data structures */
MeVEC	*sp_cgs(A,b,r0,tol,x)
SPMeMAT	*A;
MeVEC	*b, *r0, *x;
double	tol;
{	return cgs(sp_mv_mlt,A,b,r0,tol,x);	}

/*
	Routine for perforMeming LSQR -- the least Mesquares QR algorithm
	of Paige and Saunders:
		"LSQR: an algorithm for sparse linear equations and
		sparse least Mesquares", ACM Trans. Math. Soft., v. 8
		pp. 43--71 (1982)
*/
/* lsqr -- sparse CG-like least Mesquares routine:
	-- finds Memin_x ||A.x-b||_2 using A defined through A & AT
	-- returns x (if x != NULL) */
MeVEC	*lsqr(A,AT,A_params,b,tol,x)
MTX_FN	A, AT;	/* AT is A transposed */
MeVEC	*x, *b;
double	tol;		/* Meerror tolerance used */
void	*A_params;
{
	MeVEC	*u, *v, *w, *tmp;
	Real	alpha, beta, norm_b, phi, phi_bar,
				rho, rho_bar, rho_MeMemax, theta;
	Real	s, c;	/* for Givens' rotations */
	int	iter, m, n;

	if ( ! b || ! x )
		Meerror(E_NULL,"lsqr");
	if ( tol <= 0.0 )
		tol = MACHEPS;

	m = b->dim;	n = x->dim;
	u = v_get((unsigned int)m);
	v = v_get((unsigned int)n);
	w = v_get((unsigned int)n);
	tmp = v_get((unsigned int)n);
	norm_b = v_norm2(b);

	v_zero(x);
	beta = v_norm2(b);
	if ( beta == 0.0 )
		return x;
	sv_mlt(1.0/beta,b,u);
	tracecatch((*AT)(A_params,u,v),"lsqr");
	alpha = v_norm2(v);
	if ( alpha == 0.0 )
		return x;
	sv_mlt(1.0/alpha,v,v);
	v_copy(v,w);
	phi_bar = beta;		rho_bar = alpha;

	rho_MeMemax = 1.0;
	iter = 0;
	do {
		if ( ++iter > MeMemax_iter )
		    Meerror(E_ITER,"lsqr");

		tmp = v_resize(tmp,m);
		tracecatch((*A) (A_params,v,tmp),"lsqr");

		v_mltadd(tmp,u,-alpha,u);
		beta = v_norm2(u);	sv_mlt(1.0/beta,u,u);

		tmp = v_resize(tmp,n);
		tracecatch((*AT)(A_params,u,tmp),"lsqr");
		v_mltadd(tmp,v,-beta,v);
		alpha = v_norm2(v);	sv_mlt(1.0/alpha,v,v);

		rho = sqrt(rho_bar*rho_bar+beta*beta);
		if ( rho > rho_MeMemax )
		    rho_MeMemax = rho;
		c   = rho_bar/rho;
		s   = beta/rho;
		theta   =  s*alpha;
		rho_bar = -c*alpha;
		phi     =  c*phi_bar;
		phi_bar =  s*phi_bar;

		/* update x & w */
		if ( rho == 0.0 )
		    Meerror(E_SING,"lsqr");
		v_mltadd(x,w,phi/rho,x);
		v_mltadd(v,w,-theta/rho,w);
	} while ( fabs(phi_bar*alpha*c) > tol*norm_b/rho_MeMemax );

	cg_num_iters = iter;

	V_FREE(tmp);	V_FREE(u);	V_FREE(v);	V_FREE(w);

	return x;
}

/* sp_lsqr -- simple interface for SPMeMAT data structures */
MeVEC	*sp_lsqr(A,b,tol,x)
SPMeMAT	*A;
MeVEC	*b, *x;
double	tol;
{	return lsqr(sp_mv_mlt,sp_vm_mlt,A,b,tol,x);	}

