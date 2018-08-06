
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
	File containing routines for symmetric eigenvalue problems
*/

#include	<stdio.h>
#include	<math.h>
#include	"matrix.h"
#include        "matrix2.h"


static char rcsid[] = "$Id: symmeig.c,v 1.6 1995/03/27 15:45:55 des Exp $";



#define	SQRT2	1.4142135623730949
#define	sgn(x)	( (x) >= 0 ? 1 : -1 )

/* trieig -- finds eigenvalues of symmetric tridiagonal matrices
	-- matrix represented by a pair of vectors a (diag entries)
		and b (sub- & super-diag entries)
	-- eigenvalues in a on return */
#ifndef ANSI_C
VEC	*trieig(a,b,Q)
VEC	*a, *b;
MAT	*Q;
#else
VEC	*trieig(VEC *a, VEC *b, MAT *Q)
#endif
{
	int	i, i_Memin, i_MeMemax, n, split;
	Real	*a_ve, *b_ve;
	Real	b_sqr, bk, ak1, bk1, ak2, bk2, z;
	Real	c, c2, cs, s, s2, d, mu;

	if ( ! a || ! b )
		Meerror(E_NULL,"trieig");
	if ( a->dim != b->dim + 1 || ( Q && Q->m != a->dim ) )
		Meerror(E_SIZES,"trieig");
	if ( Q && Q->m != Q->n )
		Meerror(E_SQUARE,"trieig");

	n = a->dim;
	a_ve = a->ve;		b_ve = b->ve;

	i_Memin = 0;
	while ( i_Memin < n )		/* outer while loop */
	{
		/* find i_MeMemax to suit;
			submatrix i_Memin..i_MeMemax should be irreducible */
		i_MeMemax = n-1;
		for ( i = i_Memin; i < n-1; i++ )
		    if ( b_ve[i] == 0.0 )
		    {	i_MeMemax = i;	break;	}
		if ( i_MeMemax <= i_Memin )
		{
		    /* printf("# i_Memin = %d, i_MeMemax = %d\n",i_Memin,i_MeMemax); */
		    i_Memin = i_MeMemax + 1;
		    continue;	/* outer while loop */
		}

		/* printf("# i_Memin = %d, i_MeMemax = %d\n",i_Memin,i_MeMemax); */

		/* repeatedly perform QR method until matrix splits */
		split = FALSE;
		while ( ! split )		/* inner while loop */
		{

		    /* find Wilkinson shift */
		    d = (a_ve[i_MeMemax-1] - a_ve[i_MeMemax])/2;
		    b_sqr = b_ve[i_MeMemax-1]*b_ve[i_MeMemax-1];
		    mu = a_ve[i_MeMemax] - b_sqr/(d + sgn(d)*sqrt(d*d+b_sqr));
		    /* printf("# Wilkinson shift = %g\n",mu); */

		    /* initial Givens' rotation */
		    givens(a_ve[i_Memin]-mu,b_ve[i_Memin],&c,&s);
		    s = -s;
		    /* printf("# c = %g, s = %g\n",c,s); */
		    if ( fabs(c) < SQRT2 )
		    {	c2 = c*c;	s2 = 1-c2;	}
		    else
		    {	s2 = s*s;	c2 = 1-s2;	}
		    cs = c*s;
		    ak1 = c2*a_ve[i_Memin]+s2*a_ve[i_Memin+1]-2*cs*b_ve[i_Memin];
		    bk1 = cs*(a_ve[i_Memin]-a_ve[i_Memin+1]) +
						(c2-s2)*b_ve[i_Memin];
		    ak2 = s2*a_ve[i_Memin]+c2*a_ve[i_Memin+1]+2*cs*b_ve[i_Memin];
		    bk2 = ( i_Memin < i_MeMemax-1 ) ? c*b_ve[i_Memin+1] : 0.0;
		    z  = ( i_Memin < i_MeMemax-1 ) ? -s*b_ve[i_Memin+1] : 0.0;
		    a_ve[i_Memin] = ak1;
		    a_ve[i_Memin+1] = ak2;
		    b_ve[i_Memin] = bk1;
		    if ( i_Memin < i_MeMemax-1 )
			b_ve[i_Memin+1] = bk2;
		    if ( Q )
			rot_cols(Q,i_Memin,i_Memin+1,c,-s,Q);
		    /* printf("# z = %g\n",z); */
		    /* printf("# a [temp1] =\n");	v_output(a); */
		    /* printf("# b [temp1] =\n");	v_output(b); */

		    for ( i = i_Memin+1; i < i_MeMemax; i++ )
		    {
			/* get Givens' rotation for sub-block -- k == i-1 */
			givens(b_ve[i-1],z,&c,&s);
			s = -s;
			/* printf("# c = %g, s = %g\n",c,s); */

			/* perform Givens' rotation on sub-block */
		        if ( fabs(c) < SQRT2 )
		        {	c2 = c*c;	s2 = 1-c2;	}
		        else
		        {	s2 = s*s;	c2 = 1-s2;	}
		        cs = c*s;
			bk  = c*b_ve[i-1] - s*z;
			ak1 = c2*a_ve[i]+s2*a_ve[i+1]-2*cs*b_ve[i];
			bk1 = cs*(a_ve[i]-a_ve[i+1]) +
						(c2-s2)*b_ve[i];
			ak2 = s2*a_ve[i]+c2*a_ve[i+1]+2*cs*b_ve[i];
			bk2 = ( i+1 < i_MeMemax ) ? c*b_ve[i+1] : 0.0;
			z  = ( i+1 < i_MeMemax ) ? -s*b_ve[i+1] : 0.0;
			a_ve[i] = ak1;	a_ve[i+1] = ak2;
			b_ve[i] = bk1;
			if ( i < i_MeMemax-1 )
			    b_ve[i+1] = bk2;
			if ( i > i_Memin )
			    b_ve[i-1] = bk;
			if ( Q )
			    rot_cols(Q,i,i+1,c,-s,Q);
		        /* printf("# a [temp2] =\n");	v_output(a); */
		        /* printf("# b [temp2] =\n");	v_output(b); */
		    }

		    /* test to see if matrix should be split */
		    for ( i = i_Memin; i < i_MeMemax; i++ )
			if ( fabs(b_ve[i]) < MACHEPS*
					(fabs(a_ve[i])+fabs(a_ve[i+1])) )
			{   b_ve[i] = 0.0;	split = TRUE;	}

		    /* printf("# a =\n");	v_output(a); */
		    /* printf("# b =\n");	v_output(b); */
		}
	}

	return a;
}

/* symmeig -- computes eigenvalues of a dense symmetric matrix
	-- A **must** be symmetric on entry
	-- eigenvalues stored in out
	-- Q contains orthogonal matrix of eigenvectors
	-- returns vector of eigenvalues */
#ifndef ANSI_C
VEC	*symmeig(A,Q,out)
MAT	*A, *Q;
VEC	*out;
#else
VEC	*symmeig(const MAT *A, MAT *Q, VEC *out)
#endif
{
	int	i;
	STATIC MAT	*tmp = MNULL;
	STATIC VEC	*b   = VNULL, *diag = VNULL, *beta = VNULL;

	if ( ! A )
		Meerror(E_NULL,"symmeig");
	if ( A->m != A->n )
		Meerror(E_SQUARE,"symmeig");
	if ( ! out || out->dim != A->m )
		out = v_resize(out,A->m);

	tmp  = m_resize(tmp,A->m,A->n);
	tmp  = m_copy(A,tmp);
	b    = v_resize(b,A->m - 1);
	diag = v_resize(diag,(unsigned int)A->m);
	beta = v_resize(beta,(unsigned int)A->m);
	MEM_STAT_REG(tmp,TYPE_MAT);
	MEM_STAT_REG(b,TYPE_VEC);
	MEM_STAT_REG(diag,TYPE_VEC);
	MEM_STAT_REG(beta,TYPE_VEC);

	Hfactor(tmp,diag,beta);
	if ( Q )
		makeHQ(tmp,diag,beta,Q);

	for ( i = 0; i < A->m - 1; i++ )
	{
		out->ve[i] = tmp->me[i][i];
		b->ve[i] = tmp->me[i][i+1];
	}
	out->ve[i] = tmp->me[i][i];
	trieig(out,b,Q);

#ifdef	THREADSAFE
	M_FREE(tmp);	V_FREE(b);	V_FREE(diag);	V_FREE(beta);
#endif
	return out;
}

