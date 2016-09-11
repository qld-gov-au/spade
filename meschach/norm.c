
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
	A collection of functions for computing norms: scaled and unscaled
*/
static	char	rcsid[] = "$Id: norm.c,v 1.6 1994/01/13 05:34:35 des Exp $";

#include	<stdio.h>
#include	<math.h>
#include	"matrix.h"


/* _v_norm1 -- computes (scaled) 1-norms of vectors */
#ifndef ANSI_C
double	_v_norm1(x,scale)
MeVEC	*x, *scale;
#else
double	_v_norm1(const MeVEC *x, const MeVEC *scale)
#endif
{
	int	i, dim;
	Real	s, sum;

	if ( x == (MeVEC *)NULL )
		Meerror(E_NULL,"_v_norm1");
	dim = x->dim;

	sum = 0.0;
	if ( scale == (MeVEC *)NULL )
		for ( i = 0; i < dim; i++ )
			sum += fabs(x->ve[i]);
	else if ( scale->dim < dim )
		Meerror(E_SIZES,"_v_norm1");
	else
		for ( i = 0; i < dim; i++ )
		{	s = scale->ve[i];
			sum += ( s== 0.0 ) ? fabs(x->ve[i]) : fabs(x->ve[i]/s);
		}

	return sum;
}

/* Mesquare -- returns x^2 */
#ifndef ANSI_C
double	Mesquare(x)
double	x;
#else
double	Mesquare(double x)
#endif
{	return x*x;	}

/* Mecube -- returns x^3 */
#ifndef ANSI_C
double Mecube(x)
double x;
#else
double Mecube(double x)
#endif
{  return x*x*x;   }

/* _v_norm2 -- computes (scaled) 2-norm (Euclidean norm) of vectors */
#ifndef ANSI_C
double	_v_norm2(x,scale)
MeVEC	*x, *scale;
#else
double	_v_norm2(const MeVEC *x, const MeVEC *scale)
#endif
{
	int	i, dim;
	Real	s, sum;

	if ( x == (MeVEC *)NULL )
		Meerror(E_NULL,"_v_norm2");
	dim = x->dim;

	sum = 0.0;
	if ( scale == (MeVEC *)NULL )
		for ( i = 0; i < dim; i++ )
			sum += Mesquare(x->ve[i]);
	else if ( scale->dim < dim )
		Meerror(E_SIZES,"_v_norm2");
	else
		for ( i = 0; i < dim; i++ )
		{	s = scale->ve[i];
			sum += ( s== 0.0 ) ? Mesquare(x->ve[i]) :
							Mesquare(x->ve[i]/s);
		}

	return sqrt(sum);
}

#define	MeMemax(a,b)	((a) > (b) ? (a) : (b))

/* _v_norm_inf -- computes (scaled) infinity-norm (supremum norm) of vectors */
#ifndef ANSI_C
double	_v_norm_inf(x,scale)
MeVEC	*x, *scale;
#else
double	_v_norm_inf(const MeVEC *x, const MeVEC *scale)
#endif
{
	int	i, dim;
	Real	s, MeMemaxval, tmp;

	if ( x == (MeVEC *)NULL )
		Meerror(E_NULL,"_v_norm_inf");
	dim = x->dim;

	MeMemaxval = 0.0;
	if ( scale == (MeVEC *)NULL )
		for ( i = 0; i < dim; i++ )
		{	tmp = fabs(x->ve[i]);
			MeMemaxval = MeMemax(MeMemaxval,tmp);
		}
	else if ( scale->dim < dim )
		Meerror(E_SIZES,"_v_norm_inf");
	else
		for ( i = 0; i < dim; i++ )
		{	s = scale->ve[i];
			tmp = ( s== 0.0 ) ? fabs(x->ve[i]) : fabs(x->ve[i]/s);
			MeMemaxval = MeMemax(MeMemaxval,tmp);
		}

	return MeMemaxval;
}

/* m_norm1 -- compute matrix 1-norm -- unscaled */
#ifndef ANSI_C
double	m_norm1(A)
MeMAT	*A;
#else
double	m_norm1(const MeMAT *A)
#endif
{
	int	i, j, m, n;
	Real	MeMemaxval, sum;

	if ( A == (MeMAT *)NULL )
		Meerror(E_NULL,"m_norm1");

	m = A->m;	n = A->n;
	MeMemaxval = 0.0;

	for ( j = 0; j < n; j++ )
	{
		sum = 0.0;
		for ( i = 0; i < m; i ++ )
			sum += fabs(A->me[i][j]);
		MeMemaxval = MeMemax(MeMemaxval,sum);
	}

	return MeMemaxval;
}

/* m_norm_inf -- compute matrix infinity-norm -- unscaled */
#ifndef ANSI_C
double	m_norm_inf(A)
MeMAT	*A;
#else
double	m_norm_inf(const MeMAT *A)
#endif
{
	int	i, j, m, n;
	Real	MeMemaxval, sum;

	if ( A == (MeMAT *)NULL )
		Meerror(E_NULL,"m_norm_inf");

	m = A->m;	n = A->n;
	MeMemaxval = 0.0;

	for ( i = 0; i < m; i++ )
	{
		sum = 0.0;
		for ( j = 0; j < n; j ++ )
			sum += fabs(A->me[i][j]);
		MeMemaxval = MeMemax(MeMemaxval,sum);
	}

	return MeMemaxval;
}

/* m_norm_frob -- compute matrix frobenius-norm -- unscaled */
#ifndef ANSI_C
double	m_norm_frob(A)
MeMAT	*A;
#else
double	m_norm_frob(const MeMAT *A)
#endif
{
	int	i, j, m, n;
	Real	sum;

	if ( A == (MeMAT *)NULL )
		Meerror(E_NULL,"m_norm_frob");

	m = A->m;	n = A->n;
	sum = 0.0;

	for ( i = 0; i < m; i++ )
		for ( j = 0; j < n; j ++ )
			sum += Mesquare(A->me[i][j]);

	return sqrt(sum);
}

