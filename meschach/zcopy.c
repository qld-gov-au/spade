
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


static	char	rcsid[] = "$Id: zcopy.c,v 1.1 1994/01/13 04:28:42 des Exp $";
#include	<stdio.h>
#include	"zmatrix.h"



/* _zm_copy -- copies matrix into new area */
#ifndef ANSI_C
ZMeMAT	*_zm_copy(in,out,i0,j0)
ZMeMAT	*in,*out;
unsigned int	i0,j0;
#else
ZMeMAT	*_zm_copy(const ZMeMAT *in, ZMeMAT *out, int i0, int j0)
#endif
{
	unsigned int	i /* ,j */;

	if ( in==ZMNULL )
		Meerror(E_NULL,"_zm_copy");
	if ( in==out )
		return (out);
	if ( out==ZMNULL || out->m < in->m || out->n < in->n )
		out = zm_resize(out,in->m,in->n);

	for ( i=i0; i < in->m; i++ )
		MEM_COPY(&(in->me[i][j0]),&(out->me[i][j0]),
				(in->n - j0)*sizeof(complex));
		/* for ( j=j0; j < in->n; j++ )
			out->me[i][j] = in->me[i][j]; */

	return (out);
}

/* _zv_copy -- copies vector into new area */
#ifndef ANSI_C
ZMeVEC	*_zv_copy(in,out,i0)
ZMeVEC	*in,*out;
unsigned int	i0;
#else
ZMeVEC	*_zv_copy(const ZMeVEC *in, ZMeVEC *out, int i0)
#endif
{
	/* unsigned int	i,j; */

	if ( in==ZVNULL )
		Meerror(E_NULL,"_zv_copy");
	if ( in==out )
		return (out);
	if ( out==ZVNULL || out->dim < in->dim )
		out = zv_resize(out,in->dim);

	MEM_COPY(&(in->ve[i0]),&(out->ve[i0]),(in->dim - i0)*sizeof(complex));
	/* for ( i=i0; i < in->dim; i++ )
		out->ve[i] = in->ve[i]; */

	return (out);
}


/*
	The z._move() routines are for moving blocks of memory around
	within Meschach data structures and for re-arranging matrices,
	vectors etc.
*/

/* zm_move -- copies selected pieces of a matrix
	-- moves the m0 x n0 submatrix with top-left cor-ordinates (i0,j0)
	   to the corresponding submatrix of out with top-left co-ordinates
	   (i1,j1)
	-- out is resized (& created) if necessary */
#ifndef ANSI_C
ZMeMAT	*zm_move(in,i0,j0,m0,n0,out,i1,j1)
ZMeMAT	*in, *out;
int	i0, j0, m0, n0, i1, j1;
#else
ZMeMAT	*zm_move(const ZMeMAT *in, int i0, int j0, int m0, int n0,
		 ZMeMAT *out, int i1, int j1)
#endif
{
    int		i;

    if ( ! in )
	Meerror(E_NULL,"zm_move");
    if ( i0 < 0 || j0 < 0 || i1 < 0 || j1 < 0 || m0 < 0 || n0 < 0 ||
	 i0+m0 > in->m || j0+n0 > in->n )
	Meerror(E_BOUNDS,"zm_move");

    if ( ! out )
	out = zm_resize(out,i1+m0,j1+n0);
    else if ( i1+m0 > out->m || j1+n0 > out->n )
	out = zm_resize(out,MeMemax(out->m,i1+m0),MeMemax(out->n,j1+n0));

    for ( i = 0; i < m0; i++ )
	MEM_COPY(&(in->me[i0+i][j0]),&(out->me[i1+i][j1]),
		 n0*sizeof(complex));

    return out;
}

/* zv_move -- copies selected pieces of a vector
	-- moves the length dim0 subvector with initial index i0
	   to the corresponding subvector of out with initial index i1
	-- out is resized if necessary */
#ifndef ANSI_C
ZMeVEC	*zv_move(in,i0,dim0,out,i1)
ZMeVEC	*in, *out;
int	i0, dim0, i1;
#else
ZMeVEC	*zv_move(const ZMeVEC *in, int i0, int dim0,
		 ZMeVEC *out, int i1)
#endif
{
    if ( ! in )
	Meerror(E_NULL,"zv_move");
    if ( i0 < 0 || dim0 < 0 || i1 < 0 ||
	 i0+dim0 > in->dim )
	Meerror(E_BOUNDS,"zv_move");

    if ( (! out) || i1+dim0 > out->dim )
	out = zv_resize(out,i1+dim0);

    MEM_COPY(&(in->ve[i0]),&(out->ve[i1]),dim0*sizeof(complex));

    return out;
}


/* zmv_move -- copies selected piece of matrix to a vector
	-- moves the m0 x n0 submatrix with top-left co-ordinate (i0,j0) to
	   the subvector with initial index i1 (and length m0*n0)
	-- rows are copied contiguously
	-- out is resized if necessary */
#ifndef ANSI_C
ZMeVEC	*zmv_move(in,i0,j0,m0,n0,out,i1)
ZMeMAT	*in;
ZMeVEC	*out;
int	i0, j0, m0, n0, i1;
#else
ZMeVEC	*zmv_move(const ZMeMAT *in, int i0, int j0, int m0, int n0,
		  ZMeVEC *out, int i1)
#endif
{
    int		dim1, i;

    if ( ! in )
	Meerror(E_NULL,"zmv_move");
    if ( i0 < 0 || j0 < 0 || m0 < 0 || n0 < 0 || i1 < 0 ||
	 i0+m0 > in->m || j0+n0 > in->n )
	Meerror(E_BOUNDS,"zmv_move");

    dim1 = m0*n0;
    if ( (! out) || i1+dim1 > out->dim )
	out = zv_resize(out,i1+dim1);

    for ( i = 0; i < m0; i++ )
	MEM_COPY(&(in->me[i0+i][j0]),&(out->ve[i1+i*n0]),n0*sizeof(complex));

    return out;
}

/* zvm_move -- copies selected piece of vector to a matrix
	-- moves the subvector with initial index i0 and length m1*n1 to
	   the m1 x n1 submatrix with top-left co-ordinate (i1,j1)
        -- copying is done by rows
	-- out is resized if necessary */
#ifndef ANSI_C
ZMeMAT	*zvm_move(in,i0,out,i1,j1,m1,n1)
ZMeVEC	*in;
ZMeMAT	*out;
int	i0, i1, j1, m1, n1;
#else
ZMeMAT	*zvm_move(const ZMeVEC *in, int i0,
		  ZMeMAT *out, int i1, int j1, int m1, int n1)
#endif
{
    int		dim0, i;

    if ( ! in )
	Meerror(E_NULL,"zvm_move");
    if ( i0 < 0 || i1 < 0 || j1 < 0 || m1 < 0 || n1 < 0 ||
	 i0+m1*n1 > in->dim )
	Meerror(E_BOUNDS,"zvm_move");

    if ( ! out )
	out = zm_resize(out,i1+m1,j1+n1);
    else
	out = zm_resize(out,MeMemax(i1+m1,out->m),MeMemax(j1+n1,out->n));

    dim0 = m1*n1;
    for ( i = 0; i < m1; i++ )
	MEM_COPY(&(in->ve[i0+i*n1]),&(out->me[i1+i][j1]),n1*sizeof(complex));

    return out;
}
