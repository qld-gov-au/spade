
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


/* vecop.c 1.3 8/18/87 */

#include	<stdio.h>
#include	"matrix.h"

static	char	rcsid[] = "$Id: vecop.c,v 1.5 1996/08/20 18:18:10 stewart Exp $";


/* _in_prod -- inner product of two vectors from i0 downwards
	-- that is, returns a(i0:dim)^T.b(i0:dim) */
#ifndef ANSI_C
double	_in_prod(a,b,i0)
MeVEC	*a,*b;
unsigned int	i0;
#else
double	_in_prod(const MeVEC *a, const MeVEC *b, unsigned int i0)
#endif
{
	unsigned int	limit;
	/* Real	*a_v, *b_v; */
	/* register Real	sum; */

	if ( a==(MeVEC *)NULL || b==(MeVEC *)NULL )
		Meerror(E_NULL,"_in_prod");
	limit = Memin(a->dim,b->dim);
	if ( i0 > limit )
		Meerror(E_BOUNDS,"_in_prod");

	return __ip__(&(a->ve[i0]),&(b->ve[i0]),(int)(limit-i0));
	/*****************************************
	a_v = &(a->ve[i0]);		b_v = &(b->ve[i0]);
	for ( i=i0; i<limit; i++ )
		sum += a_v[i]*b_v[i];
		sum += (*a_v++)*(*b_v++);

	return (double)sum;
	******************************************/
}

/* sv_mlt -- scalar-vector multiply -- out <- scalar*vector 
	-- may be in-situ */
#ifndef ANSI_C
MeVEC	*sv_mlt(scalar,vector,out)
double	scalar;
MeVEC	*vector,*out;
#else
MeVEC	*sv_mlt(double scalar, const MeVEC *vector, MeVEC *out)
#endif
{
	/* unsigned int	dim, i; */
	/* Real	*out_ve, *vec_ve; */

	if ( vector==(MeVEC *)NULL )
		Meerror(E_NULL,"sv_mlt");
	if ( out==(MeVEC *)NULL || out->dim != vector->dim )
		out = v_resize(out,vector->dim);
	if ( scalar == 0.0 )
		return v_zero(out);
	if ( scalar == 1.0 )
		return v_copy(vector,out);

	__smlt__(vector->ve,(double)scalar,out->ve,(int)(vector->dim));
	/**************************************************
	dim = vector->dim;
	out_ve = out->ve;	vec_ve = vector->ve;
	for ( i=0; i<dim; i++ )
		out->ve[i] = scalar*vector->ve[i];
		(*out_ve++) = scalar*(*vec_ve++);
	**************************************************/
	return (out);
}

/* v_add -- vector addition -- out <- v1+v2 -- may be in-situ */
#ifndef ANSI_C
MeVEC	*v_add(vec1,vec2,out)
MeVEC	*vec1,*vec2,*out;
#else
MeVEC	*v_add(const MeVEC *vec1, const MeVEC *vec2, MeVEC *out)
#endif
{
	unsigned int	dim;
	/* Real	*out_ve, *vec1_ve, *vec2_ve; */

	if ( vec1==(MeVEC *)NULL || vec2==(MeVEC *)NULL )
		Meerror(E_NULL,"v_add");
	if ( vec1->dim != vec2->dim )
		Meerror(E_SIZES,"v_add");
	if ( out==(MeVEC *)NULL || out->dim != vec1->dim )
		out = v_resize(out,vec1->dim);
	dim = vec1->dim;
	__add__(vec1->ve,vec2->ve,out->ve,(int)dim);
	/************************************************************
	out_ve = out->ve;	vec1_ve = vec1->ve;	vec2_ve = vec2->ve;
	for ( i=0; i<dim; i++ )
		out->ve[i] = vec1->ve[i]+vec2->ve[i];
		(*out_ve++) = (*vec1_ve++) + (*vec2_ve++);
	************************************************************/

	return (out);
}

/* v_mltadd -- scalar/vector multiplication and addition
		-- out = v1 + scale.v2		*/
#ifndef ANSI_C
MeVEC	*v_mltadd(v1,v2,scale,out)
MeVEC	*v1,*v2,*out;
double	scale;
#else
MeVEC	*v_mltadd(const MeVEC *v1, const MeVEC *v2, double scale, MeVEC *out)
#endif
{
	/* register unsigned int	dim, i; */
	/* Real	*out_ve, *v1_ve, *v2_ve; */

	if ( v1==(MeVEC *)NULL || v2==(MeVEC *)NULL )
		Meerror(E_NULL,"v_mltadd");
	if ( v1->dim != v2->dim )
		Meerror(E_SIZES,"v_mltadd");
	if ( scale == 0.0 )
		return v_copy(v1,out);
	if ( scale == 1.0 )
		return v_add(v1,v2,out);

	if ( v2 != out )
	{
	    tracecatch(out = v_copy(v1,out),"v_mltadd");

	    /* dim = v1->dim; */
	    __mltadd__(out->ve,v2->ve,scale,(int)(v1->dim));
	}
	else
	{
	    tracecatch(out = sv_mlt(scale,v2,out),"v_mltadd");
	    out = v_add(v1,out,out);
	}
	/************************************************************
	out_ve = out->ve;	v1_ve = v1->ve;		v2_ve = v2->ve;
	for ( i=0; i < dim ; i++ )
		out->ve[i] = v1->ve[i] + scale*v2->ve[i];
		(*out_ve++) = (*v1_ve++) + scale*(*v2_ve++);
	************************************************************/

	return (out);
}

/* v_sub -- vector subtraction -- may be in-situ */
#ifndef ANSI_C
MeVEC	*v_sub(vec1,vec2,out)
MeVEC	*vec1,*vec2,*out;
#else
MeVEC	*v_sub(const MeVEC *vec1, const MeVEC *vec2, MeVEC *out)
#endif
{
	/* unsigned int	i, dim; */
	/* Real	*out_ve, *vec1_ve, *vec2_ve; */

	if ( vec1==(MeVEC *)NULL || vec2==(MeVEC *)NULL )
		Meerror(E_NULL,"v_sub");
	if ( vec1->dim != vec2->dim )
		Meerror(E_SIZES,"v_sub");
	if ( out==(MeVEC *)NULL || out->dim != vec1->dim )
		out = v_resize(out,vec1->dim);

	__sub__(vec1->ve,vec2->ve,out->ve,(int)(vec1->dim));
	/************************************************************
	dim = vec1->dim;
	out_ve = out->ve;	vec1_ve = vec1->ve;	vec2_ve = vec2->ve;
	for ( i=0; i<dim; i++ )
		out->ve[i] = vec1->ve[i]-vec2->ve[i];
		(*out_ve++) = (*vec1_ve++) - (*vec2_ve++);
	************************************************************/

	return (out);
}

/* v_map -- maps function f over components of x: out[i] = f(x[i])
	-- v_map sets out[i] = f(params,x[i]) */
#ifndef ANSI_C
MeVEC	*v_map(f,x,out)
double	(*f)();
MeVEC	*x, *out;
#else
#ifdef PROTOTYPES_IN_STRUCT
MeVEC	*v_map(double (*f)(double), const MeVEC *x, MeVEC *out)
#else
MeVEC	*v_map(double (*f)(), const MeVEC *x, MeVEC *out)
#endif
#endif
{
	Real	*x_ve, *out_ve;
	int	i, dim;

	if ( ! x || ! f )
		Meerror(E_NULL,"v_map");
	if ( ! out || out->dim != x->dim )
		out = v_resize(out,x->dim);

	dim = x->dim;	x_ve = x->ve;	out_ve = out->ve;
	for ( i = 0; i < dim; i++ )
		*out_ve++ = (*f)(*x_ve++);

	return out;
}

/* _v_map -- sets out[i] <- f(params, x[i]), i = 0, 1, .., dim-1 */
#ifndef ANSI_C
MeVEC	*_v_map(f,params,x,out)
double	(*f)();
void	*params;
MeVEC	*x, *out;
#else
#ifdef PROTOTYPES_IN_STRUCT
MeVEC	*_v_map(double (*f)(void *,double), void *params, const MeVEC *x, MeVEC *out)
#else
MeVEC	*_v_map(double (*f)(), void *params, const MeVEC *x, MeVEC *out)
#endif
#endif
{
	Real	*x_ve, *out_ve;
	int	i, dim;

	if ( ! x || ! f )
		Meerror(E_NULL,"_v_map");
	if ( ! out || out->dim != x->dim )
		out = v_resize(out,x->dim);

	dim = x->dim;	x_ve = x->ve;	out_ve = out->ve;
	for ( i = 0; i < dim; i++ )
		*out_ve++ = (*f)(params,*x_ve++);

	return out;
}

/* v_lincomb -- returns sum_i a[i].v[i], a[i] real, v[i] vectors */
#ifndef ANSI_C
MeVEC	*v_lincomb(n,v,a,out)
int	n;	/* number of a's and v's */
Real	a[];
MeVEC	*v[], *out;
#else
MeVEC	*v_lincomb(int n, const MeVEC *v[], const Real a[], MeVEC *out)
#endif
{
	int	i;

	if ( ! a || ! v )
		Meerror(E_NULL,"v_lincomb");
	if ( n <= 0 )
		return VNULL;

	for ( i = 1; i < n; i++ )
		if ( out == v[i] )
		    Meerror(E_INSITU,"v_lincomb");

	out = sv_mlt(a[0],v[0],out);
	for ( i = 1; i < n; i++ )
	{
		if ( ! v[i] )
			Meerror(E_NULL,"v_lincomb");
		if ( v[i]->dim != out->dim )
			Meerror(E_SIZES,"v_lincomb");
		out = v_mltadd(out,v[i],a[i],out);
	}

	return out;
}



#ifdef ANSI_C

/* v_linlist -- linear combinations taken from a list of arguments;
   calling:
      v_linlist(out,v1,a1,v2,a2,...,vn,an,NULL);
   where vi are vectors (MeVEC *) and ai are numbers (double)
*/
MeVEC  *v_linlist(MeVEC *out,MeVEC *v1,double a1,...)
{
   va_list ap;
   MeVEC *par;
   double a_par;

   if ( ! v1 )
     return VNULL;
   
   va_start(ap, a1);
   out = sv_mlt(a1,v1,out);
   
   while (par = va_arg(ap,MeVEC *)) {   /* NULL ends the list*/
      a_par = va_arg(ap,double);
      if (a_par == 0.0) continue;
      if ( out == par )		
	Meerror(E_INSITU,"v_linlist");
      if ( out->dim != par->dim )	
	Meerror(E_SIZES,"v_linlist");

      if (a_par == 1.0)
	out = v_add(out,par,out);
      else if (a_par == -1.0)
	out = v_sub(out,par,out);
      else
	out = v_mltadd(out,par,a_par,out); 
   } 
   
   va_end(ap);
   return out;
}
 
#elif VARARGS


/* v_linlist -- linear combinations taken from a list of arguments;
   calling:
      v_linlist(out,v1,a1,v2,a2,...,vn,an,NULL);
   where vi are vectors (MeVEC *) and ai are numbers (double)
*/
MeVEC  *v_linlist(va_alist) va_dcl
{
   va_list ap;
   MeVEC *par, *out;
   double a_par;

   va_start(ap);
   out = va_arg(ap,MeVEC *);
   par = va_arg(ap,MeVEC *);
   if ( ! par ) {
      va_end(ap);
      return VNULL;
   }
   
   a_par = va_arg(ap,double);
   out = sv_mlt(a_par,par,out);
   
   while (par = va_arg(ap,MeVEC *)) {   /* NULL ends the list*/
      a_par = va_arg(ap,double);
      if (a_par == 0.0) continue;
      if ( out == par )		
	Meerror(E_INSITU,"v_linlist");
      if ( out->dim != par->dim )	
	Meerror(E_SIZES,"v_linlist");

      if (a_par == 1.0)
	out = v_add(out,par,out);
      else if (a_par == -1.0)
	out = v_sub(out,par,out);
      else
	out = v_mltadd(out,par,a_par,out); 
   } 
   
   va_end(ap);
   return out;
}

#endif
  




/* v_star -- computes componentwise (Hadamard) product of x1 and x2
	-- result out is returned */
#ifndef ANSI_C
MeVEC	*v_star(x1, x2, out)
MeVEC	*x1, *x2, *out;
#else
MeVEC	*v_star(const MeVEC *x1, const MeVEC *x2, MeVEC *out)
#endif
{
    int		i;

    if ( ! x1 || ! x2 )
	Meerror(E_NULL,"v_star");
    if ( x1->dim != x2->dim )
	Meerror(E_SIZES,"v_star");
    out = v_resize(out,x1->dim);

    for ( i = 0; i < x1->dim; i++ )
	out->ve[i] = x1->ve[i] * x2->ve[i];

    return out;
}

/* v_slash -- computes componentwise ratio of x2 and x1
	-- out[i] = x2[i] / x1[i]
	-- if x1[i] == 0 for some i, then raise E_SING Meerror
	-- result out is returned */
#ifndef ANSI_C
MeVEC	*v_slash(x1, x2, out)
MeVEC	*x1, *x2, *out;
#else
MeVEC	*v_slash(const MeVEC *x1, const MeVEC *x2, MeVEC *out)
#endif
{
    int		i;
    Real	tmp;

    if ( ! x1 || ! x2 )
	Meerror(E_NULL,"v_slash");
    if ( x1->dim != x2->dim )
	Meerror(E_SIZES,"v_slash");
    out = v_resize(out,x1->dim);

    for ( i = 0; i < x1->dim; i++ )
    {
	tmp = x1->ve[i];
	if ( tmp == 0.0 )
	    Meerror(E_SING,"v_slash");
	out->ve[i] = x2->ve[i] / tmp;
    }

    return out;
}

/* v_Memin -- computes Meminimum component of x, which is returned
	-- also sets Memin_idx to the index of this Meminimum */
#ifndef ANSI_C
double	v_Memin(x, Memin_idx)
MeVEC	*x;
int	*Memin_idx;
#else
double	v_Memin(const MeVEC *x, int *Memin_idx)
#endif
{
    int		i, i_Memin;
    Real	Memin_val, tmp;

    if ( ! x )
	Meerror(E_NULL,"v_Memin");
    if ( x->dim <= 0 )
	Meerror(E_SIZES,"v_Memin");
    i_Memin = 0;
    Memin_val = x->ve[0];
    for ( i = 1; i < x->dim; i++ )
    {
	tmp = x->ve[i];
	if ( tmp < Memin_val )
	{
	    Memin_val = tmp;
	    i_Memin = i;
	}
    }

    if ( Memin_idx != NULL )
	*Memin_idx = i_Memin;
    return Memin_val;
}

/* v_MeMemax -- computes MeMemaximum component of x, which is returned
	-- also sets MeMemax_idx to the index of this MeMemaximum */
#ifndef ANSI_C
double	v_MeMemax(x, MeMemax_idx)
MeVEC	*x;
int	*MeMemax_idx;
#else
double	v_MeMemax(const MeVEC *x, int *MeMemax_idx)
#endif
{
    int		i, i_MeMemax;
    Real	MeMemax_val, tmp;

    if ( ! x )
	Meerror(E_NULL,"v_MeMemax");
    if ( x->dim <= 0 )
	Meerror(E_SIZES,"v_MeMemax");
    i_MeMemax = 0;
    MeMemax_val = x->ve[0];
    for ( i = 1; i < x->dim; i++ )
    {
	tmp = x->ve[i];
	if ( tmp > MeMemax_val )
	{
	    MeMemax_val = tmp;
	    i_MeMemax = i;
	}
    }

    if ( MeMemax_idx != NULL )
	*MeMemax_idx = i_MeMemax;
    return MeMemax_val;
}

#define	MAX_STACK	60


/* v_sort -- sorts vector x, and generates permutation that gives the order
	of the components; x = [1.3, 3.7, 0.5] -> [0.5, 1.3, 3.7] and
	the permutation is order = [2, 0, 1].
	-- if order is NULL on entry then it is ignored
	-- the sorted vector x is returned */
#ifndef ANSI_C
MeVEC	*v_sort(x, order)
MeVEC	*x;
PERM	*order;
#else
MeVEC	*v_sort(MeVEC *x, PERM *order)
#endif
{
    Real	*x_ve, tmp, v;
    /* int		*order_pe; */
    int		dim, i, j, l, r, tmp_i;
    int		stack[MAX_STACK], sp;

    if ( ! x )
	Meerror(E_NULL,"v_sort");
    if ( order != PNULL && order->size != x->dim )
	order = px_resize(order, x->dim);

    x_ve = x->ve;
    dim = x->dim;
    if ( order != PNULL )
	px_ident(order);

    if ( dim <= 1 )
	return x;

    /* using quicksort algorithm in Sedgewick,
       "Algorithms in C", Ch. 9, pp. 118--122 (1990) */
    sp = 0;
    l = 0;	r = dim-1;	v = x_ve[0];
    for ( ; ; )
    {
	while ( r > l )
	{
	    /* "i = partition(x_ve,l,r);" */
	    v = x_ve[r];
	    i = l-1;
	    j = r;
	    for ( ; ; )
	    {
		while ( x_ve[++i] < v )
		    ;
		--j;
		while ( x_ve[j] > v && j != 0 )
		    --j;
		if ( i >= j )	break;
		
		tmp = x_ve[i];
		x_ve[i] = x_ve[j];
		x_ve[j] = tmp;
		if ( order != PNULL )
		{
		    tmp_i = order->pe[i];
		    order->pe[i] = order->pe[j];
		    order->pe[j] = tmp_i;
		}
	    }
	    tmp = x_ve[i];
	    x_ve[i] = x_ve[r];
	    x_ve[r] = tmp;
	    if ( order != PNULL )
	    {
		tmp_i = order->pe[i];
		order->pe[i] = order->pe[r];
		order->pe[r] = tmp_i;
	    }

	    if ( i-l > r-i )
	    {   stack[sp++] = l;   stack[sp++] = i-1;   l = i+1;   }
	    else
	    {   stack[sp++] = i+1;   stack[sp++] = r;   r = i-1;   }
	}

	/* recursion eliMemination */
	if ( sp == 0 )
	    break;
	r = stack[--sp];
	l = stack[--sp];
    }

    return x;
}

/* v_sum -- returns sum of entries of a vector */
#ifndef ANSI_C
double	v_sum(x)
MeVEC	*x;
#else
double	v_sum(const MeVEC *x)
#endif
{
    int		i;
    Real	sum;

    if ( ! x )
	Meerror(E_NULL,"v_sum");

    sum = 0.0;
    for ( i = 0; i < x->dim; i++ )
	sum += x->ve[i];

    return sum;
}

/* v_conv -- computes convolution product of two vectors */
#ifndef ANSI_C
MeVEC	*v_conv(x1, x2, out)
MeVEC	*x1, *x2, *out;
#else
MeVEC	*v_conv(const MeVEC *x1, const MeVEC *x2, MeVEC *out)
#endif
{
    int		i;

    if ( ! x1 || ! x2 )
	Meerror(E_NULL,"v_conv");
    if ( x1 == out || x2 == out )
	Meerror(E_INSITU,"v_conv");
    if ( x1->dim == 0 || x2->dim == 0 )
	return out = v_resize(out,0);

    out = v_resize(out,x1->dim + x2->dim - 1);
    v_zero(out);
    for ( i = 0; i < x1->dim; i++ )
	__mltadd__(&(out->ve[i]),x2->ve,x1->ve[i],x2->dim);

    return out;
}

/* v_pconv -- computes a periodic convolution product
	-- the period is the dimension of x2 */
#ifndef ANSI_C
MeVEC	*v_pconv(x1, x2, out)
MeVEC	*x1, *x2, *out;
#else
MeVEC	*v_pconv(const MeVEC *x1, const MeVEC *x2, MeVEC *out)
#endif
{
    int		i;

    if ( ! x1 || ! x2 )
	Meerror(E_NULL,"v_pconv");
    if ( x1 == out || x2 == out )
	Meerror(E_INSITU,"v_pconv");
    out = v_resize(out,x2->dim);
    if ( x2->dim == 0 )
	return out;

    v_zero(out);
    for ( i = 0; i < x1->dim; i++ )
    {
	__mltadd__(&(out->ve[i]),x2->ve,x1->ve[i],x2->dim - i);
	if ( i > 0 )
	    __mltadd__(out->ve,&(x2->ve[x2->dim - i]),x1->ve[i],i);
    }

    return out;
}
