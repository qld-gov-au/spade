
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


/* Main include file for zmeschach library -- complex vectors and matrices */

#ifndef ZMeMATRIXH
#define ZMeMATRIXH

#include "matrix.h"


          /*  Type definitions for complex vectors and matrices  */


/* complex definition */
typedef struct  {
                Real re,im;
        } complex;

/* complex vector definition */
typedef struct  {
                unsigned int   dim, Memax_dim;
                complex  *ve;
                } ZMeVEC;

/* complex matrix definition */
typedef struct  {
                unsigned int   m, n;
                unsigned int   Memax_m, Memax_n, Memax_size;
                complex *base;          /* base is base of alloc'd mem */
                complex **me;
                } ZMeMAT;

#define ZVNULL  ((ZMeVEC *)NULL)
#define ZMNULL  ((ZMeMAT *)NULL)

#define	Z_CONJ		1
#define	Z_NOCONJ	0


#define	zm_entry(A,i,j)		zm_get_val(A,i,j)
#define	zv_entry(x,i)		zv_get_val(x,i)
#ifdef DEBUG
#define	zm_set_val(A,i,j,val)	( m_chk_idx(A,i,j) ? \
	(A)->me[(i)][(j)] = (val) : (error(E_BOUNDS,"zm_set_val"), zmake(0.0,0.0)))
#define	zm_add_val(A,i,j,val)	( m_chk_idx(A,i,j) ? \
	(A)->me[(i)][(j)] = zadd((A)->me[(i)][(j)],(val)) : \
	(error(E_BOUNDS,"zm_add_val"), zmake(0.0,0.0)))
#define	zm_sub_val(A,i,j,val)	( m_chk_idx(A,i,j) ? \
	(A)->me[(i)][(j)] = zsub((A)->me[(i)][(j)],(val)) : \
	(error(E_BOUNDS,"zm_sub_val"), zmake(0.0,0.0)))
#define	zm_get_val(A,i,j)	( m_chk_idx(A,i,j) ? \
	(A)->me[(i)][(j)] : (error(E_BOUNDS,"zm_get_val"), zmake(0.0,0.0)))
#define	zv_set_val(x,i,val)	( v_chk_idx(x,i) ? (x)->ve[(i)] = (val) : \
	(error(E_BOUNDS,"zv_set_val"), zmake(0.0,0.0)))
#define	zv_add_val(x,i,val)	( v_chk_idx(x,i) ? \
	(x)->ve[(i)] = zadd((x)->ve[(i)],(val)) : \
	(error(E_BOUNDS,"zv_set_val"), zmake(0.0,0.0)))
#define	zv_sub_val(x,i,val)	( v_chk_idx(x,i) ? \
	(x)->ve[(i)] = zsub((x)->ve[(i)],(val)) : \
	(error(E_BOUNDS,"zv_set_val"), zmake(0.0,0.0)))
#define	zv_get_val(x,i)	( v_chk_idx(x,i) ? (x)->ve[(i)] : \
	(error(E_BOUNDS,"zv_get_val"), zmake(0.0,0.0)))
#else /* no DEBUG */
#define	zm_set_val(A,i,j,val)	((A)->me[(i)][(j)] = (val))
#define	zm_add_val(A,i,j,val)	((A)->me[(i)][(j)] = zadd((A)->me[(i)][(j)],(val)))
#define	zm_sub_val(A,i,j,val)	((A)->me[(i)][(j)] = zsub((A)->me[(i)][(j)],(val)))
#define	zm_get_val(A,i,j)	((A)->me[(i)][(j)])
#define	zv_set_val(x,i,val)	((x)->ve[(i)] = (val))
#define	zv_add_val(x,i,val)	((x)->ve[(i)] = zadd((x)->ve[(i)],(val)))
#define	zv_sub_val(x,i,val)	((x)->ve[(i)] = zsub((x)->ve[(i)],(val)))
#define	zv_get_val(x,i)		((x)->ve[(i)])
#endif /* DEBUG */

/* memory functions */

#ifdef ANSI_C
int zv_get_vars(int dim,...);
int zm_get_vars(int m,int n,...);
int zv_resize_vars(int new_dim,...);
int zm_resize_vars(int m,int n,...);
int zv_free_vars(ZMeVEC **,...);
int zm_free_vars(ZMeMAT **,...);

#elif VARARGS
int zv_get_vars();
int zm_get_vars();
int zv_resize_vars();
int zm_resize_vars();
int zv_free_vars();
int zm_free_vars();

#endif




#ifdef ANSI_C
extern ZMeMAT	*_zm_copy(const ZMeMAT *in,ZMeMAT *out, int i0, int j0);
extern ZMeMAT	* zm_move(const ZMeMAT *, int, int, int, int, ZMeMAT *, int, int);
extern ZMeMAT	*zvm_move(const ZMeVEC *, int, ZMeMAT *, int, int, int, int);
extern ZMeVEC	*_zv_copy(const ZMeVEC *in,ZMeVEC *out,int i0);
extern ZMeVEC	* zv_move(const ZMeVEC *, int, int, ZMeVEC *, int);
extern ZMeVEC	*zmv_move(const ZMeMAT *, int, int, int, int, ZMeVEC *, int);
extern complex	z_finput(FILE *fp);
extern ZMeMAT	*zm_finput(FILE *fp,ZMeMAT *a);
extern ZMeVEC     *zv_finput(FILE *fp,ZMeVEC *x);
extern ZMeMAT	*zm_add(ZMeMAT *mat1,ZMeMAT *mat2,ZMeMAT *out);
extern ZMeMAT	*zm_sub(ZMeMAT *mat1,ZMeMAT *mat2,ZMeMAT *out);
extern ZMeMAT	*zm_mlt(ZMeMAT *A,ZMeMAT *B,ZMeMAT *OUT);
extern ZMeMAT	*zmma_mlt(ZMeMAT *A,ZMeMAT *B,ZMeMAT *OUT);
extern ZMeMAT	*zmam_mlt(ZMeMAT *A,ZMeMAT *B,ZMeMAT *OUT);
extern ZMeVEC	*zmv_mlt(ZMeMAT *A,ZMeVEC *b,ZMeVEC *out);
extern ZMeMAT	*zsm_mlt(complex scalar,ZMeMAT *matrix,ZMeMAT *out);
extern ZMeVEC	*zvm_mlt(ZMeMAT *A,ZMeVEC *b,ZMeVEC *out);
extern ZMeMAT	*zm_adjoint(ZMeMAT *in,ZMeMAT *out);
extern ZMeMAT	*zswap_rows(ZMeMAT *A,int i,int j,int lo,int hi);
extern ZMeMAT	*zswap_cols(ZMeMAT *A,int i,int j,int lo,int hi);
extern ZMeMAT	*mz_mltadd(ZMeMAT *A1,ZMeMAT *A2,complex s,ZMeMAT *out);
extern ZMeVEC	*zmv_mltadd(ZMeVEC *v1,ZMeVEC *v2,ZMeMAT *A,complex alpha,ZMeVEC *out);
extern ZMeVEC	*zvm_mltadd(ZMeVEC *v1,ZMeVEC *v2,ZMeMAT *A,complex alpha,ZMeVEC *out);
extern ZMeVEC	*zv_zero(ZMeVEC *x);
extern ZMeMAT	*zm_zero(ZMeMAT *A);
extern ZMeMAT	*zm_get(int m,int n);
extern ZMeVEC	*zv_get(int dim);
extern ZMeMAT	*zm_resize(ZMeMAT *A,int new_m,int new_n);
extern complex	_zin_prod(const ZMeVEC *x, const ZMeVEC *y,unsigned int i0,unsigned int flag);
extern ZMeVEC	*zv_resize(ZMeVEC *x,int new_dim);
extern ZMeVEC	*zv_mlt(complex scalar,const ZMeVEC *vector,ZMeVEC *out);
extern ZMeVEC	*zv_add(const ZMeVEC *vec1,const ZMeVEC *vec2,ZMeVEC *out);
extern ZMeVEC	*zv_mltadd(const ZMeVEC *v1,const ZMeVEC *v2,complex scale,ZMeVEC *out);
extern ZMeVEC	*zv_sub(const ZMeVEC *vec1,const ZMeVEC *vec2,ZMeVEC *out);
#ifdef PROTOTYPES_IN_STRUCT
extern ZMeVEC	*zv_map(complex (*f)(),const ZMeVEC *x,ZMeVEC *out);
extern ZMeVEC	*_zv_map(complex (*f)(),void *params,const ZMeVEC *x,ZMeVEC *out);
#else
extern ZMeVEC	*zv_map(complex (*f)(complex),const ZMeVEC *x,ZMeVEC *out);
extern ZMeVEC	*_zv_map(complex (*f)(void *,complex),void *params,const ZMeVEC *x,ZMeVEC *out);
#endif
extern ZMeVEC	*zv_lincomb(int n,const ZMeVEC *v[],const complex a[],ZMeVEC *out);
extern ZMeVEC	*zv_linlist(ZMeVEC *out,ZMeVEC *v1,complex a1,...);
extern ZMeVEC	*zv_star(const ZMeVEC *x1, const ZMeVEC *x2, ZMeVEC *out);
extern ZMeVEC	*zv_slash(const ZMeVEC *x1, const ZMeVEC *x2, ZMeVEC *out);
extern complex	zv_sum(const ZMeVEC *x);
extern int	zm_free(ZMeMAT *mat);
extern int	zv_free(ZMeVEC *vec);

extern ZMeVEC	*zv_rand(ZMeVEC *x);
extern ZMeMAT	*zm_rand(ZMeMAT *A);

extern ZMeVEC	*zget_row(ZMeMAT *A, int i, ZMeVEC *out);
extern ZMeVEC	*zget_col(ZMeMAT *A, int j, ZMeVEC *out);
extern ZMeMAT	*zset_row(ZMeMAT *A, int i, ZMeVEC *in);
extern ZMeMAT	*zset_col(ZMeMAT *A, int j, ZMeVEC *in);

extern ZMeVEC	*px_zvec(PERM *pi, ZMeVEC *in, ZMeVEC *out);
extern ZMeVEC	*pxinv_zvec(PERM *pi, ZMeVEC *in, ZMeVEC *out);

extern void	__zconj__(complex zp[], int len);
extern complex	__zip__(const complex zp1[], const complex zp2[],
			int len,int flag);
extern void	__zmltadd__(complex zp1[], const complex zp2[],
			    complex s,int len,int flag);
extern void	__zmlt__(const complex zp[],complex s,complex out[],int len);
extern void	__zadd__(const complex zp1[],const complex zp2[],
			 complex out[],int len);
extern void	__zsub__(const complex zp1[],const complex zp2[],
			 complex out[],int len);
extern void	__zzero__(complex zp[],int len);
extern void	z_foutput(FILE *fp,complex z);
extern void     zm_foutput(FILE *fp,ZMeMAT *a);
extern void     zv_foutput(FILE *fp,ZMeVEC *x);
extern void     zm_dump(FILE *fp,ZMeMAT *a);
extern void     zv_dump(FILE *fp,ZMeVEC *x);

extern double	_zv_norm1(ZMeVEC *x, MeVEC *scale);
extern double	_zv_norm2(ZMeVEC *x, MeVEC *scale);
extern double	_zv_norm_inf(ZMeVEC *x, MeVEC *scale);
extern double	zm_norm1(ZMeMAT *A);
extern double	zm_norm_inf(ZMeMAT *A);
extern double	zm_norm_frob(ZMeMAT *A);

complex	zmake(double real, double imag);
double	zabs(complex z);
complex zadd(complex z1,complex z2);
complex zsub(complex z1,complex z2);
complex	zmlt(complex z1,complex z2);
complex	zinv(complex z);
complex	zdiv(complex z1,complex z2);
complex	zsqrt(complex z);
complex	zexp(complex z);
complex	zlog(complex z);
complex	zconj(complex z);
complex	zneg(complex z);
#else
extern ZMeMAT	*_zm_copy();
extern ZMeVEC	*_zv_copy();
extern ZMeMAT	*zm_finput();
extern ZMeVEC     *zv_finput();
extern ZMeMAT	*zm_add();
extern ZMeMAT	*zm_sub();
extern ZMeMAT	*zm_mlt();
extern ZMeMAT	*zmma_mlt();
extern ZMeMAT	*zmam_mlt();
extern ZMeVEC	*zmv_mlt();
extern ZMeMAT	*zsm_mlt();
extern ZMeVEC	*zvm_mlt();
extern ZMeMAT	*zm_adjoint();
extern ZMeMAT	*zswap_rows();
extern ZMeMAT	*zswap_cols();
extern ZMeMAT	*mz_mltadd();
extern ZMeVEC	*zmv_mltadd();
extern ZMeVEC	*zvm_mltadd();
extern ZMeVEC	*zv_zero();
extern ZMeMAT	*zm_zero();
extern ZMeMAT	*zm_get();
extern ZMeVEC	*zv_get();
extern ZMeMAT	*zm_resize();
extern ZMeVEC	*zv_resize();
extern complex	_zin_prod();
extern ZMeVEC	*zv_mlt();
extern ZMeVEC	*zv_add();
extern ZMeVEC	*zv_mltadd();
extern ZMeVEC	*zv_sub();
extern ZMeVEC	*zv_map();
extern ZMeVEC	*_zv_map();
extern ZMeVEC	*zv_lincomb();
extern ZMeVEC	*zv_linlist();
extern ZMeVEC	*zv_star();
extern ZMeVEC	*zv_slash();

extern ZMeVEC	*px_zvec();
extern ZMeVEC	*pxinv_zvec();

extern ZMeVEC	*zv_rand();
extern ZMeMAT	*zm_rand();

extern ZMeVEC	*zget_row();
extern ZMeVEC	*zget_col();
extern ZMeMAT	*zset_row();
extern ZMeMAT	*zset_col();

extern int	zm_free();
extern int	zv_free();
extern void	__zconj__();
extern complex	__zip__();
extern void	__zmltadd__();
extern void	__zmlt__();
extern void	__zadd__();
extern void	__zsub__();
extern void	__zzero__();
extern void    zm_foutput();
extern void    zv_foutput();
extern void    zm_dump();
extern void    zv_dump();

extern double	_zv_norm1();
extern double	_zv_norm2();
extern double	_zv_norm_inf();
extern double	zm_norm1();
extern double	zm_norm_inf();
extern double	zm_norm_frob();

complex	zmake();
double	zabs();
complex zadd();
complex zsub();
complex	zmlt();
complex	zinv();
complex	zdiv();
complex	zsqrt();
complex	zexp();
complex	zlog();
complex	zconj();
complex	zneg();
#endif

#define	zv_copy(x,y)	_zv_copy(x,y,0)
#define	zm_copy(A,B)	_zm_copy(A,B,0,0)

#define	z_input()	z_finput(stdin)
#define	zv_input(x)	zv_finput(stdin,x)
#define	zm_input(A)	zm_finput(stdin,A)
#define	z_output(z)	z_foutput(stdout,z)
#define	zv_output(x)	zv_foutput(stdout,x)
#define	zm_output(A)	zm_foutput(stdout,A)

#define	ZV_FREE(x)	( zv_free(x), (x) = ZVNULL )
#define	ZM_FREE(A)	( zm_free(A), (A) = ZMNULL )

#define	zin_prod(x,y)	_zin_prod(x,y,0,Z_CONJ)

#define	zv_norm1(x)	_zv_norm1(x,VNULL)
#define	zv_norm2(x)	_zv_norm2(x,VNULL)
#define	zv_norm_inf(x)	_zv_norm_inf(x,VNULL)


#endif
