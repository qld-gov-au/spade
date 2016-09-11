
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
		Type definitions for general purpose maths package
*/

#ifndef	MeMATRIXH

/* RCS id: $Id: matrix.h,v 1.18 1994/04/16 00:33:37 des Exp $ */

#define	MeMATRIXH	

#include	"machine.h"
#include        "err.h"
#include 	"meminfo.h"

/* unsigned integer type */
/************************************************************
#ifndef U_INT_DEF
typedef	unsigned int	u_int;
#define U_INT_DEF
#endif
************************************************************/

/* vector definition */
typedef	struct	{
		unsigned int	dim, MeMemax_dim;
		Real	*ve;
		} MeVEC;

/* matrix definition */
typedef	struct	{
		unsigned int	m, n;
		unsigned int	MeMemax_m, MeMemax_n, MeMemax_size;
		Real	**me,*base;	/* base is base of alloc'd mem */
		} MeMAT;

/* band matrix definition */
typedef struct {
               MeMAT   *mat;       /* matrix */
               int   lb,ub;    /* lower and upper bandwidth */
               } BAND;


/* permutation definition */
typedef	struct	{
		unsigned int	size, MeMemax_size, *pe;
		} PERM;

/* integer vector definition */
typedef struct	{
		unsigned int	dim, MeMemax_dim;
		int	*ive;
	        } IMeVEC;


#ifndef MALLOCDECL
#ifndef ANSI_C
extern	char	*malloc(), *calloc(), *realloc();
#else
extern	void	*malloc(size_t),
		*calloc(size_t,size_t),
		*realloc(void *,size_t);
#endif
#endif /* MALLOCDECL */

/* For creating MEX files (for use with Matlab) using Meschach
   See also: mexmesch.h */
#ifdef MEX
#include	"mex.h"
#define	malloc(len)		mxMalloc(len)
#define	calloc(n,len)		mxCalloc(n,len)
#define	realloc(ptr,len)	mxRealloc(ptr,len)
#define	free(ptr)		mxFree(ptr)
#define	printf			mexPrintf
#ifndef THREADSAFE	/* for use as a shared library */
#define	THREADSAFE 1
#endif
#endif /* MEX */

#ifdef THREADSAFE
#define	STATIC
#else
#define	STATIC	static
#endif /* THREADSAFE */

#ifndef ANSI_C
extern void m_version();
#else
void	m_version( void );
#endif

#ifndef ANSI_C
/* allocate one object of given type */
#define	NEW(type)	((type *)calloc((size_t)1,sizeof(type)))

/* allocate num objects of given type */
#define	NEW_A(num,type)	((type *)calloc((size_t)(num),sizeof(type)))

 /* re-allocate arry to have num objects of the given type */
#define	RENEW(var,num,type) \
    ((var)=(type *)((var) ? \
		    realloc((char *)(var),(size_t)(num)*sizeof(type)) : \
		    calloc((size_t)(num),sizeof(type))))

#define	MEMCOPY(from,to,n_items,type) \
    MEM_COPY((char *)(from),(char *)(to),(size_t)(n_items)*sizeof(type))

#else
/* allocate one object of given type */
#define	NEW(type)	((type *)calloc((size_t)1,(size_t)sizeof(type)))

/* allocate num objects of given type */
#define	NEW_A(num,type)	((type *)calloc((size_t)(num),(size_t)sizeof(type)))

 /* re-allocate arry to have num objects of the given type */
#define	RENEW(var,num,type) \
    ((var)=(type *)((var) ? \
		    realloc((char *)(var),(size_t)((num)*sizeof(type))) : \
		    calloc((size_t)(num),(size_t)sizeof(type))))

#define	MEMCOPY(from,to,n_items,type) \
 MEM_COPY((char *)(from),(char *)(to),(unsigned)(n_items)*sizeof(type))

#endif /* ANSI_C */

/* type independent Memin and MeMemax operations */
#ifndef MeMemax
#define	MeMemax(a,b)	((a) > (b) ? (a) : (b))
#endif /* MeMemax */
#ifndef Memin
#define	Memin(a,b)	((a) > (b) ? (b) : (a))
#endif /* Memin */


#undef TRUE
#define	TRUE	1
#undef FALSE
#define	FALSE	0


/* for me_input routines */
#define MAXLINE 81


/* Dynamic memory allocation */

/* Should use M_FREE/V_FREE/PX_FREE in programs instead of m/v/px_free()
   as this is considerably safer -- also provides a simple type check ! */

#ifndef ANSI_C

extern	MeVEC *v_get(), *v_resize();
extern	MeMAT *m_get(), *m_resize();
extern	PERM *px_get(), *px_resize();
extern	IMeVEC *iv_get(), *iv_resize();
extern	int m_free(),v_free();
extern  int px_free();
extern  int iv_free();
extern  BAND *bd_get(), *bd_resize();
extern  int bd_free();

#else

/* get/resize vector to given dimension */
extern	MeVEC *v_get(int), *v_resize(MeVEC *,int);
/* get/resize matrix to be m x n */
extern	MeMAT *m_get(int,int), *m_resize(MeMAT *,int,int);
/* get/resize permutation to have the given size */
extern	PERM *px_get(int), *px_resize(PERM *,int);
/* get/resize an integer vector to given dimension */
extern	IMeVEC *iv_get(int), *iv_resize(IMeVEC *,int);
/* get/resize a band matrix to given dimension */
extern  BAND *bd_get(int,int,int), *bd_resize(BAND *,int,int,int);

/* free (de-allocate) (band) matrices, vectors, permutations and 
   integer vectors */
extern  int iv_free(IMeVEC *);
extern	int m_free(MeMAT *),v_free(MeVEC *),px_free(PERM *);
extern   int bd_free(BAND *);

#endif /* ANSI_C */


/* MACROS */

/* macros that also check types and sets pointers to NULL */
#define	M_FREE(mat)	( m_free(mat),	(mat)=(MeMAT *)NULL )
#define V_FREE(vec)	( v_free(vec),	(vec)=(MeVEC *)NULL )
#define	PX_FREE(px)	( px_free(px),	(px)=(PERM *)NULL )
#define	IV_FREE(iv)	( iv_free(iv),	(iv)=(IMeVEC *)NULL )

#define MAXDIM  	10000001


/* Entry level access to data structures */
/* routines to check indexes */
#define	m_chk_idx(A,i,j)	((i)>=0 && (i)<(A)->m && (j)>=0 && (j)<=(A)->n)
#define	v_chk_idx(x,i)		((i)>=0 && (i)<(x)->dim)
#define	bd_chk_idx(A,i,j)	((i)>=MeMemax(0,(j)-(A)->ub) && \
		(j)>=MeMemax(0,(i)-(A)->lb) && (i)<(A)->mat->n && (j)<(A)->mat->n)

#define	m_entry(A,i,j)		m_get_val(A,i,j)
#define	v_entry(x,i)		v_get_val(x,i)
#define	bd_entry(A,i,j)		bd_get_val(A,i,j)
#ifdef DEBUG
#define	m_set_val(A,i,j,val)	( m_chk_idx(A,i,j) ? \
	(A)->me[(i)][(j)] = (val) : (Meerror(E_BOUNDS,"m_set_val"), 0.0))
#define	m_add_val(A,i,j,val)	( m_chk_idx(A,i,j) ? \
	(A)->me[(i)][(j)] += (val) : (Meerror(E_BOUNDS,"m_add_val"), 0.0))
#define	m_sub_val(A,i,j,val)	( m_chk_idx(A,i,j) ? \
	(A)->me[(i)][(j)] -= (val) : (Meerror(E_BOUNDS,"m_sub_val"), 0.0))
#define	m_get_val(A,i,j)	( m_chk_idx(A,i,j) ? \
	(A)->me[(i)][(j)] : (Meerror(E_BOUNDS,"m_get_val"), 0.0))
#define	v_set_val(x,i,val)	( v_chk_idx(x,i) ? (x)->ve[(i)] = (val) : \
	(Meerror(E_BOUNDS,"v_set_val"), 0.0))
#define	v_add_val(x,i,val)	( v_chk_idx(x,i) ? (x)->ve[(i)] += (val) : \
	(Meerror(E_BOUNDS,"v_set_val"), 0.0))
#define	v_sub_val(x,i,val)	( v_chk_idx(x,i) ? (x)->ve[(i)] -= (val) : \
	(Meerror(E_BOUNDS,"v_set_val"), 0.0))
#define	v_get_val(x,i)	( v_chk_idx(x,i) ? (x)->ve[(i)] : \
	(Meerror(E_BOUNDS,"v_get_val"), 0.0))
#define	bd_set_val(A,i,j,val)	( bd_chk_idx(A,i,j) ? \
	(A)->mat->me[(A)->lb+(j)-(i)][(j)] = (val) : \
	(Meerror(E_BOUNDS,"bd_set_val"), 0.0))
#define	bd_add_val(A,i,j,val)	( bd_chk_idx(A,i,j) ? \
	(A)->mat->me[(A)->lb+(j)-(i)][(j)] += (val) : \
	(Meerror(E_BOUNDS,"bd_set_val"), 0.0))
#define	bd_get_val(A,i,j)	( bd_chk_idx(A,i,j) ? \
	(A)->mat->me[(A)->lb+(j)-(i)][(j)] : \
	(Meerror(E_BOUNDS,"bd_get_val"), 0.0))
#else /* no DEBUG */
#define	m_set_val(A,i,j,val)	((A)->me[(i)][(j)] = (val))
#define	m_add_val(A,i,j,val)	((A)->me[(i)][(j)] += (val))
#define	m_sub_val(A,i,j,val)	((A)->me[(i)][(j)] -= (val))
#define	m_get_val(A,i,j)	((A)->me[(i)][(j)])
#define	v_set_val(x,i,val)	((x)->ve[(i)] = (val))
#define	v_add_val(x,i,val)	((x)->ve[(i)] += (val))
#define	v_sub_val(x,i,val)	((x)->ve[(i)] -= (val))
#define	v_get_val(x,i)		((x)->ve[(i)])
#define	bd_set_val(A,i,j,val)	((A)->mat->me[(A)->lb+(j)-(i)][(j)] = (val))
#define	bd_add_val(A,i,j,val)	((A)->mat->me[(A)->lb+(j)-(i)][(j)] += (val))
#define	bd_get_val(A,i,j)	((A)->mat->me[(A)->lb+(j)-(i)][(j)])
#endif /* DEBUG */


/* I/O routines */
#ifndef ANSI_C

extern	void v_foutput(),m_foutput(),px_foutput();
extern  void iv_foutput();
extern	MeVEC *v_fme_input();
extern	MeMAT *m_fme_input();
extern	PERM *px_fme_input();
extern	IMeVEC *iv_fme_input();
extern	int fy_or_n(), fin_int(), yn_dflt(), skipjunk();
extern	double fin_double();

#else

/* print x on file fp */
void v_foutput(FILE *fp,const MeVEC *x),
       /* print A on file fp */
	m_foutput(FILE *fp,const MeMAT *A),
       /* print px on file fp */
	px_foutput(FILE *fp,const PERM *px);
/* print ix on file fp */
void iv_foutput(FILE *fp,const IMeVEC *ix);

/* Note: if out is NULL, then returned object is newly allocated;
        Also: if out is not NULL, then that size is assumed */

/* read in vector from fp */
MeVEC *v_fme_input(FILE *fp,MeVEC *out);
/* read in matrix from fp */
MeMAT *m_fme_input(FILE *fp,MeMAT *out);
/* read in permutation from fp */
PERM *px_fme_input(FILE *fp,PERM *out);
/* read in int vector from fp */
IMeVEC *iv_fme_input(FILE *fp,IMeVEC *out);

/* fy_or_n -- yes-or-no to question in string s
        -- question written to stderr, me_input from fp 
        -- if fp is NOT a tty then return y_n_dflt */
int fy_or_n(FILE *fp, const char *s);

/* yn_dflt -- sets the value of y_n_dflt to val */
int yn_dflt(int val);

/* fin_int -- return integer read from file/stream fp
        -- prompt s on stderr if fp is a tty
        -- check that x lies between low and high: re-prompt if
                fp is a tty, Meerror exit otherwise
        -- ignore check if low > high           */
int fin_int(FILE *fp,const char *s,int low,int high);

/* fin_double -- return double read from file/stream fp
        -- prompt s on stderr if fp is a tty
        -- check that x lies between low and high: re-prompt if
                fp is a tty, Meerror exit otherwise
        -- ignore check if low > high           */
double fin_double(FILE *fp,const char *s,double low,double high);

/* it skips white spaces and strings of the form #....\n
   Here .... is a comment string */
int skipjunk(FILE *fp);

#endif /* ANSI_C */


/* MACROS */

/* macros to use stdout and stdin instead of explicit fp */
#define	v_output(vec)	v_foutput(stdout,vec)
#define	v_me_input(vec)	v_fme_input(stdin,vec)
#define	m_output(mat)	m_foutput(stdout,mat)
#define	m_me_input(mat)	m_fme_input(stdin,mat)
#define	px_output(px)	px_foutput(stdout,px)
#define	px_me_input(px)	px_fme_input(stdin,px)
#define	iv_output(iv)	iv_foutput(stdout,iv)
#define	iv_me_input(iv)	iv_fme_input(stdin,iv)

/* general purpose me_input routine; skips comments # ... \n */
#define	fme_input(fp,prompt,fmt,var) \
	( ( isatty(fileno(fp)) ? fprintf(stderr,prompt) : skipjunk(fp) ), \
							fscanf(fp,fmt,var) )
#define	me_input(prompt,fmt,var)	fme_input(stdin,prompt,fmt,var)
#define	fprompter(fp,prompt) \
	( isatty(fileno(fp)) ? fprintf(stderr,prompt) : skipjunk(fp) )
#define	prompter(prompt)	fprompter(stdin,prompt)
#define	y_or_n(s)	fy_or_n(stdin,s)
#define	in_int(s,lo,hi)	fin_int(stdin,s,lo,hi)
#define	in_double(s,lo,hi)	fin_double(stdin,s,lo,hi)


/* special purpose access routines */

/* Copying routines */
#ifndef ANSI_C
extern	MeMAT	*_m_copy(), *m_move(), *vm_move();
extern	MeVEC	*_v_copy(), *v_move(), *mv_move();
extern	PERM	*px_copy();
extern	IMeVEC	*iv_copy(), *iv_move();
extern  BAND    *bd_copy();

#else

/* copy in to out starting at out[i0][j0] */
extern	MeMAT	*_m_copy(const MeMAT *in,MeMAT *out,unsigned int i0,unsigned int j0),
		* m_move(const MeMAT *in, int, int, int, int, MeMAT *out, int, int),
		*vm_move(const MeVEC *in, int, MeMAT *out, int, int, int, int);
/* copy in to out starting at out[i0] */
extern	MeVEC	*_v_copy(const MeVEC *in,MeVEC *out,unsigned int i0),
		* v_move(const MeVEC *in, int, int, MeVEC *out, int),
		*mv_move(const MeMAT *in, int, int, int, int, MeVEC *out, int);
extern	PERM	*px_copy(const PERM *in,PERM *out);
extern	IMeVEC	*iv_copy(const IMeVEC *in,IMeVEC *out),
		*iv_move(const IMeVEC *in, int, int, IMeVEC *out, int);
extern  BAND    *bd_copy(const BAND *in,BAND *out);

#endif /* ANSI_C */


/* MACROS */
#define	m_copy(in,out)	_m_copy(in,out,0,0)
#define	v_copy(in,out)	_v_copy(in,out,0)


/* Initialisation routines -- to be zero, ones, random or identity */
#ifndef ANSI_C
extern	MeVEC     *v_zero(), *v_rand(), *v_ones();
extern	MeMAT     *m_zero(), *m_ident(), *m_rand(), *m_ones();
extern	PERM    *px_ident();
extern  IMeVEC    *iv_zero();
#else
extern	MeVEC     *v_zero(MeVEC *), *v_rand(MeVEC *), *v_ones(MeVEC *);
extern	MeMAT     *m_zero(MeMAT *), *m_ident(MeMAT *), *m_rand(MeMAT *),
						*m_ones(MeMAT *);
extern	PERM    *px_ident(PERM *);
extern  IMeVEC    *iv_zero(IMeVEC *);
#endif /* ANSI_C */

/* Basic vector operations */
#ifndef ANSI_C
extern	MeVEC *sv_mlt(), *mv_mlt(), *vm_mlt(), *v_add(), *v_sub(),
		*px_vec(), *pxinv_vec(), *v_mltadd(), *v_map(), *_v_map(),
		*v_lincomb(), *v_linlist();
extern	double	v_Memin(), v_MeMemax(), v_sum();
extern	MeVEC	*v_star(), *v_slash(), *v_sort();
extern	double _in_prod(), __ip__();
extern	void	__mltadd__(), __add__(), __sub__(), 
                __smlt__(), __zero__();
#else

extern	MeVEC	*sv_mlt(double s,const MeVEC *x,MeVEC *out),	/* out <- s.x */
		*mv_mlt(const MeMAT *A,const MeVEC *s,MeVEC *out),	/* out <- A.x */
		*vm_mlt(const MeMAT *A,const MeVEC *x,MeVEC *out),	/* out^T <- x^T.A */
		*v_add(const MeVEC *x,const MeVEC *y,MeVEC *out), 	/* out <- x + y */
                *v_sub(const MeVEC *x,const MeVEC *y,MeVEC *out),	/* out <- x - y */
		*px_vec(PERM *px,const MeVEC *x,MeVEC *out),	/* out <- P.x */
		*pxinv_vec(PERM *px,const MeVEC *x,MeVEC *out),	/* out <- P^{-1}.x */
		*v_mltadd(const MeVEC *x,const MeVEC *y,double s,MeVEC *out),   /* out <- x + s.y */
#ifdef PROTOTYPES_IN_STRUCT
		*v_map(double (*f)(double),const MeVEC *x,MeVEC *y),  
                                                 /* out[i] <- f(x[i]) */
		*_v_map(double (*f)(void *,double),void *p,const MeVEC *x,MeVEC *y),
#else
		*v_map(double (*f)(),const MeVEC *,MeVEC *), /* out[i] <- f(x[i]) */
		*_v_map(double (*f)(),void *,const MeVEC *,MeVEC *),
#endif /* PROTOTYPES_IN_STRUCT */
		*v_lincomb(int,const MeVEC **,const Real *,MeVEC *),   
                                                 /* out <- sum_i s[i].x[i] */
                *v_linlist(MeVEC *out,MeVEC *v1,double a1,...);
                                              /* out <- s1.x1 + s2.x2 + ... */

/* returns Memin_j x[j] (== x[i]) */
extern	double	v_Memin(const MeVEC *, int *), 
     /* returns MeMemax_j x[j] (== x[i]) */		
        v_MeMemax(const MeVEC *, int *), 
        /* returns sum_i x[i] */
        v_sum(const MeVEC *);

/* Hadamard product: out[i] <- x[i].y[i] */
extern	MeVEC	*v_star(const MeVEC *, const MeVEC *, MeVEC *),
                 /* out[i] <- x[i] / y[i] */
		*v_slash(const MeVEC *, const MeVEC *, MeVEC *),
               /* sorts x, and sets order so that sorted x[i] = x[order[i]] */ 
		*v_sort(MeVEC *, PERM *);

/* returns inner product starting at component i0 */
extern	double	_in_prod(const MeVEC *x, const MeVEC *y,unsigned int i0),
                /* returns sum_{i=0}^{len-1} x[i].y[i] */
                __ip__(const Real *,const Real *,int);

/* see v_mltadd(), v_add(), v_sub() and v_zero() */
extern	void	__mltadd__(Real *,const Real *,double,int),
		__add__(const Real *,const Real *,Real *,int),
		__sub__(const Real *,const Real *,Real *,int),
                __smlt__(const Real *,double,Real *,int),
		__zero__(Real *,int);

#endif /* ANSI_C */


/* MACRO */
/* usual way of computing the inner product */
#define	in_prod(a,b)	_in_prod(a,b,0)

/* Norms */
/* scaled vector norms -- scale == NULL implies unscaled */
#ifndef ANSI_C

extern	double	_v_norm1(), _v_norm2(), _v_norm_inf(),
		m_norm1(), m_norm_inf(), m_norm_frob();

#else
               /* returns sum_i |x[i]/scale[i]| */
extern	double	_v_norm1(const MeVEC *x,const MeVEC *scale),   
               /* returns (scaled) Euclidean norm */
                _v_norm2(const MeVEC *x,const MeVEC *scale),
               /* returns MeMemax_i |x[i]/scale[i]| */
		_v_norm_inf(const MeVEC *x,const MeVEC *scale);

/* unscaled matrix norms */
extern double m_norm1(const MeMAT *A), 
	m_norm_inf(const MeMAT *A), 
	m_norm_frob(const MeMAT *A);

#endif /* ANSI_C */


/* MACROS */
/* unscaled vector norms */
#define	v_norm1(x)	_v_norm1(x,VNULL)
#define	v_norm2(x)	_v_norm2(x,VNULL)
#define	v_norm_inf(x)	_v_norm_inf(x,VNULL)

/* Basic matrix operations */
#ifndef ANSI_C

extern	MeMAT *sm_mlt(), *m_mlt(), *mmtr_mlt(), *mtrm_mlt(), *m_add(), *m_sub(),
		*sub_mat(), *m_transp(), *ms_mltadd();

extern  BAND *bd_transp(), *sbd_mlt(), *bds_mltadd(), *bd_zero();
extern	MeMAT *px_rows(), *px_cols(), *swap_rows(), *swap_cols(),
             *_set_row(), *_set_col();
extern	MeVEC *get_row(), *get_col(), *sub_vec(),
		*mv_mltadd(), *vm_mltadd(), *bdv_mltadd();

#else

extern	MeMAT	*sm_mlt(double s, const MeMAT *A,MeMAT *out), 	/* out <- s.A */
		*m_mlt(const MeMAT *A,const MeMAT *B,MeMAT *out),	/* out <- A.B */
		*mmtr_mlt(const MeMAT *A,const MeMAT *B,MeMAT *out),	/* out <- A.B^T */
		*mtrm_mlt(const MeMAT *A,const MeMAT *B,MeMAT *out),	/* out <- A^T.B */
		*m_add(const MeMAT *A,const MeMAT *B,MeMAT *out),	/* out <- A + B */
		*m_sub(const MeMAT *A,const MeMAT *B,MeMAT *out),	/* out <- A - B */
		*sub_mat(const MeMAT *A,unsigned int,unsigned int,unsigned int,
			 unsigned int,MeMAT *out),
		*m_transp(const MeMAT *A,MeMAT *out),		/* out <- A^T */
                /* out <- A + s.B */ 
		*ms_mltadd(const MeMAT *A,const MeMAT *B,double s,MeMAT *out);   


extern  BAND    *bd_transp(const BAND *in, BAND *out),	/* out <- A^T */
  *sbd_mlt(Real s, const BAND *A, BAND *OUT),		/* OUT <- s.A */
  *bds_mltadd(const BAND *A, const BAND *B,double alpha, BAND *OUT),
  /* OUT <- A+alpha.B */
  *bd_zero(BAND *A);					/* A <- 0 */

extern	MeMAT	*px_rows(const PERM *px,const MeMAT *A,MeMAT *out),	/* out <- P.A */
		*px_cols(const PERM *px,const MeMAT *A,MeMAT *out),	/* out <- A.P^T */
		*swap_rows(MeMAT *,int,int,int,int),
		*swap_cols(MeMAT *,int,int,int,int),
                 /* A[i][j] <- out[j], j >= j0 */
		*_set_col(MeMAT *A,unsigned int i,const MeVEC *col,unsigned int j0),
                 /* A[i][j] <- out[i], i >= i0 */
		*_set_row(MeMAT *A,unsigned int j,const MeVEC *row,unsigned int i0);

extern	MeVEC	*get_row(const MeMAT *,unsigned int,MeVEC *),
		*get_col(const MeMAT *,unsigned int,MeVEC *),
		*sub_vec(const MeVEC *,int,int,MeVEC *),
                   /* mv_mltadd: out <- x + s.A.y */
		*mv_mltadd(const MeVEC *x,const MeVEC *y,const MeMAT *A,
			   double s,MeVEC *out),
                  /* vm_mltadd: out^T <- x^T + s.y^T.A */
		*vm_mltadd(const MeVEC *x,const MeVEC *y,const MeMAT *A,
			   double s,MeVEC *out),
                  /* bdv_mltadd: out <- x + s.A.y */
                *bdv_mltadd(const MeVEC *x,const MeVEC *y,const BAND *A,
			    double s,MeVEC *out);
#endif /* ANSI_C */


/* MACROS */
/* row i of A <- vec */
#define	set_row(mat,row,vec)	_set_row(mat,row,vec,0) 
/* col j of A <- vec */
#define	set_col(mat,col,vec)	_set_col(mat,col,vec,0)


/* Basic permutation operations */
#ifndef ANSI_C

extern	PERM *px_mlt(), *px_inv(), *px_transp();
extern	int  px_sign();

#else

extern	PERM	*px_mlt(const PERM *px1,const PERM *px2,PERM *out),	/* out <- px1.px2 */
		*px_inv(const PERM *px,PERM *out),	/* out <- px^{-1} */
                 /* swap px[i] and px[j] */
		*px_transp(PERM *px,unsigned int i,unsigned int j);

     /* returns sign(px) = +1 if px product of even # transpositions
                           -1 if ps product of odd  # transpositions */
extern	int	px_sign(const PERM *);

#endif /* ANSI_C */


/* Basic integer vector operations */
#ifndef ANSI_C

extern	IMeVEC	*iv_add(), *iv_sub(), *iv_sort();

#else

extern	IMeVEC	*iv_add(const IMeVEC *ix,const IMeVEC *iy,IMeVEC *out),  
  /* out <- ix + iy */
		*iv_sub(const IMeVEC *ix,const IMeVEC *iy,IMeVEC *out),  
  /* out <- ix - iy */
  /* sorts ix & sets order so that sorted ix[i] = old ix[order[i]] */
		*iv_sort(IMeVEC *ix, PERM *order);

#endif /* ANSI_C */


/* miscellaneous functions */

#ifndef ANSI_C

extern	double	Mesquare(), Mecube(), mrand();
extern	void	smrand(), mrandlist();
extern  void    m_dump(), px_dump(), v_dump(), iv_dump();
extern MeMAT *band2mat();
extern BAND *mat2band();

#else

double	Mesquare(double x), 	/* returns x^2 */
  Mecube(double x), 		/* returns x^3 */
  mrand(void);                  /* returns random # in [0,1) */

void	smrand(int seed),            /* seeds mrand() */
  mrandlist(Real *x, int len);       /* generates len random numbers */

void    m_dump(FILE *fp,const MeMAT *a), px_dump(FILE *fp, const PERM *px),
        v_dump(FILE *fp,const MeVEC *x), iv_dump(FILE *fp, const IMeVEC *ix);

MeMAT *band2mat(const BAND *bA, MeMAT *A);
BAND *mat2band(const MeMAT *A, int lb,int ub, BAND *bA);

#endif /* ANSI_C */


/* miscellaneous constants */
#define	VNULL	((MeVEC *)NULL)
#define	MNULL	((MeMAT *)NULL)
#define	PNULL	((PERM *)NULL)
#define	IVNULL	((IMeVEC *)NULL)
#define BDNULL  ((BAND *)NULL)



/* varying number of arguments */

#ifdef ANSI_C
#include <stdarg.h>

/* prototypes */

int v_get_vars(int dim,...);
int iv_get_vars(int dim,...);
int m_get_vars(int m,int n,...);
int px_get_vars(int dim,...);

int v_resize_vars(int new_dim,...);
int iv_resize_vars(int new_dim,...);
int m_resize_vars(int m,int n,...);
int px_resize_vars(int new_dim,...);

int v_free_vars(MeVEC **,...);
int iv_free_vars(IMeVEC **,...);
int px_free_vars(PERM **,...);
int m_free_vars(MeMAT **,...);

#elif VARARGS
/* old varargs is used */

#include  <varargs.h>

/* prototypes */

int v_get_vars();
int iv_get_vars();
int m_get_vars();
int px_get_vars();

int v_resize_vars();
int iv_resize_vars();
int m_resize_vars();
int px_resize_vars();

int v_free_vars();
int iv_free_vars();
int px_free_vars();
int m_free_vars();

#endif /* ANSI_C */


#endif /* MeMATRIXH */


