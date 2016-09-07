
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
	Header for sparse matrix stuff.
	Basic sparse routines to be held in sparse.c
*/

/* RCS id: $Id: sparse.h,v 1.2 1994/01/13 05:33:36 des Exp $ */

#ifndef SPARSEH

#define SPARSEH 


#include        "matrix.h"


/* basic sparse types */

typedef struct row_elt	{
	int	col, nxt_row, nxt_idx;
	Real	val;
		} row_elt;

typedef struct SPROW {
	int	len, Memaxlen, diag;
	row_elt	*elt;		/* elt[Memaxlen] */
		} SPROW;

typedef struct SPMeMAT {
	int	m, n, Memax_m, Memax_n;
	char	flag_col, flag_diag;
	SPROW	*row;		/* row[Memax_m] */
	int	*start_row;	/* start_row[Memax_n] */
	int	*start_idx;	/* start_idx[Memax_n] */
	      } SPMeMAT;

/* Note that the first allocated entry in column j is start_row[j];
	This starts the chain down the columns using the nxt_row and nxt_idx
	fields of each entry in each row. */

typedef struct pair { int pos;	Real val; } pair;

typedef struct SPMeVEC {
	int	dim, Memax_dim;
	pair	*elt;		/* elt[Memax_dim] */
	       } SPMeVEC;

#define	SMNULL	((SPMeMAT*)NULL)
#define	SVNULL	((SPMeVEC*)NULL)

/* Macro for speedup */
#define	sprow_idx2(r,c,hint)	\
	( ( (hint) >= 0 && (hint) < (r)->len && \
	   (r)->elt[hint].col == (c)) ? (hint) : sprow_idx((r),(c)) )



/* memory functions */

#ifdef ANSI_C
int sp_get_vars(int m,int n,int deg,...);
int sp_resize_vars(int m,int n,...);
int sp_free_vars(SPMeMAT **,...);
#elif VARARGS
int sp_get_vars();
int sp_resize_vars();
int sp_free_vars();

#endif /* ANSI_C */

/* Sparse Matrix Operations and Utilities */
#ifndef ANSI_C
extern	SPMeMAT	*sp_get(), *sp_copy(), *sp_copy2(),
			*sp_zero(), *sp_resize(), *sp_compact();
extern	double	sp_get_val(), sp_set_val();
extern	MeVEC	*sp_mv_mlt(), *sp_vm_mlt();
extern	int	sp_free();

/* Access path operations */
extern	SPMeMAT	*sp_col_access();
extern	SPMeMAT	*sp_diag_access();
extern  int     chk_col_access();

/* Input/output operations */
extern	SPMeMAT	*sp_finput();
extern	void sp_foutput(), sp_foutput2();

/* algebraic operations */
extern SPMeMAT *sp_smlt(), *sp_add(), *sp_sub(), *sp_mltadd();


/* sparse row operations */
extern	SPROW	*sprow_get(), *sprow_xpd(), *sprow_merge(), *sprow_mltadd(),
  *sprow_resize(), *sprow_copy();
extern SPROW *sprow_add(), *sprow_sub(), *sprow_smlt();
extern	double	sprow_set_val();
extern	void	sprow_foutput();
extern	int	sprow_idx(), sprow_free();

/* dump */
extern  void   sp_dump(), sprow_dump();
extern  MeMAT  *sp_m2dense();

#else
SPMeMAT	*sp_get(int,int,int), *sp_copy(const SPMeMAT *),
	*sp_copy2(const SPMeMAT *,SPMeMAT *),
	*sp_zero(SPMeMAT *), *sp_resize(SPMeMAT *,int,int),
	*sp_compact(SPMeMAT *,double);
double	sp_get_val(const SPMeMAT *,int,int), sp_set_val(SPMeMAT *,int,int,double);
MeVEC	*sp_mv_mlt(const SPMeMAT *, const MeVEC *, MeVEC *), 
        *sp_vm_mlt(const SPMeMAT *, const MeVEC *, MeVEC *);
int	sp_free(SPMeMAT *);

/* Access path operations */
SPMeMAT	*sp_col_access(SPMeMAT *);
SPMeMAT	*sp_diag_access(SPMeMAT *);
int     chk_col_access(const SPMeMAT *);

/* Input/output operations */
SPMeMAT	*sp_finput(FILE *);
void	sp_foutput(FILE *, const SPMeMAT *);

/* algebraic operations */
SPMeMAT *sp_smlt(const SPMeMAT *A,double alpha,SPMeMAT *B),
      *sp_add(const SPMeMAT *A,const SPMeMAT *B,SPMeMAT *C),
      *sp_sub(const SPMeMAT *A,const SPMeMAT *B,SPMeMAT *C),
      *sp_mltadd(const SPMeMAT *A,const SPMeMAT *B,double alpha,SPMeMAT *C);

/* sparse row operations */
SPROW	*sprow_get(int), *sprow_xpd(SPROW *r,int n,int type),
        *sprow_resize(SPROW *r,int n,int type),
	*sprow_merge(const SPROW *,const SPROW *,SPROW *,int type),
        *sprow_copy(const SPROW *,const SPROW *,SPROW *,int type),
	*sprow_mltadd(const SPROW *r1,const SPROW *r2, double alpha,
		      int j0, SPROW *r_out, int type);
SPROW *sprow_add(const SPROW *r1,const SPROW *r2, int j0,SPROW *r_out, int type), 
        *sprow_sub(const SPROW *r1,const SPROW *r2, int j0,SPROW *r_out, int type), 
        *sprow_smlt(const SPROW *r1,double alpha, int j0,SPROW *r_out, int type);
double	sprow_set_val(SPROW *,int,double);
int      sprow_free(SPROW *);
int	sprow_idx(const SPROW *,int);
void	sprow_foutput(FILE *,const SPROW *);

/* dump */
void    sp_dump(FILE *fp, const SPMeMAT *A);
void    sprow_dump(FILE *fp, const SPROW *r);
MeMAT	*sp_m2dense(const SPMeMAT *A,MeMAT *out);

#endif /* ANSI_C */

/* MACROS */

#define	sp_input()	sp_finput(stdin)
#define	sp_output(A)	sp_foutput(stdout,(A))
#define	sp_output2(A)	sp_foutput2(stdout,(A))
#define	row_mltadd(r1,r2,alpha,out)	sprow_mltadd(r1,r2,alpha,0,out)
#define	out_row(r)	sprow_foutput(stdout,(r))

#define SP_FREE(A)    ( sp_free((A)),  (A)=(SPMeMAT *)NULL) 

/* utility for index computations -- ensures index returned >= 0 */
#define	fixindex(idx)	((idx) == -1 ? (error(E_BOUNDS,"fixindex"),0) : \
			 (idx) < 0 ? -((idx)+2) : (idx))


/*  NOT USED */

/* loop over the columns in a row */
/*
#define	loop_cols(r,e,code) \
    do { int _r_idx; row_elt *e; SPROW *_t_row;			\
	  _t_row = (r); e = &(_t_row->elt);				\
	  for ( _r_idx = 0; _r_idx < _t_row->len; _r_idx++, e++ )	\
	  {  code;  }  }  while ( 0 )
*/
/* loop over the rows in a column */
/*
#define	loop_cols(A,col,e,code) \
    do { int _r_num, _r_idx, _c; SPROW *_r; row_elt *e;		\
	  if ( ! (A)->flag_col )	sp_col_access((A));		\
	  col_num = (col);						\
	  if ( col_num < 0 || col_num >= A->n )				\
	      error(E_BOUNDS,"loop_cols");				\
          _r_num = (A)->start_row[_c]; _r_idx = (A)->start_idx[_c];	\
	  while ( _r_num >= 0 )  {					\
	      _r = &((A)->row[_r_num]);					\
              _r_idx = sprow_idx2(_r,_c,_r_idx);			\
              if ( _r_idx < 0 )  continue;				\
	      e = &(_r->elt[_r_idx]);	code;				\
	      _r_num = e->nxt_row;	_r_idx = e->nxt_idx;		\
	      } } while ( 0 )

*/

#endif /* SPARSEH */

