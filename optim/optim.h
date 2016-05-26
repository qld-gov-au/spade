#ifndef SPADE_OPTIM_H
#define SPADE_OPTIM_H

VEC * bfgs(

	   VEC * (*model)(VEC *,struct DATA *,VEC *,double *),
	   VEC *x,
	   struct DATA *data

	   );

int cstep(

	  double *stx,
	  double *fx,
	  double *dx,
	  double *sty,
	  double *fy,
	  double *dy,
	  double *stp,
	  double fp,
	  double dp,
	  int *brackt,
	  double stpmin,
	  double stpmax

	  );

int cvsrch(VEC *(*f)(VEC *,struct DATA *,VEC *,double *),VEC *,double,VEC *,VEC *,double,double,double,double,double,double,int,struct DATA *); // More-Thuente line search taken from code by Nocedal and Dianne O'Leary
#endif
