#ifndef SPADE_OPTIM_H
#define SPADE_OPTIM_H

#include "../parameters.h"

VEC * bfgs(

	   VEC * (*model)(VEC *,Data *,VEC *,double *,Parameters *),
	   VEC *x,
	   Data *data,
     Parameters * parameters,
     OptimControl optim

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

int cvsrch(VEC *(*f)(VEC *,Data *,VEC *,double *,Parameters *),VEC *,double,VEC *,VEC *,double,double,double,double,double,double,int,Data *,Parameters *,double *fv); // More-Thuente line search taken from code by Nocedal and Dianne O'Leary
#endif
