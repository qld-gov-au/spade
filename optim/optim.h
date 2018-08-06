#ifndef SPADE_OPTIM_H
#define SPADE_OPTIM_H

#include "../parameters.h"

VEC * bfgs(

	   VEC * (*model)(VEC *,VEC *,Real *,Parameters *),
	   VEC *x,
     Parameters * parameters,
     OptimControl optim

	   );

int cstep(

	  Real *stx,
	  Real *fx,
	  Real *dx,
	  Real *sty,
	  Real *fy,
	  Real *dy,
	  Real *stp,
	  Real fp,
	  Real dp,
	  int *brackt,
	  Real stpmin,
	  Real stpmax

	  );

int cvsrch(VEC *(*f)(VEC *,VEC *,Real *,Parameters *),VEC *,Real,VEC *,VEC *,Real,Real,Real,Real,Real,Real,int,Parameters *,Real *fv); // More-Thuente line search taken from code by Nocedal and Dianne O'Leary
#endif
