#ifndef SPADE_OBJFNS_H
#define SPADE_OBJFNS_H

#include "../meschach/matrix.h"
#include "../common.h"
#include "../parameters.h"

double K(

     Parameters * parameters,
		 Data *data,
		 Solve_Core_Args *core_args
		 	
		 );

double K_dr(

  Parameters * parameters,
  Data *data
  );

double G(

	 MAT *p,
	 MAT *x,
	 MAT *u,
	 Data *data,
	 double iota
	
	 );

double G_ni(

	    MAT *p,
	    MAT *x,
	    MAT *u,
	    Data *data,
	    double iota

	    );

#endif
