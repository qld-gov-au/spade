#ifndef SPADE_SOLVE_H
#define SPADE_SOLVE_H

#include "../../meschach/matrix.h"

VEC *initial(

	     VEC *theta,
	     VEC *x,
	     VEC *u

	     );

void solve(

		 VEC *theta,
		 MAT *x,
		 MAT *u,
		 MAT *xhh,
		 MAT *xh,
		 MAT *xn,
		 MAT *uh,
		 MAT *un,
		 VEC *Ui,
		 VEC *Uh,
		 VEC *Uhh,
		 IVEC *idxi,
		 VEC *eff,
		 double k,		 
		 int S

		 );
#endif