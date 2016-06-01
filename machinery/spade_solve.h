#ifndef SPADE_SOLVE_H
#define SPADE_SOLVE_H

#include "../meschach/matrix.h"

VEC *initial(

	     VEC *theta,
	     VEC *x,
	     VEC *u

	     );

void solve(

		 VEC *theta,		 
		 VEC *eff,
		 double k,		 
		 int S,
		 Solve_Core_Args *core_args  // modified

		 );
#endif