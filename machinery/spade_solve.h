#ifndef SPADE_SOLVE_H
#define SPADE_SOLVE_H

#include "../meschach/matrix.h"
#include "../common.h"
#include "../parameters.h"

VEC *initial(

	     Parameters *parameters,
	     VEC *x,
	     VEC *u

	     );

void solve(

		 Parameters *parameters,
		 VEC *eff,
		 double k,		 
		 int S,
		 Solve_Core_Args *core_args  // modified

		 );
#endif