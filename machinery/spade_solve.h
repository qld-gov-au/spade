#ifndef SPADE_SOLVE_H
#define SPADE_SOLVE_H

#include "../meschach/matrix.h"
#include "../common.h"
#include "../parameters.h"

MeVEC *initial(

	     Parameters *parameters,
	     MeVEC *x,
	     MeVEC *u

	     );

void solve(

		 Parameters *parameters,
		 MeVEC *eff,
		 Real k,		 
		 int S,
     int Y,
		 Solve_Core_Args *core_args  // modified

		 );
#endif