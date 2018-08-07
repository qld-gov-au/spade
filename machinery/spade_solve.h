// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

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
		 Real k,		 
		 int S,
     int Y,
		 Solve_Core_Args *core_args  // modified

		 );
#endif
