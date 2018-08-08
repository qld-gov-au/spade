// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#ifndef SPADE_ZSTAR_H
#define SPADE_ZSTAR_H

#include "../meschach/matrix.h"

Real zstar(

	      VEC *eff,
	      Real b,
	      Real g,
	      Real k,
	      Real i,
	      Real t,
	      Real x,
	      Real U,
	      Real r,
        int Y

	     );

#endif
