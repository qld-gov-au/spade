// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#ifndef SPADE_UTIL_H
#define SPADE_UTIL_H

#include <math.h>
#include "../meschach/matrix.h"


VEC *idxremove(

	       VEC *zn,
	       VEC *z,
	       int idx

	       );

Real idxselect(

		 Real omega,
		 VEC *xn

		 );

Real get_bw(

	      VEC *dt

	      );
#endif
