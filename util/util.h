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