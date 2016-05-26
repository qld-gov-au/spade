#ifndef SPADE_UTIL_H
#define SPADE_UTIL_H

#include <math.h>
#include "../meschach/matrix.h"


VEC *idxremove(

	       VEC *zn,
	       VEC *z,
	       int idx

	       );

double idxselect(

		 double omega,
		 VEC *xn

		 );

double get_bw(

	      VEC *dt

	      );
#endif