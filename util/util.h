#ifndef SPADE_UTIL_H
#define SPADE_UTIL_H

#include <math.h>
#include "../meschach/matrix.h"


MeVEC *idxremove(

	       MeVEC *zn,
	       MeVEC *z,
	       int idx

	       );

Real idxselect(

		 Real omega,
		 MeVEC *xn

		 );

Real get_bw(

	      MeVEC *dt

	      );
#endif