#ifndef SPADE_INITIAL_H
#define SPADE_INITIAL_H

#include "../meschach/matrix.h"
#include "../common.h"

VEC * calc_alpha(

			Real a,
			Real k,
			Real w,
			Real bt,
			Real f

			);

VEC * VMGMM_eq(

		  VEC *x,
		  Data *d,
		  VEC *grad,
		  Real *f

		  );

#endif