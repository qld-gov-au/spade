#ifndef SPADE_INITIAL_H
#define SPADE_INITIAL_H

#include "../meschach/matrix.h"
#include "../common.h"

MeVEC * calc_alpha(

			Real a,
			Real k,
			Real w,
			Real bt,
			Real f

			);

MeVEC * VMGMM_eq(

		  MeVEC *x,
		  Data *d,
		  MeVEC *grad,
		  Real *f

		  );

#endif