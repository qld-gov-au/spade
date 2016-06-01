#ifndef SPADE_INITIAL_H
#define SPADE_INITIAL_H

#include "../meschach/matrix.h"
#include "../common.h"

VEC * calc_alpha(

			double a,
			double k,
			double w,
			double bt,
			double f

			);

VEC * VMGMM_eq(

		  VEC *x,
		  Data *d,
		  VEC *grad,
		  double *f

		  );

#endif