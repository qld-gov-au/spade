#ifndef SPADE_VMGMM_H
#define SPADE_VMGMM_H

#include "../meschach/matrix.h"

VEC *VMGMM(

		  VEC *theta,
		  struct DATA *dataptr,
		  VEC *grad,
		  double *f

		  );
#endif
