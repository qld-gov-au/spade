#ifndef SPADE_MACHINERY_H
#define SPADE_MACHINERY_H

#include "../meschach/matrix.h"
#include "../parameters.h"

VEC *VMGMM(

		  VEC *,
		  Data *,
		  VEC *,
		  Real *,
      Parameters * parameters

		  );

void check_derivative(Parameters * parameters, Grad_Args * args, Data * d, Real * f, Solve_Core_Args core_args);
#endif
