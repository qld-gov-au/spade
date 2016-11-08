#ifndef SPADE_MACHINERY_H
#define SPADE_MACHINERY_H

#include "../meschach/matrix.h"
#include "../parameters.h"
#include "../common.h"


VEC *_VMGMM(

		  VEC *,
		  VEC *,
		  Real *,
		  Parameters * parameters

		  );


VEC *VMGMM(

		  VEC *,
		  Data *,
		  VEC *,
		  Real *,
      Parameters * parameters

		  );

void check_derivative2(Parameters * parameters, Real * f);
void check_derivative(Parameters * parameters, Grad_Args * args, Data * d, Real * f, Solve_Core_Args core_args);
void check_derivative_no_meschach(Parameters * parameters, Grad_Args_No_MESCHACH * args, Data * d, Real * f, Solve_Core_Args core_args);
#endif
