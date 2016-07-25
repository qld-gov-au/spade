#ifndef SPADE_MATHPROP_H
#define SPADE_MATHPROP_H

#include "../meschach/matrix.h"
#include "../parameters.h"

Real sple(

	  const int nk,
	  const Real xval,
	  const VEC *knots,
	  const VEC *coef

	  );



VEC *numgrad(

	     Real (*model)(VEC *,void *),
	     void *stuff,
	     VEC *par,
	     Real epsilon

	     );

Real ConditionNumber(Parameters *,Data *);

 #endif
