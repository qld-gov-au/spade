#ifndef SPADE_MATHPROP_H
#define SPADE_MATHPROP_H

#include "../meschach/matrix.h"
#include "../parameters.h"

VEC *numgrad(

	     Real (*model)(VEC *,void *),
	     void *stuff,
	     VEC *par,
	     Real epsilon

	     );

Real ConditionNumber(Parameters *,Data *);

 #endif
