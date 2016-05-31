#ifndef SPADE_MATHPROP_H
#define SPADE_MATHPROP_H

#include "../meschach/matrix.h"

VEC *numgrad(

	     double (*model)(VEC *,void *),
	     void *stuff,
	     VEC *par,
	     double epsilon

	     );

double ConditionNumber(VEC *,struct DATA *);

 #endif
