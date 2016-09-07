#ifndef SPADE_MeMATHPROP_H
#define SPADE_MeMATHPROP_H

#include "../meschach/matrix.h"
#include "../parameters.h"

Real sple(

	  const int nk,
	  const Real xval,
	  const MeVEC *knots,
	  const MeVEC *coef

	  );



MeVEC *numgrad(

	     Real (*model)(MeVEC *,void *),
	     void *stuff,
	     MeVEC *par,
	     Real epsilon

	     );

Real ConditionNumber(Parameters *,Data *);

 #endif
