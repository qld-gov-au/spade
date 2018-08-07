// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#include <math.h>
#include "../../common.h"

Real b(

	 const Real a,
	 const Real x

	 )
{ /* birth function */
  return a*(A1*x + A2*pow(x,2.));
}
