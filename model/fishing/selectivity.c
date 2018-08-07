// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#include <math.h>
#include "../../common.h"
#include "../../meschach/matrix.h"

Real s(Real x) 
{ 
  if (x<58)
    return 0;
  else if (x <=60)
    {
      Real m = exp(-pow(60-phi*iota1,2.)/(2*iota2*pow(phi,2.)))/2;
      return m*(x-58);
    }
  else
    return exp(-pow(x-phi*iota1,2.)/(2*iota2*pow(phi,2.)));
}
