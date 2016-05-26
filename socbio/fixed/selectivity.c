﻿#include <math.h>
#include "../../common.h"

double s(double x) 
{ 
  if (x<58)
    return 0;
  else if (x <=60)
    {
      double m = exp(-pow(60-phi*iota1,2.)/(2*iota2*pow(phi,2.)))/2;
      return m*(x-58);
    }
  else
    return exp(-pow(x-phi*iota1,2.)/(2*iota2*pow(phi,2.)));
}
