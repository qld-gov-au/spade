#include "../meschach/matrix.h"
#include "fishing/selectivity.h"
#include "fishing/effort.h"

double death(

	      VEC *eff,
	      double b,
	      double g,
	      double i,
	      double t,
	      double x,
	      double U,
	      double r

	     )
{ /* death function: beta + gamma U +  s(x) f(t) */

  return b + g*U + s(x)*i*e(eff,r,t);
  
}
