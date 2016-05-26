#include "../../meschach/matrix.h"
#include "../../socbio/fixed/selectivity.h"
#include "../../socbio/variable/effort.h"

double zstar(

	      VEC *eff,
	      double b,
	      double g,
	      double k,
	      double i,
	      double t,
	      double x,
	      double U,
	      double r

	     )
{ /* death function: beta + gamma U + s(x)f(t) - kappa */

  return b + g*U + s(x)*i*e(eff,r,t) - k;

}