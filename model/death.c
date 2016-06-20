#include "../meschach/matrix.h"
#include "fishing/selectivity.h"
#include "fishing/effort.h"

Real death(

	      VEC *eff,
	      Real b,
	      Real g,
	      Real i,
	      Real t,
	      Real x,
	      Real U,
	      Real r,
        int Y

	     )
{ /* death function: beta + gamma U +  s(x) f(t) */

  return b + g*U + s(x)*i*e(eff,r,t,Y);
  
}
