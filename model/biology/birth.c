#include <math.h>
#include "../../common.h"

Real b(

	 const Real a,
	 const Real x

	 )
{ /* birth function */
  return a*(A1*x + A2*pow(x,2.));
}

Real _b(

	 const Real a1,
	 const Real a2,
	 const Real x

	 )
{ /* birth function */
  return a1*x + a2*x*x;
}
