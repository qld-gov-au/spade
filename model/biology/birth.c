#include <math.h>
#include "../../common.h"

Real b(

	 const Real a,
	 const Real x

	 )
{ /* birth function */
  return a*(A1*x + A2*pow(x,2.));
}
