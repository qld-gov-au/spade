#include <math.h>
#include "../../common.h"

double b(

	 const double a,
	 const double x

	 )
{ /* birth function */
  return a*(A1*x + A2*pow(x,2.));
}
