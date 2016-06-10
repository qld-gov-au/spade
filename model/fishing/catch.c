#include <math.h>
#include "../../meschach/matrix.h"

Real c(

	  VEC * ca,
	  Real r,
	  Real t

	  )
{

  if (t<0)
    return 0;
  else
    {
      #if REAL == DOUBLE
          int idx = floor((t + (r/2) - DBL_EPSILON)/r);
      #elif REAL == FLOAT
          int idx = floor((t + (r/2) - FLT_EPSILON)/r);
      #elif REAL == LONGDOUBLE
          int idx = floor((t + (r/2) - LDBL_EPSILON)/r);
      #endif

      return ca->ve[idx];
    }
}