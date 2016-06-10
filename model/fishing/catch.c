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
      int idx = floor((t + (r/2) - 1e-12)/r); // better to use double epsilon?
      return ca->ve[idx];
    }
}