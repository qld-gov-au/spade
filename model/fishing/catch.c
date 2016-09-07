#include <math.h>
#include "../../meschach/matrix.h"

Real c(

	  MeVEC * ca,
	  Real r,
	  Real t

	  )
{

  if (t<0)
    return 0;
  else
    {
      int idx = floor((t + (r/2) - MACHEPS)/r);
      return ca->ve[idx];
    }
}