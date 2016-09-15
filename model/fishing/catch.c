#include <math.h>
#include "../../meschach/matrix.h"
#include "../../common.h"

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


Real _c(

	  Real * ca,
	  Real t

	  )
{

  int idx = floor((t + (k/2) - MACHEPS)/k);
  return ca[idx];
    
}

