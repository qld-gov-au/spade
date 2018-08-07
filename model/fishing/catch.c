// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

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
      int idx = floor((t + (r/2) - MACHEPS)/r);
      return ca->ve[idx];
    }
}
