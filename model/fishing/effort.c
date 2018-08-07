// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#include <math.h>
#include "../../meschach/matrix.h"
#include "../../common.h"

Real e(

	  VEC * ef,
	  Real r,
	  Real t,
    int Y

	  )
{

  Real cek;
  if (QUARTER)
    cek = r/4;
  else    
    cek = r/2;

  if (t<0)
    {      
      Real ept5 = ef->ve[(int)floor((.5 + (cek/2) - MACHEPS)/cek)];		
      Real m = ept5/Y;	

      return m*t+ ept5;
    }
  else
    {
      int idx = floor((t + (cek/2) - MACHEPS)/cek);
      return ef->ve[idx];
    }
}
