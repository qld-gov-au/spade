﻿#include <math.h>
#include "../../meschach/matrix.h"
#include "../../common.h"

Real e(

	  MeVEC * ef,
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


Real _e(

	  Real * ef,
	  Real r,
	  Real t

	  )
{

  Real cek;

  cek = r/2;

  Real tt = t + 24;
    
  int idx = floor((tt + (cek/2) - MACHEPS)/cek);
  return ef[idx];
    
}
