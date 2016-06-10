#include <math.h>
#include "../../meschach/matrix.h"

Real e(

	  VEC * ef,
	  Real r,
	  Real t

	  )
{

  if (t<0)
    {
      Real cek = r/4;
      Real ept5 = ef->ve[(int)floor((.5 + (cek/2) - 1e-12)/cek)];
      Real m = ept5/24;

      return m*t+ ept5;
    }
  else
    {
      Real cek = r/4;
      int idx = floor((t + (cek/2) - 1e-12)/cek); // better to use Real epsilon?
      return ef->ve[idx];
    }
}