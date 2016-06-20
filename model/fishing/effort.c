#include <math.h>
#include "../../meschach/matrix.h"

Real e(

	  VEC * ef,
	  Real r,
	  Real t,
    int Y

	  )
{

  if (t<0)
    {
      Real cek = r/4;
      Real ept5 = ef->ve[(int)floor((.5 + (cek/2) - MACHEPS)/cek)];
      Real m = ept5/Y;

      return m*t+ ept5;
    }
  else
    {
      Real cek = r/4;
      int idx = floor((t + (cek/2) - MACHEPS)/cek);
      return ef->ve[idx];
    }
}