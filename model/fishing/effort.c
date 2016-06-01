#include <math.h>
#include "../../meschach/matrix.h"

double e(

	  VEC * ef,
	  double r,
	  double t

	  )
{

  if (t<0)
    {
      double cek = r/4;
      double ept5 = ef->ve[(int)floor((.5 + (cek/2) - 1e-12)/cek)];
      double m = ept5/24;

      return m*t+ ept5;
    }
  else
    {
      double cek = r/4;
      int idx = floor((t + (cek/2) - 1e-12)/cek); // better to use double epsilon?
      return ef->ve[idx];
    }
}