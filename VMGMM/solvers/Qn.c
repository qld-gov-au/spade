#include "../../meschach/matrix.h"

double Qn(

	  VEC * x,
	  VEC * u,
	  int n

	  )
{ /* Q for quadrature. Integrates u over x. */

  //if (x->dim != u->dim)
  //  error(E_SIZES,"Q");

  double rt = 0.;

  for (int i=0;i<n-1;i++) 
    rt = rt + .5 * (x->ve[i+1] - x->ve[i]) * (u->ve[i] + u->ve[i+1]);
  
  return rt;

}
