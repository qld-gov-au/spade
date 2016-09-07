#include "../meschach/matrix.h"

Real Q(

	 MeVEC * x,
	 MeVEC * u

	 )
{ /* Q for quadrature. Integrates u over x. */

  if (x->dim != u->dim)
    error(E_SIZES,"Q");

  Real rt = 0.;

  for (int i=0;i<x->dim-1;i++) 
    rt = rt + .5 * (x->ve[i+1] - x->ve[i]) * (u->ve[i] + u->ve[i+1]);
   
  
  return rt;

}