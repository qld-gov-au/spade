#include "../meschach/matrix.h"

VEC *numgrad(

	     double (*model)(VEC *,void *),
	     void *stuff,
	     VEC *par,
	     double epsilon

	     )
{

  double f0 = (*model)(par,stuff);

  VEC *ei = v_get(par->dim);
  VEC *d = v_get(par->dim);
  VEC *fg = v_get(par->dim);

  for (int i=0;i<par->dim;i++)
    {
      ei->ve[i] = 1;
      sv_mlt(epsilon,ei,d);
      double f1 = (*model)(v_add(par,d,VNULL),stuff);
      fg->ve[i] = (f1 - f0) / d->ve[i];
      ei->ve[i] = 0;
    }

  return fg;

}