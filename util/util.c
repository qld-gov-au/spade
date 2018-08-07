// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#include <math.h>
#include "../meschach/matrix.h"

VEC *idxremove(

	       VEC *zn,
	       VEC *z,
	       int idx

	       )

{

  for (int i=0;i<idx;i++)
    z->ve[i]=zn->ve[i];

  for (int i=idx;i<z->dim;i++)
    z->ve[i]=zn->ve[i+1];

  return z;
}

Real idxselect(

		 Real omega,
		 VEC *xn

		 )
{
  int idx = -1;
  Real val = omega;
  for (int i=1;i<xn->dim-1;i++)
    {
      Real dif = xn->ve[i+1]-xn->ve[i-1];
      if (dif < val)
	{
	  val = dif;
	  idx = i;
	}
    }

  return idx;

}

Real get_bw(

	      VEC *dt

	      )
{ /* Bandwidth. using Silverman's rule of thumb. */

  Real sm = v_sum(dt);
  Real mn = sm/(float)dt->dim;
  Real sd = 0;
  
  for (int i=0;i<dt->dim;i++)
    sd = sd + pow(dt->ve[i] - mn,2.);
  sd = sd / (float)dt->dim;
  sd = sqrt(sd);
  Real hi = sd;
  
  PERM *order = px_get(dt->dim);
  v_sort(dt,order);
  
  int fstQi = floor(dt->dim/4);
  Real fstQ = dt->ve[fstQi];
  int thrQi = floor(3*dt->dim/4);
  Real thrQ = dt->ve[thrQi];
  Real IQR = thrQ - fstQ;
  
  IQR = IQR/1.34;

  Real lo = (hi < IQR) ? hi : IQR;
  Real bw = 0.9 * lo * pow((Real)dt->dim,-0.2);

  PX_FREE(order);

  return bw;

}
