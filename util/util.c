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

double idxselect(

		 double omega,
		 VEC *xn

		 )
{
  int idx = -1;
  double val = omega;
  for (int i=1;i<xn->dim-1;i++)
    {
      double dif = xn->ve[i+1]-xn->ve[i-1];
      if (dif < val)
	{
	  val = dif;
	  idx = i;
	}
    }

  return idx;

}

double get_bw(

	      VEC *dt

	      )
{ /* Bandwidth. using Silverman's rule of thumb. */

  double sm = v_sum(dt);
  double mn = sm/(float)dt->dim;
  double sd = 0;
  
  for (int i=0;i<dt->dim;i++)
    sd = sd + pow(dt->ve[i] - mn,2.);
  sd = sd / (float)dt->dim;
  sd = sqrt(sd);
  double hi = sd;
  
  PERM *order = px_get(dt->dim);
  v_sort(dt,order);
  
  int fstQi = floor(dt->dim/4);
  double fstQ = dt->ve[fstQi];
  int thrQi = floor(3*dt->dim/4);
  double thrQ = dt->ve[thrQi];
  double IQR = thrQ - fstQ;
  
  IQR = IQR/1.34;

  double lo = (hi < IQR) ? hi : IQR;
  double bw = 0.9 * lo * pow((double)dt->dim,-0.2);

  PX_FREE(order);

  return bw;

}