#include "../../../meschach/matrix.h"
#include "../../../socbio/variable/birth.h"


/* Quadrature plus implicit.. */
void Q2_kappa(

	      double a,
	      double a2,
	      double k,
	      double w,
	      VEC *x,
	      VEC *u,
	      VEC *p

	      )
{

  // this function has not been updated for new birth function

  if (x->dim != p->dim)
    {
      printf("%d %d\n",x->dim,p->dim);
      error(E_SIZES,"Q2");
    }

  double rt = x->ve[1] * b(a,x->ve[1])*p->ve[1] - x->ve[0]*b(a,x->ve[1])*p->ve[1]; 

  for (int j=1;j<x->dim-1;j++) 
    rt = rt + (b(a,x->ve[j])*p->ve[j] + b(a,x->ve[j+1])*p->ve[j+1]) * (x->ve[j+1]-x->ve[j]);

  p->ve[0] = (rt - 2*w*u->ve[0]) / (2*k*w + x->ve[0]*b(a,x->ve[0]) - x->ve[1]*b(a,x->ve[0])); 

}
