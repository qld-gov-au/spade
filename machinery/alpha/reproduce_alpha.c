// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#include "../../meschach/matrix.h"
#include "../../common.h"

/* Quadrature plus implicit.. */
void Q2_alpha(

	       Real a,
	       Real k,
	       Real w,
	       VEC *x,
	       VEC *u,
	       VEC *p

	       )
{

  if (x->dim != u->dim)
    {
      printf("%d %d\n",x->dim,u->dim);
      error(E_SIZES,"Q2");
    }

  Real x0 = x->ve[0];
  Real x1 = x->ve[1];

  Real u0 = u->ve[0];
  Real u1 = u->ve[1];

  Real p0 = p->ve[0];
  Real p1 = p->ve[1];

  Real rt = A1*x0*u0*x1 + A2*x0*x0*u0*x1 + a*A1*x1*p1*x1 + a*A2*x1*x1*p1*x1 + A1*x1*u1*x1 + A2*x1*x1*u1*x1 - A1*x0*u0*x0 - A2*x0*x0*u0*x0 - a*A1*x1*p1*x0 - a*A2*x1*x1*p1*x0 - A1*x1*u1*x0 - A2*x1*x1*u1*x0;

  for (int j=1;j<x->dim-1;j++) 
    {

      x0 = x->ve[j];
      x1 = x->ve[j+1];

      u0 = u->ve[j];
      u1 = u->ve[j+1];

      p0 = p->ve[j];
      p1 = p->ve[j+1];

      rt += ( a*A1*x0*p0 + a*A2*x0*x0*p0 + A1*x0*u0 + A2*x0*x0*u0 + a*A1*x1*p1 + a*A2*x1*x1*p1 + A1*x1*u1 + A2*x1*x1*u1 ) * (x1-x0);

    }

  x0 = x->ve[0];
  x1 = x->ve[1];

  p->ve[0] =  rt / (2*k*w + a*A1*x0*x0 + a*A2*x0*x0*x0 - a*A1*x0*x1 - a*A2*x0*x0*x1); 

}
