// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#include "../meschach/matrix.h"
#include "../model/biology/birth.h"

/* Quadrature plus implicit.. */
void Q2(

	Real a,
	Real k,
	Real w,
	VEC *x,
	VEC *u
  
	)
{

  if (x->dim != u->dim)
    {
      printf("%d %d\n",x->dim,u->dim);
      error(E_SIZES,"Q2");
    }

  Real rt = x->ve[1] * b(a,x->ve[1])*u->ve[1] - x->ve[0]*b(a,x->ve[1])*u->ve[1]; 

  for (int j=1;j<x->dim-1;j++) 
    rt = rt + (b(a,x->ve[j])*u->ve[j] + b(a,x->ve[j+1])*u->ve[j+1]) * (x->ve[j+1]-x->ve[j]);

  u->ve[0] = rt / (2*k*w + x->ve[0]*b(a,x->ve[0]) - x->ve[1]*b(a,x->ve[0]));

}
