// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#include "../meschach/matrix.h"
#include "../model/death.h"

Real zstar(

	      VEC *eff,
	      Real b,
	      Real g,
	      Real kappa,
	      Real i,
	      Real t,
	      Real x,
	      Real U,
	      Real r,
        int Y

	     )
{ // equation X in Y

  return death(eff,b,g,i,t,x,U,r,Y) - kappa;

}
