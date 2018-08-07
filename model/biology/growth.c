// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#include "../../meschach/matrix.h"

Real g(

	 const Real kappa,
	 const Real w,
	 const Real x

	 )
{ /* von-Bertalanffy growth */
  return kappa*(w - x);
}
