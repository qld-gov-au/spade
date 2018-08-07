// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#ifndef SPADE_OPTIM_H
#define SPADE_OPTIM_H

#include "../parameters.h"

VEC * bfgs(

	   VEC * (*model)(VEC *,Data *,VEC *,Real *,Parameters *),
	   VEC *x,
	   Data *data,
     Parameters * parameters,
     OptimControl optim

	   );

#endif
