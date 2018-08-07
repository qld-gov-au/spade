// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#ifndef SPADE_OBJFNS_H
#define SPADE_OBJFNS_H

#include "../meschach/matrix.h"
#include "../common.h"
#include "../parameters.h"

Real K(

     Parameters * parameters,
		 Data *data,
		 Solve_Core_Args *core_args
		 	
		 );

Real K_dr(

  Parameters * parameters,
  Data *data
  );

Real G(

	 MAT *p,
	 MAT *x,
	 MAT *u,
	 Data *data,
	 Real iota
	
	 );

Real G_ni(

	    MAT *p,
	    MAT *x,
	    MAT *u,
	    Data *data,
	    Real iota

	    );

#endif
