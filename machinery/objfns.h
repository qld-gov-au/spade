﻿#ifndef SPADE_OBJFNS_H
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

	 MeMAT *p,
	 MeMAT *x,
	 MeMAT *u,
	 Data *data,
	 Real iota
	
	 );

Real G_ni(

	    MeMAT *p,
	    MeMAT *x,
	    MeMAT *u,
	    Data *data,
	    Real iota

	    );

#endif
