#ifndef SPADE_OBJFNS_H
#define SPADE_OBJFNS_H

#include "../meschach/matrix.h"
#include "../common.h"
#include "../parameters.h"


void Kveryfast(

	   Parameters * parameters

	   );


void Kfast(

	   Parameters * parameters

	   );

Real K_no(

	  void *args	  
	  );

/*Real K_no2(

	  Parameters * parameters
	  
	  );
*/

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
