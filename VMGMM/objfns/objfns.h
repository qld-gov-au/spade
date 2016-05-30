#ifndef SPADE_OBJFNS_H
#define SPADE_OBJFNS_H

#include "../../meschach/matrix.h"

double H(

		 MAT *x,
		 MAT *u,
		 struct DATA *data,
		 double iota
	
		 );

double G(

	 MAT *p,
	 MAT *x,
	 MAT *u,
	 struct DATA *data,
	 double iota
	
	 );

double G_ni(

	    MAT *p,
	    MAT *x,
	    MAT *u,
	    struct DATA *data,
	    double iota

	    );

#endif
