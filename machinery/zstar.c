#include "../meschach/matrix.h"
#include "../model/death.h"

double zstar(

	      VEC *eff,
	      double b,
	      double g,
	      double kappa,
	      double i,
	      double t,
	      double x,
	      double U,
	      double r

	     )
{ // equation X in Y

  return death(eff,b,g,i,t,x,U,r) - kappa;

}
