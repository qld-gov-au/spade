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
	      Real r

	     )
{ // equation X in Y

  return death(eff,b,g,i,t,x,U,r) - kappa;

}
