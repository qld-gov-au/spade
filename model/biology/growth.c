#include "../../meschach/matrix.h"

Real g(

	 const Real kappa,
	 const Real w,
	 const Real x

	 )
{ /* von-Bertalanffy growth */
  return kappa*(w - x);
}
