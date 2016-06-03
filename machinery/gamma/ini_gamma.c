#include <math.h>
#include "../../meschach/matrix.h"
#include "../../common.h"

VEC *ini_gamma(

	       Parameters *parameters,
	       VEC *x,
	       VEC *p

		)
{

  x->ve[x->dim-1] -= 1e-5;

  double a = parameters->alpha.value;
  double b = parameters->beta.value;
  double g = parameters->gamma.value;
  double k = parameters->kappa.value;
  double w = parameters->omega.value;

  double zeta = sqrt( 81*k*k*w*w*pow(a*A1+2*a*A2*w,2.) - 12*k*pow(a*A1*w+k,3.) );
  double eta = 9*a*A1*k*k*w + 18*a*A2*k*k*w*w + k*zeta;
  double Z = pow(eta,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*k*(a*A1*w+k) / pow(eta,1./3);

  for (int j=0;j<x->dim;j++)
    p->ve[j] = -pow(1-x->ve[j]/w,Z/k - 2.) * ((Z-b-k)/(1e-7*pow(g,2.)*Z)) * (a*A1 + 2*a*A2*k*w/(Z+k));

}