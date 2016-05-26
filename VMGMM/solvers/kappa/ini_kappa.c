#include <math.h>
#include "../../../meschach/matrix.h"
#include "../../../common.h"

VEC *ini_kappa(

	       VEC *theta,
	       VEC *x,
	       VEC *p

	       )
{

  // not updated for new birth function yet

  x->ve[x->dim-1] -= 1e-5;

  double a1 = theta->ve[0];
  double a2 = theta->ve[1];
  double b = theta->ve[2];
  double g = theta->ve[3];
  double k = theta->ve[4];
  double w = theta->ve[5];

  double zeta = sqrt( 81*k*k*w*w*pow(a1+2*a2*w,2.) - 12*k*pow(a1*w+k,3.) );
  double in = 9*a1*k*k*w + 18*a2*k*k*w*w + k*zeta;
  double Z = pow(in,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*k*(a1*w+k) / pow(in,1./3);

  double Zk =  6*k*(a1*w*zeta+2*a2*w*w*zeta-(2*k+a1*w)*pow(k+a1*w,2.)+9*k*w*w*pow(a1+2*a2*w,2.) ) / (zeta*pow(9*a1*k*k*w+18*a2*k*k*w*w+k*zeta,2./3)) * ( 1./(pow(2.,1/3.)*pow(3.,2/3.)) - pow(2/3.,1/3.)*k*(k+a1*w) / pow(9*a1*k*k*w+18*a2*k*k*w*w+k*zeta,2./3) ) + pow(2/3.,1/3.)*(2*k+a1*w) / pow(9*a1*k*k*w+18*a2*k*k*w*w+k*zeta,1./3) ; // check order of precedence

  for (int j=0;j<x->dim;j++)
    p->ve[j] = pow(1-x->ve[j]/w,Z/k - 2.) * ( Zk/(g*Z) *( (a1 + (2*a2*k*w/(Z+k))) * ( (b+k)/Z + log(1-x->ve[j]/w) * (Z-b-k)/k ) - 2*a2*k*w*(Z-b-k)/pow(Z+k,2.) ) + ((Z-b-k)/(g*Z*(Z+k))) * (2*a2*w + 2*a2*k*w/(Z+k)) - (a1 + 2*a2*k*w/(Z+k)) * (1/(g*Z) + log(1-x->ve[j]/w) * (Z-b-k)/(g*k*k)) ); 

}