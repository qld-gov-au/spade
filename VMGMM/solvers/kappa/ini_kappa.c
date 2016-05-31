#include <math.h>
#include "../../../meschach/matrix.h"
#include "../../../common.h"

VEC *ini_kappa(

	       VEC *theta,
	       VEC *x,
	       VEC *p

	       )
{

  x->ve[x->dim-1] -= 1e-5;

  double a = theta->ve[0];
  double b = theta->ve[1];
  double g = theta->ve[2];
  double k = theta->ve[3];
  double w = theta->ve[4];

  double zeta = sqrt( 81*k*k*w*w*pow(a*A1+2*a*A2*w,2.) - 12*k*pow(a*A1*w+k,3.) );
  double eta = 9*a*A1*k*k*w + 18*a*A2*k*k*w*w + k*zeta;
  double Z = pow(eta,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*k*(a*A1*w+k) / pow(eta,1./3);

  double Zk =  6*k*(a*A1*w*zeta+2*a*A2*w*w*zeta-(2*k+a*A1*w)*pow(k+a*A1*w,2.)+9*k*w*w*pow(a*A1+2*a*A2*w,2.) ) / (zeta*pow(9*a*A1*k*k*w+18*a*A2*k*k*w*w+k*zeta,2./3)) * ( 1./(pow(2.,1/3.)*pow(3.,2/3.)) - pow(2/3.,1/3.)*k*(k+a*A1*w) / pow(9*a*A1*k*k*w+18*a*A2*k*k*w*w+k*zeta,2./3) ) + pow(2/3.,1/3.)*(2*k+a*A1*w) / pow(9*a*A1*k*k*w+18*a*A2*k*k*w*w+k*zeta,1./3) ; // check order of precedence

  for (int j=0;j<x->dim;j++)
    p->ve[j] = pow(1-x->ve[j]/w,Z/k - 2.) * ( Zk/(g*Z) *( (a*A1 + (2*a*A2*k*w/(Z+k))) * ( (b+k)/Z + log(1-x->ve[j]/w) * (Z-b-k)/k ) - 2*a*A2*k*w*(Z-b-k)/pow(Z+k,2.) ) + ((Z-b-k)/(g*Z*(Z+k))) * (2*a*A2*w + 2*a*A2*k*w/(Z+k)) - (a*A1 + 2*a*A2*k*w/(Z+k)) * (1/(g*Z) + log(1-x->ve[j]/w) * (Z-b-k)/(g*k*k)) ); 

}
