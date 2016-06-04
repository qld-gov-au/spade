﻿#include <math.h>
#include "../../meschach/matrix.h"
#include "../../common.h"

VEC *ini_omega(

	       Parameters *parameters,
	       VEC *x,
	       VEC *p

		)
{

  x->ve[x->dim-1] -= 1e-5;

  double a = parameters->alpha.value;
  double b = parameters->beta.value;
  double g = parameters->gamma.value*1e-7;
  double k = parameters->kappa.value;
  double w = parameters->omega.value;

  double zeta = sqrt( 81*k*k*w*w*pow(a*A1+2*a*A2*w,2.) - 12*k*pow(a*A1*w+k,3.) );
  double in = 9*a*A1*k*k*w + 18*a*A2*k*k*w*w + k*zeta;
  double Z = pow(in,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*k*(a*a*A1*w+k) / pow(in,1./3);

  double Zw =  (3*k*k*(a*A1*zeta+4*a*A2*w*zeta-2*a*A1*pow(k+a*A1*w,2.)+18*a*A2*k*w*w*(a*A1+2*a*A2*w)+9*k*w*pow(a*A1+2*a*A2*w,2.) ) / (zeta*pow(9*a*A1*k*k*w+18*a*A2*k*k*w*w+k*zeta,2./3)) ) * ( 1./(pow(2.,1/3.)*pow(3.,2/3.)) - pow(2/3.,1/3.)*k*(k+a*A1*w) / pow(9*a*A1*k*k*w+18*a*A2*k*k*w*w+k*zeta,2./3) ) + pow(2/3.,1/3.)*a*A1*k / pow(9*a*A1*k*k*w+18*a*A2*k*k*w*w+k*zeta,1./3) ;

  for (int j=0;j<x->dim;j++)
    p->ve[j] = pow(1-x->ve[j]/w,Z/k - 2.) * ( (Z-b-k)/(g*Z) *( 2*k*w/(Z+k) - (2*a*A2*k*w/pow(Z+k,2.))*Zw + (a*A1 + 2*a*A2*k*w/(Z+k))*x->ve[j]*(Z-2*k)/(k*w*(w-x->ve[j]))) + (Zw/(g*Z)) * (a*A1+2*a*A2*k*w/(Z+k))*( (b+k)/Z + log(1-x->ve[j]/w) * (Z-b-k)/k ) );

}