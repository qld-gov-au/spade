#include <math.h>
#include "../../common.h"
#include "../../parameters.h"
#include "../objfns.h"
#include "../zstar.h"
#include "../Q.h"
#include "../reproduce.h"
#include "../../util/util.h"
#include "../../model/fishing/selectivity.h"
#include "../../model/fishing/effort.h"
#include "../../model/fishing/catch.h"
#include "../../model/biology/birth.h"
#include "../../model/biology/weight.h"


void grad_alpha1_fast(void *args)
{

  Parameters *parameters = (Parameters *)args;
  
  // only valid for models where fish are size zero at birth
       
  Real a1 = parameters->alpha1.value;
  Real a2 = parameters->alpha2.value;
  Real bb = parameters->beta.value;
  Real gg = parameters->gamma.value*1e-7;
  Real kk = parameters->kappa.value;
  Real ww = parameters->omega.value;
  Real ii = parameters->iota.value*1e-3;
  
  Real * restrict x = (Real *)calloc(d.J+1,sizeof(Real));
  Real * restrict u = (Real *)calloc(d.J+1,sizeof(Real));
  Real * restrict p = (Real *)calloc(d.J+1,sizeof(Real));
  Real * restrict r = (Real *)calloc(d.J+2,sizeof(Real));
  Real * restrict l = (Real *)calloc(d.J+2,sizeof(Real));
  Real * restrict o = (Real *)calloc(d.J+2,sizeof(Real));
    
  /* 
     initialize x
  */

  for (int j=0;j<=d.J;j++) 
    x[j] = h*j;
  x[d.J] = ww - 1e-12;
  
  /*
     ok, now initialize u and p
  */
  
  // prelims
  Real zeta = sqrt( 81*kk*kk*ww*ww*pow(a1+2*a2*ww,2.) - 12*kk*pow(a1*ww+kk,3.) );
  Real eta = 9*a1*kk*kk*ww + 18*a2*kk*kk*ww*ww + kk*zeta;
  Real Z = pow(eta,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*kk*(a1*ww+kk) / pow(eta,1./3);

  Real ubar = (Z - bb - kk) / gg; 
  Real vbar = (kk*ww*ubar) / (bb+gg*ubar+kk);
  Real wbar = (2*kk*ww*vbar) / (bb+gg*ubar+2*kk);  
  
  // set
  for (int j=0;j<=d.J;j++) 
    u[j] = (a1*vbar+a2*wbar)*pow(ww-x[j],(bb+gg*ubar)/kk-1) / (kk*pow(ww,(bb+gg*ubar)/kk));

  Real Za1 = ((kk*kk*ww*(3*zeta-6*pow(kk+a1*ww,2.)+27*kk*ww*(a1+2*a2*ww)) ) / (zeta*pow(9*a1*pow(kk,2.)*ww + 18*a2*kk*kk*ww*ww + kk*zeta,2./3.)) ) * ( (1/(pow(2.,1./3.)*pow(3.,2./3.))) - ( pow(2./3,1./3.)*kk*(kk+a1*ww)) / (pow(9*a1*kk*kk*ww + 18*a2*kk*kk*ww*ww + kk*zeta,2./3.)) ) + ( pow(2./3,1./3.)*kk*ww ) / pow(9*a1*kk*kk*ww + 18*a2*kk*kk*ww*ww + kk*zeta,1./3.);

  for (int j=0;j<=d.J;j++)
    {
      p[j] = pow(1-x[j]/ww,Z/kk - 2.) * ( (Z-bb-kk)/(gg*Z) * (1 - 2*a2*kk*ww/pow(Z+kk,2.) * Za1) + (Za1/(gg*Z))*(a1 + 2*a2*kk*ww/(Z+kk))*((bb+kk)/Z + log(1 - x[j]/ww)* (Z-bb-kk)/kk) );
      r[j] = 1e-3*w(x[j])*s(x[j])*u[j]*ii*_e(d.eff,0);
      l[j] = 1e-3*w(x[j])*s(x[j])*p[j]*ii*_e(d.eff,0);
    }
  

  Real ff = 0;
  
  // integrate sequentially
  Real U = 0;
  Real C = 0;
  Real P = 0;
  Real L = 0;
  
  for (int j=0;j<d.J;j++)
    {      
      U += .5 * (x[j+1] - x[j]) * (u[j+1] + u[j]);
      P += .5 * (x[j+1] - x[j]) * (p[j+1] + p[j]);
      C += .5 * (x[j+1] - x[j]) * (r[j+1] + r[j]);
      L += .5 * (x[j+1] - x[j]) * (l[j+1] + l[j]);
    }
  
  for (int j=0;j<=d.J;j++)
    o[j] = 2*(r[j] - _c(d.cat,0) * (d.p[0][j] + (1-d.Qp[0]) * r[j]/C))*(l[j] - _c(d.cat,0)*(1 - d.Qp[0]) * ( l[j] * C - r[j] * L ) / (C*C)); 
  
  for (int j=0;j<d.J;j++)
    ff += k * .5 * (x[j+1] - x[j]) * (o[j+1] + o[j]);
  
  Real * restrict xh = (Real *) calloc(d.J+1,sizeof(Real));  
  Real * restrict uh = (Real *) calloc(d.J+1,sizeof(Real));  
  Real * restrict ph = (Real *) calloc(d.J+1,sizeof(Real));  

  Real * restrict xn = (Real *) calloc(d.J+2,sizeof(Real));  
  Real * restrict un = (Real *) calloc(d.J+2,sizeof(Real));  
  Real * restrict pn = (Real *) calloc(d.J+2,sizeof(Real));  
  
  for (int i=1;i<=d.I;i++)
  {
  
      Real t = k*(i-1);
      Real th = k*(i-.5);

      // calculate uh
      for (int j=0;j<=d.J;j++)
	{	  
	  xh[j] = x[j] + k/2 * kk*(ww - x[j]);
	  uh[j] = u[j] * exp( -k/2 * (bb + gg*U + s(x[j])* ii * _e(d.eff,t) - kk) );
	}

      xh[d.J] = ww - 1e-12;
      
      // reproduce uh
      Real uh_0 = xh[0] * _b(a1,a2,xh[0])*uh[0];

      for (int j=0;j<d.J;j++)
	uh_0 += (_b(a1,a2,xh[j])*uh[j] + _b(a1,a2,xh[j+1])*uh[j+1]) * (xh[j+1]-xh[j]);

      uh_0 /= (2*kk*ww - xh[0]*_b(a1,a2,0));

      // integrate uh
      Real Uh = .5 * xh[0] * (uh_0 + uh[0]);
      for (int j=0;j<d.J;j++)
	Uh += .5 * (xh[j+1] - xh[j]) * (uh[j+1] + uh[j]);

      // calculate ph
      for (int j=0;j<=d.J;j++)
	ph[j] = p[j]*exp(-(k/2)*(bb + gg*U + s(x[j]) * ii * _e(d.eff,t) - kk)) - exp(-(k/2)*(bb + gg*U + s(x[j]) * ii * _e(d.eff,t) - kk))*(k/2)*gg*P*u[j];

      // reproduce ph
      Real ph_0 = a1*xh[0]*ph[0]*xh[0] + a2*xh[0]*xh[0]*ph[0]*xh[0] + xh[0]*uh[0]*xh[0];

      for (int j=0;j<d.J;j++)
	ph_0 += ( a1*xh[j]*ph[j] + a2*xh[j]*xh[j]*ph[j] + xh[j]*uh[j] + a1*xh[j+1]*ph[j+1] + a2*xh[j+1]*xh[j+1]*ph[j+1] + xh[j+1]*uh[j+1]) * (xh[j+1]-xh[j]);

      ph_0 /= 2*kk*ww;
      
      // integrate ph
      Real Ph = .5 * xh[0] * (ph_0 + ph[0]);
      for (int j=0;j<d.J;j++)
	Ph += .5 * (xh[j+1] - xh[j]) * (ph[j+1] + ph[j]);

      // calculate pn      
      for (int j=d.J+1;j>0;j--)
	{

	  xn[j] = x[j-1] + k * kk*(ww - xh[j-1]);
	  un[j] = u[j-1] * exp ( -k * (bb + gg*Uh + s(xh[j-1]) * ii * _e(d.eff,th) - kk));
	  Real tmp1 = exp(-k*(bb+gg*Uh+s(xh[j-1])*ii*_e(d.eff,th)-kk));
	  Real tmp2 = exp((k/2)*(bb+gg*Uh+s(xh[j-1])*ii*_e(d.eff,th)-kk));
	  pn[j] = p[j-1]*tmp1 - tmp1 * k*gg*Ph*uh[j-1] * tmp2;
	}

      xn[d.J+1] = ww - 1e-12;
      
      xn[0] = 0;

      // reproduce un
      un[0] = xn[1] * _b(a1,a2,xn[1])*un[1];

      for (int j=1;j<=d.J;j++)
	un[0] += (_b(a1,a2,xn[j])*un[j] + _b(a1,a2,xn[j+1])*un[j+1]) * (xn[j+1]-xn[j]);

      un[0] /= (2*kk*ww - xn[1]*_b(a1,a2,0));      
      
      // reproduce pn
      pn[0] = a1*xn[1]*pn[1]*xn[1] + a2*xn[1]*xn[1]*pn[1]*xn[1] + xn[1]*un[1]*xn[1];

      for (int j=1;j<=d.J;j++)
	pn[0] += ( a1*xn[j]*pn[j] + a2*xn[j]*xn[j]*pn[j] + xn[j]*un[j] + a1*xn[j+1]*pn[j+1] + a2*xn[j+1]*xn[j+1]*pn[j+1] + xn[j+1]*un[j+1]) * (xn[j+1]-xn[j]);

      pn[0] /= 2*kk*ww;

      // r and l for objective function
      for (int j=0;j<=d.J+1;j++)
	{
	  r[j] = 1e-3*w(xn[j])*s(xn[j])*un[j]*ii*_e(d.eff,k*i);
	  l[j] = 1e-3*w(xn[j])*s(xn[j])*pn[j]*ii*_e(d.eff,k*i);
	}

      U = 0;
      C = 0;
      P = 0;
      L = 0;

      // integrate
      for (int j=0;j<=d.J;j++)
	{
	  U += .5 * (xn[j+1] - xn[j]) * (un[j+1] + un[j]);
	  P += .5 * (xn[j+1] - xn[j]) * (pn[j+1] + pn[j]);
	  C += .5 * (xn[j+1] - xn[j]) * (r[j+1] + r[j]);
	  L += .5 * (xn[j+1] - xn[j]) * (l[j+1] + l[j]);
	}

      // objective function over x
      for (int j=0;j<=d.J+1;j++)
	o[j] = 2*(r[j] - _c(d.cat,k*i) * (d.p[i][j] + (1-d.Qp[i]) * r[j]/C))*(l[j] - _c(d.cat,k*i)*(1 - d.Qp[i]) * ( l[j] * C - r[j] * L ) / (C*C));

      // integrate to get final objective function value
      for (int j=0;j<=d.J;j++)
	ff += k * .5 * (xn[j+1] - xn[j]) * (o[j+1] + o[j]);

      // remove node
      for (int j=0;j<idx[i-1];j++)
	x[j] = xn[j];
      for (int j=idx[i-1];j<=d.J;j++)
	x[j] = xn[j+1];      

      for (int j=0;j<idx[i-1];j++)
	u[j] = un[j];
      for (int j=idx[i-1];j<=d.J;j++)
	u[j] = un[j+1];      

      for (int j=0;j<idx[i-1];j++)
	p[j] = pn[j];
      for (int j=idx[i-1];j<=d.J;j++)
	p[j] = pn[j+1];            
      
  }

  parameters->alpha1.gradient = ff;

  free(x);
  free(u);
  free(p);
  free(r);
  free(l);
  free(o);
  free(xh);
  free(uh);
  free(ph);
  free(xn);
  free(un);
  free(pn);
  
}


void grad_alpha1(void *args)
{

  Parameters *parameters = (Parameters *)args;
  
  // only valid for models where fish are size zero at birth
       
  Real a1 = parameters->alpha1.value;
  Real a2 = parameters->alpha2.value;
  Real bb = parameters->beta.value;
  Real gg = parameters->gamma.value*1e-7;
  Real kk = parameters->kappa.value;
  Real ww = parameters->omega.value;
  Real ii = parameters->iota.value*1e-3;
  
  Real * restrict x = (Real *)calloc(d.J,sizeof(Real));
  Real * restrict u = (Real *)calloc(d.J,sizeof(Real));
  Real * restrict p = (Real *)calloc(d.J,sizeof(Real));
  Real * restrict r = (Real *)calloc(d.J,sizeof(Real));
  Real * restrict l = (Real *)calloc(d.J,sizeof(Real));
  Real * restrict o = (Real *)calloc(d.J,sizeof(Real));
  //Real * restrict m = (Real *)calloc(d.J,sizeof(Real));
    
  /* 
     initialize x
  */
  int J = (d.J+1) - (d.I+1);

  // 'active': first J values in first 'row'
  for (int j=0;j<J;j++) 
    x[j] = h*j;

  // 'neutral': all the rest (J+1.. I+J)
  for (int j=J;j<=d.J;j++)
    x[j] = ww-1e-9;  // !! should be ww but that causes a divide by zero in equation 34: 
  //   for (int j=0;j<core_x->n;j++)
  //    p[j] = pow(1-x[j]/ww,Z/kk - 2.) * ( (Z-bb-kk)/(gg*Z) * (A1 + 2*A2*kk*ww/(Z+kk) - 2*aa*A2*kk*ww/pow(Z+kk,2.) * Za) + (Za/(gg*Z))*(aa*A1 + 2*aa*A2*kk*ww/(Z+kk))*((bb+kk)/Z + log(1 - x[j]/ww)* (Z-bb-kk)/kk) );

  /*
     ok, now initialize u and p
  */
  
  // prelims
  Real zeta = sqrt( 81*kk*kk*ww*ww*pow(a1+2*a2*ww,2.) - 12*kk*pow(a1*ww+kk,3.) );
  Real eta = 9*a1*kk*kk*ww + 18*a2*kk*kk*ww*ww + kk*zeta;
  Real Z = pow(eta,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*kk*(a1*ww+kk) / pow(eta,1./3);

  Real ubar = (Z - bb - kk) / gg; 
  Real vbar = (kk*ww*ubar) / (bb+gg*ubar+kk);
  Real wbar = (2*kk*ww*vbar) / (bb+gg*ubar+2*kk);  
  
  // set
  for (int j=0;j<=d.J;j++) 
    u[j] = (a1*vbar+a2*wbar)*pow(ww-x[j],(bb+gg*ubar)/kk-1) / (kk*pow(ww,(bb+gg*ubar)/kk));

  Real Za1 = ((kk*kk*ww*(3*zeta-6*pow(kk+a1*ww,2.)+27*kk*ww*(a1+2*a2*ww)) ) / (zeta*pow(9*a1*pow(kk,2.)*ww + 18*a2*kk*kk*ww*ww + kk*zeta,2./3.)) ) * ( (1/(pow(2.,1./3.)*pow(3.,2./3.))) - ( pow(2./3,1./3.)*kk*(kk+a1*ww)) / (pow(9*a1*kk*kk*ww + 18*a2*kk*kk*ww*ww + kk*zeta,2./3.)) ) + ( pow(2./3,1./3.)*kk*ww ) / pow(9*a1*kk*kk*ww + 18*a2*kk*kk*ww*ww + kk*zeta,1./3.);

  for (int j=0;j<=d.J;j++)
    {
      p[j] = pow(1-x[j]/ww,Z/kk - 2.) * ( (Z-bb-kk)/(gg*Z) * (1 - 2*a2*kk*ww/pow(Z+kk,2.) * Za1) + (Za1/(gg*Z))*(a1 + 2*a2*kk*ww/(Z+kk))*((bb+kk)/Z + log(1 - x[j]/ww)* (Z-bb-kk)/kk) );

      r[j] = 1e-3*w(x[j])*s(x[j])*u[j]*ii*_e(d.eff,0);
      l[j] = 1e-3*w(x[j])*s(x[j])*p[j]*ii*_e(d.eff,0);
    }
  

  Real ff = 0;
  
  // integrate sequentially
  Real U = 0;
  Real C = 0;
  Real P = 0;
  Real L = 0;
  
  for (int j=0;j<d.J;j++)
    {      
      U += .5 * (x[j+1] - x[j]) * (u[j+1] + u[j]);
      P += .5 * (x[j+1] - x[j]) * (p[j+1] + p[j]);
      C += .5 * (x[j+1] - x[j]) * (r[j+1] + r[j]);
      L += .5 * (x[j+1] - x[j]) * (l[j+1] + l[j]);
    }
  
  for (int j=0;j<=d.J;j++)
    o[j] = 2*(r[j] - _c(d.cat,0) * (d.p[0][j] + (1-d.Qp[0]) * r[j]/C))*(l[j] - _c(d.cat,0)*(1 - d.Qp[0]) * ( l[j] * C - r[j] * L ) / (C*C)); 

  //for (int j=0;j<=d.J;j++)
  //m[j] = pow(_c(d.cat,k,k*(0-d.N)) * (d.p[0][j] + (1-d.Qp[0]) * r[j]/C) - r[j],2.0);                      

  
  for (int j=0;j<d.J;j++)
    ff += k * .5 * (x[j+1] - x[j]) * (o[j+1] + o[j]);
  
  Real * restrict xh = (Real *) calloc(d.J,sizeof(Real));  
  Real * restrict uh = (Real *) calloc(d.J,sizeof(Real));  
  Real * restrict ph = (Real *) calloc(d.J,sizeof(Real));  
  
  for (int i=1;i<=d.I;i++)
  {
    //int i=1;
  
      Real t = k*(i-1);
      Real th = k*(i-.5);

      // calculate uh
      for (int j=0;j<=d.J;j++)
	{	  
	  xh[j] = x[j] + k/2 * kk*(ww - x[j]);
	  uh[j] = u[j] * exp( -k/2 * (bb + gg*U + s(x[j])* ii * _e(d.eff,t) - kk) );
	}

      // reproduce uh
      Real uh_0 = xh[0] * _b(a1,a2,xh[0])*uh[0];

      for (int j=0;j<d.J;j++)
	uh_0 += (_b(a1,a2,xh[j])*uh[j] + _b(a1,a2,xh[j+1])*uh[j+1]) * (xh[j+1]-xh[j]);

      uh_0 /= (2*kk*ww - xh[0]*_b(a1,a2,0));

      // integrate uh
      Real Uh = .5 * xh[0] * (uh_0 + uh[0]);
      for (int j=0;j<d.J;j++)
	Uh += .5 * (xh[j+1] - xh[j]) * (uh[j+1] + uh[j]);

      // calculate ph
      for (int j=0;j<=d.J;j++)
	ph[j] = p[j]*exp(-(k/2)*(bb + gg*U + s(x[j]) * ii * _e(d.eff,t) - kk)) - exp(-(k/2)*(bb + gg*U + s(x[j]) * ii * _e(d.eff,t) - kk))*(k/2)*gg*P*u[j];

      // reproduce ph
      Real ph_0 = a1*xh[0]*ph[0]*xh[0] + a2*xh[0]*xh[0]*ph[0]*xh[0] + xh[0]*uh[0]*xh[0];

      for (int j=0;j<d.J;j++)
	ph_0 += ( a1*xh[j]*ph[j] + a2*xh[j]*xh[j]*ph[j] + xh[j]*uh[j] + a1*xh[j+1]*ph[j+1] + a2*xh[j+1]*xh[j+1]*ph[j+1] + xh[j+1]*uh[j+1]) * (xh[j+1]-xh[j]);

      ph_0 /= 2*kk*ww;

      /*
  printf("\n");
  printf("%lf %lf\n",0.0,ph_0);
  for (int j=0;j<10;j++)
    printf("%lf %lf\n",xh[j],ph[j]);
  printf("e\n");

  Real delt=0.00000001;

  a1 += delt;
  
  zeta = sqrt( 81*kk*kk*ww*ww*pow(a1+2*a2*ww,2.) - 12*kk*pow(a1*ww+kk,3.) );
  eta = 9*a1*kk*kk*ww + 18*a2*kk*kk*ww*ww + kk*zeta;
  Z = pow(eta,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*kk*(a1*ww+kk) / pow(eta,1./3);

  ubar = (Z - bb - kk) / gg; 
  vbar = (kk*ww*ubar) / (bb+gg*ubar+kk);
  wbar = (2*kk*ww*vbar) / (bb+gg*ubar+2*kk);  

  Real * restrict u2 = (Real *)calloc(d.J,sizeof(Real));
  Real * restrict r2 = (Real *)calloc(d.J,sizeof(Real));
  
  // set
  for (int j=0;j<=d.J;j++)
    {
      u2[j] = (a1*vbar+a2*wbar)*pow(ww-x[j],(bb+gg*ubar)/kk-1) / (kk*pow(ww,(bb+gg*ubar)/kk));
      r2[j] = w(x[j])*s(x[j])*u2[j]*ii*_e(d.eff,k,k*(0-d.N));
    }
      
  // integrate sequentially
  Real U2 = 0;
  
  for (int j=0;j<d.J;j++)
    U2 += .5 * (x[j+1] - x[j]) * (u2[j+1] + u2[j]);
  
  Real * restrict uh2 = (Real *) calloc(d.J,sizeof(Real));  
  Real * restrict ph2 = (Real *) calloc(d.J,sizeof(Real));  
  
      // calculate uh2
      for (int j=0;j<=d.J;j++)
	{	  
	  uh2[j] = u2[j] * exp( -k/2 * (bb + gg*U2 + s(x[j])* ii * _e(d.eff,k,t) - kk) );
	}

      // reproduce uh
      Real uh_02 = xh[0] * _b(a1,a2,xh[0])*uh2[0];

      for (int j=0;j<d.J;j++)
	uh_02 += (_b(a1,a2,xh[j])*uh2[j] + _b(a1,a2,xh[j+1])*uh2[j+1]) * (xh[j+1]-xh[j]);

      uh_02 /= (2*kk*ww - xh[0]*_b(a1,a2,0));

      printf("%lf %lf\n",0.0,(uh_02 - uh_0)/delt);
      for (int j=0;j<10;j++)
	printf("%lf %lf\n",xh[j],(uh2[j]-uh[j])/delt);
      exit(1);
      
      */

      
      // integrate ph
      Real Ph = .5 * xh[0] * (ph_0 + ph[0]);
      for (int j=0;j<d.J;j++)
	Ph += .5 * (xh[j+1] - xh[j]) * (ph[j+1] + ph[j]);

      // calculate p      
      for (int j=d.J;j>0;j--)
	{

	  x[j] = x[j-1] + k * kk*(ww - xh[j-1]);
	  u[j] = u[j-1] * exp ( -k * (bb + gg*Uh + s(xh[j-1]) * ii * _e(d.eff,th) - kk));
	  //Real xxh = xh[j-1];
	  //Real uUh = Uh;
	  //Real zstr = (bb+gg*Uh+s(xh[j-1])*ii*_e(d.eff,th));
	  //Real bbb = k*gg*Ph*uh[j-1];
	  Real tmp1 = exp(-k*(bb+gg*Uh+s(xh[j-1])*ii*_e(d.eff,th)-kk));
	  Real tmp2 = exp((k/2)*(bb+gg*Uh+s(xh[j-1])*ii*_e(d.eff,th)-kk));
	  p[j] = p[j-1]*tmp1 - tmp1 * k*gg*Ph*uh[j-1] * tmp2;
	}
      
      x[0] = 0;

      // reproduce u
      u[0] = x[1] * _b(a1,a2,x[1])*u[1];

      for (int j=1;j<d.J;j++)
	u[0] += (_b(a1,a2,x[j])*u[j] + _b(a1,a2,x[j+1])*u[j+1]) * (x[j+1]-x[j]);

      u[0] /= (2*kk*ww - x[1]*_b(a1,a2,0));      
      
      // reproduce p
      p[0] = a1*x[1]*p[1]*x[1] + a2*x[1]*x[1]*p[1]*x[1] + x[1]*u[1]*x[1];

      for (int j=1;j<d.J;j++)
	p[0] += ( a1*x[j]*p[j] + a2*x[j]*x[j]*p[j] + x[j]*u[j] + a1*x[j+1]*p[j+1] + a2*x[j+1]*x[j+1]*p[j+1] + x[j+1]*u[j+1]) * (x[j+1]-x[j]);

      p[0] /= 2*kk*ww;

      for (int j=0;j<=d.J;j++)
	{
	  r[j] = 1e-3*w(x[j])*s(x[j])*u[j]*ii*_e(d.eff,k*i);
	  l[j] = 1e-3*w(x[j])*s(x[j])*p[j]*ii*_e(d.eff,k*i);
	}

      U = 0;
      C = 0;
      P = 0;
      L = 0;
  
      for (int j=0;j<d.J;j++)
	{
	  U += .5 * (x[j+1] - x[j]) * (u[j+1] + u[j]);
	  P += .5 * (x[j+1] - x[j]) * (p[j+1] + p[j]);
	  C += .5 * (x[j+1] - x[j]) * (r[j+1] + r[j]);
	  L += .5 * (x[j+1] - x[j]) * (l[j+1] + l[j]);
	}
  
      for (int j=0;j<=d.J;j++)
	o[j] = 2*(r[j] - _c(d.cat,k*i) * (d.p[i][j] + (1-d.Qp[i]) * r[j]/C))*(l[j] - _c(d.cat,k*i)*(1 - d.Qp[i]) * ( l[j] * C - r[j] * L ) / (C*C));

      for (int j=0;j<d.J;j++)
	ff += k * .5 * (x[j+1] - x[j]) * (o[j+1] + o[j]);

      //for (int j=0;j<=d.J;j++)
      //m[j] = pow(_c(d.cat,d.k,d.k*(i-d.N)) * (d.p[i][j] + (1-d.Qp[i]) * r[j]/C) - r[j],2.0);
      //for (int j=0;j<=d.J;j++)
      //ff2 += .5 * (x[j+1] - x[j]) * (m[j+1] + m[j]);
  
  }

  parameters->alpha1.gradient = ff;

  free(x);
  free(u);
  free(p);
  free(r);
  free(l);
  free(o);
  free(xh);
  free(uh);
  free(ph);
  
}
