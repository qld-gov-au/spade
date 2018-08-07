// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#include <math.h>
#include "../meschach/matrix.h"
#include "../common.h"
#include "../parameters.h"
#include "../model/biology/growth.h"
#include "../model/fishing/effort.h"
#include "../model/biology/birth.h"
#include "../model/fishing/selectivity.h"

void solve_clean(

		 Parameters * parameters,
		 VEC *eff,
		 Real k,		 
		 int S,
		 int Y,
		 Solve_Core_Args *core_args

		 )
{

  MAT *core_x = core_args->x;
  MAT *core_u = core_args->u; 
  
  Real aa = parameters->alpha.value;
  Real bb = parameters->beta.value;
  Real gg = parameters->gamma.value*1e-7;
  Real kk = parameters->kappa.value;
  Real ww = parameters->omega.value;
  Real ii = parameters->iota.value*1e-3;
  
  Real * restrict x = (Real *)calloc(core_x->n-1,sizeof(Real));
  Real * restrict u = (Real *)calloc(core_x->n-1,sizeof(Real));
    
  /* 
     initialize x
  */
  int J = core_x->n - core_x->m;

  // 'active': first J values in first 'row'
  for (int j=0;j<J;j++) {
    x[j] = h*j;
  }

  // 'neutral': all the rest (J+1.. I+J)
  for (int j=J;j<core_x->n;j++)
    x[j] = ww;

  /*
     ok, now initialize u
  */
  
  // prelims
  Real zeta = sqrt( 81*kk*kk*ww*ww*pow(aa*A1+2*aa*A2*ww,2.) - 12*kk*pow(aa*A1*ww+kk,3.) );
  Real eta = 9*aa*A1*kk*kk*ww + 18*aa*A2*kk*kk*ww*ww + kk*zeta;
  Real Z = pow(eta,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*kk*(aa*A1*ww+kk) / pow(eta,1./3);

  Real ubar = (Z - bb - kk) / gg; 
  Real vbar = (kk*ww*ubar) / (bb+gg*ubar+kk);
  Real wbar = (2*kk*ww*vbar) / (bb+gg*ubar+2*kk);  
  
  // set
  for (int j=0;j<core_x->n;j++) 
    u[j] = (aa*A1*vbar+aa*A2*wbar)*pow(ww-x[j],(bb+gg*ubar)/kk-1) / (kk*pow(ww,(bb+gg*ubar)/kk));

  for (int j=0;j<core_x->n;j++)     // load
    {
      core_x->me[0][j] = x[j];
      core_u->me[0][j] = u[j];
    }
  
  // integrate sequentially
  Real U = 0;
  for (int j=0;j<core_x->n;j++) 
    U += .5 * (x[j+1] - x[j]) * (u[j+1] + u[j]);
  
  Real * restrict xh = (Real *) calloc(core_x->n,sizeof(Real));  
  Real * restrict uh = (Real *) calloc(core_x->n,sizeof(Real));  

  for (int i=1;i<core_x->m;i++)
    {
  
      Real t = k*(i-S-1);
      Real th = k*(i-S-.5);
	 
      for (int j=0;j<core_x->n;j++)
	{	  
	  xh[j] = x[j] + k/2 * kk*(ww - x[j]);
	  uh[j] = u[j] * exp( -k/2 * (bb + gg*U + s(x[j])* ii * e(eff,k,t,Y) - kk) );
	}
      
      // only valid for models where fish are size zero at birth
      Real uh_0 = xh[0] * b(aa,xh[0])*uh[0];

      for (int j=0;j<core_x->n;j++)
	uh_0 += (b(aa,xh[j])*uh[j] + b(aa,xh[j+1])*uh[j+1]) * (xh[j+1]-xh[j]);

      uh_0 /= (2*kk*ww - xh[0]*b(aa,0));
            
      Real Uh = .5 * xh[0] * (uh_0 + uh[0]);
      for (int j=0;j<core_x->n-1;j++)
	Uh += .5 * (xh[j+1] - xh[j]) * (uh[j+1] + uh[j]);
      
      for (int j=core_x->n-1;j>0;j--) 
	{
	  x[j] = x[j-1] + k * kk*(ww-xh[j-1]);      
	  u[j] = u[j-1] * exp ( -k * (bb + gg*Uh + s(xh[j-1]) * ii * e(eff,k,th,Y) - kk));
	}      

      x[0] = 0;
      u[0] = x[1] * b(aa,x[1])*u[1];
      
      for (int j=1;j<core_x->n-1;j++)
	u[0] += (b(aa,x[j])*u[j] + b(aa,x[j+1])*u[j+1]) * (x[j+1]-x[j]);

      u[0] /= (2*kk*ww - x[1]*b(aa,0));

      for (int j=0;j<core_x->n;j++)     // load
	{
	  core_x->me[i][j] = x[j];
	  core_u->me[i][j] = u[j];
	}
      
      U = 0;
      for (int j=0;j<core_x->n-1;j++)
	U += .5 * (x[j+1] - x[j]) * (u[j+1] + u[j]);
      
    }
  
  free(x);
  free(u);
  free(xh);
  free(uh);

}
