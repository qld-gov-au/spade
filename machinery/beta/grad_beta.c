// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#include <math.h>
#include "../../meschach/matrix.h"
#include "../../common.h"
#include "../../parameters.h"
#include "../objfns.h"
#include "../zstar.h"
#include "../Q.h"
#include "../reproduce.h"
#include "ini_beta.h"
#include "../../util/util.h"
#include "../../model/fishing/selectivity.h"
#include "../../model/fishing/effort.h"
#include "../../model/biology/birth.h"

void grad_beta_clean(void* args)
{

  // only valid for models where fish are size zero at birth
  
  Grad_Args_No_MESCHACH * grad_args = (Grad_Args_No_MESCHACH *)args;

  Data *d = (*grad_args).d;

  int I,J;
  I = d->I;
  J = d->J;

  MAT *core_x = m_get(I+1,J+1);
  MAT *core_u = m_get(I+1,J+1);
  MAT *core_p = m_get(I+1,J+1);
  
  VEC *eff = (*grad_args).eff;
  Real k = (*grad_args).k;
  int S = (*grad_args).S;
  Parameters *parameters = (*grad_args).parameters;
     
  Real aa = parameters->alpha.value;
  Real bb = parameters->beta.value;
  Real gg = parameters->gamma.value*1e-7;
  Real kk = parameters->kappa.value;
  Real ww = parameters->omega.value;
  Real ii = parameters->iota.value*1e-3;
  
  Real * restrict x = (Real *)calloc(core_x->n-1,sizeof(Real));
  Real * restrict u = (Real *)calloc(core_x->n-1,sizeof(Real));
  Real * restrict p = (Real *)calloc(core_x->n-1,sizeof(Real));
    
  /* 
     initialize x
  */
  J = core_x->n - core_x->m;

  // 'active': first J values in first 'row'
  for (int j=0;j<J;j++) 
    x[j] = h*j;

  // 'neutral': all the rest (J+1.. I+J)
  for (int j=J;j<core_x->n;j++)
    x[j] = ww-1e-9;

  /*
     ok, now initialize u and p
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

  Real Za = ( kk*kk*ww*(3*A1*zeta+6*A2*ww*zeta-6*A1*pow(kk+aa*A1*ww,2.)+27*kk*ww*(A1+2*A2*ww)*(aa*A1+2*aa*A2*ww)) / ( zeta*pow(eta,2/3.) ) ) * ( (1/ (pow(2,1/3.) * pow(3,2/3.) )) - ( pow(2/3.,1/3.)*kk*(kk+aa*A1*ww) / pow(eta,2/3.) ) ) + ( pow(2/3.,1/3.)*A1*kk*ww / pow(eta,1/3.) );

  for (int j=0;j<core_x->n;j++)
    p[j] = -pow(1-x[j]/ww,Z/kk-2) * (1/(gg*Z)) * (aa*A1 + 2*aa*A2*kk*ww/(Z+kk));
  
  for (int j=0;j<core_x->n;j++)     // load
    {
      core_x->me[0][j] = x[j];
      core_u->me[0][j] = u[j];
      core_p->me[0][j] = p[j];
    }
  
  // integrate sequentially
  Real U = 0;
  for (int j=0;j<core_x->n;j++) 
    U += .5 * (x[j+1] - x[j]) * (u[j+1] + u[j]);

  Real P = 0;
  for (int j=0;j<core_x->n;j++) 
    P += .5 * (x[j+1] - x[j]) * (p[j+1] + p[j]);
  
  Real * restrict xh = (Real *) calloc(core_x->n,sizeof(Real));  
  Real * restrict uh = (Real *) calloc(core_x->n,sizeof(Real));  
  Real * restrict ph = (Real *) calloc(core_x->n,sizeof(Real));  

  for (int i=1;i<core_x->m;i++)
    {
  
      Real t = k*(i-1);
      Real th = k*(i-.5);

      // calculate uh
      for (int j=0;j<core_x->n;j++)
	{	  
	  xh[j] = x[j] + k/2 * kk*(ww - x[j]);
	  uh[j] = u[j] * exp( -k/2 * (bb + gg*U + s(x[j])* ii * e(eff,k,t,d->Y) - kk) );
	}

      // reproduce uh
      Real uh_0 = xh[0] * b(aa,xh[0])*uh[0];

      for (int j=0;j<core_x->n-1;j++)
	uh_0 += (b(aa,xh[j])*uh[j] + b(aa,xh[j+1])*uh[j+1]) * (xh[j+1]-xh[j]);

      uh_0 /= (2*kk*ww - xh[0]*b(aa,0));

      // integrate uh
      Real Uh = .5 * xh[0] * (uh_0 + uh[0]);
      for (int j=0;j<core_x->n-1;j++)
	Uh += .5 * (xh[j+1] - xh[j]) * (uh[j+1] + uh[j]);

      // calculate ph
      for (int j=0;j<core_x->n;j++)
	ph[j] = p[j]*exp(-(k/2)*(bb + gg*U + s(x[j]) * ii * e(eff,k,t,d->Y) - kk)) - exp(-(k/2)*(bb + gg*U + s(x[j]) * ii * e(eff,k,t,d->Y) - kk))*(k/2)*(1+gg*P)*u[j];

      // reproduce ph
      Real ph_0 = xh[0] * b(aa,xh[0])*ph[0]; 

      for (int j=0;j<core_x->n-1;j++)
	ph_0 += (b(aa,xh[j])*ph[j] + b(aa,xh[j+1])*ph[j+1]) * (xh[j+1]-xh[j]);

      ph_0 /= (2*kk*ww - xh[0]*b(aa,0));
                  
      // integrate ph
      Real Ph = .5 * xh[0] * (ph_0 + ph[0]);
      for (int j=0;j<core_x->n-1;j++)
	Ph += .5 * (xh[j+1] - xh[j]) * (ph[j+1] + ph[j]);

      // calculate p      
      for (int j=core_x->n-1;j>0;j--)
	{

	  x[j] = x[j-1] + k * kk*(ww - xh[j-1]);
	  u[j] = u[j-1] * exp ( -k * (bb + gg*Uh + s(xh[j]) * ii * e(eff,k,th,d->Y) - kk));		  
	  
	  Real tmp1 = exp(-k*(bb+gg*Uh+s(xh[j-1])*ii*e(eff,k,th,d->Y)-kk));
	  Real tmp2 = exp((k/2)*(bb+gg*Uh+s(xh[j-1])*ii*e(eff,k,th,d->Y)-kk));
	  p[j] = p[j-1]*tmp1 - tmp1 * k*(1+gg*Ph)*uh[j-1] * tmp2;
	}
      
      x[0] = 0;

      // reproduce u
      u[0] = x[1] * b(aa,x[1])*u[1];

      for (int j=1;j<core_x->n-1;j++)
	u[0] += (b(aa,x[j])*u[j] + b(aa,x[j+1])*u[j+1]) * (x[j+1]-x[j]);

      u[0] /= (2*kk*ww - x[1]*b(aa,0));      
      
      // reproduce p
      p[0] = x[1] * b(aa,x[1])*p[1];

      for (int j=1;j<core_x->n-1;j++)
	p[0] += (b(aa,x[j])*p[j] + b(aa,x[j+1])*p[j+1]) * (x[j+1]-x[j]);

      p[0] /= (2*kk*ww - x[0]*b(aa,0));
      
      // integrate u
      U = .5 * x[0] * (u[0] + u[1]);
      for (int j=1;j<core_x->n-1;j++)
	U += .5 * (x[j+1] - x[j]) * (u[j+1] + u[j]);
      
      // integrate p
      P = .5 * x[1] * (p[0] + p[1]);
      for (int j=1;j<core_x->n-1;j++)
	P += .5 * (x[j+1] - x[j]) * (p[j+1] + p[j]);
      
      for (int j=0;j<core_x->n;j++)     // load
	{
	  core_x->me[i][j] = x[j];
	  core_u->me[i][j] = u[j];
	  core_p->me[i][j] = p[j];
	}

    }
  
  parameters->alpha.gradient = G_ni(core_p, core_x, core_u, d, parameters->iota.value);

  free(x);
  free(u);
  free(p);
  free(xh);
  free(uh);
  free(ph);
  
  M_FREE(core_p);
  M_FREE(core_x);
  M_FREE(core_u);

}

void grad_beta(void* args)
{
  Grad_Args * grad_args = (Grad_Args *)args;

  Data *d = (*grad_args).d;
  MAT *x = (*grad_args).core_args->x;
  MAT *u = (*grad_args).core_args->u;
  MAT *xhh;
  if (QUARTER)
    xhh = (*grad_args).core_args->xhh;
  MAT *xh = (*grad_args).core_args->xh;
  MAT *xn = (*grad_args).core_args->xn;
  MAT *uh = (*grad_args).core_args->uh;
  VEC *Ui = (*grad_args).core_args->Ui;
  VEC *Uh = (*grad_args).core_args->Uh;
  VEC *Uhh;
  if (QUARTER)
    Uhh = (*grad_args).core_args->Uhh;
  IVEC *idxi = (*grad_args).core_args->idxi;
  VEC *eff = (*grad_args).eff;
  Real k = (*grad_args).k;
  int S = (*grad_args).S; 
  Parameters *parameters = (*grad_args).parameters;
  MAT *p = m_get(x->m,x->n);

  VEC *xt; VEC *ut; VEC *pt;
  VEC *xhht; VEC *xht; VEC *xnt; 
  VEC *uht; VEC *ph; VEC *pn;

  int J;
  if (!SGNM)
    J = x->n - x->m;
  else
    J = x->n - 1;
  
  xt = v_get(J+1); ut = v_get(J+1); pt = v_get(J+1);

  if (SGNM)
    {        
      xhht = v_get(J+2); xht = v_get(J+2); xnt = v_get(J+2);
      uht = v_get(J+2); ph = v_get(J+2); pn = v_get(J+2);
    }
  else
    {        
      if (QUARTER)
	xhht = v_get(J+1);
      xht = v_get(J+1); xnt = v_get(J+1);
      uht = v_get(J+1); ph = v_get(J+1); pn = v_get(J+1);
    }

  Real aa = parameters->alpha.value;
  Real bb = parameters->beta.value;
  Real gg = parameters->gamma.value*1e-7;
  Real kk = parameters->kappa.value;
  Real ww = parameters->omega.value;
  Real ii = parameters->iota.value*1e-3;

  VEC *Pi;
  Pi = v_get(x->m);

  get_row(x,0,xt);

  if (!SGNM) 
    xt = v_resize(xt,J+1);

  ini_beta(parameters,xt,pt);
  set_row(p,0,pt);

  Pi->ve[0] = Q(xt,pt);
  
  for (int i=1;i<x->m;i++)
    { 

      Real t = k*(i-S-1);
      Real th = k*(i-S-.5);
      Real thh = k*(i-S-.75);

      get_row(x,i-1,xt);
      get_row(p,i-1,pt);
      get_row(xh,i-1,xht);
      if (QUARTER)
	get_row(xhh,i-1,xhht);
      get_row(uh,i-1,uht);
      get_row(u,i-1,ut);

      if (SGNM)
	{
	  get_row(xn,i-1,xnt);
	}
      else
	{
	  get_row(x,i,xnt);
	}	  
      
      int terminator;
      if(SGNM) 
        {
          terminator = J;
        } 
      else 
        {
          terminator = J+i-1;
          xt = v_resize(xt,terminator+1);
          ut = v_resize(ut,terminator+1);
          pt = v_resize(pt,terminator+1);
	  if (QUARTER)
	    xhht = v_resize(xhht,terminator+2);
          xht = v_resize(xht,terminator+2);
          uht = v_resize(uht,terminator+2);
          xnt = v_resize(xnt,terminator+2);
          ph = v_resize(ph,terminator+2);
          pn = v_resize(pn,terminator+2);
	}

      for (int j=0;j<=terminator;j++)
        {
          ph->ve[j+1] = pt->ve[j]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j],Ui->ve[i-1],k,d->Y)) - exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j],Ui->ve[i-1],k,d->Y))*(k/2)*(1+gg*Pi->ve[i-1])*ut->ve[j];
	}
   
      Q2(aa,kk,ww,xht,ph);
      Real Ph = Q(xht,ph);

      if (QUARTER)
	for (int j=0;j<=terminator;j++)
	  {
	    Real b = k*(1+gg*Ph)*uht->ve[j+1];
	    pn->ve[j+1] = pt->ve[j]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j+1],Uh->ve[i-1],k,d->Y)) - b*exp((k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j+1],Uhh->ve[i-1],k,d->Y)-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j+1],Uh->ve[i-1],k,d->Y));
	  }
      else
	for (int j=0;j<=terminator;j++)
	  {
	    Real b = k*(1+gg*Ph)*uht->ve[j+1];
	    pn->ve[j+1] = pt->ve[j]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j+1],Uh->ve[i-1],k,d->Y)) - b*exp((k/2)*zstar(eff,bb,gg,kk,ii,th,xht->ve[j+1],Uh->ve[i-1],k,d->Y)-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j+1],Uh->ve[i-1],k,d->Y));
	  }
	
      Q2(aa,kk,ww,xnt,pn);
      Pi->ve[i] = Q(xnt,pn);

      if(SGNM) 
        {
          idxremove(pn,pt,idxi->ive[i-1]); 
          set_row(p,i,pt);
        } 
      else 
        {
          pn = v_resize(pn,p->n);
          set_row(p,i,pn);
        }

    }
   
  V_FREE(xt); 
  V_FREE(xht);
  V_FREE(xnt); 
  V_FREE(ut); 
  V_FREE(uht);
  V_FREE(pt);
  V_FREE(ph);
  V_FREE(pn);
  V_FREE(Pi);
  if (QUARTER)
    V_FREE(xhht);
  
  parameters->beta.gradient = G_ni(p, x, u, d, parameters->iota.value);
  
  M_FREE(p);
  
}
