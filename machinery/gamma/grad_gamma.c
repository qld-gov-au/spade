#include <math.h>
#include "../../meschach/matrix.h"
#include "../../common.h"
#include "../../parameters.h"
#include "../objfns.h"
#include "../zstar.h"
#include "../Q.h"
#include "../reproduce.h"
#include "ini_gamma.h"
#include "../../util/util.h"

void grad_gamma(void* args)
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
  if (SGNM)
    J = x->n - 1;
  else
    J = x->n - x->m;
  
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

  ini_gamma(parameters,xt,pt);
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
	  ph->ve[j+1] = pt->ve[j]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j],Ui->ve[i-1],k,d->Y)) - exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j],Ui->ve[i-1],k,d->Y))*(k/2)*(1e-7*Ui->ve[i-1]+gg*Pi->ve[i-1])*ut->ve[j];
	}
      
      Q2(aa,kk,ww,xht,ph);
      Real Ph = Q(xht,ph);

      if (QUARTER)
	for (int j=0;j<=terminator;j++)
	  {
	    Real b = k*(1e-7*Uh->ve[i-1]+gg*Ph)*uht->ve[j+1];
	    pn->ve[j+1] = pt->ve[j]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j+1],Uh->ve[i-1],k,d->Y)) - b*exp((k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j+1],Uhh->ve[i-1],k,d->Y)-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j+1],Uh->ve[i-1],k,d->Y));
	  }
      else
	for (int j=0;j<=terminator;j++)
	  {
	    Real b = k*(1e-7*Uh->ve[i-1]+gg*Ph)*uht->ve[j+1];
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
  
  parameters->gamma.gradient = G_ni(p, x, u, d, parameters->iota.value);

  M_FREE(p);

}
