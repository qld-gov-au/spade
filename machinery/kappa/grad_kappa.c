#include <math.h>
#include "../../meschach/matrix.h"
#include "../../common.h"
#include "../../parameters.h"
#include "../objfns.h"
#include "../zstar.h"
#include "../Q.h"
#include "../reproduce.h"
#include "reproduce_kappa.h"
#include "ini_kappa.h"
#include "../../util/util.h"

void grad_kappa(void* args)
{

  Grad_Args * grad_args = (Grad_Args *)args;
  Data *d = (*grad_args).d;
  MAT *x = (*grad_args).core_args->x;
  MAT *u = (*grad_args).core_args->u;
  MAT *xhh = (*grad_args).core_args->xhh;
  MAT *xh = (*grad_args).core_args->xh;
  MAT *xn = (*grad_args).core_args->xn;
  MAT *uh = (*grad_args).core_args->uh;
  MAT *un = (*grad_args).core_args->un;
  VEC *Ui = (*grad_args).core_args->Ui;
  VEC *Uh = (*grad_args).core_args->Uh;
  VEC *Uhh = (*grad_args).core_args->Uhh;
  IVEC *idxi = (*grad_args).core_args->idxi;
  VEC *eff = (*grad_args).eff;
  Real k = (*grad_args).k;
  int S = (*grad_args).S; 
  Parameters *parameters = (*grad_args).parameters;
  MAT *p = m_get(x->m,x->n);

  VEC *xt; VEC *xht; VEC *xnt;
  VEC *ut; VEC *uht; VEC *pt; VEC *unt;
  VEC *xhht; VEC *ph; VEC *pn;
  
  int J;
  if (BIGMATRICES)
    J = x->n - x->m;
  else
    J = x->n - 1;
  
  xt = v_get(J+1); ut = v_get(J+1); pt = v_get(J+1);

  if (BIGMATRICES)
    {        
      xnt = v_get(J+1);  unt = v_get(J+1);
      xht = v_get(J+1);  uht = v_get(J+1);
      xhht = v_get(J+1); ph = v_get(J+1);
      pn = v_get(J+1);
    }
  else
    {        
      xnt = v_get(J+2);  unt = v_get(J+2);
      xht = v_get(J+2);  uht = v_get(J+2);
      xhht = v_get(J+2); ph = v_get(J+2);
      pn = v_get(J+2);
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

  if (BIGMATRICES) 
    {
      xt->dim = J+1;
      pt->dim = J+1;
    }

  ini_kappa(parameters,xt,pt);
  set_row(p,0,pt);

  Pi->ve[0] = Q(get_row(x,0,xt),get_row(p,0,pt));


  if (BIGMATRICES)
    pt = v_resize(pt,p->n);


  for (int i=1;i<x->m;i++)
    { 

      Real t = k*(i-S-1);
      Real th = k*(i-S-.5);
      Real thh = k*(i-S-.75);

      get_row(x,i-1,xt);
      get_row(p,i-1,pt);
      get_row(xh,i-1,xht);
      get_row(xhh,i-1,xhht);
      get_row(uh,i-1,uht);
      get_row(xn,i-1,xnt);
      get_row(u,i-1,ut);

      int terminator;
      if(BIGMATRICES) 
        {
          terminator = J+i-1;
          xt = v_resize(xt,terminator+1);
          ut = v_resize(ut,terminator+1);
          pt = v_resize(pt,terminator+1);        
          xhht = v_resize(xhht,terminator+2);
          xht = v_resize(xht,terminator+2);
          uht = v_resize(uht,terminator+2);
          xnt = v_resize(xnt,terminator+2);
          unt = v_resize(unt,terminator+2);
          ph = v_resize(ph,terminator+2);
          pn = v_resize(pn,terminator+2);
        } 
      else 
        {
          terminator = J;
        }

      int j=0;
      ph->ve[j+1] = pt->ve[j]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j],Ui->ve[i-1],k,d->Y)) - exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j],Ui->ve[i-1],k,d->Y))*(k/2)*( (gg*Pi->ve[i-1]-1)*ut->ve[j] + (ww - xt->ve[j])*(ut->ve[j+1]-ut->ve[j])/(xt->ve[j+1]-xt->ve[j]) );

      for (int j=1;j<=terminator - 1;j++) {

        Real z_star_result = zstar(eff,bb,gg,kk,ii,t,xt->ve[j],Ui->ve[i-1],k,d->Y);
        Real z_star_result_2 = zstar(eff,bb,gg,kk,ii,t,xt->ve[j],Ui->ve[i-1],k,d->Y);
        Real expression1 = pt->ve[j];
        Real expression2 = exp(-(k/2)*z_star_result);
        Real expression3 = exp(-(k/2)*z_star_result_2);
        Real expression4 = (gg*Pi->ve[i-1]-1);
        Real expression5 = ut->ve[j];
        Real expression6 = (ww - xt->ve[j]);

        Real expression7 = (ut->ve[j + 1]-ut->ve[j]); // todo: Review this
        Real expression8 = (xt->ve[j + 1]-xt->ve[j]);
        Real expression9 = (ut->ve[j]-ut->ve[j-1]);
        Real expression10 = (xt->ve[j]-xt->ve[j-1]) ;


	      ph->ve[j + 1] = expression1*expression2 - expression3*(k/2)*(
          expression4 *
          expression5 +
          expression6 *
          .5*
          (
            expression7/
            expression8 +
            expression9/
            expression10
          )
        );
      }
      j=terminator;
      ph->ve[j+1] = pt->ve[j]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j],Ui->ve[i-1],k,d->Y)) - exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j],Ui->ve[i-1],k,d->Y))*(k/2)*( (gg*Pi->ve[i-1]-1)*ut->ve[j] );

      Q2_kappa(aa,kk,ww,xht,uht,ph);

      Real Ph = Q(xht,ph);

      for (int j=0;j<=terminator - 1;j++)
	     {
          Real b = k*( (gg*Ph-1)*uht->ve[j+1] + (ww - xht->ve[j+1])*.5*( (uht->ve[j+1]-uht->ve[j])/(xht->ve[j+1]-xht->ve[j]) + (uht->ve[j+2]-uht->ve[j+1])/(xht->ve[j+2]-xht->ve[j+1]) ) );
          pn->ve[j+1] = pt->ve[j]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j+1],Uh->ve[i-1],k,d->Y)) - b*exp((k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j+1],Uhh->ve[i-1],k,d->Y)-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j+1],Uh->ve[i-1],k,d->Y));
   	     }
 
      j = terminator;
      Real b = k*(gg*Ph-1)*uht->ve[j+1];
      pn->ve[j+1] = pt->ve[j]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j+1],Uh->ve[i-1],k,d->Y)) - b*exp((k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j+1],Uhh->ve[i-1],k,d->Y)-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j+1],Uh->ve[i-1],k,d->Y));

      get_row(un,i-1,unt);
      Q2_kappa(aa,kk,ww,xnt,unt,pn);
      Pi->ve[i] = Q(xnt,pn);

      if(BIGMATRICES) 
        {
          pn = v_resize(pn,p->n);
          set_row(p,i,pn);
        } 
      else 
        {
          idxremove(pn,pt,idxi->ive[i-1]); 
          set_row(p,i,pt);
        }

    }

  V_FREE(xt);
  V_FREE(xht);
  V_FREE(xnt); 
  V_FREE(ut); 
  V_FREE(uht);
  V_FREE(pt);
  V_FREE(unt);
  V_FREE(ph);
  V_FREE(pn);
  V_FREE(Pi);
  V_FREE(xhht);

  parameters->kappa.gradient = G_ni(p, x, u, d, parameters->iota.value);

  M_FREE(p);

}
