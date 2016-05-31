﻿#include <math.h>
#include "../../../meschach/matrix.h"
#include "../../../common.h"
#include "ini_kappa.h"
#include "../zstar.h"
#include "../Q.h"
#include "../reproduce.h"
#include "reproduce_kappa.h"
#include "../../../util/util.h"
#include "../../objfns/objfns.h"

void solve_p_kappa(void* args)
{

  Solve_Args * solve_args = (Solve_Args *)args;

  VEC* theta = (*solve_args).theta;
  struct DATA *dataptr = (*solve_args).dataptr;
  VEC *grad = (*solve_args).grad;
  MAT *p = (*solve_args).p;
  MAT *x = (*solve_args).x;
  MAT *u = (*solve_args).u;
  MAT *xhh = (*solve_args).xhh;
  MAT *xh = (*solve_args).xh;
  MAT *xn = (*solve_args).xn;
  MAT *uh = (*solve_args).uh;
  MAT *un = (*solve_args).un;
  VEC *Ui = (*solve_args).Ui;
  VEC *Uh = (*solve_args).Uh;
  VEC *Uhh = (*solve_args).Uhh;
  IVEC *idxi = (*solve_args).idxi;
  VEC *eff = (*solve_args).eff;
  double k = (*solve_args).k;
  int S = (*solve_args).S;

  VEC *xt; VEC *xht; VEC *xnt; 
  VEC *ut; VEC *uht; VEC *pt; VEC *unt;
  xt = v_get(x->n);
  xht = v_get(x->n+1);
  xnt = v_get(x->n+1);
  ut = v_get(x->n);
  uht = v_get(x->n+1);
  unt = v_get(x->n+1);
  pt = v_get(x->n);

  VEC *xhht; 
  xhht = v_get(x->n+1);

  VEC *ph; VEC *pn;
  ph = v_get(x->n+1);
  pn = v_get(x->n+1);

  double aa = theta->ve[0];
  double bb = theta->ve[1];
  double gg = theta->ve[2]*1e-7;
  double kk = theta->ve[4];
  double ww = omega;
  double ii = theta->ve[3]*1e-3;

  VEC *Pi;
  Pi = v_get(x->m);

  VEC *thextra = v_get(5);

  thextra->ve[0] = theta->ve[0];
  thextra->ve[1] = theta->ve[1];
  thextra->ve[2] = theta->ve[2]*1e-7;
  thextra->ve[3] = theta->ve[4];
  thextra->ve[4] = omega;

  get_row(x,0,xt);
  ini_kappa(thextra,xt,pt);
  set_row(p,0,pt);
  
  V_FREE(thextra);

  Pi->ve[0] = Q(get_row(x,0,xt),get_row(p,0,pt));

  for (int i=1;i<x->m;i++)
    { 

      double t = k*(i-S-1);
      double th = k*(i-S-.5);
      double thh = k*(i-S-.75);

      get_row(x,i-1,xt);
      get_row(p,i-1,pt);
      get_row(xh,i-1,xht);
      get_row(xhh,i-1,xhht);
      get_row(uh,i-1,uht);
      get_row(xn,i-1,xnt);

      int j=1;
      ph->ve[j] = pt->ve[j-1]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k)) - exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k))*(k/2)*( (gg*Pi->ve[i-1]-1)*ut->ve[j-1] + (ww - xt->ve[j-1])*(ut->ve[j]-ut->ve[j-1])/(xt->ve[j]-xt->ve[j-1]) );

      for (int j=2;j<=x->n;j++)
	ph->ve[j] = pt->ve[j-1]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k)) - exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k))*(k/2)*( (gg*Pi->ve[i-1]-1)*ut->ve[j-1] + (ww - xt->ve[j-1])*.5*( (ut->ve[j]-ut->ve[j-1])/(xt->ve[j]-xt->ve[j-1]) + (ut->ve[j-1]-ut->ve[j-2])/(xt->ve[j-1]-xt->ve[j-2]) ) );
	
      Q2_kappa(aa,kk,ww,xht,uht,ph);
      double Ph = Q(xht,ph);

      for (int j=1;j<x->n;j++)
	{
	  double b = k*( (gg*Ph-1)*uht->ve[j] + (ww - xht->ve[j])*.5*( (uht->ve[j]-uht->ve[j-1])/(xht->ve[j]-xht->ve[j-1]) + (uht->ve[j+1]-uht->ve[j])/(xht->ve[j+1]-xht->ve[j]) ) );
	  pn->ve[j] = pt->ve[j-1]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k)) - b*exp((k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j],Uhh->ve[i-1],k)-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k));
	}

      j= x->n;
      double b = k*(gg*Ph-1)*uht->ve[j];
      pn->ve[j] = pt->ve[j-1]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k)) - b*exp((k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j],Uhh->ve[i-1],k)-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k));

      get_row(un,i-1,unt);
      Q2_kappa(aa,kk,ww,xnt,unt,pn);
      Pi->ve[i] = Q(xnt,pn);

      idxremove(pn,pt,idxi->ive[i-1]); 
      set_row(p,i,pt);

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

  grad->ve[4] = G_ni(p, x, u, dataptr, theta->ve[3]);

}
