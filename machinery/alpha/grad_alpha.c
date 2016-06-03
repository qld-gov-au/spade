#include <math.h>
#include "../../meschach/matrix.h"
#include "../../common.h"
#include "../objfns.h"
#include "../zstar.h"
#include "../Q.h"
#include "../reproduce.h"
#include "reproduce_alpha.h"
#include "ini_alpha.h"
#include "../../util/util.h"

void grad_alpha(void* args)
{
  Grad_Args * grad_args = (Grad_Args *)args;

  Data *d = (*grad_args).d;
  VEC *grad = (*grad_args).g;
  MAT *p = (*grad_args).p;
  MAT *x = ((*grad_args).core_args)->x;
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
  double k = (*grad_args).k;
  int S = (*grad_args).S;
  Parameters *parameters = (*grad_args).parameters;

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
 
  double aa = parameters->alpha.value;
  double bb = parameters->beta.value;
  double gg = parameters->gamma.value*1e-7;
  double kk = parameters->kappa.value;
  double ww = parameters->omega.value;
  double ii = parameters->iota.value*1e-3;

  VEC *Pi;
  Pi = v_get(x->m);

  get_row(x,0,xt);
  ini_alpha(parameters,xt,pt);
  set_row(p,0,pt);
  
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

      for (int j=1;j<=x->n;j++)
	ph->ve[j] = pt->ve[j-1]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k)) - exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k))*(k/2)*gg*Pi->ve[i-1]*ut->ve[j-1];
	 
      Q2_alpha(aa,kk,ww,xht,uht,ph);
      double Ph = Q(xht,ph);

      for (int j=1;j<=x->n;j++)
	{

	  double b = k*gg*Ph*uht->ve[j];
	  pn->ve[j] = pt->ve[j-1]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k)) - b*exp((k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j],Uhh->ve[i-1],k)-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k));

	}

      get_row(un,i-1,unt);
      Q2_alpha(aa,kk,ww,xnt,unt,pn);
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

  grad->ve[parameters->alpha.index] = G_ni(p, x, u, d, parameters->iota.value);

}
