#include <math.h>
#include "../../meschach/matrix.h"
#include "../../common.h"
#include "../objfns.h"
#include "../zstar.h"
#include "../Q.h"
#include "../reproduce.h"
#include "ini_beta.h"
#include "../../util/util.h"

void grad_beta(void* args)
{
  Grad_Args * grad_args = (Grad_Args *)args;

  Data *d = (*grad_args).d;
  VEC *grad = (*grad_args).g;
  MAT *p = (*grad_args).p;
  MAT *x = (*grad_args).core_args->x;
  MAT *u = (*grad_args).core_args->u;
  MAT *xhh = (*grad_args).core_args->xhh;
  MAT *xh = (*grad_args).core_args->xh;
  MAT *xn = (*grad_args).core_args->xn;
  MAT *uh = (*grad_args).core_args->uh;
  VEC *Ui = (*grad_args).core_args->Ui;
  VEC *Uh = (*grad_args).core_args->Uh;
  VEC *Uhh = (*grad_args).core_args->Uhh;
  IVEC *idxi = (*grad_args).core_args->idxi;
  VEC *eff = (*grad_args).eff;
  double k = (*grad_args).k;
  int S = (*grad_args).S; 
  Parameters *parameters = (*grad_args).parameters;

  VEC *xt; VEC *xht; VEC *xnt;
  VEC *ut; VEC *uht; VEC *pt;
  xt = v_get(x->n);
  xht = v_get(x->n+1);
  xnt = v_get(x->n+1);
  ut = v_get(x->n);
  uht = v_get(x->n+1);
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
  ini_beta(parameters,xt,pt);
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
	ph->ve[j] = pt->ve[j-1]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k)) - exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k))*(k/2)*(1+gg*Pi->ve[i-1])*ut->ve[j-1];
	 
      Q2(aa,kk,ww,xht,ph);
      double Ph = Q(xht,ph);

      for (int j=1;j<=x->n;j++)
	{
	  double b = k*(1+gg*Ph)*uht->ve[j];
	  pn->ve[j] = pt->ve[j-1]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k)) - b*exp((k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j],Uhh->ve[i-1],k)-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k));
	}

      Q2(aa,kk,ww,xnt,pn);
      Pi->ve[i] = Q(xnt,pn);

      idxremove(pn,pt,idxi->ive[i-1]); 
      set_row(p,i,pt);

    }
   
  /*  printf("\n");
  for (int j=0;j<300;j++)
    {
      printf("%f %f\n",x->me[x->m-1][j],p->me[x->m-1][j]); //s(x->me[tmi.I+2][j])*w(x->me[tmi.I+2][j])*tmi.dp.iota*e(2*k)*); //-tmi.I+4][i]);
      //printf("%f %f\n",x->me[x->m-1][i],s(x->me[x->m-1][i])*w(x->me[x->m-1][i])*e(k*(x->m-9-(tmi.I+1)))*u->me[x->m-1][i] + s(x->me[x->m-1][i])*w(x->me[x->m-1][i])*tmi.dp.iota*e(k*(x->m-9-(tmi.I+1)))*p->me[x->m-1][i]);
    }
  printf("e\n");
  exit(1);
  */

  /*
  VEC *xxx = v_get(x->m);
  for (int j=0;j<x->m;j++)
    xxx->ve[j] = x->me[x->m-1][j];

  VEC *yyy = v_get(x->m);
  for (int j=0;j<x->m;j++)
    yyy->ve[j] = s(x->me[x->m-59][j])*w(x->me[x->m-59][j])*e(k*(x->m-59-(tmi.I+1)))*u->me[x->m-59][j] + s(x->me[x->m-59][j])*w(x->me[x->m-59][j])*iota*e(k*(x->m-59-(tmi.I+1)))*p->me[x->m-59][j];

  printf("%f\n",Q(xxx,yyy));

  exit(1);  
  */

  V_FREE(xt); 
  V_FREE(xht);
  V_FREE(xnt); 
  V_FREE(ut); 
  V_FREE(uht);
  V_FREE(pt);
  V_FREE(ph);
  V_FREE(pn);
  V_FREE(Pi);
  V_FREE(xhht);
  
  grad->ve[1] = G_ni(p, x, u, d, parameters->iota.value);
}
