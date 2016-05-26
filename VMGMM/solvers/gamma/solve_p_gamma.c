#include <math.h>
#include "../../../meschach/matrix.h"
#include "../../../common.h"
#include "../../objfns/objfns.h"
#include "../zstar.h"
#include "../Q.h"
#include "../reproduce.h"
#include "ini_gamma.h"
#include "../../../util/util.h"

void solve_p_gamma(void* args)
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
  //MAT *un = (*solve_args).un;
  VEC *Ui = (*solve_args).Ui;
  VEC *Uh = (*solve_args).Uh;
  VEC *Uhh = (*solve_args).Uhh;
  IVEC *idxi = (*solve_args).idxi;
  VEC *eff = (*solve_args).eff;
  double k = (*solve_args).k;
  int S = (*solve_args).S;

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

  double aa = theta->ve[0];
  double bb = theta->ve[1];
  double gg = theta->ve[2]*1e-7;
  double kk = kappa;
  double ww = omega;
  double ii = theta->ve[3]*1e-3;
 
  VEC *Pi;
  Pi = v_get(x->m);

  VEC * thextra = v_get(5);

  thextra->ve[0] = theta->ve[0];
  thextra->ve[1] = theta->ve[1];
  thextra->ve[2] = theta->ve[2];
  thextra->ve[3] = kappa;
  thextra->ve[4] = omega;

  get_row(x,0,xt);
  ini_gamma(thextra,xt,pt);
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

      for (int j=1;j<=x->n;j++)
	ph->ve[j] = pt->ve[j-1]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k)) - exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k))*(k/2)*(1e-7*Ui->ve[i-1]+gg*Pi->ve[i-1])*ut->ve[j-1];
	 
      Q2(aa,kk,ww,xht,ph);
      double Ph = Q(xht,ph);

      for (int j=1;j<=x->n;j++)
	{
	  double b = k*(1e-7*Uh->ve[i-1]+gg*Ph)*uht->ve[j];
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
  
  grad->ve[2] = G_ni(p, x, u, dataptr, theta->ve[3]);

}