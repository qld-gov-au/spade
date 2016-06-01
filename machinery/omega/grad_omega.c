#include <math.h>
#include "../../meschach/matrix.h"
#include "../../common.h"
#include "../objfns.h"
#include "../zstar.h"
#include "../Q.h"
#include "../reproduce.h"
#include "reproduce_omega.h"
#include "ini_omega.h"
#include "../../util/util.h"

void grad_omega(void* args)
{
  Grad_Args * grad_args = (Grad_Args *)args;

  VEC* theta = (*grad_args).theta;
  Data *d = (*grad_args).d;
  VEC *grad = (*grad_args).g;
  MAT *p = (*grad_args).p;
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
  double k = (*grad_args).k;
  int S = (*grad_args).S;
  
  // this function has not been updated for new birth function

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
 
  double a1 = theta->ve[0];
  double a2 = 0;
  double bb = theta->ve[1];
  double gg = theta->ve[2]*1e-7;
  double kk = theta->ve[3];
  double ww = theta->ve[4];
  double ii = theta->ve[5];

  VEC *Pi;
  Pi = v_get(x->m);

  get_row(x,0,xt);
  ini_omega(theta,xt,pt);
  set_row(p,0,pt);
  
  Pi->ve[0] = Q(get_row(x,0,xt),get_row(p,0,pt));

  if (PLOTSOLVP) 
    {

      FILE *p2 = fopen("plot.txt","w");

      for (int j=0;j<x->n;j++) 
	fprintf(p2,"%f %f\n",xt->ve[j],pt->ve[j]);

      fclose(p2);

      system("./plo1-nr > plotp_w_ini.pdf");

    }

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
      ph->ve[j] = pt->ve[j-1]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k)) - exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k))*(k/2)*( gg*Pi->ve[i-1]*ut->ve[j-1] + kk*(ut->ve[j]-ut->ve[j-1])/(xt->ve[j]-xt->ve[j-1]) );

      for (int j=2;j<=x->n;j++)
	ph->ve[j] = pt->ve[j-1]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k)) - exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k))*(k/2)*( gg*Pi->ve[i-1]*ut->ve[j-1] + kk*.5*( (ut->ve[j]-ut->ve[j-1])/(xt->ve[j]-xt->ve[j-1]) + (ut->ve[j-1]-ut->ve[j-2])/(xt->ve[j-1]-xt->ve[j-2]) ) );
	
      Q2_omega(a1,a2,kk,ww,xht,uht,ph);
      double Ph = Q(xht,ph);

      for (int j=1;j<x->n;j++)
	{
	  double b = k*( gg*Ph*uht->ve[j] + kk*.5*( (uht->ve[j]-uht->ve[j-1])/(xht->ve[j]-xht->ve[j-1]) + (uht->ve[j+1]-uht->ve[j])/(xht->ve[j+1]-xht->ve[j]) ) );
	  pn->ve[j] = pt->ve[j-1]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k)) - b*exp((k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j],Uhh->ve[i-1],k)-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k));
	}

      j= x->n;
      double b = k*gg*Ph*uht->ve[j];
      pn->ve[j] = pt->ve[j-1]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k)) - b*exp((k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j],Uhh->ve[i-1],k)-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k));

      get_row(un,i-1,unt);
      Q2_omega(a1,a2,kk,ww,xnt,unt,pn);
      Pi->ve[i] = Q(xnt,pn);

      idxremove(pn,pt,idxi->ive[i-1]); 
      set_row(p,i,pt);

    }

  if (PLOTSOLVP) 
    {

      FILE *p2 = fopen("plot.txt","w");

      int i=100;
      for (int j=0;j<x->n;j++)
	fprintf(p2,"%f %f\n",x->me[i][j],p->me[i][j]);

      fclose(p2);

      system("./plo1_nr > plotp_w_100.pdf");

      p2 = fopen("plot.txt","w");

      i=1000;
      for (int j=0;j<x->n;j++)
	fprintf(p2,"%f %f\n",x->me[i][j],p->me[i][j]);

      fclose(p2);

      system("./plo1_nr > plotp_w_1000.pdf");

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
  V_FREE(xhht);

}