#include <math.h>
#include "../meschach/matrix.h"
#include "../common.h"
#include "../model/biology/growth.h"
#include "zstar.h"
#include "Q.h"
#include "reproduce.h"
#include "../util/util.h"

VEC *initial(

	     Parameters *parameters,
	     VEC *x,
	     VEC *u

	     )

{

  x->ve[x->dim-1] -= 1e-5;

  double a = parameters->alpha.value;
  double b = parameters->beta.value;
  double g = parameters->gamma.value*1e-7;
  double k = parameters->kappa.value;
  double w = parameters->omega.value;

  double zeta = sqrt( 81*k*k*w*w*pow(a*A1+2*a*A2*w,2.) - 12*k*pow(a*A1*w+k,3.) );
  double eta = 9*a*A1*k*k*w + 18*a*A2*k*k*w*w + k*zeta;
  double Z = pow(eta,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*k*(a*A1*w+k) / pow(eta,1./3);

  double ubar = (Z - b - k) / g; 
  double vbar = (k*w*ubar) / (b+g*ubar+k);
  double wbar = (2*k*w*vbar) / (b+g*ubar+2*k);  

  for (int j=0;j<x->dim;j++) 
    u->ve[j] = (a*A1*vbar+a*A2*wbar)*pow(w-x->ve[j],(b+g*ubar)/k-1) / (k*pow(w,(b+g*ubar)/k));

  return u;
}

void solve(

     Parameters * parameters,
		 VEC *eff,
		 double k,		 
		 int S,
		 Solve_Core_Args *core_args  // modified

		 )
{

  VEC *xnt; VEC *unt; 
  VEC *xht; VEC *uht;
  VEC *xhht; VEC *uhh;
  VEC *xt; VEC *ut;

  MAT *x = core_args->x;
  MAT *u = core_args->u; 
  MAT *xhh = core_args->xhh;
  MAT *xh = core_args->xh;
  MAT *xn = core_args->xn;
  MAT *uh = core_args->uh;
  MAT *un = core_args->un;  
  VEC *Ui = core_args->Ui;
  VEC *Uh = core_args->Uh;
  VEC *Uhh = core_args->Uhh;  
  IVEC *idxi = core_args->idxi;

  xnt = v_get(x->n+1);  unt = v_get(x->n+1);
  xht = v_get(x->n+1);  uht = v_get(x->n+1);
  xhht = v_get(x->n+1); uhh = v_get(x->n+1);
  xt = v_get(x->n);     ut = v_get(x->n);

  for (int j=1;j<x->n;j++) 
    core_args->x->me[0][j] = h*j;

  set_row(core_args->u,0,initial(parameters,get_row(core_args->x,0,xt),ut));
 
  Ui->ve[0] = Q(xt,ut);

  double aa = parameters->alpha.value;
  double bb = parameters->beta.value;
  double gg = parameters->gamma.value*1e-7;
  double kk = parameters->kappa.value;
  double ww = parameters->omega.value;
  double ii = parameters->iota.value*1e-3;

  //  printf("\n");

  for (int i=1;i<x->m;i++)
    {

      double t = k*(i-S-1);
      double th = k*(i-S-.5);
      double thh = k*(i-S-.75);

      get_row(core_args->x,i-1,xt);
      get_row(core_args->u,i-1,ut);

      for (int j=1;j<=x->n;j++)
	{
	  xhht->ve[j] = xt->ve[j-1] + (k/4)*g(kk,ww,xt->ve[j-1]);
	  uhh->ve[j] = ut->ve[j-1]*exp(-(k/4)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],core_args->Ui->ve[i-1],k));
	}

      Q2(aa,kk,ww,xhht,uhh);
      core_args->Uhh->ve[i-1] = Q(xhht,uhh);
      set_row(core_args->xhh,i-1,xhht);
      
      for (int j=1;j<=x->n;j++)
	{
	  xht->ve[j] = xt->ve[j-1] + (k/2)*g(kk,ww,xhht->ve[j]);
	  uht->ve[j] = ut->ve[j-1]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j],core_args->Uhh->ve[i-1],k));
	}

      Q2(aa,kk,ww,xht,uht);
      core_args->Uh->ve[i-1] = Q(xht,uht);
      set_row(core_args->xh,i-1,xht);
      set_row(core_args->uh,i-1,uht);

      for (int j=1;j<=x->n;j++)
	{
	  xnt->ve[j] = xt->ve[j-1] + k*g(kk,ww,xht->ve[j]);
	  unt->ve[j] = ut->ve[j-1]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],core_args->Uh->ve[i-1],k));
	}

      Q2(aa,kk,ww,xnt,unt);
      core_args->Ui->ve[i] = Q(xnt,unt);
      set_row(xn,i-1,xnt);
      set_row(un,i-1,unt);
            
      int idx = idxselect(ww,xnt);
      core_args->idxi->ive[i-1] = idx;

      set_row(core_args->x,i,idxremove(xnt,xt,idx));
      set_row(core_args->u,i,idxremove(unt,ut,idx));

      //      printf("%f\n",Ui->ve[i]);

    }

  //exit(1);

  V_FREE(xt);
  V_FREE(xnt); 
  V_FREE(unt); 
  V_FREE(xht); 
  V_FREE(uht);
  V_FREE(xhht); 
  V_FREE(uhh);
  V_FREE(ut);

}
