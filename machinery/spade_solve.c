// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#include <math.h>
#include "../meschach/matrix.h"
#include "../common.h"
#include "../parameters.h"
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

  Real a = parameters->alpha.value;
  Real b = parameters->beta.value;
  Real g = parameters->gamma.value*1e-7;
  Real k = parameters->kappa.value;
  Real w = parameters->omega.value;

  Real zeta = sqrt( 81*k*k*w*w*pow(a*A1+2*a*A2*w,2.) - 12*k*pow(a*A1*w+k,3.) );
  Real eta = 9*a*A1*k*k*w + 18*a*A2*k*k*w*w + k*zeta;
  Real Z = pow(eta,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*k*(a*A1*w+k) / pow(eta,1./3);

  Real ubar = (Z - b - k) / g; 
  Real vbar = (k*w*ubar) / (b+g*ubar+k);
  Real wbar = (2*k*w*vbar) / (b+g*ubar+2*k);  

  for (int j=0;j<x->dim;j++) 
    u->ve[j] = (a*A1*vbar+a*A2*wbar)*pow(w-x->ve[j],(b+g*ubar)/k-1) / (k*pow(w,(b+g*ubar)/k));

  return u;
}

void solve(

     Parameters * parameters,
		 VEC *eff,
		 Real k,		 
		 int S,
     int Y,
		 Solve_Core_Args *core_args  // modified

		 )
{

  VEC *xnt; VEC *unt; 
  VEC *xht; VEC *uht;
  VEC *xhht; VEC *uhh;
  VEC *xt; VEC *ut;

  MAT *x = core_args->x;
  MAT *u = core_args->u; 
  MAT *xhh;
  if (QUARTER)
    xhh = core_args->xhh;
  MAT *xh = core_args->xh;
  MAT *xn = core_args->xn;
  MAT *uh = core_args->uh;
  MAT *un = core_args->un;  
  VEC *Ui = core_args->Ui;
  VEC *Uh = core_args->Uh;
  VEC *Uhh = core_args->Uhh;  
  IVEC *idxi = core_args->idxi;

  int J; 
  int bigJ;
  if (!SGNM){
    bigJ = x->n - 1;
    J = x->n - x->m;

  } else {
    J = x->n - 1;
   
  }

  if (!SGNM)
    {        
    xt = v_get(bigJ+1); ut = v_get(bigJ+1);
      xnt = v_get(bigJ+1);  unt = v_get(bigJ+1);
      xht = v_get(bigJ+1);  uht = v_get(bigJ+1);
      if (QUARTER)	{
	xhht = v_get(bigJ+1);
	uhh = v_get(bigJ+1);
      }
    }
  else
    {       
      xnt = v_get(J+2);  unt = v_get(J+2);
      xht = v_get(J+2);  uht = v_get(J+2);
      if (QUARTER) {
	xhht = v_get(J+2);
	uhh = v_get(J+2);
      }
      xt = v_get(J+1); ut = v_get(J+1); 
    }
                           
  for (int j=0;j<=J;j++) 
    core_args->x->me[0][j] = h*j;


  get_row(core_args->x,0,xt);

  if (!SGNM) 
    {
      v_resize(xt,J+1);
      v_resize(ut,J+1);
    }


  initial(parameters,xt,ut);

  set_row(core_args->u,0,ut);

  Ui->ve[0] = Q(xt,ut);

  if (!SGNM)
    {
      v_resize(ut,u->n);
    }
      
 

  Real aa = parameters->alpha.value;
  Real bb = parameters->beta.value;
  Real gg = parameters->gamma.value*1e-7;
  Real kk = parameters->kappa.value;
  Real ww = parameters->omega.value;
  Real ii = parameters->iota.value*1e-3;

  //printf("\n");
  //for (int i=0;i<eff->dim;i++)
  //  printf("%lf\n",eff->ve[i]);
  //exit(1);
  
  for (int i=1;i<x->m;i++)
    {

      Real t = k*(i-S-1);
      Real th = k*(i-S-.5);
      Real thh;

      if (QUARTER)
	thh = k*(i-S-.75);
      
      int terminator;
      if(!SGNM) 
        {
          terminator = J+i-1;
          xt = v_resize(xt,terminator+1);
          ut = v_resize(ut,terminator+1);
	  if (QUARTER){
	    xhht = v_resize(xhht,terminator+2);
	    uhh = v_resize(uhh,terminator+2);
	  }
          xht = v_resize(xht,terminator+2);
          uht = v_resize(uht,terminator+2);
          xnt = v_resize(xnt,terminator+2);
          unt = v_resize(unt,terminator+2);
        } 
      else 
        {
          terminator = J;
        }

      get_row(x,i-1,xt);
      get_row(u,i-1,ut);

      if (QUARTER)
	{
	
	  for (int j=0;j<=terminator;j++)
	    {
	      xhht->ve[j+1] = xt->ve[j] + (k/4)*g(kk,ww,xt->ve[j]);
	      uhh->ve[j+1] = ut->ve[j]*exp(-(k/4)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j],Ui->ve[i-1],k,Y));
	    }      

	  Q2(aa,kk,ww,xhht,uhh);
	  Uhh->ve[i-1] = Q(xhht,uhh);
            
	  set_row(xhh,i-1,xhht);
	      
	  for (int j=0;j<=terminator;j++)
	    {
	      xht->ve[j+1] = xt->ve[j] + (k/2)*g(kk,ww,xhht->ve[j+1]);
	      uht->ve[j+1] = ut->ve[j]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j+1],Uhh->ve[i-1],k,Y));
	    }
	}

      else

	{

	  for (int j=0;j<=terminator;j++)
	    {
	      xht->ve[j+1] = xt->ve[j] + (k/2)*g(kk,ww,xt->ve[j]);
	      uht->ve[j+1] = ut->ve[j]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j],Ui->ve[i-1],k,Y));
	    }
	  
	}		
      
      Q2(aa,kk,ww,xht,uht);
      Uh->ve[i-1] = Q(xht,uht);
         
      set_row(xh,i-1,xht);
      set_row(uh,i-1,uht);

      for (int j=0;j<=terminator;j++)
	    {
	      xnt->ve[j+1] = xt->ve[j] + k*g(kk,ww,xht->ve[j+1]);
	      unt->ve[j+1] = ut->ve[j]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j+1],Uh->ve[i-1],k,Y));
        }
        
      Q2(aa,kk,ww,xnt,unt);

      Ui->ve[i] = Q(xnt,unt);      

      //printf("%lf\n",Ui->ve[i]);
      
      if(!SGNM) 
        {
          xnt = v_resize(xnt,xn->n);
          unt = v_resize(unt,un->n);

          set_row(x,i,xnt); // todo: test this.
          set_row(u,i,unt);
        }
      else
        {
          set_row(xn,i-1,xnt);
          set_row(un,i-1,unt);

          int idx = idxselect(ww,xnt);
          idxi->ive[i-1] = idx;
        
          xt = idxremove(xnt,xt,idx);
          ut = idxremove(unt,ut,idx);

          set_row(x,i,xt);
          set_row(u,i,ut);
        }
    }
  //  exit(1);

  V_FREE(xt);
  V_FREE(xnt); 
  V_FREE(unt); 
  V_FREE(xht); 
  V_FREE(uht);
  if (QUARTER) {
    V_FREE(xhht);
    V_FREE(uhh);
  }
  V_FREE(ut);

}
