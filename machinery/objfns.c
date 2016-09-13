#include <math.h>
#include "../meschach/matrix.h"
#include "../common.h"
#include "../parameters.h"
#include "spade_solve.h"
#include "spade_solve_clean.h"
#include "../model/fishing/selectivity.h"
#include "../model/biology/weight.h"
#include "../model/fishing/effort.h"
#include "../model/fishing/catch.h"
#include "../model/biology/birth.h"
#include "Q.h"
#include "../util/util.h"

Real K_no(

	  Parameters *parameters

	  )
{

  Real a1 = parameters->alpha1.value;
  Real a2 = parameters->alpha2.value;
  Real bb = parameters->beta.value;
  Real gg = parameters->gamma.value*1e-7;
  Real kk = parameters->kappa.value;
  Real ww = parameters->omega.value;
  Real ii = parameters->iota.value*1e-3;
  
  Real * restrict x = (Real *)calloc(d.J,sizeof(Real));
  Real * restrict u = (Real *)calloc(d.J,sizeof(Real));
  Real * restrict r = (Real *)calloc(d.J,sizeof(Real));
  Real * restrict o = (Real *)calloc(d.J,sizeof(Real));
  
  /* 
     initialize x
  */
  int J = (d.J+1) - (d.I+1);
  
  // 'active': first J values in first 'row'
  for (int j=0;j<J;j++) {
    x[j] = h*j;
  }

  // 'neutral': all the rest (J+1.. I+J)
  for (int j=J;j<=d.J;j++)
    x[j] = ww;

  /*
     ok, now initialize u
  */
  
  // prelims
  Real zeta = sqrt( 81*kk*kk*ww*ww*pow(a1+2*a2*ww,2.) - 12*kk*pow(a1*ww+kk,3.) );
  Real eta = 9*a1*kk*kk*ww + 18*a2*kk*kk*ww*ww + kk*zeta;
  Real Z = pow(eta,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*kk*(a1*ww+kk) / pow(eta,1./3);

  Real ubar = (Z - bb - kk) / gg; 
  Real vbar = (kk*ww*ubar) / (bb+gg*ubar+kk);
  Real wbar = (2*kk*ww*vbar) / (bb+gg*ubar+2*kk);  
  
  // set
  for (int j=0;j<=d.J;j++)
    {
      u[j] = (a1*vbar+a2*wbar)*pow(ww-x[j],(bb+gg*ubar)/kk-1) / (kk*pow(ww,(bb+gg*ubar)/kk));
      r[j] = w(x[j])*s(x[j])*u[j]*ii*_e(d.eff,d.k,d.k*(0-d.N));
    }

  Real ff = 0;
  
  Real U = 0;
  Real C = 0;
  for (int j=0;j<d.J;j++)
    {
      U += .5 * (x[j+1] - x[j]) * (u[j+1] + u[j]);
      C += .5 * (x[j+1] - x[j]) * (r[j+1] + r[j]);
    } 
  
  for (int j=0;j<=d.J;j++) 
    o[j] = pow(_c(d.cat,d.k,d.k*(0-d.N)) * (d.p[0][j] + (1-d.Qp[0]) * r[j]/C) - r[j],2.0);                      
  for (int j=0;j<=d.J;j++)   
    ff += .5 * (x[j+1] - x[j]) * (o[j+1] + o[j]);
  
  Real * restrict xh = (Real *) calloc(d.J,sizeof(Real));  
  Real * restrict uh = (Real *) calloc(d.J,sizeof(Real));  

  for (int i=1;i<=d.I;i++)
    {
  
      Real t = d.k*(i-d.N-1);
      Real th = d.k*(i-d.N-.5);
	 
      for (int j=0;j<=d.J;j++)
	{	  
	  xh[j] = x[j] + d.k/2 * kk*(ww - x[j]);
	  uh[j] = u[j] * exp( -d.k/2 * (bb + gg*U + s(x[j])* ii * _e(d.eff,d.k,t) - kk) );
	}
      
      // only valid for models where fish are size zero at birth
      Real uh_0 = xh[0] * _b(a1,a2,xh[0])*uh[0];

      for (int j=0;j<=d.J;j++)
	uh_0 += (_b(a1,a2,xh[j])*uh[j] + _b(a1,a2,xh[j+1])*uh[j+1]) * (xh[j+1]-xh[j]);

      uh_0 /= (2*kk*ww - xh[0]*_b(a1,a2,0));
            
      Real Uh = .5 * xh[0] * (uh_0 + uh[0]);
      for (int j=0;j<d.J;j++)
	Uh += .5 * (xh[j+1] - xh[j]) * (uh[j+1] + uh[j]);
      
      for (int j=d.J;j>0;j--) 
	{
	  x[j] = x[j-1] + d.k * kk*(ww-xh[j-1]);      
	  u[j] = u[j-1] * exp ( -d.k * (bb + gg*Uh + s(xh[j-1]) * ii * _e(d.eff,d.k,th) - kk));
	}      

      x[0] = 0;
      u[0] = x[1] * _b(a1,a2,x[1])*u[1];
      
      for (int j=1;j<d.J;j++)
	u[0] += (_b(a1,a2,x[j])*u[j] + _b(a1,a2,x[j+1])*u[j+1]) * (x[j+1]-x[j]);

      u[0] /= (2*kk*ww - x[1]*_b(a1,a2,0));

      for (int j=0;j<=d.J;j++)
	r[j] = w(x[j])*s(x[j])*u[j]*ii*_e(d.eff,d.k,d.k*(i-d.N));
            
      U = 0;
      C = 0;
      for (int j=0;j<d.J;j++)
	{
	  U += .5 * (x[j+1] - x[j]) * (u[j+1] + u[j]);
	  C += .5 * (x[j+1] - x[j]) * (r[j+1] + r[j]);
	}
      
      for (int j=0;j<=d.J;j++)
	o[j] = pow(_c(d.cat,d.k,d.k*(i-d.N)) * (d.p[i][j] + (1-d.Qp[i]) * r[j]/C) - r[j],2.0);
      
      for (int j=0;j<=d.J;j++)
	ff += .5 * (x[j+1] - x[j]) * (o[j+1] + o[j]);
      
    }
  
  free(o);
  free(r);
  free(x);
  free(u);
  free(xh);
  free(uh);

  return ff;
	  
}



Real K_no2(

	  Parameters *parameters

	  )
{

  Real aa = 1; //parameters->alpha.value;
  Real bb = parameters->beta.value;
  Real gg = parameters->gamma.value*1e-7;
  Real kk = parameters->kappa.value;
  Real ww = parameters->omega.value;
  Real ii = parameters->iota.value*1e-3;
  
  Real * restrict x = (Real *)calloc(d.J,sizeof(Real));
  Real * restrict u = (Real *)calloc(d.J,sizeof(Real));
  Real * restrict r = (Real *)calloc(d.J,sizeof(Real));
  Real * restrict o = (Real *)calloc(d.J,sizeof(Real));
  
   
  //   initialize x
  
  int J = (d.J+1) - (d.I+1);
  
  // 'active': first J values in first 'row'
  for (int j=0;j<J;j++) {
    x[j] = h*j;
  }

  // 'neutral': all the rest (J+1.. I+J)
  for (int j=J;j<=d.J;j++)
    x[j] = ww;

  
  //   ok, now initialize u
  
  
  // prelims
  Real zeta = sqrt( 81*kk*kk*ww*ww*pow(aa*A1+2*aa*A2*ww,2.) - 12*kk*pow(aa*A1*ww+kk,3.) );
  Real eta = 9*aa*A1*kk*kk*ww + 18*aa*A2*kk*kk*ww*ww + kk*zeta;
  Real Z = pow(eta,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*kk*(aa*A1*ww+kk) / pow(eta,1./3);

  Real ubar = (Z - bb - kk) / gg; 
  Real vbar = (kk*ww*ubar) / (bb+gg*ubar+kk);
  Real wbar = (2*kk*ww*vbar) / (bb+gg*ubar+2*kk);  
  
  // set
  for (int j=0;j<=d.J;j++)
    {
      u[j] = (aa*A1*vbar+aa*A2*wbar)*pow(ww-x[j],(bb+gg*ubar)/kk-1) / (kk*pow(ww,(bb+gg*ubar)/kk));
      r[j] = w(x[j])*s(x[j])*u[j]*ii*_e(d.eff,d.k,d.k*(0-d.N));
    }

  Real ff = 0;
  
  Real U = 0;
  Real C = 0;
  for (int j=0;j<d.J;j++)
    {
      U += .5 * (x[j+1] - x[j]) * (u[j+1] + u[j]);
      C += .5 * (x[j+1] - x[j]) * (r[j+1] + r[j]);
    } 
  
  for (int j=0;j<=d.J;j++) 
    o[j] = pow(_c(d.cat,d.k,d.k*(0-d.N)) * (d.p[0][j] + (1-d.Qp[0]) * r[j]/C) - r[j],2.0);                      
  for (int j=0;j<=d.J;j++)   
    ff += .5 * (x[j+1] - x[j]) * (o[j+1] + o[j]);

  //  printf("\n%lf\n",ff);
  
  Real * restrict xh = (Real *) calloc(d.J,sizeof(Real));  
  Real * restrict uh = (Real *) calloc(d.J,sizeof(Real));  

  for (int i=1;i<=d.I;i++)
    {
  
      Real t = d.k*(i-d.N-1);
      Real th = d.k*(i-d.N-.5);
	 
      for (int j=0;j<d.J;j++)
	{	  
	  xh[j] = x[j] + d.k/2 * kk*(ww - x[j]);
	  uh[j] = u[j] * exp( -d.k/2 * (bb + gg*U + s(x[j])* ii * _e(d.eff,d.k,t) - kk) );
	}
      
      // only valid for models where fish are size zero at birth
      Real uh_0 = xh[0] * b(aa,xh[0])*uh[0];

      for (int j=0;j<d.J-1;j++)
	uh_0 += (b(aa,xh[j])*uh[j] + b(aa,xh[j+1])*uh[j+1]) * (xh[j+1]-xh[j]);

      uh_0 /= (2*kk*ww - xh[0]*b(aa,0));
            
      Real Uh = .5 * xh[0] * (uh_0 + uh[0]);
      for (int j=0;j<d.J-1;j++)
	Uh += .5 * (xh[j+1] - xh[j]) * (uh[j+1] + uh[j]);
      
      for (int j=d.J;j>0;j--) 
	{
	  x[j] = x[j-1] + d.k * kk*(ww-xh[j-1]);      
	  u[j] = u[j-1] * exp ( -d.k * (bb + gg*Uh + s(xh[j-1]) * ii * _e(d.eff,d.k,th) - kk));
	}      

      x[0] = 0;
      u[0] = x[1] * b(aa,x[1])*u[1];
      
      for (int j=1;j<d.J;j++)
	u[0] += (b(aa,x[j])*u[j] + b(aa,x[j+1])*u[j+1]) * (x[j+1]-x[j]);

      u[0] /= (2*kk*ww - x[1]*b(aa,0));

      for (int j=0;j<=d.J;j++)
	r[j] = w(x[j])*s(x[j])*u[j]*ii*_e(d.eff,d.k,d.k*(i-d.N));
            
      U = 0;
      C = 0;
      for (int j=0;j<d.J;j++)
	{
	  U += .5 * (x[j+1] - x[j]) * (u[j+1] + u[j]);
	  C += .5 * (x[j+1] - x[j]) * (r[j+1] + r[j]);
	}
      
      for (int j=0;j<=d.J;j++)
	o[j] = pow(_c(d.cat,d.k,d.k*(i-d.N)) * (d.p[i][j] + (1-d.Qp[i]) * r[j]/C) - r[j],2.0);
      
      for (int j=0;j<=d.J;j++)
	ff += .5 * (x[j+1] - x[j]) * (o[j+1] + o[j]);

      //      printf("%lf\n",ff);
      
    }

  // exit(1);
  
  free(o);
  free(r);
  free(x);
  free(u);
  free(xh);
  free(uh);

  return ff;
	  
}


Real K_dr(

  Parameters * parameters,
  Data *d
  
  )
{

  Solve_Core_Args core_args;
  int I = d->I;
  int J = d->J;
    
  core_args.x = m_get(I,J);
  core_args.u = m_get(I,J);
  core_args.xh = m_get(I,J+1);
  core_args.uh = m_get(I,J+1);
  core_args.xn = m_get(I,J+1);
  core_args.xhh = m_get(I,J+1);
  core_args.un = m_get(I,J+1);
  core_args.Ui = v_get(I);
  core_args.Uh = v_get(I);
  core_args.Uhh = v_get(I);
  core_args.idxi = iv_get(I-1);   
  
  solve(parameters,d->eff,d->k,d->S,d->Y,&core_args);

  MeMAT *x = core_args.x;
  MeMAT *u = core_args.u;

  Real iota = parameters->iota.value;

  iota *= 1e-3;
  int S = d->S;
  Real k = d->k;

  int lfi=0;

  MeVEC *ht = v_get(x->m);
  MeVEC *tt = v_get(x->m);

  for (int i=S;i<x->m;i++)
    {

      MeVEC *xt = v_get(x->n);
      get_row(x,i,xt);      
      MeVEC *ut = v_get(x->n);
      get_row(u,i,ut);      
      MeVEC *v = v_get(x->n);

      for (int j=0;j<x->n;j++)
  v->ve[j] = iota*e(d->eff,k,k*(i-S),d->Y)*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];

      if(lfi < d->n && d->t_id[lfi]==i) 
  {
      
    MeVEC *dt = v_get(d->t_sz[lfi]);

    for (int j=0;j<dt->dim;j++)
      dt->ve[j] = d->lf[lfi][j];

    Real bw = get_bw(dt);

    MeVEC *l = v_get(xt->dim);

    for (int j=0;j<xt->dim;j++)
      for (int jj=0;jj<dt->dim;jj++)
        l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

    Real al = 1e3*c(d->cat,k,k*(i - S)) / Q(xt,l); 

    MeVEC *ld = v_get(x->n);

    for (int j=0;j<xt->dim;j++)
      ld->ve[j] = pow(v->ve[j] - al*l->ve[j],2.);

    ht->ve[i] = Q(xt,ld);

    lfi += 1;

    V_FREE(dt);
    V_FREE(l);
    V_FREE(ld);

  } 
      else 
  {
    ht->ve[i] = pow(Q(xt,v)-1e3*c(d->cat,k,k*(i - S)),2.);
  }

      tt->ve[i] = k*(i-S);

      V_FREE(xt);
      V_FREE(ut);
      V_FREE(v);

    }

  Real blah = Q(tt,ht);

  V_FREE(ht);
  V_FREE(tt);

  M_FREE(core_args.x);
  M_FREE(core_args.u);
  M_FREE(core_args.xh);
  M_FREE(core_args.uh);
  M_FREE(core_args.xn);
  M_FREE(core_args.xhh);
  M_FREE(core_args.un);
  V_FREE(core_args.Ui);
  V_FREE(core_args.Uh);
  V_FREE(core_args.Uhh);
  IV_FREE(core_args.idxi);

  return blah; 

}

	  


Real K(

       Parameters *parameters,
       Data *d,
       Solve_Core_Args *core_args
		 	
       )
{


  if (MESCHACH)
    solve(parameters,d->eff,d->k,d->S,d->Y,core_args);
  else
    solve_clean(parameters,d->eff,d->k,d->S,d->Y,core_args);
  
  MeMAT *x = core_args->x;
  MeMAT *u = core_args->u;

  Real iota = parameters->iota.value;

  iota *= 1e-3;
  int S = d->S;
  Real k = d->k;

  int lfi=0;

  MeVEC *ht = v_get(x->m);
  MeVEC *tt = v_get(x->m);

  int J; 
  int bigJ;
  if (!SGNM){
    bigJ = x->n - 1;
    J = x->n - x->m;
  } else {
    J = x->n - 1;   
  }
  
  for (int i=S;i<x->m;i++)
    {

      int terminator;
      if(!SGNM)
	terminator = J+i;
      else
	terminator = J;

      MeVEC *xt = v_get(terminator+1);
      get_row(x,i,xt);      
      MeVEC *ut = v_get(terminator+1);
      get_row(u,i,ut);
      
      if (!SGNM)
	{
	  xt = v_resize(xt,terminator+1);
	  ut = v_resize(ut,terminator+1);
	}
      
      MeVEC *v = v_get(terminator+1);

      for (int j=0;j<=terminator;j++)
	v->ve[j] = iota*e(d->eff,k,k*(i-S),d->Y)*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];

      if(lfi < d->n && d->t_id[lfi]==i) 
	{
      
 	  MeVEC *dt = v_get(d->t_sz[lfi]);

	  for (int j=0;j<dt->dim;j++)
	    dt->ve[j] = d->lf[lfi][j];

	  Real bw = get_bw(dt);

	  MeVEC *l = v_get(terminator+1);

	  for (int j=0;j<=terminator;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  Real al = 1e3*c(d->cat,k,k*(i - S)) / Q(xt,l); 

	  MeVEC *ld = v_get(terminator+1);

	  //printf("\n");
	  for (int j=0;j<=terminator;j++){
	    
	    ld->ve[j] = pow(v->ve[j] - al*l->ve[j],2.);
		    //ld->ve[j] = fabs(v->ve[j] - al*l->ve[j]);
	    //printf("%Lf %Lf\n",xt->ve[j],ld->ve[j]);
	  }
	  //exit(1);
	  
	  ht->ve[i] = Q(xt,ld);

	  lfi += 1;

	  V_FREE(dt);
	  V_FREE(l);
	  V_FREE(ld);

	} 
      else 
	{
	  ht->ve[i] = pow(Q(xt,v)-1e3*c(d->cat,k,k*(i - S)),2.);
		  //ht->ve[i] = fabs(Q(xt,v)-1e3*c(d->cat,k,k*(i - S)));
	}

      tt->ve[i] = k*(i-S);

      V_FREE(xt);
      V_FREE(ut);
      V_FREE(v);

    }

  Real blah = Q(tt,ht);

  V_FREE(ht);
  V_FREE(tt);

  return blah; 

}

Real G(

	 MeMAT *p,
	 MeMAT *x,
	 MeMAT *u,
	 Data *data,
	 Real iota
	
	 )
{

  iota *=1e-3;
  int lfi=0;

  int S = data->S;
  Real k = data->k;

  MeVEC *ht = v_get(x->m);
  MeVEC *tt = v_get(x->m);

  int J; 
  int bigJ;
  if (!SGNM){
    bigJ = x->n - 1;
    J = x->n - x->m;

  } else {
    J = x->n - 1;
   
  }

  
  for (int i=S;i<x->m;i++)
    {

      int terminator;
      if(!SGNM)
	terminator = J+i;
      else
	terminator = J;

      MeVEC *xt = v_get(terminator+1);
      get_row(x,i,xt);      
      MeVEC *ut = v_get(terminator+1);
      get_row(u,i,ut);      
      MeVEC *pt = v_get(terminator+1);
      get_row(p,i,pt); 

      if (!SGNM)
	{
	  xt = v_resize(xt,terminator+1);
	  ut = v_resize(ut,terminator+1);
	  pt = v_resize(pt,terminator+1);
	}
      
      MeVEC *v = v_get(terminator+1);
      MeVEC *pv = v_get(terminator+1);

      for (int j=0;j<=terminator;j++) 
	{
	  pv->ve[j] = iota*e(data->eff,k,k*(i-S),data->Y)*s(xt->ve[j])*w(xt->ve[j])*pt->ve[j] + 1e-3*e(data->eff,k,k*(i-S),data->Y)*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];
	  v->ve[j] = iota*e(data->eff,k,k*(i-S),data->Y)*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];
	}

      //      if(data->t_id[lfi]==i) 
      if(lfi < data->n && data->t_id[lfi]==i) 
	{
      
	  MeVEC *dt = v_get(data->t_sz[lfi]);

	  for (int j=0;j<dt->dim;j++)
	    dt->ve[j] = data->lf[lfi][j];

	  Real bw = get_bw(dt);

	  MeVEC *l = v_get(terminator+1);

	  for (int j=0;j<=terminator;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  Real al = 1e3*c(data->cat,k,k*(i - S)) / Q(xt,l);

	  MeVEC *ld = v_get(terminator+1);

	  for (int j=0;j<=terminator;j++) {
	    ld->ve[j] = 2*(v->ve[j] - al*l->ve[j])*pv->ve[j];

	    /*
	    if (v->ve[j] < al*l->ve[j])
	      ld->ve[j] = -pv->ve[j];
	    else
	      ld->ve[j] = pv->ve[j];
	    */
	  }
	  //	    ld->ve[j] = pv->ve[j]*(v->ve[j] - al*l->ve[j]) / fabs(v->ve[j] - al*l->ve[j]);

	  ht->ve[i] = Q(xt,ld);

	  lfi += 1;

	  V_FREE(dt);
	  V_FREE(l);
	  V_FREE(ld);

	} 
      else 
	{ 
          //printf("%d\n",i);
          //printf("%f\n",Q(xt,pt));
          //printf("%f\n",Q(xt,v));
          //printf("%f\n",c(k*(i-tmi.I)));

	  ht->ve[i] = 2*(Q(xt,v)-1e3*c(data->cat,k,k*(i-S)))*Q(xt,pv);

	  /*
	  if (Q(xt,v)<1e3*c(data->cat,k,k*(i-S)))
	    ht->ve[i] = -Q(xt,pv);//*(Q(xt,v)-1e3*c(k*(i - tmi.I)))/fabs(Q(xt,v)-1e3*c(k*(i - tmi.I)));
	  else
	  ht->ve[i]=Q(xt,pv);*/
	}

      tt->ve[i] = k*(i-S);

      V_FREE(xt);
      V_FREE(ut);
      V_FREE(pt);
      V_FREE(v);
      V_FREE(pv);

    }

  Real blah = Q(tt,ht);

  V_FREE(ht);
  V_FREE(tt);

  return blah;

}


/*Real G_ni_for_condition_number(

	    MeMAT *p,
	    MeMAT *x,
	    MeMAT *u,
	    Data *data,
	    Real iota

	    )
{

  int S = data->S;
  Real k = data->k;

  int lfi=0;

  MeVEC *ht = v_get(x->m);
  MeVEC *tt = v_get(x->m);

  for (int i=S;i<x->m;i++)
    {

      MeVEC *xt = v_get(x->n);
      get_row(x,i,xt);      
      MeVEC *ut = v_get(x->n);
      get_row(u,i,ut);      
      MeVEC *pt = v_get(x->n);
      get_row(p,i,pt); 

      MeVEC *v = v_get(x->n);
      MeVEC *pv = v_get(x->n);

      for (int j=0;j<x->n;j++) 
	{
	  pv->ve[j] = iota*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*pt->ve[j];
	  v->ve[j] = iota*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];
	}

      if(data->t_id[lfi]==i) 
	{
      
	  MeVEC *l = v_get(xt->dim);
	  for (int j=0;j<xt->dim;j++)
	    l->ve[j] = ut->ve[j]; 

	  for (int j=0;j<xt->dim;j++)
	    if (v->ve[j] < al*l->ve[j])
	      ld->ve[j] = -pv->ve[j];
	    else
	      ld->ve[j] = pv->ve[j];

	  ht->ve[i] = Q(xt,ld);

          if (lfi<data->n)
	    lfi += 1;

	} 
      else 
	{ 
          //printf("%d\n",i);
          //printf("%f\n",Q(xt,pt));
          //printf("%f\n",Q(xt,v));
          //printf("%f\n",c(k*(i-tmi.I)));
        
	  if (Q(xt,v)<1e3*c(data->cat,k,k*(i-S)))
	    ht->ve[i] = -Q(xt,pv);//*(Q(xt,v)-1e3*c(k*(i - tmi.I)))/fabs(Q(xt,v)-1e3*c(k*(i - tmi.I)));
	  else
	    ht->ve[i]=Q(xt,pv);
	}

      tt->ve[i] = k*(i-S);

    }


  if (BLAH) {

    FILE *p1 = fopen("plot.txt","w");

    for (int i=S;i<x->m;i++)
      fprintf(p1,"%f %f\n",tt->ve[i],ht->ve[i]);

    fclose(p1);

    char buffer[100];

    sprintf(buffer,"./plo1 > plotht.pdf");

    system(buffer);

  }

  return Q(tt,ht);

}
*/

Real G_ni(

	    MeMAT *p,
	    MeMAT *x,
	    MeMAT *u,
	    Data *data,
	    Real iota

	    )
{

  iota *=1e-3;
  int S = data->S;
  Real k = data->k;

  int lfi=0;

  MeVEC *ht = v_get(x->m);
  MeVEC *tt = v_get(x->m);

  int J; 
  int bigJ;
  
  if (!SGNM){
    bigJ = x->n - 1;
    J = x->n - x->m;

  } else {
    J = x->n - 1;
   
  }

  for (int i=S;i<x->m;i++)
    {

      int terminator;
      if(!SGNM)
	terminator = J+i;
      else
	terminator = J;

      MeVEC *xt = v_get(terminator+1);
      get_row(x,i,xt);      
      MeVEC *ut = v_get(terminator+1);
      get_row(u,i,ut);      
      MeVEC *pt = v_get(terminator+1);
      get_row(p,i,pt); 

      if (!SGNM)
	{
	  xt = v_resize(xt,terminator+1);
	  ut = v_resize(ut,terminator+1);
	  pt = v_resize(pt,terminator+1);
	}
      
      MeVEC *v = v_get(terminator+1);
      MeVEC *pv = v_get(terminator+1);
      
      for (int j=0;j<=terminator;j++) 
	{
	  pv->ve[j] = iota*e(data->eff,k,k*(i-S),data->Y)*s(xt->ve[j])*w(xt->ve[j])*pt->ve[j];
	  v->ve[j] = iota*e(data->eff,k,k*(i-S),data->Y)*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];
	}

      if(lfi < data->n && data->t_id[lfi]==i) 
	{
      
	  MeVEC *dt = v_get(data->t_sz[lfi]);

	  for (int j=0;j<dt->dim;j++)
	    dt->ve[j] = data->lf[lfi][j];

	  Real bw = get_bw(dt);

	  MeVEC *l = v_get(terminator+1);

	  for (int j=0;j<=terminator;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  Real al = 1e3*c(data->cat,k,k*(i - S)) / Q(xt,l);

	  MeVEC *ld = v_get(terminator+1);

	  //printf("\n");
	  for (int j=0;j<=terminator;j++)
	    {
	    
	      ld->ve[j] = 2*(v->ve[j] - al*l->ve[j])*pv->ve[j];
	    //printf("%Lf %Lf\n",xt->ve[j],pv->ve[j]);
	  
	  //exit(1);
	  /*
	      //	      for (int j=0;j<xt->dim;j++)
	      if (v->ve[j] < al*l->ve[j])
		ld->ve[j] = -pv->ve[j];
	      else
	      ld->ve[j] = pv->ve[j];*/
	    }

	  // Ld->ve[j] = pv->ve[j]*(v->ve[j] - al*l->ve[j]) / fabs(v->ve[j] - al*l->ve[j]);

	  ht->ve[i] = Q(xt,ld);

	  lfi += 1;

	  V_FREE(dt);
	  V_FREE(l);
	  V_FREE(ld);

	} 
      else 
	{ 
          //printf("%d\n",i);
          //printf("%f\n",Q(xt,pt));
          //printf("%f\n",Q(xt,v));
          //printf("%f\n",c(k*(i-tmi.I)));

	  ht->ve[i] = 2*(Q(xt,v)-1e3*c(data->cat,k,k*(i-S)))*Q(xt,pv);
	  /*        
	  if (Q(xt,v)<1e3*c(data->cat,k,k*(i-S)))
	    ht->ve[i] = -Q(xt,pv);//*(Q(xt,v)-1e3*c(k*(i - tmi.I)))/fabs(Q(xt,v)-1e3*c(k*(i - tmi.I)));
	  else
	  ht->ve[i] = Q(xt,pv);*/
	}

      tt->ve[i] = k*(i-S);

      V_FREE(xt);
      V_FREE(ut);
      V_FREE(pt);
      V_FREE(v);
      V_FREE(pv);

    }

  Real blah = Q(tt,ht);

  V_FREE(ht);
  V_FREE(tt);

  return blah;

}

/*
Real newK(
	  
     Parameters *parameters,
		 Data *d,
		 Solve_Core_Args *core_args
		 
	
		 )
{


  solve(parameters,d->eff,d->k,d->S,d->Y,core_args);

  MeMAT *x = core_args->x;
  MeMAT *u = core_args->u;

  Real iota = parameters->iota.value;

  iota *= 1e-3;
  int S = d->S;
  Real k = d->k;

  int lfi=0;

  MeVEC *ht = v_get(x->m);
  MeVEC *tt = v_get(x->m);

  int J; 
  int bigJ;
  if (BIGMeMATRICES){
    bigJ = x->n - 1;
    J = x->n - x->m;

  } else {
    J = x->n - 1;
   
  }
  
  for (int i=S;i<x->m;i++)
    {

      int terminator;
      if(BIGMeMATRICES)
	terminator = J+i;
      else
	terminator = J;

      MeVEC *xt = v_get(terminator+1);
      get_row(x,i,xt);      
      MeVEC *ut = v_get(terminator+1);
      get_row(u,i,ut);
      
      if (BIGMeMATRICES)
	{
	  xt = v_resize(xt,terminator+1);
	  ut = v_resize(ut,terminator+1);
	}
      
      MeVEC *v = v_get(terminator+1);

      for (int j=0;j<=terminator;j++)
	v->ve[j] = iota*e(d->eff,k,k*(i-S),d->Y)*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];

      if(lfi < d->n && d->t_id[lfi]==i) 
	{
      
 	  MeVEC *dt = v_get(d->t_sz[lfi]);

	  for (int j=0;j<dt->dim;j++)
	    dt->ve[j] = d->lf[lfi][j];

	  Real bw = get_bw(dt);

	  MeVEC *l = v_get(terminator+1);

	  for (int j=0;j<=terminator;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  Real al = 1e3*c(d->cat,k,k*(i - S)) / Q(xt,l); 

	  MeVEC *ld = v_get(terminator+1);

	  //printf("\n");
	  for (int j=0;j<=terminator;j++){
	    
	    ld->ve[j] = pow(v->ve[j] - al*l->ve[j],2.);
		    //ld->ve[j] = fabs(v->ve[j] - al*l->ve[j]);
	    //printf("%Lf %Lf\n",xt->ve[j],ld->ve[j]);
	  }
	  //exit(1);
	  
	  ht->ve[i] = Q(xt,ld);

	  lfi += 1;

	  V_FREE(dt);
	  V_FREE(l);
	  V_FREE(ld);

	}
*/
