#include <math.h>
#include "plot.h"
#include "../meschach/matrix.h"
#include "../common.h"
#include "../parameters.h"
#include "../model/fishing/selectivity.h"
#include "../model/fishing/effort.h"
#include "../model/fishing/catch.h"
#include "../model/biology/weight.h"
#include "../model/biology/birth.h"

void plot_fast(

	  Parameters *parameters,
	  char *label

	  )
{

  
  Real a1 = parameters->alpha1.value;
  Real a2 = parameters->alpha2.value;
  Real bb = parameters->beta.value;
  Real gg = parameters->gamma.value*1e-7;
  Real kk = parameters->kappa.value;
  Real ww = parameters->omega.value;
  Real ii = parameters->iota.value*1e-3;
  
  Real * restrict x = (Real *)calloc(d.J+1,sizeof(Real));
  Real * restrict u = (Real *)calloc(d.J+1,sizeof(Real));
  Real * restrict r = (Real *)calloc(d.J+2,sizeof(Real));
  Real * restrict o = (Real *)calloc(d.J+2,sizeof(Real));
  
  /* 
     initialize x
  */
  for (int j=0;j<=d.J;j++) {
    x[j] = h*j;
  }

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
      r[j] = 1e-3*w(x[j])*s(x[j])*u[j]*ii*_e(d.eff,0);
    }

  Real ff = 0;
  
  Real U = 0;
  Real C = 0;
  
  for (int j=0;j<d.J;j++)
    {
      U += .5 * (x[j+1] - x[j]) * (u[j+1] + u[j]);
      C += .5 * (x[j+1] - x[j]) * (r[j+1] + r[j]);
    } 

  Real * restrict robs = (Real *)calloc(d.J+2,sizeof(Real));
  
  for (int j=0;j<=d.J;j++)
    {
      robs[j] = _c(d.cat,0) * (d.p[0][j] + (1-d.Qp[0]) * r[j]/C);
      o[j] = pow(robs[j]- r[j],2.0);
    }
  
  Real R = 0;

  for (int j=0;j<d.J;j++)
    R += .5 * (x[j+1] - x[j]) * (robs[j+1] + robs[j]);      
  
  Real * restrict xh = (Real *) calloc(d.J+1,sizeof(Real));  
  Real * restrict uh = (Real *) calloc(d.J+1,sizeof(Real));  

  Real * restrict xn = (Real *)calloc(d.J+2,sizeof(Real));
  Real * restrict un = (Real *)calloc(d.J+2,sizeof(Real));

  
  FILE *sp1 = fopen("plot1.txt","w");
  FILE *sp2 = fopen("plot2.txt","w");
  FILE *sp3 = fopen("plot.txt","w");
  FILE *sp4 = fopen("plote.txt","w");
  FILE *sp5; 
  FILE *sp6;

  char sbufferblah[100];
  
  fprintf(sp1,"%lf %lf\n",0.0,C);
  fprintf(sp2,"%lf %lf\n",0.0,R);
  fprintf(sp3,"%lf %lf\n",0.0,U);
  fprintf(sp4,"%lf %lf\n",0.0,ii*_e(d.eff,0.0));

  for (int i=1;i<=d.I;i++)
    {

      sprintf(sbufferblah,"plotlfa_%d.txt",i);
      sp5 = fopen(sbufferblah,"w");

      sprintf(sbufferblah,"plotlfb_%d.txt",i);
      sp6 = fopen(sbufferblah,"w");
      		  
      Real t = k*(i-1);
      Real th = k*(i-.5);
	 
      for (int j=0;j<=d.J;j++)
	{	  
	  xh[j] = x[j] + k/2 * kk*(ww - x[j]);
	  uh[j] = u[j] * exp( -k/2 * (bb + gg*U + s(x[j])* ii * _e(d.eff,t) - kk) );
	}
      
      // only valid for models where fish are size zero at birth
      Real uh_0 = xh[0] * _b(a1,a2,xh[0])*uh[0];

      for (int j=0;j<=d.J;j++)
	uh_0 += (_b(a1,a2,xh[j])*uh[j] + _b(a1,a2,xh[j+1])*uh[j+1]) * (xh[j+1]-xh[j]);

      uh_0 /= (2*kk*ww - xh[0]*_b(a1,a2,0));
            
      Real Uh = .5 * xh[0] * (uh_0 + uh[0]);
      for (int j=0;j<d.J;j++)
	Uh += .5 * (xh[j+1] - xh[j]) * (uh[j+1] + uh[j]);
      
      for (int j=d.J+1;j>0;j--) 
	{
	  xn[j] = x[j-1] + k * kk*(ww-xh[j-1]);      
	  un[j] = u[j-1] * exp ( -k * (bb + gg*Uh + s(xh[j-1]) * ii * _e(d.eff,th) - kk));
	}      

      xn[0] = 0;
      un[0] = xn[1] * _b(a1,a2,xn[1])*un[1];
      
      for (int j=1;j<=d.J;j++)
	un[0] += (_b(a1,a2,xn[j])*un[j] + _b(a1,a2,xn[j+1])*un[j+1]) * (xn[j+1]-xn[j]);

      un[0] /= (2*kk*ww - xn[1]*_b(a1,a2,0));

      for (int j=0;j<=d.J+1;j++)
	{	  
	  r[j] = 1e-3*w(xn[j])*s(xn[j])*un[j]*ii*_e(d.eff,k*i);
	  fprintf(sp5,"%lf %lf\n",xn[j],r[j]);
	}
            	        
      U = 0;
      C = 0;
      for (int j=0;j<=d.J;j++)
	{
	  U += .5 * (xn[j+1] - xn[j]) * (un[j+1] + un[j]);
	  C += .5 * (xn[j+1] - xn[j]) * (r[j+1] + r[j]);
	}
      
      for (int j=0;j<=d.J+1;j++)
	{

	  robs[j] = _c(d.cat,k*i) * (d.p[i][j] + (1-d.Qp[i]) * r[j]/C);	  
	  o[j] = pow(robs[j] - r[j],2.0);

	  fprintf(sp6,"%lf %lf\n",xn[j],robs[j]);
	  
	}
      
      fprintf(sp1,"%lf %lf\n",k*i,C);

      R = 0;
      for (int j=0;j<=d.J;j++)	
	R += .5 * (xn[j+1] - xn[j]) * (robs[j+1] + robs[j]);

      fprintf(sp2,"%lf %lf\n",k*i,R); 
      
      for (int j=0;j<=d.J;j++)
	ff += k * .5 * (xn[j+1] - xn[j]) * (o[j+1] + o[j]);

      fprintf(sp3,"%lf %lf\n",k*i,U);
      fprintf(sp4,"%lf %lf\n",k*i,ii*_e(d.eff,k*i));

      // remove
      for (int j=0;j<idx[i-1];j++)
	x[j] = xn[j];
      for (int j=idx[i-1];j<=d.J;j++)
	x[j] = xn[j+1];      

      for (int j=0;j<idx[i-1];j++)
	u[j] = un[j];
      for (int j=idx[i-1];j<=d.J;j++)
	u[j] = un[j+1];

      fclose(sp5);
      fclose(sp6);
      
    }

  fclose(sp1);
  fclose(sp2);
  fclose(sp3);
  fclose(sp4);
  
  char sbuffer[100];
  sprintf(sbuffer,"./plo > plotc_%s.pdf",label);
  system(sbuffer);

  char sbuffer2[100];
  sprintf(sbuffer2,"./plo1 > plotu_%s.pdf",label);
  system(sbuffer2);

  char sbuffer3[100];
  sprintf(sbuffer3,"./ploe > plote_%s.pdf",label);
  system(sbuffer3);

  
  char sbuffer4[100];
  
  for (int i=1;i<=d.I;i++)
    {

      sprintf(sbufferblah,"mv plotlfa_%d.txt plot1.txt",i);
      system(sbufferblah);
      
      sprintf(sbufferblah,"mv plotlfb_%d.txt plot2.txt",i);
      system(sbufferblah);

      sprintf(sbuffer4,"./plo > plotlf_%d_%s.pdf",i,label);
      system(sbuffer4);

    }
  
  
  free(o);
  free(r);
  free(x);
  free(u);
  free(xh);
  free(uh);
  free(xn);
  free(un);

}

void plot(

	  Parameters *parameters,
	  char *label

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
      r[j] = 1e-3*w(x[j])*s(x[j])*u[j]*ii*_e(d.eff,0);
    }

  Real ff = 0;
  
  Real U = 0;
  Real C = 0;
  
  for (int j=0;j<d.J;j++)
    {
      U += .5 * (x[j+1] - x[j]) * (u[j+1] + u[j]);
      C += .5 * (x[j+1] - x[j]) * (r[j+1] + r[j]);
    } 

  Real * restrict robs = (Real *)calloc(d.J,sizeof(Real));
  
  for (int j=0;j<=d.J;j++)
    {
      robs[j] = _c(d.cat,0) * (d.p[0][j] + (1-d.Qp[0]) * r[j]/C);
      o[j] = pow(robs[j]- r[j],2.0);
    }
  
  Real R = 0;

  for (int j=0;j<=d.J;j++)
    R += .5 * (x[j+1] - x[j]) * (robs[j+1] + robs[j]);      
  
  Real * restrict xh = (Real *) calloc(d.J,sizeof(Real));  
  Real * restrict uh = (Real *) calloc(d.J,sizeof(Real));  

  FILE *sp1 = fopen("plot1.txt","w");
  FILE *sp2 = fopen("plot2.txt","w");
  FILE *sp3 = fopen("plot.txt","w");
  FILE *sp4 = fopen("plote.txt","w");

  fprintf(sp1,"%lf %lf\n",0.0,C);
  fprintf(sp2,"%lf %lf\n",0.0,R);
  fprintf(sp3,"%lf %lf\n",0.0,U);
  fprintf(sp4,"%lf %lf\n",0.0,ii*_e(d.eff,0.0));
  
  for (int i=1;i<=d.I;i++)
    {
  
      Real t = k*(i-1);
      Real th = k*(i-.5);
	 
      for (int j=0;j<=d.J;j++)
	{	  
	  xh[j] = x[j] + k/2 * kk*(ww - x[j]);
	  uh[j] = u[j] * exp( -k/2 * (bb + gg*U + s(x[j])* ii * _e(d.eff,t) - kk) );
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
	  x[j] = x[j-1] + k * kk*(ww-xh[j-1]);      
	  u[j] = u[j-1] * exp ( -k * (bb + gg*Uh + s(xh[j-1]) * ii * _e(d.eff,th) - kk));
	}      

      x[0] = 0;
      u[0] = x[1] * _b(a1,a2,x[1])*u[1];
      
      for (int j=1;j<d.J;j++)
	u[0] += (_b(a1,a2,x[j])*u[j] + _b(a1,a2,x[j+1])*u[j+1]) * (x[j+1]-x[j]);

      u[0] /= (2*kk*ww - x[1]*_b(a1,a2,0));

      for (int j=0;j<=d.J;j++)
	r[j] = 1e-3*w(x[j])*s(x[j])*u[j]*ii*_e(d.eff,k*i);
      
      U = 0;
      C = 0;
      for (int j=0;j<d.J;j++)
	{
	  U += .5 * (x[j+1] - x[j]) * (u[j+1] + u[j]);
	  C += .5 * (x[j+1] - x[j]) * (r[j+1] + r[j]);
	}
      
      for (int j=0;j<=d.J;j++)
	{

	  robs[j] = _c(d.cat,k*i) * (d.p[i][j] + (1-d.Qp[i]) * r[j]/C);
	  
	  o[j] = pow(robs[j] - r[j],2.0);
	}
      
      fprintf(sp1,"%lf %lf\n",k*i,C);

      R = 0;
      for (int j=0;j<d.J;j++)	
	R += .5 * (x[j+1] - x[j]) * (robs[j+1] + robs[j]);

      fprintf(sp2,"%lf %lf\n",k*i,R); 
      
      for (int j=0;j<=d.J;j++)
	ff += k * .5 * (x[j+1] - x[j]) * (o[j+1] + o[j]);

      fprintf(sp3,"%lf %lf\n",k*i,U);
      fprintf(sp4,"%lf %lf\n",k*i,ii*_e(d.eff,k*i));
      
    }

  fclose(sp1);
  fclose(sp2);
  fclose(sp3);
  fclose(sp4);
  
  char sbuffer[100];
  sprintf(sbuffer,"./plo > plotc_%s.pdf",label);
  system(sbuffer);

  char sbuffer2[100];
  sprintf(sbuffer2,"./plo1 > plotu_%s.pdf",label);
  system(sbuffer2);

  char sbuffer3[100];
  sprintf(sbuffer3,"./ploe > plote_%s.pdf",label);
  system(sbuffer3);

  free(o);
  free(r);
  free(x);
  free(u);
  free(xh);
  free(uh);

}

/*
void plot_old( 

	  Parameters *parameters,
	  Data *d,
	  char * label

	   )
{

  Solve_Core_Args core_args;
  
  int I,J;
  I = d->I;
  J = d->J;

  core_args.x = m_get(I+1,J+1);
  core_args.u = m_get(I+1,J+1);

  // todo: Review this. Condition has been flipped to maintain same functionality as last commit
  // in non-bigmatrices mode but may not be correct for bigmatrices mode. Affects xh in particular
  // as it was of dimension 401 here but was of dimension 402 in working copy.
  if (SGNM)
    {
       core_args.xh = m_get(I+1,J+2);
       core_args.uh = m_get(I+1,J+2);
       core_args.xn = m_get(I+1,J+2);
       core_args.xhh = m_get(I+1,J+2);
       core_args.un = m_get(I+1,J+2);
     }
  else
    {
       core_args.xh = m_get(I+1,J+1);
       core_args.uh = m_get(I+1,J+1);
       core_args.xn = m_get(I+1,J+1);
       core_args.xhh = m_get(I+1,J+1);
       core_args.un = m_get(I+1,J+1);
    }

  core_args.Ui = v_get(I+1);
  core_args.Uh = v_get(I+1);
  core_args.Uhh = v_get(I+1);
  core_args.idxi = iv_get(I);

  solve(parameters,d->eff,d->k,d->S,d->Y,&core_args);

  MAT *x = core_args.x;
  MAT *u = core_args.u;

  M_FREE(core_args.xh);
  M_FREE(core_args.uh);
  M_FREE(core_args.xn);
  M_FREE(core_args.xhh);
  M_FREE(core_args.un);
  V_FREE(core_args.Ui);
  V_FREE(core_args.Uh);
  V_FREE(core_args.Uhh);
  IV_FREE(core_args.idxi);

  Real iota = parameters->iota.value;
  iota *= 1e-3;

  int bigJ;
  if (!SGNM){
    bigJ = x->n - 1;
    J = x->n - x->m;

  } else {
    J = x->n - 1;   
  }
  
  FILE *sp1 = fopen("plot1.txt","w");

  for (int i=0;i<x->m;i++)
    {

      int terminator;
      if(!SGNM)
	terminator = J+i;
      else
	terminator = J;

      VEC *ctt = v_get(terminator+1);
      VEC *xt = v_get(terminator+1);      
      
      xt = get_row(x,i,xt);

      if (!SGNM)
	xt = v_resize(xt,terminator+1);
      
      for (int j=0;j<=terminator;j++)
	ctt->ve[j] = s(x->me[i][j])*iota*e(d->eff,d->k,d->k*(i-d->S),d->Y)*w(x->me[i][j])*u->me[i][j];
      fprintf(sp1,"%lf %lf\n",d->k*(i-d->S),Q(xt,ctt));

      V_FREE(ctt);
      V_FREE(xt);

    }

  fclose(sp1);

  FILE *sp2 = fopen("plot2.txt","w");

  for (int i=0;i<x->m;i++) 
    fprintf(sp2,"%lf %lf\n",d->k*(i-d->S),1e3*c(d->cat,d->k,d->k*(i-d->S)));

  fclose(sp2);

  char sbuffer[100];
  sprintf(sbuffer,"./plo > plotc_%s.pdf",label);
  system(sbuffer); 
  //free(sbuffer);

  int S = d->S;
  Real k = d->k;

  int lfi=0;

  for (int i=S;i<x->m;i++)
    {

      int terminator;
      if(!SGNM)
	terminator = J+i;
      else
	terminator = J;


      VEC *xt = v_get(terminator+1);
      get_row(x,i,xt);      
      VEC *ut = v_get(terminator+1);
      get_row(u,i,ut);
      
      if (!SGNM)
	{
	  xt = v_resize(xt,terminator+1);
	  ut = v_resize(ut,terminator+1);
	}
      
      VEC *v = v_get(terminator+1);

      for (int j=0;j<=terminator;j++)
	v->ve[j] = iota*e(d->eff,k,k*(i-S),d->Y)*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];

      if(lfi < d->n && d->t_id[lfi]==i) 
	{
      
	  VEC *dt = v_get(d->t_sz[lfi]);

	  for (int j=0;j<dt->dim;j++)
	    dt->ve[j] = d->lf[lfi][j];

	  Real bw = get_bw(dt);

	  VEC *l = v_get(terminator+1);

	  for (int j=0;j<=terminator;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  Real al = 1e3*c(d->cat,k,k*(i - S)) / Q(xt,l); 

	  FILE *p1 = fopen("plot1.txt","w");

	  for (int j=0;j<=terminator;j++)
	    fprintf(p1,"%lf %lf\n",xt->ve[j],al*l->ve[j]);

	  fclose(p1);

	  FILE *p2 = fopen("plot2.txt","w");

	  for (int j=0;j<=terminator;j++)
	    fprintf(p2,"%lf %lf\n",xt->ve[j],v->ve[j]);

	  fclose(p2);

	  char buffer[100];

	  sprintf(buffer,"./plo > plotl%.3lf_%s.pdf",k*(d->t_id[lfi] - S),label);	

	  system(buffer);

	  //free(buffer);

	  V_FREE(dt);

	  lfi += 1;

	}

      V_FREE(xt);
      V_FREE(ut);
      V_FREE(v);

    }

  M_FREE(core_args.x);
  M_FREE(core_args.u);

}
*/
