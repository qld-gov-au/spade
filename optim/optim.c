// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#include <math.h>
#include "../meschach/matrix.h"
#include "../common.h"
#include "../parameters.h"
#include "optim.h"
#include "more-thuente-oleary.h"

VEC * bfgs(

	   VEC * (*model)(VEC *,Data *,VEC *,Real *,Parameters *),
	   VEC *x,
	   Data *data,
	   Parameters * parameters,
     OptimControl opt

	   )   
{

  int n = x->dim;
  int nfev = 0;

  Real f;
  Real fv;

  VEC *s = v_get(n);
  VEC *y = v_get(n);
  VEC *g = v_get(n);
  VEC *p = v_get(n);
  VEC *u = v_get(n);

  VEC *oldx = v_get(n);
  VEC *oldg = v_get(n);

  MAT *B = m_get(n,n);

  m_ident(B);

  g = (*model)(x,data,g,&f,parameters);
  //exit(0);
  nfev += 1;
/*
  Real dx = 1e-7;
  x->ve[3] += dx;
  parameters->parameter[4]->value = x->ve[3];

  Real fnew;
  g = (*model)(x,data,g,&fnew,parameters);

  Real dy = fnew-f;
  Real ng = dy / dx;
  printf("%g %g\n",ng,g->ve[3]);
  
  exit(0);
*/
  fv = f;

  int count = 0;

  while (1)
    {

    #if REAL == DOUBLE
         printf("infnorm: %g fv: %lf a: %lf b: %lf g: %lf i: %lf k: %lf w: %lf\n",v_norm_inf(g),fv,parameters->alpha.value,parameters->beta.value,parameters->gamma.value,parameters->iota.value,parameters->kappa.value,parameters->omega.value);
    #elif REAL == FLOAT
         printf("infnorm: %g fv: %f a: %f b: %f g: %f i: %f k: %f w: %f\n",v_norm_inf(g),fv,parameters->alpha.value,parameters->beta.value,parameters->gamma.value,parameters->iota.value,parameters->kappa.value,parameters->omega.value);
    #elif REAL == LONGDOUBLE
         printf("infnorm: %g fv: %Lf a: %Lf b: %Lf g: %Lf i: %Lf k: %Lf w: %Lf\n",v_norm_inf(g),fv,parameters->alpha.value,parameters->beta.value,parameters->gamma.value,parameters->iota.value,parameters->kappa.value,parameters->omega.value);
    #endif

      if (v_norm_inf(g) < 1)
	break;

      mv_mlt(B,g,p);      
      sv_mlt(-1./v_norm2(p),p,p);

      //VEC * p_copy = v_get(n);
      //v_copy(p,p_copy);
      /*
      if (count==0) {      
      
      
        printf("\n");
        //for (int i=-20;i<=-3;i++) {       
        for (int i=-10;i<=10;i++) {       
        
          Real stp = (Real)i*0.0001; //exp((Real)i);
          //Real stp = (Real)i*0.0001; //exp((Real)i);

          Real nfv;
          
          VEC *newx = v_get(n);
          VEC *ptmp = v_get(n);
          ptmp = sv_mlt(stp,p,ptmp);
          v_add(x,ptmp,newx);         
          
          //v_copy(x,newx);
          //newx->ve[4] = newx->ve[4] + stp;
          
          //v_output(newx);
          V_FREE(ptmp);
          
          // apply changes to theta (x) back to parameters
          int iTheta = 0;
          for(int i = 0; i < parameters->count; i++) {
            if(parameters->parameter[i]->active == TRUE) {
              parameters->parameter[i]->value = newx->ve[iTheta++];
            }
          }

          //newx->ve[3] += stp;
          //parameters->kappa.value = newx->ve[3];
          //parameter[4]->value = newx->ve[3];

          g = (*model)(newx,data,g,&nfv,parameters);
          
          //Real dg = in_prod(g,p);
          //printf("dg: %f\n",dg);
          //printf("%g %g %g %g %g\n",stp,nfv,dg,newx->ve[4],g->ve[4]);
          printf("%Lf %.15Lf\n",stp,nfv); // newx->ve[3],g->ve[3]);
          //printf("%g %g\n",newx->ve[4],g->ve[4]); //,dg,newx->ve[4],g->ve[4]);

              
          V_FREE(newx);
             
        }

        exit(0);
        
      }*/
                

      v_copy(x,oldx);
      v_copy(g,oldg);

      // More-Thuente line search
      nfev += cvsrch(model,x,f,g,p,opt.stp,opt.ftol,opt.gtol,opt.xtol,opt.stpmin,opt.stpmax,opt.maxfev,data,parameters,&fv);
    
      v_sub(x,oldx,s);
      v_sub(g,oldg,y);

      // update inverse hessian based on julia code
      Real sy = in_prod(s,y);

      if (sy==0)
	{
	  printf("cannot find a stepsize that reduces the function along the descent direction\n");
	  exit(1);
	}

      mv_mlt(B,y,u);

      Real yBy = in_prod(u,y);
      Real c1 = (sy + yBy) / (sy*sy);
      Real c2 = 1/sy;

      // not using meschach for this - overcomplicates it.
      for (int i=0;i<n;i++)
	    for (int j=0;j<n;j++)
	      B->me[i][j] += c1 * s->ve[i] * s->ve[j] - c2 * ( u->ve[i] * s->ve[j] + u->ve[j] * s->ve[i] );

    //break;
    
      count++;
    
    }

  (*model)(x,data,g,&f,parameters);
  #if REAL == DOUBLE
    // lf
    printf("number function evals: %d, function value: %lf\n",nfev,f);
  #elif REAL == FLOAT
    // f
    printf("number function evals: %d, function value: %f\n",nfev,f);
  #elif REAL == LONGDOUBLE
    // Lf
    printf("number function evals: %d, function value: %Lf\n",nfev,f);
  #endif

  v_output(x);

  V_FREE(oldx);
  V_FREE(oldg);
  V_FREE(s);
  V_FREE(y);
  V_FREE(g);
  M_FREE(B);
  V_FREE(u);
  V_FREE(p);

  return(x);

}


