#include "../meschach/matrix.h"
#include "../meschach/matrix2.h"
#include "../common.h"
#include "../parameters.h"
#include "../machinery/spade_solve.h"
#include "../machinery/alpha/grad_alpha.h"
#include "../machinery/beta/grad_beta.h"
#include "../machinery/gamma/grad_gamma.h"
#include "../machinery/iota/grad_iota.h"
#include "../machinery/kappa/grad_kappa.h"
#include "../machinery/omega/grad_omega.h"
#include "../machinery/objfns.h"

VEC *numgrad(

	     Real (*model)(VEC *,void *),
	     void *stuff,
	     VEC *par,
	     Real epsilon

	     )
{

  Real f0 = (*model)(par,stuff);

  VEC *ei = v_get(par->dim);
  VEC *d = v_get(par->dim);
  VEC *fg = v_get(par->dim);

  for (int i=0;i<par->dim;i++)
    {
      ei->ve[i] = 1;
      sv_mlt(epsilon,ei,d);
      Real f1 = (*model)(v_add(par,d,VNULL),stuff);
      fg->ve[i] = (f1 - f0) / d->ve[i];
      ei->ve[i] = 0;
    }

  return fg;

}


Real ConditionNumber(

		       Parameters *parameters,
		       Data *d

		       )
{

  Solve_Core_Args core_args;
  
  int I,J;
  I = d->I;
  J = d->J;
  
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

  Real fv = K(parameters,d,&core_args);
  
  // get the active parameters and run their grad functions
  Grad_Args args;
  args.d = d;
  args.eff = d->eff;
  args.k = d->k;
  args.S = d->S;
  args.core_args = &core_args;
  args.parameters = parameters;
  
  Real epsilon = 1e-6;
  
  // Determine the number of parameters which are active
  int activeParameterCount = 0;
  for(int i = 0; i < parameters->count; i++) {
    if(parameters->parameter[i]->active == TRUE) {
      activeParameterCount++;
    }
  }
  int n = activeParameterCount;
  
  MAT *tmp1 = m_get(n,n);
  MAT *tmp2 = m_get(n,n);

  int iTheta=0;
  
  for(int i = 0; i < parameters->count; i++) {
    
    if(parameters->parameter[i]->active == TRUE) {
    
      parameters->parameter[i]->value += epsilon;

      int jTheta=0;  
      for (int j = 0; j < parameters->count; j++) {
        if(parameters->parameter[j]->active == TRUE) {
        
          parameters->parameter[j]->grad((void *) &(args));         
          tmp1->me[i][j] = parameters->parameter[j]->gradient;
          jTheta++;
          
        }
      }
      
      parameters->parameter[i]->value -= 2*epsilon;
      jTheta=0;

      for (int j = 0; j < parameters->count; j++) {
        if(parameters->parameter[j]->active == TRUE) {
        
          parameters->parameter[j]->grad((void *) &(args));         
          tmp2->me[i][j] = parameters->parameter[j]->gradient;
          jTheta++;
          
        }
      }
      
      parameters->parameter[i]->value += epsilon;
  
      iTheta++;
  
    }
  }

  MAT *hessin = m_get(n,n);
  
  for (int i=0;i<n;i++)
    for (int j=0;j<n;j++)
      hessin->me[i][j] = (tmp1->me[i][j] - tmp2->me[i][j])/(4*epsilon) + (tmp1->me[j][i] - tmp2->me[j][i])/(4*epsilon);

  VEC * evals = v_get(n);
  MAT *Q = m_get(n,n);

  evals = symmeig(hessin,Q,evals);

  v_output(evals);
  
  M_FREE(hessin);
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

  return evals->ve[0]/evals->ve[n-1];

}



/*
void output_plots(

		  VEC *p,
		  Data *d,
		  char *label

		  )
{

  int I = d->I+1;
  int J = d->J+1;
  MAT *x = m_get(I,J);
  MAT *u = m_get(I,J);
  MAT *xh = m_get(I,J+1);
  MAT *uh = m_get(I,J+1);
  MAT *xn = m_get(I,J+1);
  MAT *xhh = m_get(I,J+1);
  MAT *un = m_get(I,J+1);
  VEC *Ui = v_get(I);
  VEC *Uh = v_get(I);
  VEC *Uhh = v_get(I);
  IVEC *idxi = iv_get(I-1);

  solve(theta,x,u,xhh,xh,xn,uh,un,Ui,Uh,Uhh,idxi,d->eff,d->k,d->S);

  //  VEC *ctt = v_get(x->n);
  //  VEC *xt = v_get(x->n);

  /*
  FILE *p1 = fopen("plot1.txt","w");

  for (int i=0;i<x->m;i++)
    {

      xt = get_row(x,i,xt);
      for (int j=0;j<x->n;j++)
	ctt->ve[j] = s(x->me[i][j])*theta->ve[6]*e(d->eff,d->k,d->k*(i-d->S))*w(x->me[i][j])*u->me[i][j];
      fprintf(p1,"%f %f\n",d->k*(i-d->S),Q(xt,ctt)/1e3);

    }

  fclose(p1);

  FILE *p2 = fopen("plot2.txt","w");

  for (int i=0;i<x->m;i++) 
    fprintf(p2,"%f %f\n",d->k*(i-d->S),c(d->cat,d->k,d->k*(i-d->S)));

  fclose(p2);

  system("./plo > plotc.pdf");

  exit(1);
  */
