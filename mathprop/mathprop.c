#include "../meschach/matrix.h"
#include "../meschach/matrix2.h"
#include "../common.h"
#include "../machinery/spade_solve.h"
#include "../machinery/alpha/grad_alpha.h"
#include "../machinery/beta/grad_beta.h"
#include "../machinery/gamma/grad_gamma.h"
#include "../machinery/iota/grad_iota.h"
#include "../machinery/kappa/grad_kappa.h"
#include "../machinery/omega/grad_omega.h"
#include "../machinery/objfns.h"



VEC *numgrad(

	     double (*model)(VEC *,void *),
	     void *stuff,
	     VEC *par,
	     double epsilon

	     )
{

  double f0 = (*model)(par,stuff);

  VEC *ei = v_get(par->dim);
  VEC *d = v_get(par->dim);
  VEC *fg = v_get(par->dim);

  for (int i=0;i<par->dim;i++)
    {
      ei->ve[i] = 1;
      sv_mlt(epsilon,ei,d);
      double f1 = (*model)(v_add(par,d,VNULL),stuff);
      fg->ve[i] = (f1 - f0) / d->ve[i];
      ei->ve[i] = 0;
    }

  return fg;

}

/*
double ConditionNumber(

		       VEC * theta,
		       Data *dataptr

		       )
{

  int n = theta->dim;

  VEC * grad = v_get(n);

  int I = dataptr->I+1;
  int J = dataptr->J+1;
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

  solve(theta,x,u,xhh,xh,xn,uh,un,Ui,Uh,Uhh,idxi,dataptr->eff,dataptr->k,dataptr->S);
  double f = H(x,u,dataptr,theta->ve[3]);

  // p alpha
  MAT *p_a = m_get(x->m,x->n);
  Solve_Args grad_alpha_args;
  grad_alpha_args.dataptr = dataptr;
  grad_alpha_args.grad = grad;
  grad_alpha_args.p = p_a;
  grad_alpha_args.x = x;
  grad_alpha_args.u = u;
  grad_alpha_args.xhh = xhh;
  grad_alpha_args.xh = xh;
  grad_alpha_args.xn = xn;
  grad_alpha_args.uh = uh;
  grad_alpha_args.un = un;
  grad_alpha_args.Ui = Ui;
  grad_alpha_args.Uh = Uh;
  grad_alpha_args.Uhh = Uhh;
  grad_alpha_args.idxi = idxi;
  grad_alpha_args.eff = dataptr->eff;
  grad_alpha_args.k = dataptr->k;
  grad_alpha_args.S = dataptr->S;

  // p beta
  MAT *p_b = m_get(x->m,x->n);
  Solve_Args grad_beta_args;
  grad_beta_args.theta = theta;
  grad_beta_args.dataptr = dataptr;
  grad_beta_args.grad = grad;
  grad_beta_args.p = p_b;
  grad_beta_args.x = x;
  grad_beta_args.u = u;
  grad_beta_args.xhh = xhh;
  grad_beta_args.xh = xh;
  grad_beta_args.xn = xn;
  grad_beta_args.uh = uh;
  grad_beta_args.Ui = Ui;
  grad_beta_args.Uh = Uh;
  grad_beta_args.Uhh = Uhh;
  grad_beta_args.idxi = idxi;
  grad_beta_args.eff = dataptr->eff;
  grad_beta_args.k = dataptr->k;
  grad_beta_args.S = dataptr->S;

  // p gamma
  MAT *p_g = m_get(x->m,x->n);
  Solve_Args grad_gamma_args;
  grad_gamma_args.theta = theta;
  grad_gamma_args.dataptr = dataptr;
  grad_gamma_args.grad = grad;
  grad_gamma_args.p = p_g;
  grad_gamma_args.x = x;
  grad_gamma_args.u = u;
  grad_gamma_args.xhh = xhh;
  grad_gamma_args.xh = xh;
  grad_gamma_args.xn = xn;
  grad_gamma_args.uh = uh;
  grad_gamma_args.Ui = Ui;
  grad_gamma_args.Uh = Uh;
  grad_gamma_args.Uhh = Uhh;
  grad_gamma_args.idxi = idxi;
  grad_gamma_args.eff = dataptr->eff;
  grad_gamma_args.k = dataptr->k;
  grad_gamma_args.S = dataptr->S;
 
  // p iota
  MAT *p_i = m_get(x->m,x->n);
  Solve_Args grad_iota_args;
  grad_iota_args.theta = theta;
  grad_iota_args.dataptr = dataptr;
  grad_iota_args.grad = grad;
  grad_iota_args.p = p_i;
  grad_iota_args.x = x;
  grad_iota_args.u = u;
  grad_iota_args.xhh = xhh;
  grad_iota_args.xh = xh;
  grad_iota_args.xn = xn;
  grad_iota_args.uh = uh;
  grad_iota_args.Ui = Ui;
  grad_iota_args.Uh = Uh;
  grad_iota_args.Uhh = Uhh;
  grad_iota_args.idxi = idxi;
  grad_iota_args.eff = dataptr->eff;
  grad_iota_args.k = dataptr->k;
  grad_iota_args.S = dataptr->S;

  // p kappa
  MAT *p_k = m_get(x->m,x->n);
  Solve_Args grad_kappa_args;
  grad_kappa_args.theta = theta;
  grad_kappa_args.dataptr = dataptr;
  grad_kappa_args.grad = grad;
  grad_kappa_args.p = p_k;
  grad_kappa_args.x = x;
  grad_kappa_args.u = u;
  grad_kappa_args.xhh = xhh;
  grad_kappa_args.xh = xh;
  grad_kappa_args.xn = xn;
  grad_kappa_args.uh = uh;
  grad_kappa_args.un = un;
  grad_kappa_args.Ui = Ui;
  grad_kappa_args.Uh = Uh;
  grad_kappa_args.Uhh = Uhh;
  grad_kappa_args.idxi = idxi;
  grad_kappa_args.eff = dataptr->eff;
  grad_kappa_args.k = dataptr->k;
  grad_kappa_args.S = dataptr->S;

  VEC *theta_save = v_get(n);
  
  theta_save->ve[0]=theta->ve[0];
  theta_save->ve[1]=theta->ve[1];
  theta_save->ve[2]=theta->ve[2];
  theta_save->ve[3]=theta->ve[3];
  theta_save->ve[4]=theta->ve[4];

  MAT *hessin = m_get(n,n);

  double epsilon = 1e-7;

  VEC *delta = v_get(n);
  
  MAT *gp = m_get(n,n);
  MAT *gn = m_get(n,n);
  
  for (int i=0;i<n;i++)
    {

      VEC *delta = v_get(n);
      
      delta->ve[i] = epsilon;
      v_add(theta_save,delta,theta);

      grad_alpha_args.theta = theta;
      grad_beta_args.theta = theta;
      grad_gamma_args.theta = theta;
      grad_iota_args.theta = theta;
      grad_kappa_args.theta = theta;

      grad_alpha((void*)&grad_alpha_args);
      grad_beta((void*)&grad_beta_args);
      grad_gamma((void*)&grad_gamma_args);
      grad_iota((void*)&grad_iota_args);
      grad_kappa((void*)&grad_kappa_args);
  
      vm_move(grad,0,gp,i,0,1,n);
       
      delta->ve[i] = -epsilon;
      v_add(theta_save,delta,theta);

      grad_alpha_args.theta = theta;
      grad_beta_args.theta = theta;
      grad_gamma_args.theta = theta;
      grad_iota_args.theta = theta;
      grad_kappa_args.theta = theta;

      grad_alpha((void*)&grad_alpha_args);
      grad_beta((void*)&grad_beta_args);
      grad_gamma((void*)&grad_gamma_args);
      grad_iota((void*)&grad_iota_args);
      grad_kappa((void*)&grad_kappa_args);
 
      vm_move(grad,0,gn,i,0,1,n);  
            
      V_FREE(delta);
     }
  
  for (int i=0;i<n;i++)
    for (int j=0;j<n;j++)
      hessin->me[i][j] = (gp->me[i][j] - gn->me[i][j])/(4*epsilon) + (gp->me[j][i] - gn->me[j][i])/(4*epsilon);

  VEC * evals = v_get(n);
  MAT *Q = m_get(n,n);

  evals = symmeig(hessin,Q,evals);

  v_output(evals);

  M_FREE(hessin);
  M_FREE(p_a);
  M_FREE(p_b);
  M_FREE(p_g);
  M_FREE(p_i);
  M_FREE(x);
  M_FREE(u);
  M_FREE(xh);
  M_FREE(uh);
  M_FREE(xn);
  M_FREE(xhh);
  M_FREE(un);
  V_FREE(Ui);
  V_FREE(Uh);
  V_FREE(Uhh);
  IV_FREE(idxi);

  return evals->ve[4]/evals->ve[0];

}*/



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
