#include "../meschach/matrix.h"

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
		       struct DATA *dataptr

		       )
{

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
  *f = H(x,u,dataptr,theta->ve[3]);

  // p alpha
  MAT *p_a = m_get(x->m,x->n);
  Solve_Args solve_p_alpha_args;
  solve_p_alpha_args.dataptr = dataptr;
  solve_p_alpha_args.grad = grad;
  solve_p_alpha_args.p = p_a;
  solve_p_alpha_args.x = x;
  solve_p_alpha_args.u = u;
  solve_p_alpha_args.xhh = xhh;
  solve_p_alpha_args.xh = xh;
  solve_p_alpha_args.xn = xn;
  solve_p_alpha_args.uh = uh;
  solve_p_alpha_args.un = un;
  solve_p_alpha_args.Ui = Ui;
  solve_p_alpha_args.Uh = Uh;
  solve_p_alpha_args.Uhh = Uhh;
  solve_p_alpha_args.idxi = idxi;
  solve_p_alpha_args.eff = dataptr->eff;
  solve_p_alpha_args.k = dataptr->k;
  solve_p_alpha_args.S = dataptr->S;

  // p beta
  MAT *p_b = m_get(x->m,x->n);
  Solve_Args solve_p_beta_args;
  solve_p_beta_args.theta = theta;
  solve_p_beta_args.dataptr = dataptr;
  solve_p_beta_args.grad = grad;
  solve_p_beta_args.p = p_b;
  solve_p_beta_args.x = x;
  solve_p_beta_args.u = u;
  solve_p_beta_args.xhh = xhh;
  solve_p_beta_args.xh = xh;
  solve_p_beta_args.xn = xn;
  solve_p_beta_args.uh = uh;
  solve_p_beta_args.Ui = Ui;
  solve_p_beta_args.Uh = Uh;
  solve_p_beta_args.Uhh = Uhh;
  solve_p_beta_args.idxi = idxi;
  solve_p_beta_args.eff = dataptr->eff;
  solve_p_beta_args.k = dataptr->k;
  solve_p_beta_args.S = dataptr->S;

  // p gamma
  MAT *p_g = m_get(x->m,x->n);
  Solve_Args solve_p_gamma_args;
  solve_p_gamma_args.theta = theta;
  solve_p_gamma_args.dataptr = dataptr;
  solve_p_gamma_args.grad = grad;
  solve_p_gamma_args.p = p_g;
  solve_p_gamma_args.x = x;
  solve_p_gamma_args.u = u;
  solve_p_gamma_args.xhh = xhh;
  solve_p_gamma_args.xh = xh;
  solve_p_gamma_args.xn = xn;
  solve_p_gamma_args.uh = uh;
  solve_p_gamma_args.Ui = Ui;
  solve_p_gamma_args.Uh = Uh;
  solve_p_gamma_args.Uhh = Uhh;
  solve_p_gamma_args.idxi = idxi;
  solve_p_gamma_args.eff = dataptr->eff;
  solve_p_gamma_args.k = dataptr->k;
  solve_p_gamma_args.S = dataptr->S;
 
  // p iota
  MAT *p_i = m_get(x->m,x->n);
  Solve_Args solve_p_iota_args;
  solve_p_iota_args.theta = theta;
  solve_p_iota_args.dataptr = dataptr;
  solve_p_iota_args.grad = grad;
  solve_p_iota_args.p = p_i;
  solve_p_iota_args.x = x;
  solve_p_iota_args.u = u;
  solve_p_iota_args.xhh = xhh;
  solve_p_iota_args.xh = xh;
  solve_p_iota_args.xn = xn;
  solve_p_iota_args.uh = uh;
  solve_p_iota_args.Ui = Ui;
  solve_p_iota_args.Uh = Uh;
  solve_p_iota_args.Uhh = Uhh;
  solve_p_iota_args.idxi = idxi;
  solve_p_iota_args.eff = dataptr->eff;
  solve_p_iota_args.k = dataptr->k;
  solve_p_iota_args.S = dataptr->S;

  VEC *theta_save = v_get(theta->dim);
  theta_save->ve[0]=theta->ve[0];
  theta_save->ve[1]=theta->ve[1];
  theta_save->ve[2]=theta->ve[2];
  theta_save->ve[3]=theta->ve[3];

  MAT *H = m_get(4,4);

  double epsilon = 1e-7;

  for (int j=1;j<=theta->dim;j++)
    {

      VEC *delta = v_get(theta->dim);
      delta->ve[j] = epsilon;

      v_add(theta_save,delta,theta);

      solve_p_alpha_args.theta = theta;
      solve_p_beta_args.theta = theta;
      solve_p_gamma_args.theta = theta;
      solve_p_iota_args.theta = theta;

      solve_p_alpha((void*)&solve_p_alpha_args);
      solve_p_beta((void*)&solve_p_beta_args);
      solve_p_gamma((void*)&solve_p_gamma_args);
      solve_p_iota((void*)&solve_p_iota_args);

      VEC *g1 = v_get(4);

      g1 = grad;

      delta->ve[j] = -epsilon;

      v_add(theta_save,delta,theta);

      solve_p_alpha_args.theta = theta;
      solve_p_beta_args.theta = theta;
      solve_p_gamma_args.theta = theta;
      solve_p_iota_args.theta = theta;

      solve_p_alpha((void*)&solve_p_alpha_args);
      solve_p_beta((void*)&solve_p_beta_args);
      solve_p_gamma((void*)&solve_p_gamma_args);
      solve_p_iota((void*)&solve_p_iota_args);

      VEC *g2 = v_get(4);
      g2 = grad;

      g1 = grad;

      H->me[i][j];

      }

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

  return cn;

}
*/


/*
void output_plots(

		  VEC *p,
		  struct DATA *d,
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
