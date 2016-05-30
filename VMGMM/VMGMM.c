#include <pthread.h>
#include "../meschach/matrix.h"
#include "../common.h"
#include "solvers/spade_solve.h"
#include "solvers/alpha/solve_p_alpha.h"
#include "solvers/beta/solve_p_beta.h"
#include "solvers/gamma/solve_p_gamma.h"
#include "solvers/iota/solve_p_iota.h"
#include "solvers/kappa/solve_p_kappa.h"
#include "solvers/omega/solve_p_omega.h"
#include "objfns/objfns.h"
#include "../common.h"

VEC *VMGMM(

		  VEC *theta,
		  struct DATA *dataptr,
		  VEC *grad,
		  double *f

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

  //printf("%g ",*f);

  MAT *p_a = m_get(x->m,x->n);
  // p alpha
  Solve_Args solve_p_alpha_args;
  solve_p_alpha_args.theta = theta;
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
  //solve_p_beta_args.un = un;
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
  //solve_p_gamma_args.un = un;
  solve_p_gamma_args.Ui = Ui;
  solve_p_gamma_args.Uh = Uh;
  solve_p_gamma_args.Uhh = Uhh;
  solve_p_gamma_args.idxi = idxi;
  solve_p_gamma_args.eff = dataptr->eff;
  solve_p_gamma_args.k = dataptr->k;
  solve_p_gamma_args.S = dataptr->S;
 
  /*
  MAT *p_k = m_get(x->m,x->n);
  solve_p_kappa (theta,p_k ,x,u,xhh,xh,xn,uh,un,Ui,Uh,Uhh,idxi,dataptr->eff,dataptr->k,dataptr->S);
  grad->ve[4] = G_ni(p_k,x,u,dataptr,theta->ve[6]);

  MAT *p_w = m_get(x->m,x->n);
  solve_p_omega (theta,p_w ,x,u,xhh,xh,xn,uh,un,Ui,Uh,Uhh,idxi,dataptr->eff,dataptr->k,dataptr->S);
  grad->ve[5] = G_ni(p_w,x,u,dataptr,theta->ve[6]);
  */

  // p iota
  MAT *p_i = m_get(x->m,x->n);
  Solve_Args solve_p_iota_args;
  solve_p_iota_args.theta = theta;
  solve_p_iota_args.dataptr = dataptr;
  solve_p_iota_args.grad = grad;
  solve_p_iota_args.grad = grad;
  solve_p_iota_args.dataptr = dataptr;
  solve_p_iota_args.p = p_i;
  solve_p_iota_args.x = x;
  solve_p_iota_args.u = u;
  solve_p_iota_args.xhh = xhh;
  solve_p_iota_args.xh = xh;
  solve_p_iota_args.xn = xn;
  solve_p_iota_args.uh = uh;
  //solve_p_iota_args.un = un;
  solve_p_iota_args.Ui = Ui;
  solve_p_iota_args.Uh = Uh;
  solve_p_iota_args.Uhh = Uhh;
  solve_p_iota_args.idxi = idxi;
  solve_p_iota_args.eff = dataptr->eff;
  solve_p_iota_args.k = dataptr->k;
  solve_p_iota_args.S = dataptr->S;

  if (PTH)
    {

      pthread_t solve_p_alpha_thread;
      if (pthread_create(&solve_p_alpha_thread, NULL, (void *)solve_p_alpha, (void*)&solve_p_alpha_args)) {
	fprintf(stderr, "Error creating thread - solve_p_alpha_thread");
	exit(0);
      }

      pthread_t solve_p_beta_thread;
      if (pthread_create(&solve_p_beta_thread, NULL, (void *)solve_p_beta, (void*)&solve_p_beta_args)) {
	fprintf(stderr, "Error creating thread - solve_p_beta_thread");
	exit(0);
      }

      pthread_t solve_p_gamma_thread;
      if (pthread_create(&solve_p_gamma_thread, NULL, (void *)solve_p_gamma, (void*)&solve_p_gamma_args)) {
    fprintf(stderr, "Error creating thread - solve_p_gamma_thread");
    exit(0);
      }

      pthread_t solve_p_iota_thread;
      if (pthread_create(&solve_p_iota_thread, NULL, (void *)solve_p_iota, (void*)&solve_p_iota_args)) { 
	fprintf(stderr, "Error creating thread - solve_p_iota_thread");
	exit(0);
      }

      // await all threads

      pthread_join(solve_p_alpha_thread, NULL);
      pthread_join(solve_p_beta_thread, NULL);
      pthread_join(solve_p_gamma_thread, NULL);
      pthread_join(solve_p_iota_thread, NULL);

    }
  else
    {
      solve_p_alpha((void*)&solve_p_alpha_args);
      solve_p_beta((void*)&solve_p_beta_args);
      solve_p_gamma((void*)&solve_p_gamma_args);
      solve_p_iota((void*)&solve_p_iota_args);
    }  

  // p alpha debugging
  //  printf("%g ", grad->ve[0]);

  // p beta debugging
  //printf("%g ", grad->ve[1]);

  // p gamma debugging
  //printf("%g ", grad->ve[2]);

  // p iota debugging
  //printf("%g ", grad->ve[3]);

  //printf("%g ",theta->ve[0]);  printf("%g ",theta->ve[1]);  printf("%g ",theta->ve[2]);  printf("%g\n",theta->ve[3]);

  M_FREE(p_a);
  M_FREE(p_b);
  M_FREE(p_g);
  //M_FREE(p_k);
  //M_FREE(p_w);
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
  //V_FREE(ctt);
  //V_FREE(xt);

  return grad;

}
