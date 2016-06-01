#include <pthread.h>
#include "../meschach/matrix.h"
#include "../common.h"
#include "spade_solve.h"
#include "alpha/grad_alpha.h"
#include "beta/grad_beta.h"
#include "gamma/grad_gamma.h"
#include "iota/grad_iota.h"
#include "kappa/grad_kappa.h"
#include "omega/grad_omega.h"
#include "objfns.h"

VEC *VMGMM(

		  VEC *theta,
		  Data *d,
		  VEC *g,
		  double *f

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

  *f = K(theta,d,&core_args);  

  MAT *p_a = m_get(I,J);  
  Grad_Args grad_alpha_args;    
  grad_alpha_args.theta = theta;
  grad_alpha_args.d = d;
  grad_alpha_args.g = g;
  grad_alpha_args.p = p_a;
  grad_alpha_args.eff = d->eff;
  grad_alpha_args.k = d->k;
  grad_alpha_args.S = d->S;   
  grad_alpha_args.core_args = &core_args;
  
  MAT *p_b = m_get(I,J);
  Grad_Args grad_beta_args;
  grad_beta_args.theta = theta;
  grad_beta_args.d = d;
  grad_beta_args.g = g;
  grad_beta_args.p = p_b;
  grad_beta_args.eff = d->eff;
  grad_beta_args.k = d->k;
  grad_beta_args.S = d->S;
  grad_beta_args.core_args = &core_args;

  MAT *p_g = m_get(I,J);
  Grad_Args grad_gamma_args;
  grad_gamma_args.theta = theta;
  grad_gamma_args.d = d;
  grad_gamma_args.g = g;
  grad_gamma_args.p = p_g;
  grad_gamma_args.eff = d->eff;
  grad_gamma_args.k = d->k;
  grad_gamma_args.S = d->S;
  grad_gamma_args.core_args = &core_args;   

  MAT *p_i = m_get(I,J);
  Grad_Args grad_iota_args;
  grad_iota_args.theta = theta;
  grad_iota_args.d = d;
  grad_iota_args.g = g;
  grad_iota_args.p = p_i;
  grad_iota_args.eff = d->eff;
  grad_iota_args.k = d->k;
  grad_iota_args.S = d->S;
  grad_iota_args.core_args = &core_args;   

  MAT *p_k = m_get(I,J);
  Grad_Args grad_kappa_args;
  grad_kappa_args.theta = theta;
  grad_kappa_args.d = d;
  grad_kappa_args.g = g;
  grad_kappa_args.p = p_k;
  grad_kappa_args.eff = d->eff;
  grad_kappa_args.k = d->k;
  grad_kappa_args.S = d->S;
  grad_kappa_args.core_args = &core_args;
  
  if (PTH)
    {

      pthread_t grad_alpha_thread;
      if (pthread_create(&grad_alpha_thread, NULL, (void *)grad_alpha, (void*)&grad_alpha_args)) {
	fprintf(stderr, "Error creating thread - grad_alpha_thread");
	exit(0);
      }

      pthread_t grad_beta_thread;
      if (pthread_create(&grad_beta_thread, NULL, (void *)grad_beta, (void*)&grad_beta_args)) {
	fprintf(stderr, "Error creating thread - grad_beta_thread");
	exit(0);
      }

      pthread_t grad_gamma_thread;
      if (pthread_create(&grad_gamma_thread, NULL, (void *)grad_gamma, (void*)&grad_gamma_args)) {
    fprintf(stderr, "Error creating thread - grad_gamma_thread");
    exit(0);
      }

      pthread_t grad_iota_thread;
      if (pthread_create(&grad_iota_thread, NULL, (void *)grad_iota, (void*)&grad_iota_args)) { 
	fprintf(stderr, "Error creating thread - grad_iota_thread");
	exit(0);
      }

      pthread_t grad_kappa_thread;
      if (pthread_create(&grad_kappa_thread, NULL, (void *)grad_kappa, (void*)&grad_kappa_args)) { 
	fprintf(stderr, "Error creating thread - grad_kappa_thread");
	exit(0);
      }

      // await all threads

      pthread_join(grad_alpha_thread, NULL);
      pthread_join(grad_beta_thread, NULL);
      pthread_join(grad_gamma_thread, NULL);
      pthread_join(grad_iota_thread, NULL);
      pthread_join(grad_kappa_thread, NULL);

    }
  else
    {
      grad_alpha((void*)&grad_alpha_args);
      grad_beta((void*)&grad_beta_args);
      grad_gamma((void*)&grad_gamma_args);
      grad_iota((void*)&grad_iota_args);
      grad_kappa((void*)&grad_kappa_args);
    }  


  //printf("%g %g",grad->ve[4],ng);
  //exit(1);

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
  M_FREE(p_k);
  M_FREE(p_i);

  return g;

}
