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
		  double *f,
      Parameters * parameters

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

  *f = K(parameters,d,&core_args);

  // get the active parameters and run their grad functions
  MAT * matrices[theta->dim];
  Grad_Args args[theta->dim];
  pthread_t threads[theta->dim];

  int iTheta = 0;
  for(int i = 0; i < parameters->count; i++) {
    if(parameters->parameter[i]->active == TRUE) {
      matrices[iTheta] = m_get(I, J);
      args[iTheta].d = d;
      args[iTheta].g = g;
      args[iTheta].p = matrices[iTheta];
      args[iTheta].eff = d->eff;
      args[iTheta].k = d->k;
      args[iTheta].S = d->S;
      args[iTheta].core_args = &core_args;
      args[iTheta].parameters = parameters;

      iTheta++;
    }
  }

  // multi-threaded mode
  if(PTH) {
    // launch a thread for each gradient function
    int iTheta = 0;
    for(int i = 0; i < parameters->count; i++) {
      if(parameters->parameter[i]->active == TRUE) {
        if(pthread_create(&(threads[iTheta]), NULL, (void *) parameters->parameter[i]->grad, (void *) &(args[iTheta]))) {
          fprintf(stderr, "Error creating thread");
          exit(0);
        }
        iTheta++;
      }
    }

    // await all threads
    for(int i = 0; i < theta-> dim; i++) {
      pthread_join(threads[i], NULL);
    }
  }

  // single-threaded mode
  else {
    for(int i = 0; i < parameters->count; i++) {
      if(parameters->parameter[i]->active == TRUE) {
        parameters->parameter[i]->grad((void *) &(args[iTheta]));
      }
    }
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

  for(int i = 0; i < theta->dim; i++) {
    M_FREE(matrices[i]);
  }

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

  return g;

}
