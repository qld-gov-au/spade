#include <pthread.h>
#include "../meschach/matrix.h"
#include "../common.h"
#include "../parameters.h"
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
		  Real *f,
      Parameters * parameters

		  )
{

  Solve_Core_Args core_args;
  
  int I,J;
  I = d->I;
  J = d->J;

  core_args.x = m_get(I+1,J+1);
  core_args.u = m_get(I+1,J+1);

  // todo: Review this. Condition has been flipped to maintain same functionality as last commit
  // in SGNM mode but may not be correct for non-SGNM mode. Affects xh in particular
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

  *f = K(parameters,d,&core_args);

  // get the active parameters and run their grad functions
  Grad_Args args[theta->dim];
  pthread_t threads[theta->dim];

  int iTheta = 0;
  for(int i = 0; i < parameters->count; i++) {
    if(parameters->parameter[i]->active == TRUE) {
      args[iTheta].d = d;
      args[iTheta].eff = d->eff;
      args[iTheta].k = d->k;
      args[iTheta].S = d->S;
      args[iTheta].core_args = &core_args;
      args[iTheta].parameters = parameters;

      iTheta++;
    }
  }

  parameters->parameter[4]->grad((void *) &(args[4]));

  printf("analytic: %Lf\n",parameters->parameter[4]->gradient);

  Real par_save = parameters->parameter[4]->value;
  
  for (int i=-4;i>=-20;i--) {

    parameters->parameter[4]->value = par_save + exp((Real)i);

    Real dY = K(parameters,d,&core_args) - *f;

    Real dX = exp((Real)i);
    
    printf("%Lf %Lf\n",dX,dY/dX); 

  }

  exit(1);
  
  
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
    int iTheta = 0;
    for(int i = 0; i < parameters->count; i++) {
      if(parameters->parameter[i]->active == TRUE) {
        parameters->parameter[i]->grad((void *) &(args[iTheta]));
      }
    }
  }

  // load parameters[i].gradient back into g
  iTheta=0;
  for(int i = 0; i < parameters->count; i++) 
    if(parameters->parameter[i]->active == TRUE) {
      g->ve[iTheta] = parameters->parameter[i]->gradient;
      iTheta++;
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
