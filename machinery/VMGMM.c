// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#include <pthread.h>
#include <math.h>
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
#include <math.h>
#include "VMGMM.h"

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

  if (MESCHACH)
    {
      
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
       if (QUARTER)
	 core_args.xhh = m_get(I+1,J+1);
       core_args.un = m_get(I+1,J+1);
    }

  core_args.Ui = v_get(I+1);
  core_args.Uh = v_get(I+1);
  if (QUARTER)
    core_args.Uhh = v_get(I+1);
  core_args.idxi = iv_get(I);
    }
  
  *f = K(parameters,d,&core_args);
  
  // get the active parameters and run their grad functions
  Grad_Args args[theta->dim];
  Grad_Args_No_MESCHACH args_nm[theta->dim];
  
  pthread_t threads[theta->dim];

  int iTheta = 0;
  for(int i = 0; i < parameters->count; i++) {
    if(parameters->parameter[i]->active == TRUE) {

      if (MESCHACH)
	{
      args[iTheta].d = d;
      args[iTheta].eff = d->eff;
      args[iTheta].k = d->k;
      args[iTheta].S = d->S;
      args[iTheta].core_args = &core_args;
      args[iTheta].parameters = parameters;
	}
      else
	{
      args_nm[iTheta].d = d;
      args_nm[iTheta].eff = d->eff;
      args_nm[iTheta].k = d->k;
      args_nm[iTheta].S = d->S;
      args_nm[iTheta].parameters = parameters;
	}
      
      iTheta++;
    }
  }


  if (interactive_mode_requested)
    {
      if(MESCHACH){
        check_derivative(parameters, args, d, f, core_args);
      }
      else {
        check_derivative_no_meschach(parameters, args_nm, d, f, core_args);
      }

      interactive_mode_requested = 0;
    }

 
  // multi-threaded mode
  if(PTH) {
    // launch a thread for each gradient function
    int iTheta = 0;
    for(int i = 0; i < parameters->count; i++) {
      if(parameters->parameter[i]->active == TRUE) {
	if (MESCHACH)
	  {
	    if(pthread_create(&(threads[iTheta]), NULL, (void *) parameters->parameter[i]->grad, (void *) &(args[iTheta])))
	      {
		fprintf(stderr, "Error creating thread");
		exit(0);
	      }
	  }
	else
	  {
	    if(pthread_create(&(threads[iTheta]), NULL, (void *) parameters->parameter[i]->grad, (void *) &(args_nm[iTheta])))
	      {
		fprintf(stderr, "Error creating thread");
		exit(0);
	      }
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
	if (MESCHACH)
	  parameters->parameter[i]->grad((void *) &(args[iTheta]));
	else
	  parameters->parameter[i]->grad((void *) &(args_nm[iTheta]));
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

  if (MESCHACH)
    {
  M_FREE(core_args.xh);
  M_FREE(core_args.uh);
  M_FREE(core_args.xn);
  if (QUARTER)
    {
      M_FREE(core_args.xhh);
      V_FREE(core_args.Uhh);
    }
  M_FREE(core_args.un);
  V_FREE(core_args.Ui);
  V_FREE(core_args.Uh);
  IV_FREE(core_args.idxi);
    }
  
  return g;

}

void check_derivative(Parameters * parameters, Grad_Args * args, Data * d, Real * f, Solve_Core_Args core_args) {
  // Prompt user for the lower bound, upper bound, and parameter name.
  char buf[16];
  int upper_bound;
  int lower_bound;
  char requested_parameter[16];
  printf("\nDerivative Checker\n");

  printf("Enter upper bound: ");
  fgets(buf, sizeof(buf), stdin);
  if(sscanf(buf, "%d", &upper_bound) != 1) {
    printf("Invalid upper bound\n");
    return;
  }

  printf("Enter lower bound: ");
  fgets(buf, sizeof(buf), stdin);
  if(sscanf(buf, "%d", &lower_bound) != 1) {
    printf("Invalid lower bound\n");
    return;
  }

  printf("Enter parameter name (e.g. alpha, must be an active parameter): ");
  fgets(buf, sizeof(buf), stdin);
  if(sscanf(buf, "%s", requested_parameter) != 1) {
    printf("Invalid parameter name\n");
    return;
  }

  // Ensure lower and upper bound are valid
  if(lower_bound > upper_bound) {
    printf("Lower bound must be less than or equal to upper bound.\n");
    return;
  }

  printf("\n");

  // Attempt to find the parameter and args corresponding to the
  // parameter name requested by the user. If the parameter is
  // found, found_parameter will be set to 1.
  Parameter* parameter;
  Grad_Args arg;
  int found_parameter = 0;
  int iTheta = 0;
  for(int i = 0; i < parameters->count; i++) {
    if(parameters->parameter[i]->active == TRUE) {
      // Check if the requested parameter name matches this parameter's name
      if(strcmp(requested_parameter, parameters->parameter[i]->name) == 0) {
        parameter = parameters->parameter[i];
        arg = args[iTheta];
        found_parameter = 1;
      }
      iTheta++;
    }
  }

  // If the user entered an invalid parameter name or a parameter
  // which is inactive, notify the user and return to normal
  // program flow.
  if(found_parameter == 0) {
    printf("Could not find an active parameter named %s\n", requested_parameter);
    return;
  }


  parameter->grad((void *) &arg);

  #if REAL == DOUBLE
    // lf
    printf("analytic: %lf\n",parameter->gradient);
  #elif REAL == FLOAT
    // f
    printf("analytic: %f\n",parameter->gradient);
  #elif REAL == LONGDOUBLE
    // Lf
    printf("analytic: %Lf\n",parameter->gradient);
  #endif

  Real par_save = parameter->value;

  for (int i=upper_bound;i>=lower_bound;i--) {

    Real dX = exp((Real)i);

    parameter->value = par_save + dX;

    Real dY = K(parameters,d,&core_args) - *f;

  #if REAL == DOUBLE
    // lf
    printf("%lf %lf\n",dX,dY/dX);
  #elif REAL == FLOAT
    // f
    printf("%f %f\n",dX,dY/dX);
  #elif REAL == LONGDOUBLE
    // Lf
    printf("%Lf %Lf\n",dX,dY/dX);
  #endif
  }
}

void check_derivative_no_meschach(Parameters * parameters, Grad_Args_No_MESCHACH * args, Data * d, Real * f, Solve_Core_Args core_args) {
  // Prompt user for the lower bound, upper bound, and parameter name.
  char buf[16];
  int upper_bound;
  int lower_bound;
  char requested_parameter[16];
  printf("\nDerivative Checker\n");

  printf("Enter upper bound: ");
  fgets(buf, sizeof(buf), stdin);
  if(sscanf(buf, "%d", &upper_bound) != 1) {
    printf("Invalid upper bound\n");
    return;
  }

  printf("Enter lower bound: ");
  fgets(buf, sizeof(buf), stdin);
  if(sscanf(buf, "%d", &lower_bound) != 1) {
    printf("Invalid lower bound\n");
    return;
  }

  printf("Enter parameter name (e.g. alpha, must be an active parameter): ");
  fgets(buf, sizeof(buf), stdin);
  if(sscanf(buf, "%s", requested_parameter) != 1) {
    printf("Invalid parameter name\n");
    return;
  }

  // Ensure lower and upper bound are valid
  if(lower_bound > upper_bound) {
    printf("Lower bound must be less than or equal to upper bound.\n");
    return;
  }

  printf("\n");

  // Attempt to find the parameter and args corresponding to the
  // parameter name requested by the user. If the parameter is
  // found, found_parameter will be set to 1.
  Parameter* parameter;
  Grad_Args_No_MESCHACH arg;
  int found_parameter = 0;
  int iTheta = 0;
  for(int i = 0; i < parameters->count; i++) {
    if(parameters->parameter[i]->active == TRUE) {
      // Check if the requested parameter name matches this parameter's name
      if(strcmp(requested_parameter, parameters->parameter[i]->name) == 0) {
        parameter = parameters->parameter[i];
        arg = args[iTheta];
        found_parameter = 1;
      }
      iTheta++;
    }
  }

  // If the user entered an invalid parameter name or a parameter
  // which is inactive, notify the user and return to normal
  // program flow.
  if(found_parameter == 0) {
    printf("Could not find an active parameter named %s\n", requested_parameter);
    return;
  }


  parameter->grad((void *) &arg);

  #if REAL == DOUBLE
    // lf
    printf("analytic: %lf\n",parameter->gradient);
  #elif REAL == FLOAT
    // f
    printf("analytic: %f\n",parameter->gradient);
  #elif REAL == LONGDOUBLE
    // Lf
    printf("analytic: %Lf\n",parameter->gradient);
  #endif

  Real par_save = parameter->value;

  for (int i=upper_bound;i>=lower_bound;i--) {

    Real dX = exp((Real)i);

    parameter->value = par_save + dX;

    Real dY = K(parameters,d,&core_args) - *f;

  #if REAL == DOUBLE
    // lf
    printf("%lf %lf\n",dX,dY/dX);
  #elif REAL == FLOAT
    // f
    printf("%f %f\n",dX,dY/dX);
  #elif REAL == LONGDOUBLE
    // Lf
    printf("%Lf %Lf\n",dX,dY/dX);
  #endif
  }
}
