// SPADE: Stock assessment using PArtial Differential Equations
// Copyright 2016 State of Queensland

// This file is part of SPADE

// SPADE is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// SPADE is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include <fenv.h>
#include <math.h>
#include <signal.h>
#include "spade.h"
#include "common.h"
#include "arg.h"
#include "parameters.h"
#include "machinery/VMGMM.h"
#include "machinery/objfns.h"
#include "optim/optim.h"
#include "plotting/plot.h"
#include "mathprop/mathprop.h"
#include "machinery/alpha/grad_alpha.h"
#include "machinery/beta/grad_beta.h"
#include "machinery/gamma/grad_gamma.h"
#include "machinery/iota/grad_iota.h"
#include "machinery/kappa/grad_kappa.h"
#include "machinery/omega/grad_omega.h"
#include "util/util.h"

int feenableexcept(int);

void print_usage() {
    printf(
      "SPADE: Stock assessment using PArtial Differential Equations\n"
      "Copyright 2016 State of Queensland\n\n"
      "SPADE comes with ABSOLUTELY NO WARRANTY and you are welcome \n"
      "to redistribute it under certain conditions; further details below.\n"
      "\n"
      "Usage:\n"
      "  spade -fn <file> -alpha <a> -beta <b> -gamma <g> -iota <i> -kappa <k>\n"
      "        -omega <w>\n"
      "\n"
      "Options:\n"
      "  -fn <file>\n"
      "      Specifies the common name of the input data files. SPADE will attempt\n"
      "      to read the files <file>-ce.dat and <file>-lf.dat from disk. This\n"
      "      option is required.\n"
      "\n"
      "  -minfish <minfish>      Default: 250\n"
      "\n"
      "  -J <J>                  Default: 400\n"
      "\n"
      "  -timestep <timestep>    Default: 0.025\n"
      "      The model timestep interval represented as a fraction of one year.\n"
      "\n"
      "  -warmup-ratio <ratio>   Default: 1\n"
      "      The SPADE model consists of two stages: a warmup stage followed by\n"
      "      a model stage. This option specifies the number of warmup steps as\n"
      "      a function of the number of model steps. For example, a warmup ratio\n"
      "      of 1 specified the same number of warmup steps as model steps. A warmup\n"
      "      ratio of 0.5 specifies half the number of warmup steps as compared to\n"
      "      model steps. Any numeric value greater than or equal to zero is an\n"
      "      acceptable warmup ratio.\n"
      "\n"
      "Parameters:\n"
      "  To disable a parameter suffix the parameter name with '-disabled'.\n"
      "  Example: 'spade -alpha-disabled <a> ...'\n"
      "\n"
      "SPADE is free software: you can redistribute it and/or modify\n"
      "it under the terms of the GNU Lesser General Public License as published by\n"
      "the Free Software Foundation, either version 3 of the License, or\n"
      "(at your option) any later version.\n\n"
      "SPADE is distributed in the hope that it will be useful,\n"
      "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
      "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
      "GNU Lesser General Public License for more details.\n\n"
      "You should have received a copy of the GNU Lesser General Public License\n"
      "along with this program.  If not, see <https://www.gnu.org/licenses/>.\n"
      "This product includes code derived from software written by Jorge More and \n"
      "developed by the University of Chicago, as Operator of Argonne National Laboratory.\n"
      "This product includes software (Meschach) developed by David E. Stewart and Zbigniew Leyk.\n"

  );
}

int main(int argc, char *argv[])
{
  feenableexcept(FE_DIVBYZERO); 
  feenableexcept(FE_INVALID); 
  feenableexcept(FE_OVERFLOW);
  signal(SIGINT, request_interactive_mode);

  // To test interactive mode in the debug environment,
  // uncomment this line
  //request_interactive_mode(1);

  int N;
  int minfish;
  Real k;
  Real warmup_ratio;

  if (!NEWOBJ) {
    // Read model-related command line args and set defaults if
    // arguments have not been provided
    if(arg_read_int("minfish", &minfish, argc, argv) == FALSE) {
      minfish = 250;
    }
  }

  if(arg_read_int("J", &J, argc, argv) == FALSE) {
    J = 400;
  }

  if(arg_read_real("timestep", &k, argc, argv) == FALSE) {
    k = 0.025;
  }

  if(arg_read_real("warmup-ratio", &warmup_ratio, argc, argv) == FALSE) {
    // By default we have an equal number of warmup and model steps
    warmup_ratio = 1;
  }

  if(warmup_ratio < 0) {
    print_usage();
    exit(EXIT_FAILURE);
  }

  // Read and parse data files
  char * data_file_name;
  if(arg_read_string("fn", &data_file_name, argc, argv) == FALSE) {
    print_usage();
    exit(EXIT_FAILURE);
  }

  Data data;
  NewData newdata; 
  
  if (!NEWOBJ) {

      data_read_ce(data_file_name, &data, &N, k);
      data_read_lf(data_file_name, &data, N, k, minfish);
  }
  else
    {
      data_read_ce_new(data_file_name,&newdata);
      data_read_lf_new(data_file_name,&newdata);
      
    }

		       
  // The model consists of two stages - a warmup stage followed by the model stage.
  // I (total number of time steps) = warmup_steps + N (number of time steps for model)
  //   * warmup_steps = N (number of time steps for model) * warmup_ratio
  //       For example a warmup_ratio of 1 would imply an equal number of warmup steps and model steps
  //   * warmup_ratio = number of warmup steps as a fraction of number of model steps
  //   * Number of time steps for model is determined by k (time step as fraction of a year) and  data->Y (number of years of data)

  int warmup_steps = floor(N * warmup_ratio);
  data.I = warmup_steps + N;
  data.S = warmup_steps;

  if(!SGNM) {
    data.J = J + data.I;
  } else {
    data.J = J;
  }

  data.k = k;

  // Read optim options
  OptimControl optim;
  optim_control_read("control.optim", &optim);

  // Configure each parameter. This must be updated when a
  // new parameter is created.
  Parameters parameters= {
   .alpha = { .name = "alpha", .grad = &grad_alpha },
   .beta = { .name = "beta", .grad = &grad_beta },
   .gamma = { .name = "gamma", .grad = &grad_gamma },
   .iota = { .name = "iota", .grad = &grad_iota },
   .kappa = { .name = "kappa", .grad = &grad_kappa },
   .omega = { .name = "omega", .grad = &grad_omega }
  };

  if (!MESCHACH)
    {

  parameters.alpha.grad = &grad_alpha_clean;
  parameters.beta.grad = &grad_beta_clean;
  parameters.gamma.grad = &grad_gamma_clean;
  parameters.iota.grad = &grad_iota_clean;
  parameters.kappa.grad = &grad_kappa_clean;
  parameters.omega.grad = &grad_omega_clean;

    }
  
  // Map all parameters to the parameter array. This must
  // be updated when a new parameter is created.
  parameters.parameter[0] = &parameters.alpha;
  parameters.parameter[1] = &parameters.beta;
  parameters.parameter[2] = &parameters.gamma;
  parameters.parameter[3] = &parameters.iota;
  parameters.parameter[4] = &parameters.kappa;
  parameters.parameter[5] = &parameters.omega;
  parameters.count = PARAMETER_COUNT;

  // Read all parameter values from the command line.
  if(!parameters_read(&parameters, argc, argv)) {
    print_usage();
    exit(EXIT_FAILURE);
  }

  VEC *theta = parameters_to_vec(&parameters);

  h = parameters.omega.value / J;

  char labbuffer[10];
  sprintf(labbuffer,"before");
  //plot(&parameters,&data,labbuffer);
  
  theta = bfgs(VMGMM,theta,&data,&parameters,optim);

  char labbuffer2[10];
  sprintf(labbuffer2,"after");
  //plot(&parameters,&data,labbuffer2);
  
  V_FREE(theta);
  V_FREE(data.cat);
  V_FREE(data.eff);
  free(data.t_id);
  free(data.t_sz);

  for (int i=0;i<data.n;i++)
    free(data.lf[i]);
  free(data.lf);

  return(0);
}

