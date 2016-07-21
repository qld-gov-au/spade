// SPADE
// Stock assessment using PArtial Differential Equations
// Alex Campbell 'ghostofsandy' 2015 - 2016

// Testing:

// quick test:
//time ./spade -fn karumba .09 .09 1.3 .07
//number function evals: 42, function value: 790362.592819
//Vector: dim: 4
//   0.100412718    0.121894385     1.10724432    0.079491708
//real	1m1.321s
//user	2m30.212s
//sys	0m0.024s

// longer test:
// ./spade -fn karumba .8 .4 1.3 .1
// should return:
//  number function evals: 352, function value: 790362.592819
//  Vector: dim: 4
//     0.100412718    0.121894385     1.10724432    0.079

#include <fenv.h>
#include <math.h>
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

int feenableexcept(int);

void print_usage() {
    printf(
      "SPADE: Stock assessment using PArtial Differential Equations.\n"
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
  );
}

int main(int argc, char *argv[])
{
  feenableexcept(FE_DIVBYZERO); 
  feenableexcept(FE_INVALID); 
  feenableexcept(FE_OVERFLOW);

  int N;
  int minfish;
  Real k;
  Real warmup_ratio;

  // Read model-related command line args and set defaults if
  // arguments have not been provided
  if(arg_read_int("minfish", &minfish, argc, argv) == FALSE) {
    minfish = 250;
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
  data_read_ce(data_file_name, &data, &N, k);
  data_read_lf(data_file_name, &data, N, k, minfish);

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
  Parameters parameters = {
   .alpha = { .name = "alpha", .grad = &grad_alpha },
   .beta = { .name = "beta", .grad = &grad_beta },
   .gamma = { .name = "gamma", .grad = &grad_gamma },
   .iota = { .name = "iota", .grad = &grad_iota },
   .kappa = { .name = "kappa", .grad = &grad_kappa },
   .omega = { .name = "omega", .grad = &grad_omega }
  };

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

  theta = bfgs(VMGMM,theta,&data,&parameters,optim);

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
