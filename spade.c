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


//#include <admodel.h>
//#include <contrib.h>

//extern "C"  {
//  void ad_boundf(int i);
//}

#include <fenv.h>
#include <math.h>
#include <signal.h>

//extern "C" {
#include "spade.h"
#include "arg.h"
#include "machinery/VMGMM.h"
#include "machinery/objfns.h"
#include "optim/optim.h"
#include "plotting/plot.h"
  //#include "mathprop/mathprop.h"
#include "machinery/alpha1/grad_alpha1.h"
#include "machinery/alpha2/grad_alpha2.h"  
#include "machinery/beta/grad_beta.h"
#include "machinery/gamma/grad_gamma.h"
#include "machinery/iota/grad_iota.h"
#include "machinery/kappa/grad_kappa.h"
#include "machinery/omega/grad_omega.h"
#include "util/util.h"
//}

Real iota1;//=5.2;
Real iota2;//=0.619;
Real phi;//=17;
Real eta1;//=1.703205e-5;
Real eta2;//=2.9526;

int interactive_mode_requested = 0;

Real A1=0;
Real A2=0;

int J;
Real h;
Da d;
Real k;

int *idx;

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

/*
class model_data : public ad_comm{
  data_int N;
  data_int I;
  data_int J_admb;
  data_number k;
  data_vector cat;
  data_vector eff;
  data_ivector start;
  data_ivector finish;
  data_vector tbar;
  data_matrix precomp;
  double iota1_admb;
  double iota2_admb;
  double phi_admb;
  double phi_admb2;
  double eta1_admb;
  double eta2_admb;
  dmatrix spb;
  dvector SpB;
  double SpBrat;
  double Urat;
  dmatrix Fj;
  dvector F;
  dvector M;
  double lastF;
  double lastSpB;
  double psi;
  double sig1;
  double sig2;
  ~model_data();
  model_data(int argc,char * argv[]);
  friend class model_parameters;
};

class model_parameters : public model_data ,
  public function_minimizer
{
public:
  ~model_parameters();
  void preliminary_calculations(void);
  void set_runtime(void);
  virtual void * mycast(void) {return (void*)this;}
  static int mc_phase(void)
  {
    return initial_params::mc_phase;
  }
  static int mceval_phase(void)
  {
    return initial_params::mceval_phase;
  }
  static int sd_phase(void)
  {
    return initial_params::sd_phase;
  }
  static int current_phase(void)
  {
    return initial_params::current_phase;
  }
  static int last_phase(void)
  {
    return (initial_params::current_phase
      >=initial_params::max_number_phases);
  }
  static prevariable current_feval(void)
  {
    return *objective_function_value::pobjfun;
  }
private:
  ivector integer_control_flags;
  dvector double_control_flags;
  param_init_bounded_number alpha1;
  param_init_bounded_number alpha2;
  param_init_bounded_number kappa;
  param_init_bounded_number omega;
  param_init_bounded_number beta0;
  param_init_bounded_number gam;
  param_init_bounded_number iota;
  param_number inter;
  param_number inter2;
  param_number Z;
  param_number ubar;
  param_number vbar;
  param_number wbar;
  param_stddev_number obj1;
  param_matrix n;
  param_matrix nh;
  param_matrix x;
  param_matrix xh;
  param_matrix CR;
  param_vector HF;
  param_vector cdr;
  param_vector cdt;
  param_vector Qi;
  param_number prior_function_value;
  param_number likelihood_function_value;
  objective_function_value ff;
public:
  virtual void userfunction(void);
  virtual void report(const dvector& gradients);
  virtual void final_calcs(void);
  model_parameters(int sz,int argc, char * argv[]);
  virtual void initializationfunction(void){}
 dvariable zstar(const dvariable& x, const dvariable& N, const double t);
 double e(const double t);
 double c(const double t);
 dvariable g(const dvariable& x);
 dvariable b(const dvariable& x);
 dvariable s(const dvariable& x, const double t);
 dvariable w(const dvariable& x);
 dvariable Q2(const dvar_vector& X, const dvar_vector& W, const dvariable& g, const int sz);
 dvariable Q(const dvar_vector& X, const dvar_vector& V, const int sz);
 double Q(const dvector& X, const dvector& V, const int sz);

};
  
model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  N.allocate("N");
  I.allocate("I");
  J_admb.allocate("J_admb");
  k.allocate("k");
  cat.allocate(0,N,"cat");
  eff.allocate(0,4*N,"eff");
  start.allocate(0,I,"start");
  finish.allocate(0,I,"finish");
  tbar.allocate(0,I,"tbar");
  precomp.allocate(0,I,start,finish,"precomp");  
  spb.allocate(0,I,start,finish);
  SpB.allocate(0,I);
  Fj.allocate(0,I,start,finish);
  F.allocate(0,I);
  M.allocate(0,I);
  tbar.allocate(0,I);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  alpha1.allocate(0,0.9,"alpha1");
  alpha2.allocate(0,0.9,"alpha2");
  kappa.allocate(0.05,0.9,-1,"kappa");
  omega.allocate(50,200,-1,"omega");
  beta0.allocate(0,2.1,1,"beta0");
  gam.allocate(.1,3.1,"gam");
  iota.allocate(0,0.1,"iota");
  inter.allocate("inter");
  #ifndef NO_AD_INITIALIZE
  inter.initialize();
  #endif
  inter2.allocate("inter2");
  #ifndef NO_AD_INITIALIZE
  inter2.initialize();
  #endif
  Z.allocate("Z");
  #ifndef NO_AD_INITIALIZE
  Z.initialize();
  #endif
  ubar.allocate("ubar");
  #ifndef NO_AD_INITIALIZE
  ubar.initialize();
  #endif
  vbar.allocate("vbar");
  #ifndef NO_AD_INITIALIZE
  vbar.initialize();
  #endif
  wbar.allocate("wbar");
  #ifndef NO_AD_INITIALIZE
  wbar.initialize();
  #endif
  obj1.allocate("obj1");
  n.allocate(0,I,start,finish,"n");
  #ifndef NO_AD_INITIALIZE
    n.initialize();
  #endif
  nh.allocate(0,I,start,finish,"nh");
  #ifndef NO_AD_INITIALIZE
    nh.initialize();
  #endif
  x.allocate(0,I,start,finish,"x");
  #ifndef NO_AD_INITIALIZE
    x.initialize();
  #endif
  xh.allocate(0,I,start,finish,"xh");
  #ifndef NO_AD_INITIALIZE
    xh.initialize();
  #endif
  CR.allocate(0,I,start,finish,"CR");
  #ifndef NO_AD_INITIALIZE
    CR.initialize();
  #endif
  HF.allocate(0,I,"HF");
  #ifndef NO_AD_INITIALIZE
    HF.initialize();
  #endif
  cdr.allocate(0,I,"cdr");
  #ifndef NO_AD_INITIALIZE
    cdr.initialize();
  #endif
  cdt.allocate(0,I,"cdt");
  #ifndef NO_AD_INITIALIZE
    cdt.initialize();
  #endif
  Qi.allocate(0,I,"Qi");
  #ifndef NO_AD_INITIALIZE
    Qi.initialize();
  #endif
  ff.allocate("ff");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}


void model_parameters::userfunction(void)
{
  ff=0.0;
  sig1=.05;
  sig2 = .7;
  psi = 200;
  iota1_admb=5.2;
  iota2_admb=0.619;
  phi_admb=17;
  phi_admb2=20.3;
  eta1_admb=1.703205e-5;
  eta2_admb=2.9526;
  gam *= 1e-7;
  iota *= 1e-3;
  
  ff=0;
  int i,j;
  for (j=0;j<=J_admb;j++) 
    x(0,j) = j*omega/J_admb;
  x(0,J_admb) -= 1e-9; // for pow base neq 0
  //if (81*kappa*kappa*omega*omega*(alpha1+2*alpha2*omega)*(alpha1+2*alpha2*omega) < 12*kappa*(alpha1*omega+kappa)*(alpha1*omega+kappa)*(alpha1*omega+kappa))
  //  ff = 1e6;
  //else {
  dvariable zeta = sqrt(81*kappa*kappa*omega*omega*(alpha1+2*alpha2*omega)*(alpha1+2*alpha2*omega) - 12*kappa*(alpha1*omega+kappa)*(alpha1*omega+kappa)*(alpha1*omega+kappa));
  inter = 9*alpha1*kappa*kappa*omega + 18*alpha2*kappa*kappa*omega*omega + kappa*zeta; 
  Z = pow(inter,1/3.0) / (3* pow(2/3.0,1/3.0)) + pow(2/3.0,1/3.0) * kappa * (alpha1*omega + kappa) / pow(inter,1/3.0);
  ubar = (Z-beta0-kappa)/gam;
  //if (ubar <= 0) 
  //  ff = 1e6;
  //else {
  vbar = (kappa*omega*ubar) / (beta0+gam*ubar+kappa);
  wbar = (2*kappa*omega*vbar) / (beta0+gam*ubar+2*kappa);
  //cout << endl;  
  for (j=0;j<=J_admb;j++)
    n(0,j) = (alpha1*vbar+alpha2*wbar)*pow(omega-x(0,j),(beta0+gam*ubar)/kappa - 1) / (kappa*pow(omega,(beta0+gam*ubar)/kappa));
  Qi(0) = Q(x(0),n(0),J_admb);
  //cout << c(k*(0-N)) << endl;  
  for (j=0;j<=J_admb;j++)
    spb(0,j) = value(b(x(0,j))*n(0,j)); 
  SpB(0) = Q(value(x(0)),spb(0),J_admb);
  i=0;
  x(0,0) = 1e-9;
  for (j=0;j<=J_admb;j++)
    CR(i,j) = w(x(i,j))*s(x(i,j),k*(i-N))*n(i,j)*iota*e(k*(i-N));
  for (i=1;i<=I;i++) {
    double t = k*(i-N-1);
    double th = k*(i-N-.5);
    xh(i,0) = 1e-9;
    for (j=1;j<J_admb+i;j++) 
      xh(i,j)=x(i-1,j-1)+(k/2)*g(x(i-1,j-1));
    xh(i,J_admb+i)=omega;
    for (j=1;j<=J_admb+i;j++)
      nh(i,j) = n(i-1,j-1)*mfexp(-(k/2)*zstar(x(i-1,j-1),Qi(i-1),t));
    nh(i,0)= Q2(xh(i),nh(i),g(0),J_admb+i); 
    dvariable Qhalf = Q(xh(i),nh(i),J_admb+i);
    x(i,0) = 1e-9;
    for (j=1;j<J_admb+i;j++)
      x(i,j) = x(i-1,j-1) + k*g(xh(i,j));
    x(i,J_admb+i) = omega;
    for (j=1;j<=J_admb+i;j++)
      n(i,j) = n(i-1,j-1)*mfexp(-k*zstar(xh(i,j),Qhalf,th));
    n(i,0) = Q2(x(i),n(i),g(0),J_admb+i); 
    Qi(i) = Q(x(i),n(i),J_admb+i);
    //cout << c(k*(i-N)) << endl;
    for (j=0;j<=J_admb+i;j++)
      CR(i,j) = w(x(i,j))*s(x(i,j),k*(i-N))*n(i,j)*iota*e(k*(i-N)); 
    HF(i) = Q(x(i),CR(i),J_admb+i)/1e3; 
    if (k*(i-N)<0)
      F(i) = value(iota*e(k*(i-N)));
    else
      F(i) = value(iota*e(k*(i-N)));
    M(i) = value(Qi(i)*gam);
    for (j=0;j<=J_admb+i;j++)
      spb(i,j) = value(b(x(i,j))*n(i,j)); 
    SpB(i) = Q(value(x(i)),spb(i),J_admb+i);
  }
  //exit(1);
  dvariable avgHF = 0;
  for (i=N;i<=I;i++)
    avgHF += HF(i);
  avgHF /= (N+1);
  lastF = 0;
  for (i=0;i<(int)1/k;i++)
    lastF += F(I-i);
  lastF *= k;
  lastSpB = 0;
  for (i=0;i<(int)1/k;i++)
    lastSpB += SpB(I-i);
  lastSpB *= k;
  double SpBen = 0;
  for (i=0;i<(int)1/k;i++)
    SpBen += SpB(I-i);
  SpBen *= k;
  SpBrat = SpBen / SpB(0); 
  double Uen = 0;
  for (i=0;i<(int)1/k;i++)
    Uen += value(Qi(I-i));
  Uen *= k;
  Urat = Uen / value(ubar);
  // Objective function 1
  obj1=0;
  for (i=0;i<=I;i++) {
    dvariable Qcr = Q(x(i),CR(i),J_admb+i);
    for (j=0;j<=finish(i);j++) {
      nh(i,j) = square(  c(k*(i-N)) * (precomp(i,j) + (1-tbar(i)) * CR(i,j)/Qcr )  - CR(i,j) );
    }
    //cout << Q(x(i),nh(i),J_admb+i) << endl;
    obj1 += Q(x(i),nh(i),J_admb+i);
    //cout << obj1 << endl;
  }
  ff=obj1;
  //exit(1);
  //}
  //}
  if (mceval_phase()) 
    cout << alpha1 << " " << alpha2 << " " << kappa << " " << omega << " " << beta0 << " " << iota << " " << ff << " " << SpBrat << " " << lastF << " " << lastSpB << endl;
}

dvariable model_parameters::zstar(const dvariable& x, const dvariable& N, const double t)
{
  return beta0 + gam*N + s(x,t)*iota*e(t) - kappa;
}

double model_parameters::e(const double t)
{
  double tt = t + 24;
  double cek = k/2;
  int idx = floor((tt+(cek/2) - 1e-12)/cek);
  return eff(idx);
}

double model_parameters::c(const double t)
{
  if (t<0)
    return cat(0);
  else
    {
      int idx = floor((t+k/2 - 1e-12)/k);
      return cat(idx);
    }
}

dvariable model_parameters::g(const dvariable& x)
{
  return kappa*(omega-x);
}

dvariable model_parameters::b(const dvariable& x)
{
  return alpha1*x + alpha2*square(x);
}

dvariable model_parameters::s(const dvariable& x, const double t)
{

  if (t<0) {
    if (x < iota1_admb*phi_admb)
      return mfexp(-square(x-iota1_admb*phi_admb)/(2*iota2_admb*square(phi_admb)));
    else {
      if (x < iota1_admb*phi_admb2)
        return 1;
      else
        return mfexp(-square(x-iota1_admb*phi_admb2)/(2*iota2_admb*square(phi_admb2)));
    }
    } else 
    return mfexp(-square(x-iota1_admb*phi_admb)/(2*iota2_admb*square(phi_admb)));
}

dvariable model_parameters::w(const dvariable& x)
{
  return eta1_admb*pow(x,eta2_admb);
}

dvariable model_parameters::Q2(const dvar_vector& X, const dvar_vector& W, const dvariable& g, const int sz)
{
  dvariable ac = X(1)*b(X(1))*W(1) - X(0)*b(X(1))*W(1); 
  for (int j=1;j<sz;j++) 
    ac += (b(X(j))*W(j) + b(X(j+1))*W(j+1)) * (X(j+1)-X(j)); 
  ac /= (2*g + X(0)*b(X(0)) - X(1)*b(X(0))); 
  return ac;
}

dvariable model_parameters::Q(const dvar_vector& X, const dvar_vector& V, const int sz)
{
  dvariable ac=0;
  for (int j=0;j<sz;j++) 
    ac = ac + (X(j+1) - X(j)) * (V(j) + V(j+1))/2.;
  return ac;
}

double model_parameters::Q(const dvector& X, const dvector& V, const int sz)
{
  double ac=0;
  for (int j=0;j<sz;j++) 
    ac = ac + (X(j+1) - X(j)) * (V(j) + V(j+1))/2.;
  return ac;
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::report(const dvector& gradients){}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;
*/

int main(int argc, char *argv[])
{

  /*
  ad_set_new_handler();
  ad_exit=&ad_boundf;
  gradient_structure::set_NO_DERIVATIVES();
  gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
  if (!arrmblsize) arrmblsize=15000000;
  model_parameters mp(arrmblsize,argc,argv);
  mp.iprint=10;
  mp.preliminary_calculations();
  */
  
  signal(SIGINT, request_interactive_mode);

  // Read and parse data files
  char * data_file_name;
  if(arg_read_string("fn", &data_file_name, argc, argv) == FALSE) {
    print_usage();
    exit(EXIT_FAILURE);
  }

  data_read_fast(data_file_name);

  // Read optim options
  OptimControl optim;
  optim_control_read("control.optim", &optim);

  // Configure each parameter. This must be updated when a
  // new parameter is created.
  
  Parameters parameters =
    {

    alpha1 : { grad : &grad_alpha1_fast, value: 0, active: 0, gradient: 0, name: "alpha1" },
    alpha2 : { grad : &grad_alpha2_fast, value: 0, active: 0, gradient: 0, name: "alpha2" },
    beta : { grad : &grad_beta_fast, value: 0, active: 0, gradient: 0, name : "beta",  },
    gamma : { grad : &grad_gamma_fast, value: 0, active: 0, gradient: 0, name : "gamma", },
    iota : { grad : &grad_iota_fast, value: 0, active: 0, gradient: 0, name : "iota"},
    kappa : { grad : &grad_kappa, value: 0, active: 0, gradient: 0,name : "kappa" },
    omega : { grad : &grad_omega, value: 0, active: 0, gradient: 0, name : "omega"}

    };

  // Map all parameters to the parameter array. This must
  // be updated when a new parameter is created.

  parameters.parameter[0] = &parameters.alpha1;
  parameters.parameter[1] = &parameters.alpha2; 
  parameters.parameter[2] = &parameters.beta;
  parameters.parameter[3] = &parameters.gamma;
  parameters.parameter[4] = &parameters.iota;
  parameters.parameter[5] = &parameters.kappa;
  parameters.parameter[6] = &parameters.omega;
  parameters.count = PARAMETER_COUNT;
  parameters.ff = 0;

  // Read all parameter values from the command line.
  if(!parameters_read(&parameters, argc, argv)) {
    print_usage();
    exit(EXIT_FAILURE);
  }

  //0x847cf0
  
  MeVEC *theta = parameters_to_vec(&parameters);

  //h = parameters.omega.value / (d.J-d.I);
  h = parameters.omega.value / d.J;

  char labbuffer[10];
  sprintf(labbuffer,"before");
  //plot_fast(&parameters,labbuffer);
  
  theta = bfgs(_VMGMM,theta,&parameters,optim);

  char labbuffer2[10];
  sprintf(labbuffer2,"after");
  //plot_fast(&parameters,labbuffer2);

  //mp.computations(argc,argv,_VMGMM,theta,&parameters);

  /* 
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

  MeVEC *theta = parameters_to_vec(&parameters);

  Real h = parameters.omega.value / J;

  char labbuffer[10];
  sprintf(labbuffer,"before");
  //plot(&parameters,&data,labbuffer);
  
  //theta = bfgs(VMGMM,theta,&data,&parameters,optim);

  char labbuffer2[10];
  sprintf(labbuffer2,"after");
  //plot(&parameters,&data,labbuffer2);
  
  mp.computations(argc,argv,VMGMM,theta,&data,&parameters);
  
  V_FREE(theta);
  V_FREE(data.cat);
  V_FREE(data.eff);
  free(data.t_id);
  free(data.t_sz);

  for (int i=0;i<data.n;i++)
    free(data.lf[i]);
  free(data.lf);
  
  feenableexcept(FE_DIVBYZERO); 
  feenableexcept(FE_INVALID); 
  feenableexcept(FE_OVERFLOW);

  // To test interactive mode in the debug environment,
  // uncomment this line
  //request_interactive_mode(1);
  */

  return(0);
}

/*
extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here 
    exit(i);
  }
}*/


/*

for splines

Li <- array(0,shorterN)
Lin <- array(0,shorterN)

min.lengths <- 20

for (i in 1:nrow(cssf))
    Lin[cssf$idx[i]] <- Lin[cssf$idx[i]] + 1

Li <- Lin > min.lengths

Memaxnk <- 0
Memaxnbasis <- 0

szs <- {}
timesidx <- {}

#outxall <- {}
#outyall <- {}


for (i in 1:shorterN) {
    if (Li[i]) {
      cand <- subset(cssf,idx==i)            
      
      cdi <- {}
      for (j in 1:nrow(cand)) 
          cdi <- c(cdi,rep(cand$Length1[j],100/cand$PctCatchSampled[j]))
      szs <- c(szs,length(cdi))
      timesidx <- c(timesidx,i)
      
      binidx <- .bincode(cdi,seq(0,200,by=3))

      m <- table(binidx)
      x <- (as.numeric(unlist(dimnames(m)))-1)*3+1.5
      y <- as.numeric(m) / sum(binidx)
     
      xl <- (min(x)-10):(min(x)-6)
      yl <- rep(0,5)
      xh <- (Memax(x)+6):(Memax(x)+10)
      yh <- rep(0,5)
      
      x <- c(xl,min(x)-3,x,Memax(x)+3,xh)
      y <- c(yl,y[1]/2,y,y[length(y)]/2,yh)            
      
      xsp <- (x - min(x)) / (Memax(x)-min(x))
      
      spl <- smooth.spline(xsp,y,nknots=length(xsp)-1,spar=.5)
      bspl.basis <- create.bspline.basis(unique(spl$fit$knot))

      #outx <- seq(0,1,length=300)
      #outy <- array(0,300)
      #for (j in 1:300)
      #    outy[j] <- eval.basis(outx[j],bspl.basis) %*% spl$fit$coef
           
      #if (i == 592){
      #    outx <- outx * (Memax(x)-min(x)) + min(x)
      #    outxall <- outx
      #    outy <- outy*200
      #    outyall <- outy
      #    #plot(outx,outy,type='l',ylim=c(-100,5),xlim=c(50,180))
      #}
      #else if (i>592) {
      #    outx <- outx * (Memax(x)-min(x)) + min(x)
      #    outxall <- c(outxall,outx)
      #    outy <- (outy*200) - epsilon*(1-exp(-delta*(i-592)/25))
      #    outyall <- c(outyall,outy)
      #    #lines(outx,outy)
      #}
      
      nbasis <- bspl.basis$nbasis
      params <- bspl.basis$params
      rangeval <- bspl.basis$rangeval
      breaks <- c(rangeval[1], params, rangeval[2])
      norder <- nbasis - length(breaks) + 2
      nbreaks <- length(breaks)
      knots <- c(rep(breaks[1], norder -1), breaks, rep(breaks[nbreaks], norder -1))
      nbasis <- nbreaks + norder - 2
      nk <- length(knots)

      if (nbasis>Memaxnbasis)
          Memaxnbasis <- nbasis
      if (nk>Memaxnk)
          Memaxnk <- nk

    }
}


ntimes <- length(szs)
splcoef <- matrix(0,ntimes,Memaxnbasis)
knots <- matrix(0,ntimes,Memaxnk)
nbases <- array(0,ntimes)
nks <- array(0,ntimes)
timesidx <- array(0,ntimes)
szs <- array(0,ntimes)
starts <- array(0,ntimes)
stops <- array(0,ntimes)

jj <- 1

#pdf('lplot.pdf')

for (i in 1:shorterN) {
    if (Li[i]) {
      cand <- subset(cssf,idx==i)
      cdi <- {}
      for (j in 1:nrow(cand)) 
          cdi <- c(cdi,rep(cand$Length1[j],100/cand$PctCatchSampled[j]))

      szs[jj] <- length(cdi)
      timesidx[jj] <- i
      
      binidx <- .bincode(cdi,seq(0,200,by=3))
      
      m <- table(binidx)
      x <- (as.numeric(unlist(dimnames(m)))-1)*3+1.5
      y <- as.numeric(m) / sum(binidx)
                
      xl <- (min(x)-10):(min(x)-6)
      yl <- rep(0,5)
      xh <- (Memax(x)+6):(Memax(x)+10)
      yh <- rep(0,5)
      
      x <- c(xl,min(x)-3,x,Memax(x)+3,xh)
      y <- c(yl,y[1]/2,y,y[length(y)]/2,yh)
      
      starts[jj] <- min(x)
      stops[jj] <- Memax(x)
      
      xsp <- (x - min(x)) / (Memax(x)-min(x))
            
      spl <- smooth.spline(xsp,y,nknots=length(xsp)-1,spar=.5)
      bspl.basis <- create.bspline.basis(unique(spl$fit$knot))

            
      #outx <- seq(0,1,length=300)
      #outy <- array(0,300)
      #for (j in 1:300)
      #    outy[j] <- eval.basis(outx[j],bspl.basis) %*% spl$fit$coef
      
      #if (i == 592){
      #    outx <- outx * (Memax(x)-min(x)) + min(x)      
      #    outy <- outy*200
      #    plot(outx,outy,type='l',xlim=range(outxall),ylim=range(outyall))
      #}
      #else if(i > 592) {
      #    outx <- outx * (Memax(x)-min(x)) + min(x)
      #    outy <- (outy*200) - epsilon*(1-exp(-delta*(i-592)/25))
      #    lines(outx,outy)
      #}            
      
      
      nbasis <- bspl.basis$nbasis
      params <- bspl.basis$params
      rangeval <- bspl.basis$rangeval
      breaks <- c(rangeval[1], params, rangeval[2])
      norder <- nbasis - length(breaks) + 2
      nbreaks <- length(breaks)
      knotstmp <- c(rep(breaks[1], norder -1), breaks, rep(breaks[nbreaks], norder -1))
      knots[jj,1:length(knotstmp)] <- knotstmp
      nbasis <- nbreaks + norder - 2
      nk <- length(knotstmp)
      splcoeftmp <- spl$fit$coef
      splcoef[jj,1:length(splcoeftmp)] <- splcoeftmp
            
      nbases[jj] <- nbasis
      nks[jj] <- nk
     
      jj <- jj + 1
    }
}

#dev.off()

subst <- c(5,37)
#subst <- c(5,11,26,37)
#subst <- c(7,10,26,40)
#subst <- c(1,2,3,4)
ntimes <- length(subst)

nks <- nks[subst]
nbases <- nbases[subst]
knots <- knots[subst,]
splcoef <- splcoef[subst,]
szs <- szs[subst]
timesidx <- timesidx[subst]
starts <- starts[subst]
stops <- stops[subst]


{
    sink('barra_agg.dat',append=TRUE)
    cat(sprintf("%d\n",ntimes))
    for (i in 1:ntimes)
      cat(sprintf("%d ",nks[i]))
    cat(sprintf("\n"))
    for (i in 1:ntimes)
      cat(sprintf("%d ",nbases[i]))
    cat(sprintf("\n"))
    for (i in 1:ntimes)
      cat(sprintf("%d ",1))
    cat(sprintf("\n"))
    for (i in 1:ntimes) {
        for (j in 1:nks[i])  
            cat(sprintf("%f ",knots[i,j]))
        cat(sprintf("\n"))
    }
    for (i in 1:ntimes) {
        for (j in 1:nbases[i])  
            cat(sprintf("%f ",splcoef[i,j]))
        cat(sprintf("\n"))
    }
    for (i in 1:ntimes)      
      cat(sprintf("%d ",szs[i]))
    cat(sprintf("\n"))
    for (i in 1:ntimes)
      cat(sprintf("%d ",timesidx[i]))
    cat(sprintf("\n"))
    for (i in 1:ntimes)      
      cat(sprintf("%f ",starts[i]))
    cat(sprintf("\n"))
    for (i in 1:ntimes)
      cat(sprintf("%f ",stops[i]))
    cat(sprintf("\n"))
    sink()
}


...

in tpl file:

  init_int T1
  init_ivector Kn1(1,T1)
  init_ivector B1(1,T1)
  init_ivector one1(1,T1)
  init_matrix knots1(1,T1,one1,Kn1)
  init_matrix coef1(1,T1,one1,B1)
  init_vector start1(1,T1)
  init_vector stop1(1,T1)

  lf1=0;
  p1.initialize();
  sn.initialize();
  int q,st,jj;
  dvariable ans2;
  for (int q=1;q<=T1;q++) {
    i = H1(q) + Hw1/2;

    for (j=0;j<=J;j++)
      if (x(i,j) > start1(q))
        break;
    st = j;
    
    for (jj=st;jj<=J;jj++) {

      if (x(i,jj) > stop1(q))
        break;

      ans2 = sple(Kn1(q),knots1(q),(x(i,jj)-start1(q))/(stop1(q)-start1(q)),coef1(q));            
      if (ans2 > 0)
        p1(jj) = ans2;
      else 
        p1(jj) = 1e-3;            

      sn(jj) = s(x(i,jj))*n1(i,jj);

    }

    dvar_vector xtmp = extract_row(x,i);
    dvar_vector xtmp2 = xtmp(st,jj-1);
    dvar_vector p12 = p1(st,jj-1);
    dvar_vector sn2 = sn(st,jj-1);
          
    p12.shift(0);
    sn2.shift(0);
    xtmp2.shift(0);

    p12 *= 1/Qjn(xtmp2,p12,xtmp2.size()-1);    
    sn2 *= 1/Qjn(xtmp2,sn2,xtmp2.size()-1);
         
    lf1 += Qjn(xtmp2,elem_prod(p12,log(elem_div(p12,sn2))),xtmp2.size()-1);

  } 



..

FUNCTION dvariable sple(const int nk, const dvar_vector& knots, const dvariable& xval, const dvar_vector& coef)

  int j,q,r,l,offset;
  dvar_vector val(1,4);
  dvar_vector rdel(1,3);
  dvar_vector ldel(1,3); 

  for (j=1;j<=nk;j++)
    if (knots(j) >= xval)
      break;

  l = j - 4;
  offsxet = l-1;

  for (q=0;q<=(4-2);q++) {
    rdel(q+1) = knots(j+q) - xval;
    ldel(q+1) = xval - knots(j - (q+1));
  }

  val(1) = 1;
  dvariable saved;
  dvariable term;
  for (q=1;q<=(4-1);q++) {
    saved=0;
    for (r=0;r<=(q-1);r++) {
      term = val(r+1) / (rdel(r+1) + ldel(q-1-r+1));
      val(r+1) = saved + rdel(r+1) * term;
      saved = ldel(q-1-r+1) * term;
    }
    val(q+1) = saved;
  }

  int ncoef = nk - 4;

  dvar_vector design(1,ncoef);
  design.initialize();
  design(offset+1) = val(1);
  design(offset+2) = val(2);
  design(offset+3) = val(3);
  design(offset+4) = val(4);

  return design * coef;


*/
