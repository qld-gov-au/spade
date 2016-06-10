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

int main(int argc, char *argv[])
{
 
  feenableexcept(FE_DIVBYZERO); 
  feenableexcept(FE_INVALID); 
  feenableexcept(FE_OVERFLOW);

  Real *ct; 
  Real *ti; 
  Real *ln;
  Real *tl;
  int Nce,Nlf;

  int N;
  int minfish;
  Real k;

  FILE *fp;
  fp = fopen("control.model","r");
  fscanf(fp,"%d\n",&minfish);
  fscanf(fp,"%d\n",&J);

  #if REAL == DOUBLE
    // lf
    fscanf(fp,"%lf\n",&k);
  #elif REAL == FLOAT
    // f
    fscanf(fp,"%f\n",&k);
  #elif REAL == LONGDOUBLE
    // Lf
    fscanf(fp,"%Lf\n",&k);
  #endif

  fclose(fp);

  Data data;

  if (argc < 2)
    {
      printf("\n");
      printf("	proper usage requires at least one filename specified\n");
      printf("		e.g. spade -ce ce.dat or ... \n");
      printf("\n");
      printf("	Other arguments:\n");
      printf("			-ce   [no default] catch effort data file\n");
      exit(0);
    }

  for (int i = 1; i < argc; i++) { 
    if (i != argc) {      // Check that we haven't finished parsing already
      if (!strcmp(argv[i], "-fn")) 
	{

	  if (i+1 == argc)
	    {
	      printf("\n no file specified\n");
	      exit(0);
	    }

	  char buffer[30];

	  sprintf(buffer,"%s-ce.dat",argv[i+1]);

	  FILE *fp1;
	  fp1 = fopen(buffer,"r");
	  fscanf(fp1,"%d",&Nce);
	  ct = (Real *) calloc(Nce,sizeof(Real));
	  ti = (Real *) calloc(Nce,sizeof(Real));

	  for (int i=0;i<Nce;i++)
    #if REAL == DOUBLE
      // lf
      fscanf(fp1,"%lf %lf", &ct[i],&ti[i]);
    #elif REAL == FLOAT
      // f
      fscanf(fp1,"%f %f", &ct[i],&ti[i]);
    #elif REAL == LONGDOUBLE
      // Lf
      fscanf(fp1,"%Lf %Lf", &ct[i],&ti[i]);
    #endif

	  VEC *vti = v_get(Nce);
	  for (int i=0;i<Nce;i++)
	    vti->ve[i] = ti[i];

	  PERM *order = px_get(vti->dim);
	  v_sort(vti,order);

	  int Y = ceil(vti->ve[Nce-1]);
	  N = Y/k;

	  PX_FREE(order);
	  V_FREE(vti);

          data.cat = v_get(N+1);
          data.eff = v_get(4*N+1);

	  Real cek = k/4;

	  for (int i=0;i<Nce;i++)
	    {
	      int idx_e = floor((ti[i]+cek/2)/cek);
	      int idx_c = floor((ti[i]+k/2)/k);
	      data.cat->ve[idx_c] += (ct[i]/k)/1e3;
	      data.eff->ve[idx_e] += 1.0/cek;
	    }
  
	  free(ct);
	  free(ti);

	  fclose(fp1);

	  /*
	  printf("\n");
	  for (int i=0;i<data.eff->dim;i++)
	    printf("%f\n",data.eff->ve[i]);
	  exit(1);
	  */

	  sprintf(buffer,"%s-lf.dat",argv[i+1]);

	  fp1 = fopen(buffer,"r");

	  fscanf(fp1,"%d",&Nlf);

	  ln = (Real *) calloc(Nlf,sizeof(Real));
	  tl = (Real *) calloc(Nlf,sizeof(Real));

	  for (int i=0;i<Nlf;i++)
      #if REAL == DOUBLE
          // lf
          fscanf(fp1,"%lf %lf", &ln[i],&tl[i]);
      #elif REAL == FLOAT
          // f
          fscanf(fp1,"%f %f", &ln[i],&tl[i]);
      #elif REAL == LONGDOUBLE
          // Lf
          fscanf(fp1,"%Lf %Lf", &ln[i],&tl[i]);
      #endif

	  fclose(fp1);

	  VEC *lfv = v_get(Nlf);
	  VEC *ilv = v_get(Nlf);
	  VEC *cnt = v_get((int)2*Y/k);

	  for (int i=0;i<Nlf;i++)
	    {
	      lfv->ve[i] = (Real)ln[i];
	      ilv->ve[i] = (int)(Y/k + floor((tl[i]+k/2)/k));
	      cnt->ve[(int)(Y/k + floor((tl[i]+k/2)/k))] += 1;
	    }

	  data.n = 0;
	  for (int i=0;i<cnt->dim;i++)	   
	    if (cnt->ve[i] > minfish)
	      data.n += 1;

          data.lf = (Real **) calloc(data.n,sizeof(Real *));
	  data.t_id = (int *) calloc(data.n,sizeof(int));
	  data.t_sz = (int *) calloc(data.n,sizeof(int));

	  int kk=0;
	  for (int i=0;i<cnt->dim;i++)
	    {
	      if (cnt->ve[i] > minfish)
		{

		  data.t_id[kk] = i;
		  data.t_sz[kk] = cnt->ve[i];		
		  data.lf[kk] = calloc(data.t_sz[kk],sizeof(Real));

		  int jj=0;
		  for (int j=0;j<Nlf;j++)
		    if (ilv->ve[j]==i)
		      {
			data.lf[kk][jj] = lfv->ve[j];
			jj += 1;
		      }

		  kk+=1;

		}
	    }

	  V_FREE(lfv);
	  V_FREE(ilv);
	  V_FREE(cnt);

	  free(ln);
	  free(tl);

	  //sscanf(argv[i+2],"%lf",&alpha);
	  //sscanf(argv[i+3],"%lf",&beta);
	  //sscanf(argv[i+4],"%lf",&gamma);
	  //sscanf(argv[i+5],"%lf",&iota);
	  //sscanf(argv[i+6],"%lf",&kappa);
	  //sscanf(argv[i+7],"%lf",&omega);

	  i += 8;

	} 
      else
	{
	  //printf("problem\n");
          //exit(1);
	}
    }
  }

  OptimControl optim;

  fp = fopen("control.optim","r");
  #if REAL == DOUBLE
    // lf
    fscanf(fp,"%lf\n",&optim.stp);
    fscanf(fp,"%lf\n",&optim.ftol);
    fscanf(fp,"%lf\n",&optim.gtol);
    fscanf(fp,"%lf\n",&optim.stpmin);
    fscanf(fp,"%lf\n",&optim.stpmax);
  #elif REAL == FLOAT
    // f
    fscanf(fp,"%f\n",&optim.stp);
    fscanf(fp,"%f\n",&optim.ftol);
    fscanf(fp,"%f\n",&optim.gtol);
    fscanf(fp,"%f\n",&optim.stpmin);
    fscanf(fp,"%f\n",&optim.stpmax);
  #elif REAL == LONGDOUBLE
    // Lf
    fscanf(fp,"%Lf\n",&optim.stp);
    fscanf(fp,"%Lf\n",&optim.ftol);
    fscanf(fp,"%Lf\n",&optim.gtol);
    fscanf(fp,"%Lf\n",&optim.stpmin);
    fscanf(fp,"%Lf\n",&optim.stpmax);
  #endif

  fscanf(fp,"%d",&optim.maxfev);
  fclose(fp);
  optim.xtol = DBL_EPSILON;

  data.I = 2*N;
  data.S = N;
  data.J = J;
  data.k = k;

  /*
  VEC *th = v_get(2);
  th->ve[0] = .44;
  th->ve[1] = .92;

  VEC *result = bfgs(VMGMM_eq,th,&data);

  //v_output(result);

  Real e_bar = v_sum(data.eff)/data.eff->dim;
  Real iota = result->ve[1]/e_bar;

  //printf("%f %f\n",e_bar,est_iota);

  Real a2 = alpha2;
  VEC *out = v_get(2);

  for (int i=0;i<10;i++) 
    { 

      out = calc_alpha2(alpha1,a2,kappa,omega,result->ve[0],result->ve[1]);
      a2 = a2 - out->ve[0]/out->ve[1];      
    }

  printf("%f %f\n",iota,a2);
  */

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
    printf("Usage: spade -fn <datafile> -alpha <alpha> -beta <beta> -gamma <gamma> -iota <iota> -kappa <kappa> -omega <omega>\n\n");
    printf("  To disable a parameter suffix the parameter name with '-disabled' e.g. '-alpha-disabled'\n");
    exit(1);
  }

  VEC *theta = parameters_to_vec(&parameters);

  h = parameters.omega.value / J;

  //Real cn = ConditionNumber(&parameters, &data);
  //printf("%f\n",cn);
  //exit(0);


/*
  Solve_Core_Args core_args;
  
  int I;
  I = data.I;
  
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

  Real fv = K(&parameters,&data,&core_args);
  Real save = parameters.kappa.value;


  printf("\n");
  for (int i=-10;i<=10;i++) {
    Real delta = (Real)i*.00001; //exp((Real)i);
    parameters.kappa.value = save + delta;
    Real nfv = K_dr(&parameters,&data);
    //Real ch = (nfv - fv) / delta;
    printf("%f %f\n",parameters.kappa.value,nfv);
    
  }
  */
  
  /* 
  // get the active parameters and run their grad functions
  MAT *p = m_get(I,J);
  Grad_Args args;
  args.d = &data;
  VEC *g = v_get(theta->dim);
  args.eff = data.eff;
  args.k = data.k;
  args.S = data.S;
  args.core_args = &core_args;
  args.parameters = &parameters;
  
  parameters.kappa.grad((void *) &(args));
  
  printf("%g\n",parameters.kappa.gradient);*/
  //exit(0);
  
  //char lab1[10]="before";
  //plot(&parameters,&data,lab1);

  theta = bfgs(VMGMM,theta,&data,&parameters,optim);

  //char lab2[10]="after";

  V_FREE(theta);
  V_FREE(data.cat);
  V_FREE(data.eff);
  free(data.t_id);
  free(data.t_sz);

  for (int i=0;i<data.n;i++)
    free(data.lf[i]);
  free(data.lf);

  return(0);

  //  VEC *gr = v_get(2);
  //Real f;
  //  gr = VMGMM_linear_eq(th,&data,gr,&f);
  // printf("%f\n",f);
  //v_output(gr);
  //exit(1);

}
