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
#include "VMGMM/VMGMM.h"
#include "optim/optim.h"

int feenableexcept(int);

double gtol = 0.9;

int main(int argc, char *argv[])
{
 
  feenableexcept(FE_DIVBYZERO); 
  feenableexcept(FE_INVALID); 
  feenableexcept(FE_OVERFLOW);

  float *ct; 
  float *ti; 
  float *ln;
  float *tl;
  int Nce,Nlf;
  int N;

  float beta,gamma,alpha,iota;

  int minfish = 250;

  J = 400;
  double k = 0.025;

  struct DATA data;

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
	  ct = (float *) calloc(Nce,sizeof(float));
	  ti = (float *) calloc(Nce,sizeof(float));

	  for (int i=0;i<Nce;i++)
	    fscanf(fp1,"%f %f", &ct[i],&ti[i]);

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

	  double cek = k/4;

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

	  ln = (float *) calloc(Nlf,sizeof(float));
	  tl = (float *) calloc(Nlf,sizeof(float));

	  for (int i=0;i<Nlf;i++)
	    fscanf(fp1,"%f %f", &ln[i],&tl[i]);

	  fclose(fp1);

	  VEC *lfv = v_get(Nlf);
	  VEC *ilv = v_get(Nlf);
	  VEC *cnt = v_get((int)2*Y/k);

	  for (int i=0;i<Nlf;i++)
	    {
	      lfv->ve[i] = (double)ln[i];
	      ilv->ve[i] = (int)(Y/k + floor((tl[i]+k/2)/k));
	      cnt->ve[(int)(Y/k + floor((tl[i]+k/2)/k))] += 1;
	    }

	  data.n = 0;
	  for (int i=0;i<cnt->dim;i++)	   
	    if (cnt->ve[i] > minfish)
	      data.n += 1;

          data.lf = (double **) calloc(data.n,sizeof(double *));
	  data.t_id = (int *) calloc(data.n,sizeof(int));
	  data.t_sz = (int *) calloc(data.n,sizeof(int));

	  int kk=0;
	  for (int i=0;i<cnt->dim;i++)
	    {
	      if (cnt->ve[i] > minfish)
		{

		  data.t_id[kk] = i;
		  data.t_sz[kk] = cnt->ve[i];		
		  data.lf[kk] = calloc(data.t_sz[kk],sizeof(double));

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

	  sscanf(argv[i+2],"%f",&alpha);
	  sscanf(argv[i+3],"%f",&beta);
	  sscanf(argv[i+4],"%f",&gamma);
	  sscanf(argv[i+5],"%f",&iota);
	  sscanf(argv[i+6],"%f",&kappa);
	  //sscanf(argv[i+6],"%f",&omega);
	  //	  kappa = .1;
	  omega = 160;

	  i += 7;

	} 
      else
	{
	  printf("problem\n");
          exit(1);
	}
    }
  }

  h = omega / J;

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

  double e_bar = v_sum(data.eff)/data.eff->dim;
  double iota = result->ve[1]/e_bar;

  //printf("%f %f\n",e_bar,est_iota);

  double a2 = alpha2;
  VEC *out = v_get(2);

  for (int i=0;i<10;i++) 
    { 

      out = calc_alpha2(alpha1,a2,kappa,omega,result->ve[0],result->ve[1]);
      a2 = a2 - out->ve[0]/out->ve[1];      
    }

  printf("%f %f\n",iota,a2);
  */

  VEC *theta = v_get(4);

  theta->ve[0] = alpha;
  theta->ve[1] = beta;
  theta->ve[2] = gamma;
  theta->ve[3] = iota;
  //theta->ve[4] = kappa;

  char lab1[10]="before";
  plot(theta,&data,lab1);

  theta = bfgs(VMGMM,theta,&data);

  char lab2[10]="after";

  plot(theta,&data,lab2);

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
  //double f;
  //  gr = VMGMM_linear_eq(th,&data,gr,&f);
  // printf("%f\n",f);
  //v_output(gr);
  //exit(1);

}
