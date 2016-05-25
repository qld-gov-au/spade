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
#include <sys/time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include "meschach/matrix.h"
#include "meschach/matrix2.h"
#include "meschach/sparse.h"

#include "spade.h"
#include "common.h"
#include "socbio/fixed/w.h"
#include "socbio/fixed/s.h"
#include "VMGMM/VMGMM.h"
#include "optim/optim.h"
#include "util/util.h"
#include "VMGMM/solvers/alpha/ini_alpha.h"
#include "VMGMM/solvers/beta/ini_beta.h"
#include "VMGMM/solvers/gamma/ini_gamma.h"
#include "VMGMM/solvers/kappa/ini_kappa.h"
#include "VMGMM/solvers/omega/ini_omega.h"
#include "VMGMM/solvers/spade_solve.h"
#include "socbio/variable/growth.h"
#include "socbio/variable/birth.h"
#include "VMGMM/objfns/objfns.h"
#include "VMGMM/solvers/alpha/solve_p_alpha.h"
#include "VMGMM/solvers/beta/solve_p_beta.h"
#include "VMGMM/solvers/gamma/solve_p_gamma.h"
#include "VMGMM/solvers/iota/solve_p_iota.h"
#include "VMGMM/solvers/kappa/solve_p_kappa.h"
#include "VMGMM/solvers/omega/solve_p_omega.h"
#include "socbio/variable/catch.h"
#include "VMGMM/solvers/Q.h"
#include "VMGMM/solvers/Qn.h"
#include "VMGMM/solvers/reproduce.h"
#include "VMGMM/solvers/alpha/reproduce_alpha.h"
#include "VMGMM/solvers/kappa/reproduce_kappa.h"
#include "VMGMM/solvers/omega/reproduce_omega.h"
#include "VMGMM/solvers/zstar.h"
#include "mathprop/mathprop.h"

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
	    if (cnt->ve[i] > 100)
	      data.n += 1;

          data.lf = (double **) calloc(data.n,sizeof(double *));
	  data.t_id = (int *) calloc(data.n,sizeof(int));
	  data.t_sz = (int *) calloc(data.n,sizeof(int));

	  int kk=0;
	  for (int i=0;i<cnt->dim;i++)
	    {
	      if (cnt->ve[i] > 100)
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
	  //sscanf(argv[i+5],"%f",&kappa);
	  //sscanf(argv[i+6],"%f",&omega);
	  kappa = .1;
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

  bfgs(VMGMM,theta,&data);

  //char lab1[10]="before";
  //output_plots(theta,VMGMM,&data,lab);

  //char lab2[10]="after";
  //output_plots(theta,VMGMM,&data,lab2);

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



















/*
double ConditionNumber(

		       VEC * theta,
		       struct DATA *dataptr

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

  // p alpha
  MAT *p_a = m_get(x->m,x->n);
  Solve_Args solve_p_alpha_args;
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
  solve_p_gamma_args.Ui = Ui;
  solve_p_gamma_args.Uh = Uh;
  solve_p_gamma_args.Uhh = Uhh;
  solve_p_gamma_args.idxi = idxi;
  solve_p_gamma_args.eff = dataptr->eff;
  solve_p_gamma_args.k = dataptr->k;
  solve_p_gamma_args.S = dataptr->S;
 
  // p iota
  MAT *p_i = m_get(x->m,x->n);
  Solve_Args solve_p_iota_args;
  solve_p_iota_args.theta = theta;
  solve_p_iota_args.dataptr = dataptr;
  solve_p_iota_args.grad = grad;
  solve_p_iota_args.p = p_i;
  solve_p_iota_args.x = x;
  solve_p_iota_args.u = u;
  solve_p_iota_args.xhh = xhh;
  solve_p_iota_args.xh = xh;
  solve_p_iota_args.xn = xn;
  solve_p_iota_args.uh = uh;
  solve_p_iota_args.Ui = Ui;
  solve_p_iota_args.Uh = Uh;
  solve_p_iota_args.Uhh = Uhh;
  solve_p_iota_args.idxi = idxi;
  solve_p_iota_args.eff = dataptr->eff;
  solve_p_iota_args.k = dataptr->k;
  solve_p_iota_args.S = dataptr->S;

  VEC *theta_save = v_get(theta->dim);
  theta_save->ve[0]=theta->ve[0];
  theta_save->ve[1]=theta->ve[1];
  theta_save->ve[2]=theta->ve[2];
  theta_save->ve[3]=theta->ve[3];

  MAT *H = m_get(4,4);

  double epsilon = 1e-7;

  for (int j=1;j<=theta->dim;j++)
    {

      VEC *delta = v_get(theta->dim);
      delta->ve[j] = epsilon;

      v_add(theta_save,delta,theta);

      solve_p_alpha_args.theta = theta;
      solve_p_beta_args.theta = theta;
      solve_p_gamma_args.theta = theta;
      solve_p_iota_args.theta = theta;

      solve_p_alpha((void*)&solve_p_alpha_args);
      solve_p_beta((void*)&solve_p_beta_args);
      solve_p_gamma((void*)&solve_p_gamma_args);
      solve_p_iota((void*)&solve_p_iota_args);

      VEC *g1 = v_get(4);

      g1 = grad;

      delta->ve[j] = -epsilon;

      v_add(theta_save,delta,theta);

      solve_p_alpha_args.theta = theta;
      solve_p_beta_args.theta = theta;
      solve_p_gamma_args.theta = theta;
      solve_p_iota_args.theta = theta;

      solve_p_alpha((void*)&solve_p_alpha_args);
      solve_p_beta((void*)&solve_p_beta_args);
      solve_p_gamma((void*)&solve_p_gamma_args);
      solve_p_iota((void*)&solve_p_iota_args);

      VEC *g2 = v_get(4);
      g2 = grad;

      g1 = grad;

      H->me[i][j];

      }

  M_FREE(p_a);
  M_FREE(p_b);
  M_FREE(p_g);
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

  return cn;

}
*/


/*
void output_plots(

		  VEC *p,
		  struct DATA *d,
		  char *label

		  )
{

  int I = d->I+1;
  int J = d->J+1;
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

  solve(theta,x,u,xhh,xh,xn,uh,un,Ui,Uh,Uhh,idxi,d->eff,d->k,d->S);

  //  VEC *ctt = v_get(x->n);
  //  VEC *xt = v_get(x->n);

  /*
  FILE *p1 = fopen("plot1.txt","w");

  for (int i=0;i<x->m;i++)
    {

      xt = get_row(x,i,xt);
      for (int j=0;j<x->n;j++)
	ctt->ve[j] = s(x->me[i][j])*theta->ve[6]*e(d->eff,d->k,d->k*(i-d->S))*w(x->me[i][j])*u->me[i][j];
      fprintf(p1,"%f %f\n",d->k*(i-d->S),Q(xt,ctt)/1e3);

    }

  fclose(p1);

  FILE *p2 = fopen("plot2.txt","w");

  for (int i=0;i<x->m;i++) 
    fprintf(p2,"%f %f\n",d->k*(i-d->S),c(d->cat,d->k,d->k*(i-d->S)));

  fclose(p2);

  system("./plo > plotc.pdf");

  exit(1);
  */
