// SPADE
// Stock assessment using PArtial Differential Equations
// Alex Campbell 'ghostofsandy' 2015 - 2016

#define PLOT 0
#define PLOTDERIV 0
#define PLOTSOLVP 0
#define GCHECK 0

#define A1 8.588e-5
#define A2 0.00144

#include <fenv.h>
#include <sys/time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "meschach/matrix.h"
#include "meschach/matrix2.h"
#include "meschach/sparse.h"

struct DATA {
  VEC *eff;
  VEC *cat;
  double **lf; 
  int n; 
  int *t_id; 
  int *t_sz;
  int I,J,S;
  double k;
};

float kappa,omega;

double iota1=5.2;
double iota2=0.619;
double phi=17;
double eta1=1.703205e-5;
double eta2=2.9526;

double w(double x) { return eta1*pow(x,eta2); }
double s(double x) 
{ 
  if (x<58)
    return 0;
  else if (x <=60)
    {
      double m = exp(-pow(60-phi*iota1,2.)/(2*iota2*pow(phi,2.)))/2;
      return m*(x-58);
    }
  else
    return exp(-pow(x-phi*iota1,2.)/(2*iota2*pow(phi,2.)));
}
double g(const double,const double,const double);
double b(const double,const double);

int mthls(VEC *(*f)(VEC *,struct DATA *,VEC *,double *),VEC *,double,VEC *,VEC *,double,double,double,double,double,double,int,struct DATA *); // More-Thuente line search taken from code by Nocedal and Dianne O'Leary
int cstep(double*,double*,double*,double*,double*,double*,double*,double,double,int*,double,double); // cstep from More-Thuente line search
VEC *bfgs(VEC * (*)(VEC *,struct DATA *,VEC *,double *),VEC *,struct DATA *);
MAT *UpdateHessian(MAT*,VEC*,VEC*);

void solve(VEC *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,VEC *,VEC *,VEC *,IVEC *,VEC *,double,int);
void solve_p_alpha(VEC *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,VEC *,VEC *,VEC *,IVEC *,VEC *,double,int);
void solve_p_beta  (VEC *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,VEC *,VEC *,VEC *,IVEC *,VEC *,double,int);
void solve_p_gamma (VEC *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,VEC *,VEC *,VEC *,IVEC *,VEC *,double,int);
void solve_p_delta (VEC *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,VEC *,VEC *,VEC *,IVEC *,VEC *,double,int);
void solve_p_kappa (VEC *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,VEC *,VEC *,VEC *,IVEC *,VEC *,double,int);
void solve_p_omega (VEC *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,VEC *,VEC *,VEC *,IVEC *,VEC *,double,int);
void solve_p_iota  (VEC *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,VEC *,VEC *,VEC *,IVEC *,VEC *,double,int);

double H(MAT *,MAT *,struct DATA *,double);
double G(MAT *,MAT *,MAT *,struct DATA *,double);
double G_ni(MAT*,MAT *,MAT *,struct DATA *,double);
double G_w(MAT*,MAT *,MAT *,struct DATA *,double);

VEC *VMGMM(VEC *,struct DATA *,VEC *,double *); // The model
VEC *VMGMM_eq(VEC *,struct DATA *,VEC *,double *); // equlibrium model
VEC *calc_alpha(double,double,double,double,double);

//double ConditionNumber(VEC *,struct DATA *);

double zstar(VEC *,double,double,double,double,double,double,double,double);
double e(VEC *,double,double);
double c(VEC *,double,double);

double Q(VEC *,VEC *);
double Qn(VEC *,VEC *,int);

void Q2(double,double,double,VEC *,VEC *);
void Q2_alpha(double,double,double,VEC *,VEC *,VEC *);
void Q2_kappa(double,double,double,double,VEC *,VEC *,VEC *);
void Q2_omega(double,double,double,double,VEC *,VEC *,VEC *);

VEC *initial(VEC *,VEC *,VEC *);
VEC *ini_alpha(VEC *,VEC *,VEC *);
VEC *ini_beta(VEC *,VEC *,VEC *);
VEC *ini_gamma(VEC *,VEC *,VEC *);
VEC *ini_omega(VEC *,VEC *,VEC *);

double idxselect(double,VEC *);
VEC *idxremove(VEC *,VEC *,int);

double get_bw(VEC *);

VEC *numgrad(double (*)(VEC *,void *),void *,VEC *,double); // numerical gradient. not used.

double h;
int J;

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

  V_FREE(theta);
  V_FREE(data.cat);
  V_FREE(data.eff);
  free(data.t_id);
  free(data.t_sz);

  for (int i=0;i<data.n;i++)
    free(data.lf[i]);;
  free(data.lf);

  return(0);

  //  VEC *gr = v_get(2);
  //double f;
  //  gr = VMGMM_linear_eq(th,&data,gr,&f);
  // printf("%f\n",f);
  //v_output(gr);
  //exit(1);

}

VEC *VMGMM(

		  VEC *theta,
		  struct DATA *dataptr,
		  VEC *grad,
		  double *f

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

  //  VEC *ctt = v_get(x->n);
  //  VEC *xt = v_get(x->n);

  /*
  FILE *p1 = fopen("plot1.txt","w");

  for (int i=0;i<x->m;i++)
    {

      xt = get_row(x,i,xt);
      for (int j=0;j<x->n;j++)
	ctt->ve[j] = s(x->me[i][j])*theta->ve[6]*e(dataptr->eff,dataptr->k,dataptr->k*(i-dataptr->S))*w(x->me[i][j])*u->me[i][j];
      fprintf(p1,"%f %f\n",dataptr->k*(i-dataptr->S),Q(xt,ctt)/1e3);

    }

  fclose(p1);

  FILE *p2 = fopen("plot2.txt","w");

  for (int i=0;i<x->m;i++) 
    fprintf(p2,"%f %f\n",dataptr->k*(i-dataptr->S),c(dataptr->cat,dataptr->k,dataptr->k*(i-dataptr->S)));

  fclose(p2);

  system("./plo > plotc.pdf");

  exit(1);
  */

  *f = H(x,u,dataptr,theta->ve[3]);

  printf("%g ",*f);

  MAT *p_a = m_get(x->m,x->n);
  solve_p_alpha(theta,p_a,x,u,xhh,xh,xn,uh,un,Ui,Uh,Uhh,idxi,dataptr->eff,dataptr->k,dataptr->S);
  grad->ve[0] = G_ni(p_a,x,u,dataptr,theta->ve[3]);

  printf("%g ",grad->ve[0]);

  MAT *p_b = m_get(x->m,x->n);
  solve_p_beta  (theta,p_b ,x,u,xhh,xh,xn,uh,Ui,Uh,Uhh,idxi,dataptr->eff,dataptr->k,dataptr->S);
  grad->ve[1]= G_ni(p_b,x,u,dataptr,theta->ve[3]);

  printf("%g ",grad->ve[1]);

  MAT *p_g = m_get(x->m,x->n);
  solve_p_gamma (theta,p_g,x,u,xhh,xh,xn,uh,Ui,Uh,Uhh,idxi,dataptr->eff,dataptr->k,dataptr->S);
  grad->ve[2] = G_ni(p_g,x,u,dataptr,theta->ve[3]);

  printf("%g ",grad->ve[2]);

  /*
  MAT *p_k = m_get(x->m,x->n);
  solve_p_kappa (theta,p_k ,x,u,xhh,xh,xn,uh,un,Ui,Uh,Uhh,idxi,dataptr->eff,dataptr->k,dataptr->S);
  grad->ve[4] = G_ni(p_k,x,u,dataptr,theta->ve[6]);

  MAT *p_w = m_get(x->m,x->n);
  solve_p_omega (theta,p_w ,x,u,xhh,xh,xn,uh,un,Ui,Uh,Uhh,idxi,dataptr->eff,dataptr->k,dataptr->S);
  grad->ve[5] = G_ni(p_w,x,u,dataptr,theta->ve[6]);
  */

  MAT *p_i = m_get(x->m,x->n);
  solve_p_iota(theta,p_i,x,u,xhh,xh,xn,uh,Ui,Uh,Uhh,idxi,dataptr->eff,dataptr->k,dataptr->S);
  grad->ve[3] = G(p_i,x,u,dataptr,theta->ve[3]);

  printf("%g ",grad->ve[3]);

  printf("%g ",theta->ve[0]);  printf("%g ",theta->ve[1]);  printf("%g ",theta->ve[2]);  printf("%g\n",theta->ve[3]);

  M_FREE(p_a);
  M_FREE(p_b);
  M_FREE(p_g);
  //M_FREE(p_k);
  //M_FREE(p_w);
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
  //V_FREE(ctt);
  //V_FREE(xt);

  return grad;

}

VEC * bfgs(

	     VEC * (*model)(VEC *,struct DATA *,VEC *,double *),
	     VEC *x,
	     struct DATA *dataptr

	       )	       
{

  VEC *oldx = v_get(x->dim);
  VEC *oldgrad = v_get(x->dim);
  //  VEC *x_new = v_get(x->dim);
  VEC *delta_x = v_get(x->dim);
  VEC *delta_grad = v_get(x->dim);
  VEC *grad = v_get(x->dim);
  double f;
  int stop,iter;

  MAT *B = m_get(x->dim,x->dim);
  stop=0; iter=0;

  m_ident(B);

  //v_output(x);
  grad = (*model)(x,dataptr,grad,&f); 
  //  exit(1);

  //printf("%f ",f);
  //v_output(grad);

  /*
  VEC *nx = v_get(x->dim);

  VEC *sd = v_get(x->dim);
  mv_mlt(HH,grad,sd);
  sv_mlt(-1.0/v_norm2(sd),sd,sd);

  v_output(sd);

  v_add(x,sv_mlt(9.99e-7,sd,VNULL),nx);

  VEC *grad2 = v_get(x->dim);
  double f2=0;

  v_output(nx);
  grad2 = (*model)(nx,dataptr,grad2,&f2);
  v_output(grad2);

  double dg = in_prod(grad2,sd);
  printf("dg: %g\n",dg);

  */

  /*
  printf("\nChecking G def d H / d beta \n\n");
  printf("\ng: %f\n",grad->ve[3]);

  int I = dataptr->I+1;
  int J = dataptr->J+1;
  MAT *xx = m_get(I,J);
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

  solve(x,xx,u,xhh,xh,xn,uh,un,Ui,Uh,Uhh,idxi,dataptr->eff,dataptr->k,dataptr->S);

  double h = H(xx,u,dataptr,x->ve[3]);

  printf("h=%g\n",h);

  double ng = 0;     
  double save = x->ve[3];
  for (int j=-2;j>-25;j--)
    {
	  double delta = exp((double)j);
	  x->ve[3] = save + delta;

	  solve(x,xx,u,xhh,xh,xn,uh,un,Ui,Uh,Uhh,idxi,dataptr->eff,dataptr->k,dataptr->S);

	  double h_d = H(xx,u,dataptr,x->ve[3]);
	  ng = (h_d-h)/delta;
	  printf("%g %g\n",delta,ng);
    }

  printf("%g %g %g\n",grad->ve[3],grad->ve[3]-ng,(grad->ve[3]-ng)/ng);

  exit(1);  
  */
  //  MAT *H = m_get(x->dim,x->dim);

  while (1)
    {
      iter=iter+1;
      if (iter == 3)
	break;

      VEC *dir = v_get(x->dim);
      mv_mlt(B,grad,dir);
      sv_mlt(-1.0/v_norm2(dir),dir,dir);

      //v_output(dir);
      double dirD = in_prod(grad,dir);
      printf("dirD: %f\n",dirD);

      if (fabs(dirD) < 1)
	return(x);
      
      v_copy(x,oldx);
      v_copy(grad,oldgrad);

      /*
      if (fabs(dirD)<1200000) 
	{

	  int I = dataptr->I+1;
	  int J = dataptr->J+1;
	  MAT *y = m_get(I,J);
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

	  solve(x,y,u,xhh,xh,xn,uh,un,Ui,Uh,Uhh,idxi,dataptr->eff,dataptr->k,dataptr->S);
	  double h = H(y,u,dataptr,x->ve[3]);

	  printf("h=%g\n",h);

	  double ng = 0;
	  double save = x->ve[3];
	  for (int j=-3;j>-15;j--)
	    {
	      double delta = exp((double)j);
	      x->ve[3] = save + delta;
	      solve(x,y,u,xhh,xh,xn,uh,un,Ui,Uh,Uhh,idxi,dataptr->eff,dataptr->k,dataptr->S);

	      double h_d = H(y,u,dataptr,x->ve[3]);
	      ng = (h_d-h)/delta;
	      printf("%g %g\n",delta,ng);
	    }

	  printf("%g %g %g\n",grad->ve[3],grad->ve[3]-ng,(grad->ve[3]-ng)/ng);
	  exit(1);
	  /*
	  printf("\n");
	  for (int i=0;i<10;i++) 
	    {

	      double step = i*1e-12;
	      VEC *newx = v_get(dir->dim);	      
	      v_add(x,sv_mlt(-step,dir,VNULL),newx);

	      int I = dataptr->I+1;
	      int J = dataptr->J+1;
	      MAT *y = m_get(I,J);
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

	      solve(newx,y,u,xhh,xh,xn,uh,un,Ui,Uh,Uhh,idxi,dataptr->eff,dataptr->k,dataptr->S);
	      double f = H(y,u,dataptr,newx->ve[3]);
	      printf("%f %f\n",step,f);
	    }
	    exit(1);*/
      //}*/

      int rt = mthls( model,x,f,grad,dir,1e-8,1e-4,0.9,DBL_EPSILON,1e-20,1e20,60,dataptr);
    
      //printf("\n%d\n",rt);
      //v_output(x);       
      //v_output(grad);

      v_sub(x,oldx,delta_x);
      v_sub(grad,oldgrad,delta_grad);

      //      grad = (*model)(x,dataptr,grad,&f);
      
      B = UpdateHessian(B,delta_x,delta_grad);

      //m_output(B);

      V_FREE(dir);


    }

  V_FREE(oldx);
  V_FREE(oldgrad);
  V_FREE(delta_x);
  V_FREE(delta_grad);
  V_FREE(grad);
  M_FREE(B);

  return(x);
}

int mthls(

	  VEC *(*fcn)(VEC *,struct DATA *,VEC *,double *),
	  VEC *x,
	  double f,
	  VEC *gr,
	  VEC *sd,
	  double stp,
	  double ftol,
	  double gtol,
	  double xtol,
	  double stpmin,
	  double stpmax,
	  int maxfev,
	  struct DATA *d

	  )
{
  
  int xtrapf = 4;
  int info = 0;
  int infoc = 1;
  double dginit = in_prod(gr,sd); // g's must be < 0 (initial gradient in search direction must be descent)
  int brackt = 0;
  int stage1 = 1;
  int nfev = 0;
  double finit = f;
  double dgtest = ftol*dginit;
  double width = stpmax - stpmin;
  double width1 = 2*width;
  VEC *wa = v_get(x->dim);
  v_copy(x,wa);

  double stx = 0;
  double fx = finit;
  double dgx = dginit;
  double sty = 0;
  double fy = finit;
  double dgy = dginit;

  while (1) 
    {

      //Set the minimum and maximum steps to correspond
      // to the present interval of uncertainty.

      //      printf("%f\n",f);

      double stmin,stmax;

      if (brackt)
	{
	  stmin = min(stx,sty);
	  stmax = max(stx,sty);
	}
      else
	{
	  stmin = stx;
	  stmax = stp + xtrapf*(stp - stx);
	}

      // Force the step to be within the bounds stpmax and stpmin.

      stp = max(stp,stpmin);
      stp = min(stp,stpmax);

      // If an unusual termination is to occur then let stp be the lowest point obtained so far.

      if ((brackt && (stp <= stmin || stp >= stmax)) || nfev >= maxfev-1 || infoc == 0 || (brackt && stmax-stmin <= xtol*stmax))
	stp = stx;

      // Evaluate the funtion and gradient at stp and compute the directional derivative

      VEC *vtmp = v_get(x->dim);
      vtmp = sv_mlt(stp,sd,vtmp);
      v_add(wa,vtmp,x);
      V_FREE(vtmp);
      gr = (*fcn)(x,d,gr,&f);
      nfev += 1;
      double dg = in_prod(gr,sd);
      double ftest1 = finit + stp*dgtest;

      // Test for convergence.

      if ((brackt && (stp <= stmin || stp >= stmax)) || infoc == 0)
	info = 6;

      if (stp == stpmax && f <= ftest1 && dg <= dgtest)
	info = 5;

      if (stp == stpmin && (f > ftest1 || dg >= dgtest))
	info = 4;

      if (nfev >= maxfev)
	info = 3;

      if (brackt && stmax - stmin <= xtol*stmax)
	info = 2;

      if (f <= ftest1 && fabs(dg) <= gtol*(-dginit))
	info = 1;

      if (info != 0){
	printf("%d\n",info);
	V_FREE(wa);
	return info;
      }

      // In the first stage we seek a step for which the modified function has a nonpositive value and nonnegative derivative.

      if (stage1 && f <= ftest1 && dg >= min(ftol,gtol)*dginit)
	stage1 = 0;

      //A modified function is used to predict the step only if we have not obtained a step for which the modified function has a nonpositive function value and nonnegative derivative, and if a lower function value has been obtained but the decrease is not sufficient.

      if (stage1 && f <= fx && f > ftest1)
	{

	  // Define the modified function and derivative values.

	  double fm = f - stp*dgtest;
	  double fxm = fx - stx*dgtest;
	  double fym = fy - sty*dgtest;
	  double dgm = dg - dgtest;
	  double dgxm = dgx - dgtest;
	  double dgym = dgy - dgtest;
 
	  // Call cstep to update the interval of uncertainty and to compute the new step.
	  infoc = cstep(&stx,&fxm,&dgxm,&sty,&fym,&dgym,&stp,fm,dgm,&brackt,stmin,stmax);

	  // Reset the function and gradient values for f.

	  fx = fxm + stx*dgtest;
	  fy = fym + sty*dgtest;
	  dgx = dgxm + dgtest;
	  dgy = dgym + dgtest;

	}
      else
	{
	  infoc = cstep(&stx,&fx,&dgx,&sty,&fy,&dgy,&stp,f,dg,&brackt,stmin,stmax);
	}

      // Force a sufficient decrease in the size of the interval of uncertainty

      if (brackt)
	{

	  if (fabs(sty-stx) >= .66*width1)
	    stp = stx + .5*(sty - stx);
	  width1 = width;
	  width = fabs(sty-stx);
	}

      //      printf("%g ",f);
    }

  //	 printf("%d %g %g %g %g %g %d\n",info,stx,stp,sty,f,dg,infoc);
  //printf("%g ",f);

}

int cstep(

	  double *stx,
	  double *fx,
	  double *dx,
	  double *sty,
	  double *fy,
	  double *dy,
	  double *stp,
	  double fp,
	  double dp,
	  int *brackt,
	  double stpmin,
	  double stpmax

	  )
{

  int info = 0;

  // Determine if the derivatives have opposite sign.

  double sgnd = dp*((*dx)/fabs((*dx)));

  // First case. A higher function value.
  // The minimum is bracketed. If the cubic step is closer to stx than the quadratic step, the cubic step is taken, else the average of the cubic and quadratic steps is taken.

  int bound;
  double theta;
  double ss;
  double gamma;
  double stpf;

  if (fp > (*fx)) 
    {

      info = 1;
      bound = 1;
      theta = 3*((*fx) - fp)/((*stp) - (*stx)) + (*dx) + dp;
      VEC *tmp = v_get(3);
      tmp->ve[0] = theta;
      tmp->ve[1] = (*dx);
      tmp->ve[2] = dp;
      ss = v_norm_inf(tmp);
      V_FREE(tmp);
      gamma = ss*sqrt(pow(theta/ss,2.) - ((*dx)/ss)*(dp/ss));
      if ((*stp) < (*stx)) 
	gamma = -gamma;

      double p = (gamma - (*dx)) + theta;
      double q = ((gamma - (*dx)) + gamma) + dp;
      double r = p/q;
      double stpc = (*stx) + r*((*stp) - (*stx));
      double stpq = (*stx) + (((*dx)/(((*fx)-fp)/((*stp)-(*stx))+(*dx)))/2)*((*stp) - (*stx));
      if (fabs(stpc-(*stx)) < fabs(stpq-(*stx)))
	stpf = stpc;
      else
	stpf = stpc + (stpq - stpc)/2;
         
      *brackt = 1;

    } 

  // Second case. A lower function value and derivatives of opposite sign. The minimum is bracketed. If the cubic step is closer to stx than the quadratic (secant) step, the cubic step is taken, else the quadratic step is taken.

  else if (sgnd < 0.0) 
    {
      info = 2;
      bound = 0;

      theta = 3*((*fx) - fp)/((*stp) - (*stx)) + (*dx) + dp;
      VEC *tmp = v_get(3);
      tmp->ve[0] = theta;
      tmp->ve[1] = (*dx);
      tmp->ve[2] = dp;
      ss = v_norm_inf(tmp);
      V_FREE(tmp);
      gamma = ss*sqrt(pow(theta/ss,2.) - ((*dx)/ss)*(dp/ss));

      if ((*stp) > (*stx)) 
	gamma = -gamma;
         
      double p = (gamma - dp) + theta;
      double q = ((gamma - dp) + gamma) + (*dx);
      double r = p/q;

      double stpc = (*stp) + r*((*stx) - (*stp));
      double stpq = (*stp) + (dp/(dp-(*dx)))*((*stx) - (*stp));
      if (fabs(stpc-(*stp)) > fabs(stpq-(*stp)))
	stpf = stpc;
      else
	stpf = stpq;
         
      *brackt = 1;
    }

  // Third case. A lower function value, derivatives of the same sign, and the magnitude of the derivative decreases. The cubic step is only used if the cubic tends to infinity in the direction of the step or if the minimum of the cubic is beyond stp. Otherwise the cubic step is defined to be either stpmin or stpmax. The quadratic (secant) step is also computed and if the minimum is bracketed then the the step closest to stx is taken, else the step farthest away is taken.

  else if (fabs(dp) < fabs(*dx)) 
    {

      info = 3;
      bound = 1;

      theta = 3*((*fx) - fp)/((*stp) - (*stx)) + (*dx) + dp;
      VEC *tmp = v_get(3);
      tmp->ve[0] = theta;
      tmp->ve[1] = (*dx);
      tmp->ve[2] = dp;
      ss = v_norm_inf(tmp);
      V_FREE(tmp);
      // The case gamma = 0 only arises if the cubic does not tend to infinity in the direction of the step.

      gamma = ss*sqrt(max(0.,pow(theta/ss,2.) - (*dx/ss)*(dp/ss)));

      if ((*stp) > (*stx)) 
	gamma = -gamma;
         
      double p = (gamma - dp) + theta;
      double q = (gamma + ((*dx) - dp)) + gamma;
      double r = p/q;

      double stpc,stpq;

      if (r < 0.0 && gamma != 0.0)
	stpc = (*stp) + r*((*stx) - (*stp));
      else if ((*stp) > (*stx))
	stpc = stpmax;
      else
	stpc = stpmin;

      stpq = (*stp) + (dp/(dp-(*dx)))*((*stx) - (*stp));
      if (*brackt)
	if (fabs((*stp)-stpc) < fabs((*stp)-stpq))
	  stpf = stpc;
	else
	  stpf = stpq;           
      else
	if (fabs((*stp)-stpc) > fabs((*stp)-stpq))
	  stpf = stpc;
	else
	  stpf = stpq;
    }

  // Fourth case. A lower function value, derivatives of the same sign, and the magnitude of the derivative does not decrease. If the minimum is not bracketed, the step is either stpmin or stpmax, else the cubic step is taken.

  else
    {
      info = 4;
      bound = 0;

      if (*brackt)
	{

	  theta = 3*(fp - (*fy))/((*sty) - *stp) + (*dy) + dp;
	  VEC *tmp = v_get(3);
	  tmp->ve[0] = theta;
	  tmp->ve[1] = (*dy);
	  tmp->ve[2] = dp;
	  ss = v_norm_inf(tmp);
	  V_FREE(tmp);

	  gamma = ss*sqrt(pow(theta/ss,2.) - ((*dy)/ss)*(dp/ss));

	  if (*stp > (*sty))
	    gamma = -gamma;
            
          double p = (gamma - dp) + theta;
          double q = ((gamma - dp) + gamma) + (*dy);
          double r = p/q;
          double stpc = *stp + r*((*sty) - *stp);
          stpf = stpc;
	}
      else if (*stp > (*stx))
	stpf = stpmax;
      else
	stpf = stpmin;
    
    }

  // Update the interval of uncertainty. This update does not depend on the new step or the case analysis above.

  if (fp > *fx)
    {
      *sty = *stp;
      *fy = fp;
      *dy = dp;
    }
  else
    {
      if (sgnd < 0.0)
	{
	  *sty = *stx;
	  *fy = *fx;
	  *dy = *dx;
	}
          
      *stx = *stp;
      *fx = fp;
      *dx = dp;
    }

  // Compute the new step and safeguard it.

  stpf = min(stpmax,stpf);
  stpf = max(stpmin,stpf);
  *stp = stpf;
  if (*brackt && bound)
    if (*sty > *stx) 
      *stp = min(*stx+.66*(*sty-*stx),*stp);
    else
      *stp = max(*stx+.66*(*sty-*stx),*stp);
         
  return info;

  // last card of subroutine cstep

}

MAT *UpdateHessian(

		   MAT* H, 
		   VEC* s,
		   VEC* y

		   )
{

  MAT *sMr = m_get(1,H->n);
  vm_move(s,0,sMr,0,0,1,H->n);
  MAT *sMc = m_get(H->m,1);
  vm_move(s,0,sMc,0,0,H->m,1);

  MAT *yMr = m_get(1,H->n);
  vm_move(y,0,yMr,0,0,1,H->n);
  MAT *yMc = m_get(H->m,1);
  vm_move(y,0,yMc,0,0,H->m,1);

  MAT *eye = m_get(H->m,H->n);
  eye = m_ident(eye);

  double rhok =  1/in_prod(y,s);

  MAT *sMc_yMr = m_get(H->n,H->n);
  sMc_yMr = m_mlt(sMc,yMr,sMc_yMr);

  MAT *rhok_sMc_yMr = m_get(H->n,H->n);
  rhok_sMc_yMr = sm_mlt(rhok,sMc_yMr,rhok_sMc_yMr);

  MAT *AA1 = m_get(H->n,H->n);
  AA1 = m_sub(eye,rhok_sMc_yMr,AA1);

  MAT *yMr_sMc = m_get(H->n,H->n);
  yMr_sMc = m_mlt(yMc,sMr,yMr_sMc);

  MAT *rhok_yMr_sMc  = m_get(H->n,H->n);
  rhok_yMr_sMc = sm_mlt(rhok,yMr_sMc,rhok_yMr_sMc);

  MAT *AA2 = m_get(H->n,H->n);
  AA2 = m_sub(eye,rhok_yMr_sMc,AA2);

  MAT *HAA2 = m_get(H->n,H->n);
  HAA2 = m_mlt(H,AA2,HAA2);

  MAT *HAA12 = m_get(H->n,H->n);
  HAA12 = m_mlt(AA1,HAA2,HAA12);

  MAT *sMc_sMr = m_get(H->n,H->n);
  sMc_sMr = m_mlt(sMc,sMr,sMc_sMr);

  sMc_sMr = sm_mlt(rhok,sMc_sMr,sMc_sMr);

  H = m_add(HAA12,sMc_sMr,H);

  M_FREE(HAA2);
  M_FREE(HAA12);
  M_FREE(sMr);
  M_FREE(sMc);
  M_FREE(yMr);
  M_FREE(yMc);
  M_FREE(eye);
  M_FREE(AA1);
  M_FREE(AA2);
  M_FREE(sMc_yMr);
  M_FREE(rhok_sMc_yMr);
  M_FREE(yMr_sMc);
  M_FREE(rhok_yMr_sMc);
  M_FREE(sMc_sMr);

  return H;

  //H =  m_add(m_mlt(AA1,m_mlt(H,AA2,MNULL),MNULL),sm_mlt(rhok,m_mlt(sMc,sMr,MNULL),MNULL),MNULL);

}

VEC *idxremove(

	       VEC *zn,
	       VEC *z,
	       int idx

	       )

{

  for (int i=0;i<idx;i++)
    z->ve[i]=zn->ve[i];

  for (int i=idx;i<z->dim;i++)
    z->ve[i]=zn->ve[i+1];

  return z;
}

double idxselect(

		 double omega,
		 VEC *xn

		 )
{
  int idx = -1;
  double val = omega;
  for (int i=1;i<xn->dim-1;i++)
    {
      double dif = xn->ve[i+1]-xn->ve[i-1];
      if (dif < val)
	{
	  val = dif;
	  idx = i;
	}
    }

  return idx;

}
  
VEC *ini_alpha(

		VEC *theta,
		VEC *x,
		VEC *p

		)
{

  x->ve[x->dim-1] -= 1e-5;

  double a = theta->ve[0];
  double b = theta->ve[1];
  double g = theta->ve[2];
  double k = theta->ve[3];
  double w = theta->ve[4];

  double zeta = sqrt( 81*k*k*w*w*pow(a*A1+2*a*A2*w,2.) - 12*k*pow(a*A1*w+k,3.) );
  double eta = 9*a*A1*k*k*w + 18*a*A2*k*k*w*w + k*zeta;
  double eta2 = 9*A1*k*k*w + 18*A2*k*k*w*w + k*zeta;
  double Z = pow(eta,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*k*(a*A1*w+k) / pow(eta,1./3);

  double Za = ( k*k*w*(3*A1*zeta+6*A2*w*zeta-6*A1*pow(k+a*A1*w,2.)+27*k*w*(A1+2*A2*w)*(a*A1+2*a*A2*w)) / ( zeta*pow(eta2,2/3.) ) ) * ( (1/ (pow(2,1/3.) * pow(3,2/3.) )) - ( pow(2/3.,1/3.)*k*(k+a*A1*w) / pow(eta2,2/3.) ) ) + ( pow(2/3.,1/3.)*A1*k*w / pow(eta,1/3.) );

  for (int j=0;j<x->dim;j++)
    p->ve[j] = pow(1-x->ve[j]/w,Z/k - 2.) * ( (Z-b-k)/(g*Z) * (A1 + 2*A2*k*w/(Z+k) - 2*a*A2*k*w/pow(Z+k,2.) * Za) + (Za/(g*Z))*(a*A1 + 2*a*A2*k*w/(Z+k))*((b+k)/Z + log(1 - x->ve[j]/w)* (Z-b-k)/k) );

}

VEC *ini_beta(

	      VEC *theta,
	      VEC *x,
	      VEC *p

	      )
{

  x->ve[x->dim-1] -= 1e-5;

  double a = theta->ve[0];
  double b = theta->ve[1];
  double g = theta->ve[2];
  double k = theta->ve[3];
  double w = theta->ve[4];

  double zeta = sqrt( 81*k*k*w*w*pow(a*A1+2*a*A2*w,2.) - 12*k*pow(a*A1*w+k,3.) );
  double eta = 9*a*A1*k*k*w + 18*a*A2*k*k*w*w + k*zeta;
  double Z = pow(eta,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*k*(a*A1*w+k) / pow(eta,1./3);

  for (int j=0;j<x->dim;j++)
    p->ve[j] = -pow(1-x->ve[j]/w,Z/k - 2.) * (1/(g*Z)) * (a*A1 + 2*a*A2*k*w/(Z+k));

}

VEC *ini_gamma(

	       VEC *theta,
	       VEC *x,
	       VEC *p

		)
{

  x->ve[x->dim-1] -= 1e-5;

  double a = theta->ve[0];
  double b = theta->ve[1];
  double g = theta->ve[2];
  double k = theta->ve[3];
  double w = theta->ve[4];

  double zeta = sqrt( 81*k*k*w*w*pow(a*A1+2*a*A2*w,2.) - 12*k*pow(a*A1*w+k,3.) );
  double eta = 9*a*A1*k*k*w + 18*a*A2*k*k*w*w + k*zeta;
  double Z = pow(eta,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*k*(a*A1*w+k) / pow(eta,1./3);

  for (int j=0;j<x->dim;j++)
    p->ve[j] = -pow(1-x->ve[j]/w,Z/k - 2.) * ((Z-b-k)/(1e-7*pow(g,2.)*Z)) * (a*A1 + 2*a*A2*k*w/(Z+k));

}

VEC *ini_kappa(

	       VEC *theta,
	       VEC *x,
	       VEC *p

	       )
{

  // not updated for new birth function yet

  x->ve[x->dim-1] -= 1e-5;

  double a1 = theta->ve[0];
  double a2 = theta->ve[1];
  double b = theta->ve[2];
  double g = theta->ve[3];
  double k = theta->ve[4];
  double w = theta->ve[5];

  double zeta = sqrt( 81*k*k*w*w*pow(a1+2*a2*w,2.) - 12*k*pow(a1*w+k,3.) );
  double in = 9*a1*k*k*w + 18*a2*k*k*w*w + k*zeta;
  double Z = pow(in,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*k*(a1*w+k) / pow(in,1./3);

  double Zk =  6*k*(a1*w*zeta+2*a2*w*w*zeta-(2*k+a1*w)*pow(k+a1*w,2.)+9*k*w*w*pow(a1+2*a2*w,2.) ) / (zeta*pow(9*a1*k*k*w+18*a2*k*k*w*w+k*zeta,2./3)) * ( 1./(pow(2.,1/3.)*pow(3.,2/3.)) - pow(2/3.,1/3.)*k*(k+a1*w) / pow(9*a1*k*k*w+18*a2*k*k*w*w+k*zeta,2./3) ) + pow(2/3.,1/3.)*(2*k+a1*w) / pow(9*a1*k*k*w+18*a2*k*k*w*w+k*zeta,1./3) ; // check order of precedence

  for (int j=0;j<x->dim;j++)
    p->ve[j] = pow(1-x->ve[j]/w,Z/k - 2.) * ( Zk/(g*Z) *( (a1 + (2*a2*k*w/(Z+k))) * ( (b+k)/Z + log(1-x->ve[j]/w) * (Z-b-k)/k ) - 2*a2*k*w*(Z-b-k)/pow(Z+k,2.) ) + ((Z-b-k)/(g*Z*(Z+k))) * (2*a2*w + 2*a2*k*w/(Z+k)) - (a1 + 2*a2*k*w/(Z+k)) * (1/(g*Z) + log(1-x->ve[j]/w) * (Z-b-k)/(g*k*k)) ); 

}

VEC *ini_omega(

	       VEC *theta,
	       VEC *x,
	       VEC *p

		)
{

  // not updated for new birth function yet

  x->ve[x->dim-1] -= 1e-5;

  double a1 = theta->ve[0];
  double a2 = theta->ve[1];
  double b = theta->ve[2];
  double g = theta->ve[3];
  double k = theta->ve[4];
  double w = theta->ve[5];

  double zeta = sqrt( 81*k*k*w*w*pow(a1+2*a2*w,2.) - 12*k*pow(a1*w+k,3.) );
  double in = 9*a1*k*k*w + 18*a2*k*k*w*w + k*zeta;
  double Z = pow(in,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*k*(a1*w+k) / pow(in,1./3);

  double Zw =  (3*k*k*(a1*zeta+4*a2*w*zeta-2*a1*pow(k+a1*w,2.)+18*a2*k*w*w*(a1+2*a2*w)+9*k*w*pow(a1+2*a2*w,2.) ) / (zeta*pow(9*a1*k*k*w+18*a2*k*k*w*w+k*zeta,2./3)) ) * ( 1./(pow(2.,1/3.)*pow(3.,2/3.)) - pow(2/3.,1/3.)*k*(k+a1*w) / pow(9*a1*k*k*w+18*a2*k*k*w*w+k*zeta,2./3) ) + pow(2/3.,1/3.)*a1*k / pow(9*a1*k*k*w+18*a2*k*k*w*w+k*zeta,1./3) ;

  for (int j=0;j<x->dim;j++)
    p->ve[j] = pow(1-x->ve[j]/w,Z/k - 2.) * ( (Z-b-k)/(g*Z) *( 2*k*w/(Z+k) - (2*a2*k*w/pow(Z+k,2.))*Zw + (a1 + 2*a2*k*w/(Z+k))*x->ve[j]*(Z-2*k)/(k*w*(w-x->ve[j]))) + (Zw/(g*Z)) * (a1+2*a2*k*w/(Z+k))*( (b+k)/Z + log(1-x->ve[j]/w) * (Z-b-k)/k ) );

}

VEC *initial(

	     VEC *theta,
	     VEC *x,
	     VEC *u

	     )

{

  x->ve[x->dim-1] -= 1e-5;

  double a = theta->ve[0];
  double b = theta->ve[1];
  double g = theta->ve[2];
  double k = theta->ve[3];
  double w = theta->ve[4];

  double zeta = sqrt( 81*k*k*w*w*pow(a*A1+2*a*A2*w,2.) - 12*k*pow(a*A1*w+k,3.) );
  double eta = 9*a*A1*k*k*w + 18*a*A2*k*k*w*w + k*zeta;
  double Z = pow(eta,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*k*(a*A1*w+k) / pow(eta,1./3);

  double ubar = (Z - b - k) / g; 
  double vbar = (k*w*ubar) / (b+g*ubar+k);
  double wbar = (2*k*w*vbar) / (b+g*ubar+2*k);  

  for (int j=0;j<x->dim;j++) 
    u->ve[j] = (a*A1*vbar+a*A2*wbar)*pow(w-x->ve[j],(b+g*ubar)/k-1) / (k*pow(w,(b+g*ubar)/k));

  return u;
}

double g(

	 const double k,
	 const double w,
	 const double x

	 )
{ /* von-Bertalanffy growth */
  return k*(w - x);
}

double b(

	 const double a,
	 const double x

	 )
{ /* birth function */
  return a*(A1*x + A2*pow(x,2.));
}
      
double H(

		 MAT *x,
		 MAT *u,
		 struct DATA *data,
		 double iota
	
		 )
{

  iota *= 1e-3;
  int S = data->S;
  double k = data->k;

  int lfi=0;

  VEC *ht = v_get(x->m);
  VEC *tt = v_get(x->m);

  for (int i=S;i<x->m;i++)
    {

      VEC *xt = v_get(x->n);
      get_row(x,i,xt);      
      VEC *ut = v_get(x->n);
      get_row(u,i,ut);      
      VEC *v = v_get(x->n);

      for (int j=0;j<x->n;j++)
	v->ve[j] = iota*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];

      if(lfi < data->n && data->t_id[lfi]==i) 
	{
      
	  VEC *dt = v_get(data->t_sz[lfi]);

	  for (int j=0;j<dt->dim;j++)
	    dt->ve[j] = data->lf[lfi][j];

	  double bw = get_bw(dt);

	  VEC *l = v_get(xt->dim);

	  for (int j=0;j<xt->dim;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  double al = 1e3*c(data->cat,k,k*(i - S)) / Q(xt,l);

	  if (PLOT) {

	    FILE *p1 = fopen("plot1.txt","w");

	    for (int j=0;j<x->n;j++)
	      fprintf(p1,"%f %f\n",xt->ve[j],al*l->ve[j]);

	    fclose(p1);

	    FILE *p2 = fopen("plot2.txt","w");

	    for (int j=0;j<x->n;j++)
	      fprintf(p2,"%f %f\n",xt->ve[j],v->ve[j]);

	    fclose(p2);

	    char buffer[100];

	    sprintf(buffer,"./plo > plotl%.3f.pdf",k*(data->t_id[lfi] - S));

	    system(buffer); 

	  }

	  VEC *ld = v_get(x->n);

	  for (int j=0;j<xt->dim;j++)
	    ld->ve[j] = pow(v->ve[j] - al*l->ve[j],2.);

	  ht->ve[i] = Q(xt,ld);

	  lfi += 1;

	  V_FREE(dt);
	  V_FREE(l);
	  V_FREE(ld);

	} 
      else 
	{
	  ht->ve[i] = pow(Q(xt,v)-1e3*c(data->cat,k,k*(i - S)),2.);
	}

      tt->ve[i] = k*(i-S);

      V_FREE(xt);
      V_FREE(ut);
      V_FREE(v);

    }

  if (PLOT) {

    FILE *p1 = fopen("plot.txt","w");

    for (int i=S;i<x->m;i++)
      fprintf(p1,"%f %f\n",tt->ve[i],ht->ve[i]);

    fclose(p1);

    char buffer[100];

    sprintf(buffer,"./plo1 > plotht.pdf");

    system(buffer);

  }

  double blah = Q(tt,ht);

  V_FREE(ht);
  V_FREE(tt);

  return blah; 

}

double G(

	 MAT *p,
	 MAT *x,
	 MAT *u,
	 struct DATA *data,
	 double iota
	
	 )
{

  iota *=1e-3;
  int lfi=0;

  int S = data->S;
  double k = data->k;

  VEC *ht = v_get(x->m);
  VEC *tt = v_get(x->m);

  for (int i=S;i<x->m;i++)
    {

      VEC *xt = v_get(x->n);
      get_row(x,i,xt);      
      VEC *ut = v_get(x->n);
      get_row(u,i,ut);      
      VEC *pt = v_get(x->n);
      get_row(p,i,pt); 

      VEC *v = v_get(x->n);
      VEC *pv = v_get(x->n);

      for (int j=0;j<x->n;j++) 
	{
	  pv->ve[j] = iota*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*pt->ve[j] + 1e-3*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];
	  v->ve[j] = iota*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];
	}

      //      if(data->t_id[lfi]==i) 
      if(lfi < data->n && data->t_id[lfi]==i) 
	{
      
	  VEC *dt = v_get(data->t_sz[lfi]);

	  for (int j=0;j<dt->dim;j++)
	    dt->ve[j] = data->lf[lfi][j];

	  double bw = get_bw(dt);

	  VEC *l = v_get(xt->dim);

	  for (int j=0;j<xt->dim;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  double al = 1e3*c(data->cat,k,k*(i - S)) / Q(xt,l);

	  if (PLOTDERIV) {

	    FILE *p1 = fopen("plot1.txt","w");

	    for (int j=0;j<x->n;j++)
	      fprintf(p1,"%f %f\n",xt->ve[j],al*l->ve[j]);

	    fclose(p1);

	    FILE *p2 = fopen("plot2.txt","w");

	    for (int j=0;j<x->n;j++)
	      fprintf(p2,"%f %f\n",xt->ve[j],pv->ve[j]);

	    fclose(p2);

	    char buffer[100];

	    sprintf(buffer,"./plo > plotl%.3f.pdf",k*(data->t_id[lfi] - S));

	    system(buffer); 

	  }

	  VEC *ld = v_get(x->n);

	  for (int j=0;j<xt->dim;j++)
	    ld->ve[j] = 2*(v->ve[j] - al*l->ve[j])*pv->ve[j];

	  /*
	    if (v->ve[j] < al*l->ve[j])
	      ld->ve[j] = -pv->ve[j];
	    else
	    ld->ve[j] = pv->ve[j];*/

	  //	    ld->ve[j] = pv->ve[j]*(v->ve[j] - al*l->ve[j]) / fabs(v->ve[j] - al*l->ve[j]);

	  ht->ve[i] = Q(xt,ld);

	  lfi += 1;

	  V_FREE(dt);
	  V_FREE(l);
	  V_FREE(ld);

	} 
      else 
	{ 
          //printf("%d\n",i);
          //printf("%f\n",Q(xt,pt));
          //printf("%f\n",Q(xt,v));
          //printf("%f\n",c(k*(i-tmi.I)));

	  ht->ve[i] = 2*(Q(xt,v)-1e3*c(data->cat,k,k*(i-S)))*Q(xt,pv);
	  /*    
	  if (Q(xt,v)<1e3*c(data->cat,k,k*(i-S)))
	    ht->ve[i] = -Q(xt,pv);//*(Q(xt,v)-1e3*c(k*(i - tmi.I)))/fabs(Q(xt,v)-1e3*c(k*(i - tmi.I)));
	  else
	  ht->ve[i]=Q(xt,pv);*/
	}

      tt->ve[i] = k*(i-S);

      V_FREE(xt);
      V_FREE(ut);
      V_FREE(pt);
      V_FREE(v);
      V_FREE(pv);

    }


  if (PLOTDERIV) {

    FILE *p1 = fopen("plot.txt","w");

    for (int i=S;i<x->m;i++)
      fprintf(p1,"%f %f\n",tt->ve[i],ht->ve[i]);

    fclose(p1);

    char buffer[100];

    sprintf(buffer,"./plo1 > plotht.pdf");

    system(buffer);

  }

  double blah = Q(tt,ht);

  V_FREE(ht);
  V_FREE(tt);

  return blah;

}


/*double G_ni_for_condition_number(

	    MAT *p,
	    MAT *x,
	    MAT *u,
	    struct DATA *data,
	    double iota

	    )
{

  int S = data->S;
  double k = data->k;

  int lfi=0;

  VEC *ht = v_get(x->m);
  VEC *tt = v_get(x->m);

  for (int i=S;i<x->m;i++)
    {

      VEC *xt = v_get(x->n);
      get_row(x,i,xt);      
      VEC *ut = v_get(x->n);
      get_row(u,i,ut);      
      VEC *pt = v_get(x->n);
      get_row(p,i,pt); 

      VEC *v = v_get(x->n);
      VEC *pv = v_get(x->n);

      for (int j=0;j<x->n;j++) 
	{
	  pv->ve[j] = iota*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*pt->ve[j];
	  v->ve[j] = iota*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];
	}

      if(data->t_id[lfi]==i) 
	{
      
	  VEC *l = v_get(xt->dim);
	  for (int j=0;j<xt->dim;j++)
	    l->ve[j] = ut->ve[j]; 

	  for (int j=0;j<xt->dim;j++)
	    if (v->ve[j] < al*l->ve[j])
	      ld->ve[j] = -pv->ve[j];
	    else
	      ld->ve[j] = pv->ve[j];

	  ht->ve[i] = Q(xt,ld);

          if (lfi<data->n)
	    lfi += 1;

	} 
      else 
	{ 
          //printf("%d\n",i);
          //printf("%f\n",Q(xt,pt));
          //printf("%f\n",Q(xt,v));
          //printf("%f\n",c(k*(i-tmi.I)));
        
	  if (Q(xt,v)<1e3*c(data->cat,k,k*(i-S)))
	    ht->ve[i] = -Q(xt,pv);//*(Q(xt,v)-1e3*c(k*(i - tmi.I)))/fabs(Q(xt,v)-1e3*c(k*(i - tmi.I)));
	  else
	    ht->ve[i]=Q(xt,pv);
	}

      tt->ve[i] = k*(i-S);

    }


  if (BLAH) {

    FILE *p1 = fopen("plot.txt","w");

    for (int i=S;i<x->m;i++)
      fprintf(p1,"%f %f\n",tt->ve[i],ht->ve[i]);

    fclose(p1);

    char buffer[100];

    sprintf(buffer,"./plo1 > plotht.pdf");

    system(buffer);

  }

  return Q(tt,ht);

}
*/

double G_ni(

	    MAT *p,
	    MAT *x,
	    MAT *u,
	    struct DATA *data,
	    double iota

	    )
{

  iota *=1e-3;
  int S = data->S;
  double k = data->k;

  int lfi=0;

  VEC *ht = v_get(x->m);
  VEC *tt = v_get(x->m);

  for (int i=S;i<x->m;i++)
    {

      VEC *xt = v_get(x->n);
      get_row(x,i,xt);      
      VEC *ut = v_get(x->n);
      get_row(u,i,ut);      
      VEC *pt = v_get(x->n);
      get_row(p,i,pt); 

      VEC *v = v_get(x->n);
      VEC *pv = v_get(x->n);

      for (int j=0;j<x->n;j++) 
	{
	  pv->ve[j] = iota*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*pt->ve[j];
	  v->ve[j] = iota*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];
	}

      if(lfi < data->n && data->t_id[lfi]==i) 
	{
      
	  VEC *dt = v_get(data->t_sz[lfi]);

	  for (int j=0;j<dt->dim;j++)
	    dt->ve[j] = data->lf[lfi][j];

	  double bw = get_bw(dt);

	  VEC *l = v_get(xt->dim);

	  for (int j=0;j<xt->dim;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  double al = 1e3*c(data->cat,k,k*(i - S)) / Q(xt,l);

	  if (PLOTDERIV) 
	    {

	      FILE *p1 = fopen("plot1.txt","w");

	      for (int j=0;j<x->n;j++)
		fprintf(p1,"%f %f\n",xt->ve[j],al*l->ve[j]);

	      fclose(p1);

	      FILE *p2 = fopen("plot2.txt","w");

	      for (int j=0;j<x->n;j++)
		fprintf(p2,"%f %f\n",xt->ve[j],pv->ve[j]);

	      fclose(p2);

	      char buffer[100];

	      sprintf(buffer,"./plo > plotl%.3f.pdf",k*(data->t_id[lfi] - S));

	      system(buffer);

	    }

	  VEC *ld = v_get(x->n);

	  for (int j=0;j<xt->dim;j++)
	    ld->ve[j] = 2*(v->ve[j] - al*l->ve[j])*pv->ve[j];

	  /*
	  for (int j=0;j<xt->dim;j++)
	    if (v->ve[j] < al*l->ve[j])
	      ld->ve[j] = -pv->ve[j];
	    else
	    ld->ve[j] = pv->ve[j];*/

	  // Ld->ve[j] = pv->ve[j]*(v->ve[j] - al*l->ve[j]) / fabs(v->ve[j] - al*l->ve[j]);

	  ht->ve[i] = Q(xt,ld);

	  lfi += 1;

	  V_FREE(dt);
	  V_FREE(l);
	  V_FREE(ld);

	} 
      else 
	{ 
          //printf("%d\n",i);
          //printf("%f\n",Q(xt,pt));
          //printf("%f\n",Q(xt,v));
          //printf("%f\n",c(k*(i-tmi.I)));

	  ht->ve[i] = 2*(Q(xt,v)-1e3*c(data->cat,k,k*(i-S)))*Q(xt,pv);
	  /*        
	  if (Q(xt,v)<1e3*c(data->cat,k,k*(i-S)))
	    ht->ve[i] = -Q(xt,pv);//*(Q(xt,v)-1e3*c(k*(i - tmi.I)))/fabs(Q(xt,v)-1e3*c(k*(i - tmi.I)));
	  else
	  ht->ve[i]=Q(xt,pv);*/
	}

      tt->ve[i] = k*(i-S);

      V_FREE(xt);
      V_FREE(ut);
      V_FREE(pt);
      V_FREE(v);
      V_FREE(pv);

    }

  if (PLOTDERIV) {

    FILE *p1 = fopen("plot.txt","w");

    for (int i=S;i<x->m;i++)
      fprintf(p1,"%f %f\n",tt->ve[i],ht->ve[i]);

    fclose(p1);

    char buffer[100];

    sprintf(buffer,"./plo1 > plotht.pdf");

    system(buffer);

  }


  double blah = Q(tt,ht);

  V_FREE(ht);
  V_FREE(tt);

  return blah;

}

void solve(

		 VEC *theta,
		 MAT *x,
		 MAT *u,
		 MAT *xhh,
		 MAT *xh,
		 MAT *xn,
		 MAT *uh,
		 MAT *un,
		 VEC *Ui,
		 VEC *Uh,
		 VEC *Uhh,
		 IVEC *idxi,
		 VEC *eff,
		 double k,		 
		 int S

		 )
{

  VEC *xnt; VEC *unt; 
  VEC *xht; VEC *uht;
  VEC *xhht; VEC *uhh;
  VEC *xt; VEC *ut;

  xnt = v_get(x->n+1);  unt = v_get(x->n+1);
  xht = v_get(x->n+1);  uht = v_get(x->n+1);
  xhht = v_get(x->n+1); uhh = v_get(x->n+1);
  xt = v_get(x->n);     ut = v_get(x->n);

  for (int j=1;j<x->n;j++) 
    x->me[0][j] = h*j;

  VEC * thextra = v_get(5);

  thextra->ve[0] = theta->ve[0];
  thextra->ve[1] = theta->ve[1];
  thextra->ve[2] = theta->ve[2]*1e-7;
  thextra->ve[3] = kappa;
  thextra->ve[4] = omega;
 
  set_row(u,0,initial(thextra,get_row(x,0,xt),ut));

  V_FREE(thextra);
 
  Ui->ve[0] = Q(xt,ut);

  double aa = theta->ve[0];
  double bb = theta->ve[1];
  double gg = theta->ve[2]*1e-7;
  double kk = kappa;
  double ww = omega;
  double ii = theta->ve[3]*1e-3;

  //  printf("\n");

  for (int i=1;i<x->m;i++)
    {

      double t = k*(i-S-1);
      double th = k*(i-S-.5);
      double thh = k*(i-S-.75);

      get_row(x,i-1,xt);
      get_row(u,i-1,ut);

      for (int j=1;j<=x->n;j++)
	{
	  xhht->ve[j] = xt->ve[j-1] + (k/4)*g(kk,ww,xt->ve[j-1]);
	  uhh->ve[j] = ut->ve[j-1]*exp(-(k/4)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k));
	}

      Q2(aa,kk,ww,xhht,uhh);
      Uhh->ve[i-1] = Q(xhht,uhh);
      set_row(xhh,i-1,xhht);
      
      for (int j=1;j<=x->n;j++)
	{
	  xht->ve[j] = xt->ve[j-1] + (k/2)*g(kk,ww,xhht->ve[j]);
	  uht->ve[j] = ut->ve[j-1]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j],Uhh->ve[i-1],k));
	}

      Q2(aa,kk,ww,xht,uht);
      Uh->ve[i-1] = Q(xht,uht);
      set_row(xh,i-1,xht);
      set_row(uh,i-1,uht);

      for (int j=1;j<=x->n;j++)
	{
	  xnt->ve[j] = xt->ve[j-1] + k*g(kk,ww,xht->ve[j]);
	  unt->ve[j] = ut->ve[j-1]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k));
	}

      Q2(aa,kk,ww,xnt,unt);
      Ui->ve[i] = Q(xnt,unt);
      set_row(xn,i-1,xnt);
      set_row(un,i-1,unt);
            
      int idx = idxselect(ww,xnt);
      idxi->ive[i-1] = idx;

      set_row(x,i,idxremove(xnt,xt,idx));
      set_row(u,i,idxremove(unt,ut,idx));

      //      printf("%f\n",Ui->ve[i]);

    }

  //exit(1);

  V_FREE(xt);
  V_FREE(xnt); 
  V_FREE(unt); 
  V_FREE(xht); 
  V_FREE(uht);
  V_FREE(xhht); 
  V_FREE(uhh);
  V_FREE(ut);

}

void solve_p_alpha(

		    VEC *theta,
		    MAT *p,
		    MAT *x,
		    MAT *u,
		    MAT *xhh,
		    MAT *xh,
		    MAT *xn,
		    MAT *uh,
		    MAT *un,
		    VEC *Ui,
		    VEC *Uh,
		    VEC *Uhh,
		    IVEC *idxi,
		    VEC *eff,
		    double k,		    
		    int S

		 )
{

  VEC *xt; VEC *xht; VEC *xnt; 
  VEC *ut; VEC *uht; VEC *pt; VEC *unt;
  xt = v_get(x->n);
  xht = v_get(x->n+1);
  xnt = v_get(x->n+1);
  ut = v_get(x->n);
  uht = v_get(x->n+1);
  unt = v_get(x->n+1);
  pt = v_get(x->n);

  VEC *xhht; 
  xhht = v_get(x->n+1);

  VEC *ph; VEC *pn;
  ph = v_get(x->n+1);
  pn = v_get(x->n+1);
 
  double aa = theta->ve[0];
  double bb = theta->ve[1];
  double gg = theta->ve[2]*1e-7;
  double kk = kappa;
  double ww = omega;
  double ii = theta->ve[3]*1e-3;

  VEC *Pi;
  Pi = v_get(x->m);

  VEC * thextra = v_get(5);

  thextra->ve[0] = theta->ve[0];
  thextra->ve[1] = theta->ve[1];
  thextra->ve[2] = theta->ve[2]*1e-7;
  thextra->ve[3] = kappa;
  thextra->ve[4] = omega;

  get_row(x,0,xt);
  ini_alpha(thextra,xt,pt);
  set_row(p,0,pt);

  V_FREE(thextra);
  
  Pi->ve[0] = Q(get_row(x,0,xt),get_row(p,0,pt));

  if (PLOTSOLVP) 
    {

      FILE *p2 = fopen("plot.txt","w");

      for (int j=0;j<x->n;j++)
	fprintf(p2,"%f %f\n",xt->ve[j],pt->ve[j]);

      fclose(p2);

      system("./plo1-nr > plotp_a_ini.pdf");

    }

  for (int i=1;i<x->m;i++)
    { 

      double t = k*(i-S-1);
      double th = k*(i-S-.5);
      double thh = k*(i-S-.75);

      get_row(x,i-1,xt);
      get_row(p,i-1,pt);
      get_row(xh,i-1,xht);
      get_row(xhh,i-1,xhht);
      get_row(uh,i-1,uht);
      get_row(xn,i-1,xnt);

      for (int j=1;j<=x->n;j++)
	ph->ve[j] = pt->ve[j-1]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k)) - exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k))*(k/2)*gg*Pi->ve[i-1]*ut->ve[j-1];
	 
      Q2_alpha(aa,kk,ww,xht,uht,ph);
      double Ph = Q(xht,ph);

      for (int j=1;j<=x->n;j++)
	{

	  double b = k*gg*Ph*uht->ve[j];
	  pn->ve[j] = pt->ve[j-1]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k)) - b*exp((k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j],Uhh->ve[i-1],k)-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k));

	}

      get_row(un,i-1,unt);
      Q2_alpha(aa,kk,ww,xnt,unt,pn);
      Pi->ve[i] = Q(xnt,pn);

      //for (int j=0;j<=x->n;j++)
      //fprintf(ft,"%f ",pn->ve[j]);
      //fprintf(ft,"\n");

      idxremove(pn,pt,idxi->ive[i-1]); 
      set_row(p,i,pt);

    }
  //  fclose(ft);
  //exit(1);   
  /*  printf("\n");
  for (int j=0;j<300;j++)
    {
      printf("%f %f\n",x->me[x->m-1][j],p->me[x->m-1][j]); //s(x->me[tmi.I+2][j])*w(x->me[tmi.I+2][j])*tmi.dp.iota*e(2*k)*); //-tmi.I+4][i]);
      //printf("%f %f\n",x->me[x->m-1][i],s(x->me[x->m-1][i])*w(x->me[x->m-1][i])*e(k*(x->m-9-(tmi.I+1)))*u->me[x->m-1][i] + s(x->me[x->m-1][i])*w(x->me[x->m-1][i])*tmi.dp.iota*e(k*(x->m-9-(tmi.I+1)))*p->me[x->m-1][i]);
    }
  printf("e\n");
  exit(1);
  */

  /*
  VEC *xxx = v_get(x->m);
  for (int j=0;j<x->m;j++)
    xxx->ve[j] = x->me[x->m-1][j];

  VEC *yyy = v_get(x->m);
  for (int j=0;j<x->m;j++)
    yyy->ve[j] = s(x->me[x->m-59][j])*w(x->me[x->m-59][j])*e(k*(x->m-59-(tmi.I+1)))*u->me[x->m-59][j] + s(x->me[x->m-59][j])*w(x->me[x->m-59][j])*iota*e(k*(x->m-59-(tmi.I+1)))*p->me[x->m-59][j];

  printf("%f\n",Q(xxx,yyy));

  exit(1);  
  */

  V_FREE(xt); 
  V_FREE(xht);
  V_FREE(xnt); 
  V_FREE(ut); 
  V_FREE(uht);
  V_FREE(pt);
  V_FREE(unt);
  V_FREE(ph);
  V_FREE(pn);
  V_FREE(Pi);
  V_FREE(xhht);


}

void solve_p_beta(

		 VEC *theta,
		 MAT *p,
		 MAT *x,
		 MAT *u,
		 MAT *xhh,
		 MAT *xh,
		 MAT *xn,
		 MAT *uh,
		 VEC *Ui,
		 VEC *Uh,
		 VEC *Uhh,
		 IVEC *idxi,
		 VEC *eff,
		 double k,
		 
		 int S

		 )
{

  VEC *xt; VEC *xht; VEC *xnt; 
  VEC *ut; VEC *uht; VEC *pt;
  xt = v_get(x->n);
  xht = v_get(x->n+1);
  xnt = v_get(x->n+1);
  ut = v_get(x->n);
  uht = v_get(x->n+1);
  pt = v_get(x->n);

  VEC *xhht; 
  xhht = v_get(x->n+1);

  VEC *ph; VEC *pn;
  ph = v_get(x->n+1);
  pn = v_get(x->n+1);

  double aa = theta->ve[0];
  double bb = theta->ve[1];
  double gg = theta->ve[2]*1e-7;
  double kk = kappa;
  double ww = omega;
  double ii = theta->ve[3]*1e-3;

  VEC *Pi;
  Pi = v_get(x->m);

  VEC * thextra = v_get(5);

  thextra->ve[0] = theta->ve[0];
  thextra->ve[1] = theta->ve[1];
  thextra->ve[2] = theta->ve[2]*1e-7;
  thextra->ve[3] = kappa;
  thextra->ve[4] = omega;

  get_row(x,0,xt);
  ini_beta(thextra,xt,pt);
  set_row(p,0,pt);

  V_FREE(thextra);
  
  Pi->ve[0] = Q(get_row(x,0,xt),get_row(p,0,pt));

  for (int i=1;i<x->m;i++)
    { 

      double t = k*(i-S-1);
      double th = k*(i-S-.5);
      double thh = k*(i-S-.75);

      get_row(x,i-1,xt);
      get_row(p,i-1,pt);
      get_row(xh,i-1,xht);
      get_row(xhh,i-1,xhht);
      get_row(uh,i-1,uht);
      get_row(xn,i-1,xnt);

      for (int j=1;j<=x->n;j++)
	ph->ve[j] = pt->ve[j-1]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k)) - exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k))*(k/2)*(1+gg*Pi->ve[i-1])*ut->ve[j-1];
	 
      Q2(aa,kk,ww,xht,ph);
      double Ph = Q(xht,ph);

      for (int j=1;j<=x->n;j++)
	{
	  double b = k*(1+gg*Ph)*uht->ve[j];
	  pn->ve[j] = pt->ve[j-1]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k)) - b*exp((k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j],Uhh->ve[i-1],k)-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k));
	}

      Q2(aa,kk,ww,xnt,pn);
      Pi->ve[i] = Q(xnt,pn);

      idxremove(pn,pt,idxi->ive[i-1]); 
      set_row(p,i,pt);

    }
   
  /*  printf("\n");
  for (int j=0;j<300;j++)
    {
      printf("%f %f\n",x->me[x->m-1][j],p->me[x->m-1][j]); //s(x->me[tmi.I+2][j])*w(x->me[tmi.I+2][j])*tmi.dp.iota*e(2*k)*); //-tmi.I+4][i]);
      //printf("%f %f\n",x->me[x->m-1][i],s(x->me[x->m-1][i])*w(x->me[x->m-1][i])*e(k*(x->m-9-(tmi.I+1)))*u->me[x->m-1][i] + s(x->me[x->m-1][i])*w(x->me[x->m-1][i])*tmi.dp.iota*e(k*(x->m-9-(tmi.I+1)))*p->me[x->m-1][i]);
    }
  printf("e\n");
  exit(1);
  */

  /*
  VEC *xxx = v_get(x->m);
  for (int j=0;j<x->m;j++)
    xxx->ve[j] = x->me[x->m-1][j];

  VEC *yyy = v_get(x->m);
  for (int j=0;j<x->m;j++)
    yyy->ve[j] = s(x->me[x->m-59][j])*w(x->me[x->m-59][j])*e(k*(x->m-59-(tmi.I+1)))*u->me[x->m-59][j] + s(x->me[x->m-59][j])*w(x->me[x->m-59][j])*iota*e(k*(x->m-59-(tmi.I+1)))*p->me[x->m-59][j];

  printf("%f\n",Q(xxx,yyy));

  exit(1);  
  */

  V_FREE(xt); 
  V_FREE(xht);
  V_FREE(xnt); 
  V_FREE(ut); 
  V_FREE(uht);
  V_FREE(pt);
  V_FREE(ph);
  V_FREE(pn);
  V_FREE(Pi);
  V_FREE(xhht);


}

void solve_p_gamma(

		 VEC *theta,
		 MAT *p,
		 MAT *x,
		 MAT *u,
		 MAT *xhh,
		 MAT *xh,
		 MAT *xn,
		 MAT *uh,
		 VEC *Ui,
		 VEC *Uh,
		 VEC *Uhh,
		 IVEC *idxi,
		 VEC *eff,
		 double k,
		 
		 int S

		 )
{

  VEC *xt; VEC *xht; VEC *xnt; 
  VEC *ut; VEC *uht; VEC *pt;
  xt = v_get(x->n);
  xht = v_get(x->n+1);
  xnt = v_get(x->n+1);
  ut = v_get(x->n);
  uht = v_get(x->n+1);
  pt = v_get(x->n);

  VEC *xhht; 
  xhht = v_get(x->n+1);

  VEC *ph; VEC *pn;
  ph = v_get(x->n+1);
  pn = v_get(x->n+1);

  double aa = theta->ve[0];
  double bb = theta->ve[1];
  double gg = theta->ve[2]*1e-7;
  double kk = kappa;
  double ww = omega;
  double ii = theta->ve[3]*1e-3;
 
  VEC *Pi;
  Pi = v_get(x->m);

  VEC * thextra = v_get(5);

  thextra->ve[0] = theta->ve[0];
  thextra->ve[1] = theta->ve[1];
  thextra->ve[2] = theta->ve[2];
  thextra->ve[3] = kappa;
  thextra->ve[4] = omega;

  get_row(x,0,xt);
  ini_gamma(thextra,xt,pt);
  set_row(p,0,pt);

  V_FREE(thextra);
  
  Pi->ve[0] = Q(get_row(x,0,xt),get_row(p,0,pt));

  for (int i=1;i<x->m;i++)
    { 

      double t = k*(i-S-1);
      double th = k*(i-S-.5);
      double thh = k*(i-S-.75);

      get_row(x,i-1,xt);
      get_row(p,i-1,pt);
      get_row(xh,i-1,xht);
      get_row(xhh,i-1,xhht);
      get_row(uh,i-1,uht);
      get_row(xn,i-1,xnt);

      for (int j=1;j<=x->n;j++)
	ph->ve[j] = pt->ve[j-1]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k)) - exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k))*(k/2)*(1e-7*Ui->ve[i-1]+gg*Pi->ve[i-1])*ut->ve[j-1];
	 
      Q2(aa,kk,ww,xht,ph);
      double Ph = Q(xht,ph);

      for (int j=1;j<=x->n;j++)
	{
	  double b = k*(1e-7*Uh->ve[i-1]+gg*Ph)*uht->ve[j];
	  pn->ve[j] = pt->ve[j-1]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k)) - b*exp((k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j],Uhh->ve[i-1],k)-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k));
	}

      Q2(aa,kk,ww,xnt,pn);
      Pi->ve[i] = Q(xnt,pn);

      idxremove(pn,pt,idxi->ive[i-1]); 
      set_row(p,i,pt);

    }
   
  /*  printf("\n");
  for (int j=0;j<300;j++)
    {
      printf("%f %f\n",x->me[x->m-1][j],p->me[x->m-1][j]); //s(x->me[tmi.I+2][j])*w(x->me[tmi.I+2][j])*tmi.dp.iota*e(2*k)*); //-tmi.I+4][i]);
      //printf("%f %f\n",x->me[x->m-1][i],s(x->me[x->m-1][i])*w(x->me[x->m-1][i])*e(k*(x->m-9-(tmi.I+1)))*u->me[x->m-1][i] + s(x->me[x->m-1][i])*w(x->me[x->m-1][i])*tmi.dp.iota*e(k*(x->m-9-(tmi.I+1)))*p->me[x->m-1][i]);
    }
  printf("e\n");
  exit(1);
  */

  /*
  VEC *xxx = v_get(x->m);
  for (int j=0;j<x->m;j++)
    xxx->ve[j] = x->me[x->m-1][j];

  VEC *yyy = v_get(x->m);
  for (int j=0;j<x->m;j++)
    yyy->ve[j] = s(x->me[x->m-59][j])*w(x->me[x->m-59][j])*e(k*(x->m-59-(tmi.I+1)))*u->me[x->m-59][j] + s(x->me[x->m-59][j])*w(x->me[x->m-59][j])*iota*e(k*(x->m-59-(tmi.I+1)))*p->me[x->m-59][j];

  printf("%f\n",Q(xxx,yyy));

  exit(1);  
  */

  V_FREE(xt); 
  V_FREE(xht);
  V_FREE(xnt); 
  V_FREE(ut); 
  V_FREE(uht);
  V_FREE(pt);
  V_FREE(ph);
  V_FREE(pn);
  V_FREE(Pi);
  V_FREE(xhht);


}

void solve_p_kappa(

		 VEC *theta,
		 MAT *p,
		 MAT *x,
		 MAT *u,
		 MAT *xhh,
		 MAT *xh,
		 MAT *xn,
		 MAT *uh,
		 MAT *un,
		 VEC *Ui,
		 VEC *Uh,
		 VEC *Uhh,
		 IVEC *idxi,
		 VEC *eff,
		 double k,
		 
		 int S

		 )
{

  // this function has not been updated for new birth function

  VEC *xt; VEC *xht; VEC *xnt; 
  VEC *ut; VEC *uht; VEC *pt; VEC *unt;
  xt = v_get(x->n);
  xht = v_get(x->n+1);
  xnt = v_get(x->n+1);
  ut = v_get(x->n);
  uht = v_get(x->n+1);
  unt = v_get(x->n+1);
  pt = v_get(x->n);

  VEC *xhht; 
  xhht = v_get(x->n+1);

  VEC *ph; VEC *pn;
  ph = v_get(x->n+1);
  pn = v_get(x->n+1);

  double a1 = theta->ve[0];
  double a2 = 0;
  double bb = theta->ve[1];
  double gg = theta->ve[2]*1e-7;
  double kk = theta->ve[3];
  double ww = theta->ve[4];
  double ii = theta->ve[5];

  VEC *Pi;
  Pi = v_get(x->m);

  get_row(x,0,xt);
  ini_kappa(theta,xt,pt);
  set_row(p,0,pt);
  
  Pi->ve[0] = Q(get_row(x,0,xt),get_row(p,0,pt));

  if (PLOTSOLVP)
    {

      FILE *p2 = fopen("plot.txt","w");

      for (int j=0;j<x->n;j++)
	fprintf(p2,"%f %f\n",xt->ve[j],pt->ve[j]);

      fclose(p2);

      system("./plo1-nr > plotp_k_ini.pdf");

    }

  for (int i=1;i<x->m;i++)
    { 

      double t = k*(i-S-1);
      double th = k*(i-S-.5);
      double thh = k*(i-S-.75);

      get_row(x,i-1,xt);
      get_row(p,i-1,pt);
      get_row(xh,i-1,xht);
      get_row(xhh,i-1,xhht);
      get_row(uh,i-1,uht);
      get_row(xn,i-1,xnt);


      int j=1;
      ph->ve[j] = pt->ve[j-1]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k)) - exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k))*(k/2)*( (gg*Pi->ve[i-1]-1)*ut->ve[j-1] + (ww - xt->ve[j-1])*(ut->ve[j]-ut->ve[j-1])/(xt->ve[j]-xt->ve[j-1]) );

      for (int j=2;j<=x->n;j++)
	ph->ve[j] = pt->ve[j-1]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k)) - exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k))*(k/2)*( (gg*Pi->ve[i-1]-1)*ut->ve[j-1] + (ww - xt->ve[j-1])*.5*( (ut->ve[j]-ut->ve[j-1])/(xt->ve[j]-xt->ve[j-1]) + (ut->ve[j-1]-ut->ve[j-2])/(xt->ve[j-1]-xt->ve[j-2]) ) );
	
      Q2_kappa(a1,a2,kk,ww,xht,uht,ph);
      double Ph = Q(xht,ph);

      for (int j=1;j<x->n;j++)
	{
	  double b = k*( (gg*Ph-1)*uht->ve[j] + (ww - xht->ve[j])*.5*( (uht->ve[j]-uht->ve[j-1])/(xht->ve[j]-xht->ve[j-1]) + (uht->ve[j+1]-uht->ve[j])/(xht->ve[j+1]-xht->ve[j]) ) );
	  pn->ve[j] = pt->ve[j-1]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k)) - b*exp((k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j],Uhh->ve[i-1],k)-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k));
	}

      j= x->n;
      double b = k*(gg*Ph-1)*uht->ve[j];
      pn->ve[j] = pt->ve[j-1]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k)) - b*exp((k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j],Uhh->ve[i-1],k)-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k));

      get_row(un,i-1,unt);
      Q2_kappa(a1,a2,kk,ww,xnt,unt,pn);
      Pi->ve[i] = Q(xnt,pn);

      idxremove(pn,pt,idxi->ive[i-1]); 
      set_row(p,i,pt);

    }

  if (PLOTSOLVP)
    {

      FILE *p2 = fopen("plot.txt","w");

      int i=100;
      for (int j=0;j<x->n;j++)
	fprintf(p2,"%f %f\n",x->me[i][j],p->me[i][j]);

      fclose(p2);

      system("./plo1_nr > plotp_k_100.pdf");

      p2 = fopen("plot.txt","w");

      i=1000;
      for (int j=0;j<x->n;j++)
	fprintf(p2,"%f %f\n",x->me[i][j],p->me[i][j]);

      fclose(p2);

      system("./plo1_nr > plotp_k_1000.pdf");

    }

  V_FREE(xt);
  V_FREE(xht);
  V_FREE(xnt); 
  V_FREE(ut); 
  V_FREE(uht);
  V_FREE(pt);
  V_FREE(unt);
  V_FREE(ph);
  V_FREE(pn);
  V_FREE(Pi);
  V_FREE(xhht);

}

void solve_p_omega(

		 VEC *theta,
		 MAT *p,
		 MAT *x,
		 MAT *u,
		 MAT *xhh,
		 MAT *xh,
		 MAT *xn,
		 MAT *uh,
		 MAT *un,
		 VEC *Ui,
		 VEC *Uh,
		 VEC *Uhh,
		 IVEC *idxi,
		 VEC *eff,
		 double k,
		 int S

		 )
{

  // this function has not been updated for new birth function

  VEC *xt; VEC *xht; VEC *xnt; 
  VEC *ut; VEC *uht; VEC *pt; VEC *unt;
  xt = v_get(x->n);
  xht = v_get(x->n+1);
  xnt = v_get(x->n+1);
  ut = v_get(x->n);
  uht = v_get(x->n+1);
  unt = v_get(x->n+1);
  pt = v_get(x->n);

  VEC *xhht; 
  xhht = v_get(x->n+1);

  VEC *ph; VEC *pn;
  ph = v_get(x->n+1);
  pn = v_get(x->n+1);
 
  double a1 = theta->ve[0];
  double a2 = 0;
  double bb = theta->ve[1];
  double gg = theta->ve[2]*1e-7;
  double kk = theta->ve[3];
  double ww = theta->ve[4];
  double ii = theta->ve[5];

  VEC *Pi;
  Pi = v_get(x->m);

  get_row(x,0,xt);
  ini_omega(theta,xt,pt);
  set_row(p,0,pt);
  
  Pi->ve[0] = Q(get_row(x,0,xt),get_row(p,0,pt));

  if (PLOTSOLVP) 
    {

      FILE *p2 = fopen("plot.txt","w");

      for (int j=0;j<x->n;j++) 
	fprintf(p2,"%f %f\n",xt->ve[j],pt->ve[j]);

      fclose(p2);

      system("./plo1-nr > plotp_w_ini.pdf");

    }

  for (int i=1;i<x->m;i++)
    { 

      double t = k*(i-S-1);
      double th = k*(i-S-.5);
      double thh = k*(i-S-.75);

      get_row(x,i-1,xt);
      get_row(p,i-1,pt);
      get_row(xh,i-1,xht);
      get_row(xhh,i-1,xhht);
      get_row(uh,i-1,uht);
      get_row(xn,i-1,xnt);


      int j=1;
      ph->ve[j] = pt->ve[j-1]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k)) - exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k))*(k/2)*( gg*Pi->ve[i-1]*ut->ve[j-1] + kk*(ut->ve[j]-ut->ve[j-1])/(xt->ve[j]-xt->ve[j-1]) );

      for (int j=2;j<=x->n;j++)
	ph->ve[j] = pt->ve[j-1]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k)) - exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k))*(k/2)*( gg*Pi->ve[i-1]*ut->ve[j-1] + kk*.5*( (ut->ve[j]-ut->ve[j-1])/(xt->ve[j]-xt->ve[j-1]) + (ut->ve[j-1]-ut->ve[j-2])/(xt->ve[j-1]-xt->ve[j-2]) ) );
	
      Q2_omega(a1,a2,kk,ww,xht,uht,ph);
      double Ph = Q(xht,ph);

      for (int j=1;j<x->n;j++)
	{
	  double b = k*( gg*Ph*uht->ve[j] + kk*.5*( (uht->ve[j]-uht->ve[j-1])/(xht->ve[j]-xht->ve[j-1]) + (uht->ve[j+1]-uht->ve[j])/(xht->ve[j+1]-xht->ve[j]) ) );
	  pn->ve[j] = pt->ve[j-1]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k)) - b*exp((k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j],Uhh->ve[i-1],k)-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k));
	}

      j= x->n;
      double b = k*gg*Ph*uht->ve[j];
      pn->ve[j] = pt->ve[j-1]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k)) - b*exp((k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j],Uhh->ve[i-1],k)-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k));

      get_row(un,i-1,unt);
      Q2_omega(a1,a2,kk,ww,xnt,unt,pn);
      Pi->ve[i] = Q(xnt,pn);

      idxremove(pn,pt,idxi->ive[i-1]); 
      set_row(p,i,pt);

    }

  if (PLOTSOLVP) 
    {

      FILE *p2 = fopen("plot.txt","w");

      int i=100;
      for (int j=0;j<x->n;j++)
	fprintf(p2,"%f %f\n",x->me[i][j],p->me[i][j]);

      fclose(p2);

      system("./plo1_nr > plotp_w_100.pdf");

      p2 = fopen("plot.txt","w");

      i=1000;
      for (int j=0;j<x->n;j++)
	fprintf(p2,"%f %f\n",x->me[i][j],p->me[i][j]);

      fclose(p2);

      system("./plo1_nr > plotp_w_1000.pdf");

    }

  V_FREE(xt); 
  V_FREE(xht);
  V_FREE(xnt); 
  V_FREE(ut); 
  V_FREE(uht);
  V_FREE(pt);
  V_FREE(ph);
  V_FREE(pn);
  V_FREE(Pi);
  V_FREE(xhht);

}

void solve_p_iota(

		  VEC *theta,
		  MAT *p,
		  MAT *x,
		  MAT *u,
		  MAT *xhh,
		  MAT *xh,
		  MAT *xn,
		  MAT *uh,
		  VEC *Ui,
		  VEC *Uh,
		  VEC *Uhh,
		  IVEC *idxi,
		  VEC *eff,
		  double k,		  
		  int S

		 )
{

  VEC *xt; VEC *xht; VEC *xnt; 
  VEC *ut; VEC *uht; VEC *pt;
  xt = v_get(x->n);
  xht = v_get(x->n+1);
  xnt = v_get(x->n+1);
  ut = v_get(x->n);
  uht = v_get(x->n+1);
  pt = v_get(x->n);

  VEC *xhht; 
  xhht = v_get(x->n+1);

  VEC *ph; VEC *pn;
  ph = v_get(x->n+1);
  pn = v_get(x->n+1);

  double aa = theta->ve[0];
  double bb = theta->ve[1];
  double gg = theta->ve[2]*1e-7;
  double kk = kappa;
  double ww = omega;
  double ii = theta->ve[3]*1e-3;
 
  VEC *Pi;
  Pi = v_get(x->m);
  
  Pi->ve[0] = Q(get_row(x,0,xt),get_row(p,0,pt));

  for (int i=1;i<x->m;i++)
    { 

      double t = k*(i-S-1);
      double th = k*(i-S-.5);
      double thh = k*(i-S-.75);

      get_row(x,i-1,xt);
      get_row(p,i-1,pt);
      get_row(xh,i-1,xht);
      get_row(xhh,i-1,xhht);
      get_row(uh,i-1,uht);
      get_row(xn,i-1,xnt);

      for (int j=1;j<=x->n;j++)
	ph->ve[j] = pt->ve[j-1]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k)) - exp(-(k/2)*zstar(eff,bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1],k))*(k/2)*(1e-3*s(xt->ve[j-1])*e(eff,k,t)+gg*Pi->ve[i-1])*ut->ve[j-1];
	 
      Q2(aa,kk,ww,xht,ph);
      double Ph = Q(xht,ph);

      for (int j=1;j<=x->n;j++)
	{
	  double b = k*(1e-3*s(xht->ve[j])*e(eff,k,th)+gg*Ph)*uht->ve[j];
	  pn->ve[j] = pt->ve[j-1]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k)) - b*exp((k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j],Uhh->ve[i-1],k)-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1],k));
	}

      Q2(aa,kk,ww,xnt,pn);
      Pi->ve[i] = Q(xnt,pn);

      idxremove(pn,pt,idxi->ive[i-1]); 
      set_row(p,i,pt);

    }
   
  /*  printf("\n");
  for (int j=0;j<300;j++)
    {
      printf("%f %f\n",x->me[x->m-1][j],p->me[x->m-1][j]); //s(x->me[tmi.I+2][j])*w(x->me[tmi.I+2][j])*tmi.dp.iota*e(2*k)*); //-tmi.I+4][i]);
      //printf("%f %f\n",x->me[x->m-1][i],s(x->me[x->m-1][i])*w(x->me[x->m-1][i])*e(k*(x->m-9-(tmi.I+1)))*u->me[x->m-1][i] + s(x->me[x->m-1][i])*w(x->me[x->m-1][i])*tmi.dp.iota*e(k*(x->m-9-(tmi.I+1)))*p->me[x->m-1][i]);
    }
  printf("e\n");
  exit(1);
  */

  /*
  VEC *xxx = v_get(x->m);
  for (int j=0;j<x->m;j++)
    xxx->ve[j] = x->me[x->m-1][j];

  VEC *yyy = v_get(x->m);
  for (int j=0;j<x->m;j++)
    yyy->ve[j] = s(x->me[x->m-59][j])*w(x->me[x->m-59][j])*e(k*(x->m-59-(tmi.I+1)))*u->me[x->m-59][j] + s(x->me[x->m-59][j])*w(x->me[x->m-59][j])*iota*e(k*(x->m-59-(tmi.I+1)))*p->me[x->m-59][j];

  printf("%f\n",Q(xxx,yyy));

  exit(1);  
  */

  V_FREE(xt); 
  V_FREE(xht);
  V_FREE(xnt); 
  V_FREE(ut); 
  V_FREE(uht);
  V_FREE(pt);
  V_FREE(ph);
  V_FREE(pn);
  V_FREE(Pi);
  V_FREE(xhht);

}

double c(

	  VEC * ca,
	  double r,
	  double t

	  )
{

  if (t<0)
    return 0;
  else
    {
      int idx = floor((t + (r/2) - 1e-12)/r); // better to use double epsilon?
      return ca->ve[idx];
    }
}

double Q(

	 VEC * x,
	 VEC * u

	 )
{ /* Q for quadrature. Integrates u over x. */

  if (x->dim != u->dim)
    error(E_SIZES,"Q");

  double rt = 0.;

  for (int i=0;i<x->dim-1;i++) 
    rt = rt + .5 * (x->ve[i+1] - x->ve[i]) * (u->ve[i] + u->ve[i+1]);
   
  
  return rt;

}

double Qn(

	  VEC * x,
	  VEC * u,
	  int n

	  )
{ /* Q for quadrature. Integrates u over x. */

  //if (x->dim != u->dim)
  //  error(E_SIZES,"Q");

  double rt = 0.;

  for (int i=0;i<n-1;i++) 
    rt = rt + .5 * (x->ve[i+1] - x->ve[i]) * (u->ve[i] + u->ve[i+1]);
  
  return rt;

}

/* Quadrature plus implicit.. */
void Q2(

	double a,
	double k,
	double w,
	VEC *x,
	VEC *u

	)
{

  if (x->dim != u->dim)
    {
      printf("%d %d\n",x->dim,u->dim);
      error(E_SIZES,"Q2");
    }

  double rt = x->ve[1] * b(a,x->ve[1])*u->ve[1] - x->ve[0]*b(a,x->ve[1])*u->ve[1]; 

  for (int j=1;j<x->dim-1;j++) 
    rt = rt + (b(a,x->ve[j])*u->ve[j] + b(a,x->ve[j+1])*u->ve[j+1]) * (x->ve[j+1]-x->ve[j]);

  u->ve[0] = rt / (2*k*w + x->ve[0]*b(a,x->ve[0]) - x->ve[1]*b(a,x->ve[0])); 

}

/* Quadrature plus implicit.. */
void Q2_alpha(

	       double a,
	       double k,
	       double w,
	       VEC *x,
	       VEC *u,
	       VEC *p

	       )
{

  if (x->dim != u->dim)
    {
      printf("%d %d\n",x->dim,u->dim);
      error(E_SIZES,"Q2");
    }

  double x0 = x->ve[0];
  double x1 = x->ve[1];

  double u0 = u->ve[0];
  double u1 = u->ve[1];

  double p0 = p->ve[0];
  double p1 = p->ve[1];

  double rt = A1*x0*u0*x1 + A2*x0*x0*u0*x1 + a*A1*x1*p1*x1 + a*A2*x1*x1*p1*x1 + A1*x1*u1*x1 + A2*x1*x1*u1*x1 - A1*x0*u0*x0 - A2*x0*x0*u0*x0 - a*A1*x1*p1*x0 - a*A2*x1*x1*p1*x0 - A1*x1*u1*x0 - A2*x1*x1*u1*x0;

  for (int j=1;j<x->dim-1;j++) 
    {

      x0 = x->ve[j];
      x1 = x->ve[j+1];

      u0 = u->ve[j];
      u1 = u->ve[j+1];

      p0 = p->ve[j];
      p1 = p->ve[j+1];

      rt += ( a*A1*x0*p0 + a*A2*x0*x0*p0 + A1*x0*u0 + A2*x0*x0*u0 + a*A1*x1*p1 + a*A2*x1*x1*p1 + A1*x1*u1 + A2*x1*x1*u1 ) * (x1-x0);

    }

  x0 = x->ve[0];
  x1 = x->ve[1];

  p->ve[0] =  rt / (2*k*w + a*A1*x0*x0 + a*A2*x0*x0*x0 - a*A1*x0*x1 - a*A2*x0*x0*x1); 

}

/* Quadrature plus implicit.. */
void Q2_kappa(

	      double a,
	      double a2,
	      double k,
	      double w,
	      VEC *x,
	      VEC *u,
	      VEC *p

	      )
{

  // this function has not been updated for new birth function

  if (x->dim != p->dim)
    {
      printf("%d %d\n",x->dim,p->dim);
      error(E_SIZES,"Q2");
    }

  double rt = x->ve[1] * b(a,x->ve[1])*p->ve[1] - x->ve[0]*b(a,x->ve[1])*p->ve[1]; 

  for (int j=1;j<x->dim-1;j++) 
    rt = rt + (b(a,x->ve[j])*p->ve[j] + b(a,x->ve[j+1])*p->ve[j+1]) * (x->ve[j+1]-x->ve[j]);

  p->ve[0] = (rt - 2*w*u->ve[0]) / (2*k*w + x->ve[0]*b(a,x->ve[0]) - x->ve[1]*b(a,x->ve[0])); 

}


void Q2_omega(

	      double a,
	      double a2,
	      double k,
	      double w,
	      VEC *x,
	      VEC *u,
	      VEC *p

	      )
{

  // this function has not been updated for new birth function

  if (x->dim != p->dim)
    {
      printf("%d %d\n",x->dim,p->dim);
      error(E_SIZES,"Q2");
    }

  double rt = x->ve[1] * b(a,x->ve[1])*p->ve[1] - x->ve[0]*b(a,x->ve[1])*p->ve[1]; 

  for (int j=1;j<x->dim-1;j++) 
    rt = rt + (b(a,x->ve[j])*p->ve[j] + b(a,x->ve[j+1])*p->ve[j+1]) * (x->ve[j+1]-x->ve[j]);

  p->ve[0] = (rt - 2*k*u->ve[0]) / (2*k*w + x->ve[0]*b(a,x->ve[0]) - x->ve[1]*b(a,x->ve[0])); 

}



double zstar(

	      VEC *eff,
	      double b,
	      double g,
	      double k,
	      double i,
	      double t,
	      double x,
	      double U,
	      double r

	     )
{ /* death function: beta + gamma U + s(x)f(t) - kappa */

  return b + g*U + s(x)*i*e(eff,r,t) - k;

}

double e(

	  VEC * ef,
	  double r,
	  double t

	  )
{

  if (t<0)
    {
      double cek = r/4;
      double ept5 = ef->ve[(int)floor((.5 + (cek/2) - 1e-12)/cek)];
      double m = ept5/24;

      return m*t+ ept5;
    }
  else
    {
      double cek = r/4;
      int idx = floor((t + (cek/2) - 1e-12)/cek); // better to use double epsilon?
      return ef->ve[idx];
    }
}

double get_bw(

	      VEC *dt

	      )
{ /* Bandwidth. using Silverman's rule of thumb. */

  double sm = v_sum(dt);
  double mn = sm/(float)dt->dim;
  double sd = 0;
  
  for (int i=0;i<dt->dim;i++)
    sd = sd + pow(dt->ve[i] - mn,2.);
  sd = sd / (float)dt->dim;
  sd = sqrt(sd);
  double hi = sd;
  
  PERM *order = px_get(dt->dim);
  v_sort(dt,order);
  
  int fstQi = floor(dt->dim/4);
  double fstQ = dt->ve[fstQi];
  int thrQi = floor(3*dt->dim/4);
  double thrQ = dt->ve[thrQi];
  double IQR = thrQ - fstQ;
  
  IQR = IQR/1.34;

  double lo = (hi < IQR) ? hi : IQR;
  double bw = 0.9 * lo * pow((double)dt->dim,-0.2);

  PX_FREE(order);

  return bw;

}

VEC *numgrad(

	     double (*model)(VEC *,void *),
	     void *stuff,
	     VEC *par,
	     double epsilon

	     )
{

  double f0 = (*model)(par,stuff);

  VEC *ei = v_get(par->dim);
  VEC *d = v_get(par->dim);
  VEC *fg = v_get(par->dim);

  for (int i=0;i<par->dim;i++)
    {
      ei->ve[i] = 1;
      sv_mlt(epsilon,ei,d);
      double f1 = (*model)(v_add(par,d,VNULL),stuff);
      fg->ve[i] = (f1 - f0) / d->ve[i];
      ei->ve[i] = 0;
    }

  return fg;

}

VEC * calc_alpha(

			double a,
			double k,
			double w,
			double bt,
			double f

			)
{

  VEC *retur = v_get(3);

  VEC *x = v_get(500);
  VEC *v = v_get(x->dim);
  VEC *t = v_get(x->dim);

  for (int j=0;j<x->dim;j++)
    x->ve[j] = j*w/500.;

  for (int j=0;j<t->dim;j++)
    t->ve[j] = (bt + f*s(x->ve[j])) / (k*(w-x->ve[j]));

  for (int j=0;j<v->dim;j++)
    v->ve[j] = a*(A1*x->ve[j] + A2*pow(x->ve[j],2.))*exp(-Qn(x,t,j+1))/(k*(w-x->ve[j]));

  double qv = Q(x,v);

  retur->ve[0] = pow(qv-1,2.);

  VEC *z = v_get(x->dim);
  VEC *q = v_get(x->dim);

  for (int j=0;j<v->dim;j++)
    z->ve[j] = pow(x->ve[j],2.)*exp(-Qn(x,t,j+1))/(k*(w-x->ve[j]));

  retur->ve[1] = 2*Q(x,z)*(qv-1); 

  return retur;

}

VEC * VMGMM_eq(

		  VEC *x,
		  struct DATA *d,
		  VEC *grad,
		  double *f

		  )
{

  VEC *xx = v_get(d->J);
  for (int j=0;j<xx->dim;j++)
    xx->ve[j] = h*j;

  int nd=0;
  for (int i=0;i<d->n;i++)
    nd += d->t_sz[i];

  VEC *dt = v_get(nd);
  int ii=0;
  for (int i=0;i<d->n;i++)
    for (int jj=0;jj<d->t_sz[i];jj++)
      {
	dt->ve[ii] = d->lf[i][jj];
	ii++;
      }

  double bw = get_bw(dt);

  VEC *l = v_get(d->J);

  for (int j=0;j<l->dim;j++)
    for (int jj=0;jj<dt->dim;jj++)
      l->ve[j] += exp( -pow((xx->ve[j] - dt->ve[jj])/bw,2.) );

  double Ql = Q(xx,l);

  for (int j=0;j<l->dim;j++)
    l->ve[j] /= Ql;

  VEC *v = v_get(xx->dim);
  VEC *w = v_get(xx->dim);

  for (int j=0;j<w->dim;j++)
    w->ve[j] = (x->ve[0] + x->ve[1]*s(xx->ve[j])) / (kappa*(omega-xx->ve[j]));

  for (int j=1;j<v->dim;j++)
    v->ve[j] = s(xx->ve[j])*exp(-Qn(xx,w,j+1))/(omega-xx->ve[j]);
      
  double B = Q(xx,v);

  VEC * vn = v_get(v->dim);

  for (int j=0;j<v->dim;j++)
    vn->ve[j] = v->ve[j]/B;
  /*
  printf("\n");
  for (int j=0;j<w->dim;j++)
    printf("%f %f\n",xx->ve[j],vn->ve[j]);
  printf("e\n\n");
  for (int j=0;j<w->dim;j++)
    printf("%f %f\n",xx->ve[j],l->ve[j]);
  exit(1);
  */
  VEC *objf = v_get(v->dim);

  for (int j=0;j<v->dim;j++)
    objf->ve[j]=l->ve[j]*log(l->ve[j]/(vn->ve[j]+1e-12) + 1e-12);
    
  *f = Q(xx,objf);

  VEC *w2 = v_get(w->dim);

  for (int j=0;j<w->dim-1;j++)
    w2->ve[j] = 1/(kappa*(omega-xx->ve[j]));

  VEC *w3 = v_get(w->dim);
  for (int j=1;j<w->dim;j++)
    w3->ve[j] = s(xx->ve[j])*exp(-Qn(xx,w,j+1))*Qn(xx,w2,j+1)/(omega-xx->ve[j]);

  double C = Q(xx,w3);

  VEC * fc = v_get(w->dim);

  for (int j=1;j<w->dim;j++)
    fc->ve[j] = (B*Qn(xx,w2,j+1)-C ) / B;
    
  VEC *w4 = v_get(w->dim);
  
  for (int j=0;j<w->dim;j++)
    w4->ve[j] =  s(xx->ve[j])/(kappa*(omega-xx->ve[j]));

  VEC *w5 = v_get(w->dim);

  for (int j=0;j<w->dim;j++)
    w5->ve[j] = s(xx->ve[j])*exp(-Qn(xx,w,j+1))*Qn(xx,w4,j+1)/(omega-xx->ve[j]);

  double D = Q(xx,w5);

  VEC * fc2 = v_get(w->dim);

  for (int j=1;j<w->dim;j++)
    fc2->ve[j] = (B*Qn(xx,w4,j+1)-D) / B;

  VEC *dGdb = v_get(w->dim);
  VEC *dGdf = v_get(w->dim);
    
  for (int j=0;j<w->dim;j++)
    {
      dGdb->ve[j]=l->ve[j]*fc->ve[j];
      dGdf->ve[j]=l->ve[j]*fc2->ve[j];
    }

  grad->ve[0] = Q(xx,dGdb);
  grad->ve[1] = Q(xx,dGdf);

  V_FREE(v);
  V_FREE(w);
  V_FREE(w2);
  V_FREE(w3);
  V_FREE(w4);
  V_FREE(w5);
  V_FREE(dGdb);
  V_FREE(dGdf);

  return grad;

}


/*double ConditionNumber(

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

  *f = H(x,u,dataptr,theta->ve[4]);
  
  MAT *p_a1 = m_get(x->m,x->n);
  solve_p_alpha1(theta,p_a1,x,u,xhh,xh,xn,uh,un,Ui,Uh,Uhh,idxi,dataptr->eff,dataptr->k,dataptr->S);
  grad->ve[0] = G_ni(p_a1,x,u,dataptr,theta->ve[4]);

  MAT *p_a2 = m_get(x->m,x->n);
  solve_p_alpha2(theta,p_a2,x,u,xhh,xh,xn,uh,un,Ui,Uh,Uhh,idxi,dataptr->eff,dataptr->k,dataptr->S);
  grad->ve[1] = G_ni(p_a2,x,u,dataptr,theta->ve[4]); 

  MAT *p_b = m_get(x->m,x->n);
  solve_p_beta  (theta,p_b,x,u,xhh,xh,xn,uh,Ui,Uh,Uhh,idxi,dataptr->eff,dataptr->k,dataptr->S);
  grad->ve[2]= G_ni(p_b,x,u,dataptr,theta->ve[4]);

  MAT *p_g = m_get(x->m,x->n);
  solve_p_gamma (theta,p_g,x,u,xhh,xh,xn,uh,Ui,Uh,Uhh,idxi,dataptr->eff,dataptr->k,dataptr->S);
  grad->ve[3] = G_ni(p_g,x,u,dataptr,theta->ve[4]);
*/
