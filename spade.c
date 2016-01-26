// SPADE
// Stock assessment using PArtial Differential Equations
// Alex Campbell 'ghostofsandy' 2015

#ifndef DEBUG
#  define DEBUG 0
#endif
#define _GNU_SOURCE
#define debug_print(...) do { if (DEBUG) fprintf(stderr, __VA_ARGS__); } while (0)

#include <fenv.h>
#include <sys/time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "meschach/matrix.h"
#include "meschach/matrix2.h"
#include "meschach/sparse.h"

struct GP { double kappa,omega; };
struct BP { double alpha1,alpha2; };
struct DP { double beta,gamma,iota; };

struct FMS { VEC *xx; VEC *obs; struct GP gp; };
struct SMS { VEC *length; VEC *age; };
struct TMS { VEC *xt; VEC *xr; VEC *tl; };
struct TMI { struct DP dp; struct GP gp; struct BP bp; int I,J; };

double iota1=5.2;
double iota2=0.619;
double phi=16.5;
double eta1=1.703205e-5;
double eta2=2.9526;

double w(double x) { return eta1*pow(x,eta2); }
double s(double x) { return exp(-pow(x-phi*iota1,2.)/(2*iota2*pow(phi,2.))); }
double g(struct GP *,const double);
double b(struct BP *,const double);
double qn1d(double);
double firstmodel(VEC *, void *);
VEC *firstmodeld(VEC *, void *, VEC *);
double secondmodel(VEC *, void *);
VEC *secondmodeld(VEC *, void *, VEC *);
double thirdmodel(VEC *, void *);
VEC *thirdmodeld(VEC *, void *, VEC *);
double themodeli(void *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,VEC *,VEC *,VEC *,IVEC *);
double themodelid(void *,MAT *,MAT *,MAT *,MAT *,MAT *,MAT *,VEC *,VEC *,VEC *,IVEC *);
double zstar(double,double,double,double,double,double,double); 
double wstar(double,double,double,double,double);
double Q(VEC *,VEC *);
double Qn(VEC *,VEC *,int);
void Q2(struct BP *,struct GP *,VEC *,VEC *);
void xhstep(struct GP *,VEC *,VEC *);
void xstep(struct GP *,VEC *,VEC *,VEC *);
void uhstep(VEC *,VEC *,VEC *,VEC *,struct BP *,struct GP *,struct DP *,double,double);
void ustep(VEC *,VEC *,VEC *,VEC *, struct BP *,struct GP *,struct DP *,double,double);
VEC *initial(struct GP *,struct BP *,struct DP *,VEC *,VEC *);
double idxselect(double,VEC *);
VEC *idxremove(VEC *,VEC *,int);
double e(double);
double c(double);
VEC *kullback(VEC *,VEC *);
VEC *kde(VEC *,VEC *,double);
VEC *regrid(VEC *,VEC *,VEC *,VEC *);
double get_bw(VEC *);
double linbin(VEC *,VEC *,VEC *,double);
VEC *get_obs(VEC *,VEC *,int);
int linesearch(VEC *,VEC *,double,double,double *,double,double,double,double,VEC *,VEC *,double *,double (*)(VEC *,void *),VEC * (*)(VEC *,void *,VEC *),void *);
MAT *UpdateHessian(MAT*, VEC*,VEC*);
double bfgsrun(double (*)(VEC *,void *),VEC * (*)(VEC *,void *,VEC *),VEC *,void *);
double pickAlphaWithinInterval(double,double,double,double,double,double,double,double);
int bracketingPhase(VEC *,VEC *,double,double,double *,double,double,double,double,VEC *,VEC *,double *,double *,double *,double *,double *,double *,double *,double (*)(VEC *,void *),VEC *(*)(VEC *,void *,VEC *),void *);
int sectioningPhase(VEC *,VEC *,VEC *,VEC *,double (*)(VEC *,void *),VEC *(*)(VEC *,void *,VEC *),void *,double *,double *,double,double,double,double,double,double,double,double,double,double,double);
void interpolatingCubic(double,double,double,double,double,double,VEC *);
double globalMinimizerOfPolyInInterval(double,double,VEC *);
double polyval(VEC *, double);
int roots(VEC *,VEC *);
VEC *numgrad(double (*)(VEC *,void *),void *,VEC *,double);

double h,k;

VEC *cat_abs;
VEC *eff_abs;
VEC *cat;
VEC *eff;
SPMAT *spobs;

int main(int argc, char *argv[])
{
 
  feenableexcept(FE_DIVBYZERO); 
  feenableexcept(FE_INVALID); // | FE_OVERFLOW);
  //feenableexcept(FE_UNDERFLOW);
  //FE_ALL_EXCEPT & ~FE_INEXACT);

  float *ct; 
  float *ti; 
  int Nce;

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
      if (!strcmp(argv[i], "-ce")) 
	{

	  if (i+1 == argc)
	    {
	      printf("\n no file specified\n");
	      exit(0);
	    }

	  FILE *fp1;
	  fp1 = fopen(argv[i+1],"r");
	  fscanf(fp1,"%d",&Nce);
	  ct = (float *) calloc(Nce,sizeof(float));
	  ti = (float *) calloc(Nce,sizeof(float));

	  for (int i=0;i<Nce;i++)
	    fscanf(fp1,"%f %f", &ct[i],&ti[i]);

	  fclose(fp1);
	  break;
	}
      else if (!strcmp(argv[i], "-testx"))
	{

	  float kappa,omega,kread;
	  int J;
	  sscanf(argv[i+1],"%f",&kappa);
	  sscanf(argv[i+2],"%f",&omega);
	  sscanf(argv[i+3],"%d",&J);
	  sscanf(argv[i+4],"%f",&kread);

	  struct GP gp;
	  gp.kappa=kappa;
	  gp.omega=omega;

	  h = omega/J;
	  k = kread;

	  VEC *xtest = v_get(J+1);
	  VEC *xtesth = v_get(J+2);
	  VEC *xtestn = v_get(J+2);

	  for (int j=1;j<xtest->dim;j++)
	    xtest->ve[j] = h*j;

	  xhstep(&gp,xtest,xtesth);
	  xstep(&gp,xtest,xtesth,xtestn);

	  for (int j=0;j<xtestn->dim;j++)
	    printf("%f ",xtestn->ve[j]);

	  exit(0);
	}
      else if (!strcmp(argv[i], "-testu"))
	{

	  float kappa,omega,beta,gamma,alpha1,alpha2,kread,iota;
	  int J;
	  sscanf(argv[i+1],"%f",&kappa);
	  sscanf(argv[i+2],"%f",&omega);
	  sscanf(argv[i+3],"%f",&beta);
	  sscanf(argv[i+4],"%f",&gamma);
	  sscanf(argv[i+5],"%f",&alpha1);
	  sscanf(argv[i+6],"%f",&alpha2);
	  sscanf(argv[i+7],"%f",&iota);
	  sscanf(argv[i+8],"%d",&J);
	  sscanf(argv[i+9],"%f",&kread);

	  struct GP gp;
	  gp.kappa=kappa;
	  gp.omega=omega;

	  struct DP dp;
	  dp.beta=beta;
	  dp.gamma=gamma;
	  dp.iota=iota;

	  struct BP bp;
	  bp.alpha1=alpha1;
	  bp.alpha2=alpha2;

	  h = omega/J;
	  k = kread;

	  cat = v_get(1);
	  eff = v_get(1);

	  VEC *xtest = v_get(J+1);
	  VEC *xtesth = v_get(J+2);
	  VEC *xtestn = v_get(J+2);
	  VEC *utest = v_get(J+1);
	  VEC *utesth = v_get(J+2);
	  VEC *utestn = v_get(J+2);

	  for (int j=1;j<xtest->dim;j++)
	    xtest->ve[j] = h*j;

	  utest = initial(&gp,&bp,&dp,xtest,utest);

	  xhstep(&gp,xtest,xtesth);
	  uhstep(xtest,xtesth,utest,utesth,&bp,&gp,&dp,Q(xtest,utest),k/2);

	  xstep(&gp,xtest,xtesth,xtestn);
	  ustep(xtesth,utest,xtestn,utestn,&bp,&gp,&dp,Q(xtesth,utesth),k);

	  for (int j=0;j<utestn->dim;j++)
	    printf("%f ",utestn->ve[j]);

	  exit(0);
	}
      else if (!strcmp(argv[i], "-testl"))
	{

	  float kappa,omega,beta,gamma,alpha1,alpha2,kread,iota;
	  int J,I;
	  sscanf(argv[i+1],"%f",&kappa);
	  sscanf(argv[i+2],"%f",&omega);
	  sscanf(argv[i+3],"%f",&beta);
	  sscanf(argv[i+4],"%f",&gamma);
	  sscanf(argv[i+5],"%f",&alpha1);
	  sscanf(argv[i+6],"%f",&alpha2);
	  sscanf(argv[i+7],"%f",&iota);
	  sscanf(argv[i+8],"%d",&J);
	  sscanf(argv[i+9],"%f",&kread);
	  sscanf(argv[i+8],"%d",&I);

	  struct GP gp;
	  gp.kappa=kappa;
	  gp.omega=omega;

	  struct DP dp;
	  dp.beta=beta;
	  dp.gamma=gamma;
	  dp.iota=iota;

	  struct BP bp;
	  bp.alpha1=alpha1;
	  bp.alpha2=alpha2;

	  h = omega/J;
	  k = kread;

	  cat = v_get(1);
	  eff = v_get(1);

	  MAT *x = m_get(I+1,J+1);
	  MAT *u = m_get(I+1,J+1);

	  VEC *Qi = v_get(I+1);
	  VEC *xn_tmp; VEC *un_tmp; VEC *xh_tmp; VEC *uh_tmp;
	  int J1 = x->n+1;
	  xn_tmp = v_get(J1); un_tmp = v_get(J1);
	  xh_tmp = v_get(J1); uh_tmp = v_get(J1);

	  //	  VEC *C = v_get(x->m);

	  VEC *tmp1; VEC *tmp2;
	  tmp1 = v_get(x->n);
	  tmp2 = v_get(x->n);

	  for (int j=1;j<x->n;j++) 
	    x->me[0][j] = h*j;
 
	  set_row(u,0,initial(&gp,&bp,&dp,get_row(x,0,tmp1),tmp2));
 
	  Qi->ve[0] = Q(get_row(x,0,tmp1),get_row(u,0,tmp2));

	  for (int i=1;i<x->m;i++)
	    {

	      xhstep(&gp,get_row(x,i-1,tmp1),xh_tmp);
	      uhstep(get_row(x,i-1,tmp1),xh_tmp,get_row(u,i-1,tmp2),uh_tmp,&bp,&gp,&dp,Qi->ve[i-1],k*(i-1));

	      double Qh = Q(xh_tmp,uh_tmp);

	      xstep(&gp,get_row(x,i-1,tmp1),xh_tmp,xn_tmp);
	      ustep(xh_tmp,get_row(u,i-1,tmp1),xn_tmp,un_tmp,&bp,&gp,&dp,Qi->ve[i-1],k*(i-.5));

	      Qi->ve[i] = Q(xn_tmp,un_tmp);
      
	      int idx = idxselect(gp.omega,xn_tmp);

	      set_row(x,i,idxremove(xn_tmp,tmp1,idx));
	      set_row(u,i,idxremove(un_tmp,tmp1,idx));

	    }

	  /*      VEC * ctm = v_get(x->n+1);
      for (int j=0;j<x->n+1;j++)
	ctm->ve[j] = s(xn_tmp->ve[j])*iota*e(k*(i-I))*w(xn_tmp->ve[j])*un_tmp->ve[j];
      C->ve[i] = Q(xn_tmp,ctm)/1e5;
      V_FREE(ctm);
	  */

	  for (int j=0;j<u->n;j++)
	    printf("%f ",u->me[I][j]);

	  exit(0);
	}
      else
	{
	  printf("\n");
	  printf("	proper usage requires at least one filename specified\n");
	  printf("		e.g. spade -ce ce.dat or ... \n");
	  printf("\n");
	  printf("	Other arguments:\n");
	  printf("			-ce   [no default] catch effort data file\n");
	  exit(0);
	}
    }

  }

  h = 173./400.0;
  k = .025;

  cat = v_get(961);
  eff = v_get(961);

  for (int i=0;i<Nce;i++)
    {
      int idxr = floor((ti[i]+k/2)/k);
      cat->ve[idxr] += (ct[i]/k)/1e3;
      eff->ve[idxr] += 1.0/k;
    }

  free(ct);
  free(ti);

  double iota = 1.1;

  qn1d(iota);

  return(0);

}  

double qn1d(

	    double iota

	    )
{

  struct GP gp;
  gp.kappa=.1;
  gp.omega=173;

  struct BP bp;
  bp.alpha1 = .1;
  bp.alpha2 = .01;

  struct DP dp;
  dp.beta = 1;
  dp.gamma = 1e-7;
  dp.iota = iota;

  struct TMI tmi;
  tmi.gp = gp;
  tmi.dp = dp;
  tmi.bp = bp;

  tmi.I = (int)24/k;
  tmi.J = 400;

  double newiota = iota;
  double eps = 1e-4;

  for (int i=0;i<200;i++)
  {

    //  printf("\n");
    // for (int j=20;j<=54;j+=2)
    // {

      MAT *x;
      MAT *u;

      int LI = 2*tmi.I + 1;
      int J = tmi.J+1;

      x = m_get(LI,J);
      u = m_get(LI,J);

      MAT *xh; MAT *uh; MAT *xn; MAT *xhh;

      xh = m_get(LI,J+1);
      uh = m_get(LI,J+1);
      xn = m_get(LI,J+1);
      xhh = m_get(LI,J+1);

      VEC *Ui = v_get(LI);
      VEC *Uh = v_get(LI);
      VEC *Uhh = v_get(LI);
      IVEC *idxi = iv_get(LI-1);

      //      double newi = j/20.0;
      iota = newiota;
      tmi.dp.iota = iota;
      
      double fa = themodeli((void *)&tmi,x,u,xhh,xh,xn,uh,Ui,Uh,Uhh,idxi);
      double fad = themodelid((void *)&tmi,x,u,xhh,xh,xn,uh,Ui,Uh,Uhh,idxi);

      //      tmi.dp.iota = newi+exp(-10.0);
      //double f1 = themodeli((void *)&tmi,x,u,xhh,xh,xn,uh,Ui,Uh,Uhh,idxi);
      //printf("%f %f\n",newi,(f1-fa)/exp(-10.0));

      /*            
      for (int j=-2;j<10;j++)
	{
	  double delta = exp(-(double)j);
	  tmi.dp.iota = 2.5+delta;
	  double ng = (f1 - fa) / (delta);

	  printf("%g %g\n",delta,ng);
	}
      */
      //      printf("%f %f %f\n",2.5,fa,fad);
      //exit(1);
      //fa = themodeli(6.8,(void *)&tmi,x,u,xh,uh,Ui,Uh,xn,un,idxi);
      //fad = themodelid(6.8,(void *)&tmi,x,u,xh,uh,Ui,Uh,xn,un,idxi);
      //printf("%f %f\n",fa,fad);
      
      if (fabs(fad) < eps)
      	break;
      else
      	newiota = iota - 4e-6 * fad;
      
      printf("%d %f %f %f %f\n",i,newiota,iota,fa,fad);
      
      //      printf("%g %g\n",newi,fad);      

      M_FREE(x);
      M_FREE(u);
      M_FREE(xh);
      M_FREE(xhh);
      M_FREE(uh);
      M_FREE(xn);
      V_FREE(Ui);
      V_FREE(Uh);
      V_FREE(Uhh);
      IV_FREE(idxi);

      //exit(1);

  }

  /*
  printf("e\n");

  printf("\n");

  for (int j=20;j<=54;j+=2)
  {

      MAT *x;
      MAT *u;

      int LI = 2*tmi.I + 1;
      int J = tmi.J+1;

      x = m_get(LI,J);
      u = m_get(LI,J);

      MAT *xh; MAT *uh; MAT *xn; MAT *xhh;

      xh = m_get(LI,J+1);
      uh = m_get(LI,J+1);
      xn = m_get(LI,J+1);
      xhh = m_get(LI,J+1);

      VEC *Ui = v_get(LI);
      VEC *Uh = v_get(LI);
      VEC *Uhh = v_get(LI);
      IVEC *idxi = iv_get(LI-1);

      double newi = j/20.0;
      //      iota = newiota;
      tmi.dp.iota = newi;
      
      double fa = themodeli((void *)&tmi,x,u,xhh,xh,xn,uh,Ui,Uh,Uhh,idxi);
      double fad = themodelid((void *)&tmi,x,u,xhh,xh,xn,uh,Ui,Uh,Uhh,idxi);

      //      tmi.dp.iota = newi+exp(-10.0);
      //double f1 = themodeli((void *)&tmi,x,u,xhh,xh,xn,uh,Ui,Uh,Uhh,idxi);
      printf("%f %f\n",newi,fad); //(f1-fa)/exp(-10.0));

      /*            
      for (int j=-2;j<10;j++)
	{
	  double delta = exp(-(double)j);
	  tmi.dp.iota = 2.5+delta;
	  double ng = (f1 - fa) / (delta);

	  printf("%g %g\n",delta,ng);
	}
      */

      //      printf("%f %f %f\n",2.5,fa,fad);
      //exit(1);

      //fa = themodeli(6.8,(void *)&tmi,x,u,xh,uh,Ui,Uh,xn,un,idxi);
      //fad = themodelid(6.8,(void *)&tmi,x,u,xh,uh,Ui,Uh,xn,un,idxi);
      
      //printf("%f %f\n",fa,fad);
      /*
      if (fabs(fad) < eps)
      	break;
      else
      	newiota = iota - 5e-4 * fad;
      
      printf("%d %f %f %f %f\n",i,newiota,iota,fa,fad);
      */
      //      printf("%g %g\n",newi,fad);      
  /*
      M_FREE(x);
      M_FREE(u);
      M_FREE(xh);
      M_FREE(xhh);
      M_FREE(uh);
      M_FREE(xn);
      V_FREE(Ui);
      V_FREE(Uh);
      V_FREE(Uhh);
      IV_FREE(idxi);

      //exit(1);


  }

       printf("e");
  */

  return iota;


}
  
double bfgsrun(

	       double (*model)(VEC *,void *),
	       VEC * (*modeld)(VEC *,void *,VEC *),
	       VEC *x,
	       void *stuff

	       )	       
{

  VEC *oldx = v_get(x->dim);
  VEC *oldgrad = v_get(x->dim);
  VEC *x_new = v_get(x->dim);
  VEC *delta_x = v_get(x->dim);
  VEC *delta_grad = v_get(x->dim);
  MAT *H = m_get(x->dim,x->dim);
  VEC *grad = v_get(x->dim);
  double f;
  int stop,iter;

  stop=0; iter=0;

  m_ident(H);

  f = (*model)(x,stuff);

  grad = (*modeld)(x,stuff,grad); 

  /*
  v_output(grad);
  VEC * ngrad = v_get(x->dim);
  //ngrad = numgrad(model,stuff,0.001);
  //v_output(ngrad);
  printf("\n");
  for (int i=0;i<100;i++) 
    {
      double eps = exp(-(double)i);
      ngrad = numgrad(model,stuff,x,eps);

      printf("%g %g\n",eps,ngrad->ve[1]);
    }      
  printf("e\n");
  exit(1);
  */

  while (stop==0)
    {
      iter=iter+1;
      if (iter == 150)
	break;

      VEC *dir = v_get(x->dim);
      mv_mlt(H,grad,dir);
      sv_mlt(-1.0,dir,dir);

      double dirD = in_prod(grad,dir);

      if (1) 
	{
	  int nzero = 0;
	  for (int i=0;i<x->dim;i++)
	    if (grad->ve[i]==0)
	      nzero++;
	  v_output(x);
          //v_output(grad);
          printf("f %f\n",f);
	  if ((nzero>=(int)(1.0*x->dim)) || (dirD > -1e-15))
	    {
	      stop=1;
	      break;
	    }
	}

      double alpha=1.0;
      if (iter ==1)
	{
	  alpha = 1.0 / v_norm_inf(grad);
	  if (alpha>1.0) alpha = 1.0;
	}

      oldx = v_copy(x,oldx);
      oldgrad = v_copy(grad,oldgrad);

      double rho=1e-2;
      double sigma = 0.1;
      double TolFun = 1e-11;
      double fminimum = -1e+10;
      double f_new = f;

      linesearch(x,dir,f,dirD,&alpha,rho,sigma,TolFun,fminimum,grad,x_new,&f_new,model,modeld,stuff);

      f = f_new;
      sv_mlt(alpha,dir,delta_x);
      v_add(x,delta_x,x);
      v_sub(grad,oldgrad,delta_grad);

      H = UpdateHessian(H,delta_x,delta_grad);

    }

}

int linesearch(

	       VEC *x,
	       VEC *dir,
	       double f,
	       double dirD,
	       double *alpha,
	       double rho,
	       double sigma,
	       double TolFun,
	       double fminimum,
	       VEC *grad,
	       VEC *x_new,
	       double *f_new,
	       double (*model)(VEC *,void *),
	       VEC *(*modeld)(VEC *,void *,VEC *),
	       void *stuff

	       )
{
  if (dirD * 0 != 0)
    return 1;

  double a,b,f_a,fPrime_a,f_b,fPrime_b;
  int flag = bracketingPhase(x,dir,f,dirD,alpha,rho,sigma,TolFun,fminimum,grad,x_new,&a,&b,&f_a,&fPrime_a,&f_b,&fPrime_b,f_new,model,modeld,stuff);

  //printf("linesearch, bracketing flag = %d\n",flag);

  if (flag==2) 
    {
      int flag2 = sectioningPhase(x,dir,grad,x_new,model,modeld,stuff,f_new,alpha,f,dirD,a,b,f_a,fPrime_a,f_b,fPrime_b,rho,sigma,TolFun);
      //printf("linesearch, sectioning flag = %d\n",flag2);
      return flag2;
    }
  else
    return flag;

    //sectioningPhase(x,dir,f,dirD,a,b,f_a,fPrime_a,f_b,rho,sigma,TolFun,grad,x_new,f_new,alpha,gp,xx,obs);

}

int bracketingPhase(

		    VEC *xInitial,
		    VEC *dir,
		    double fInitial,
		    double fPrimeInitial,
		    double *alpha,
		    double rho,
		    double sigma,
		    double TolFun,
		    double fminimum,
		    VEC *grad,
		    VEC *x_new,
		    double *a,
		    double *b,
		    double *f_a,
		    double *fPrime_a,
		    double *f_b,
		    double *fPrime_b,
		    double *f_new,
		    double (*model)(VEC *,void *),
		    VEC *(*modeld)(VEC *,void *,VEC *),
		    void *stuff

		    )

{ /* Fletcher, R. 1987 Practical methods of optimization. 2nd Edition. Page 34 */
  double tau1 = 9.0;
  double f_alpha = fInitial;
  double fPrime_alpha = fPrimeInitial;
  double alphaMax = (fminimum - fInitial)/(rho*fPrimeInitial);
  int iter=0;
  int iterMax=1000;
  double alphaPrev=0.0;

  while (iter < iterMax)
    {

      iter++;
      double fPrev = f_alpha;
      double fPrimePrev = fPrime_alpha;

      //printf("%d %f\n",iter,f_alpha);

      VEC *dirtmp = v_get(xInitial->dim);
      sv_mlt(*alpha,dir,dirtmp);
      v_add(xInitial,dirtmp,x_new);

      double alphaold = *alpha;
      int negflag = 0;
      while (x_new->ve[0] < 0 || x_new->ve[1] < 0)
	{
          alphaold = *alpha;
          *alpha = .67*(*alpha); 
	  sv_mlt(*alpha,dir,dirtmp);
	  v_add(xInitial,dirtmp,x_new);
	  negflag=1;
	}

      if (negflag)
	alphaMax = alphaold;

      f_alpha = (*model)(x_new,stuff);
      *f_new = f_alpha;
      grad = (*modeld)(x_new,stuff,grad); 
      fPrime_alpha = in_prod(grad,dir);

      if (f_alpha < fminimum)
	return 1;

      // Bracket located - case 1
      if ((f_alpha > fInitial + (*alpha)*rho*fPrimeInitial) || (f_alpha >= fPrev))
	{
	  *a = alphaPrev; *b = *alpha;
	  *f_a = fPrev; *fPrime_a = fPrimePrev;
	  *f_b = f_alpha; *fPrime_b = fPrime_alpha;
	  return 2;
	}

      if (fabs(fPrime_alpha) <= -sigma*fPrimeInitial)
	return 0;

      //Bracket located - case 2
      if (fPrime_alpha >= 0)
	{ 
	  *a = *alpha; *b = alphaPrev;
	  *f_a = f_alpha; *fPrime_a = fPrime_alpha;
	  *f_b = fPrev; *fPrime_b = fPrimePrev;
	  return 2;
	}

      if (2.0 * (*alpha) - alphaPrev < alphaMax) 
	{
	  double brcktEndPntA = 2*(*alpha) -alphaPrev;
	  double brcktEndPntB = min(alphaMax,(*alpha)+tau1*((*alpha)-alphaPrev));
	  double alphaNew = pickAlphaWithinInterval(brcktEndPntA,brcktEndPntB,alphaPrev,(*alpha),fPrev,fPrimePrev,f_alpha,fPrime_alpha);
	  alphaPrev = (*alpha);
	  (*alpha) = alphaNew;
	}
      else
	(*alpha) = alphaMax;
    }
}

int sectioningPhase(

		    VEC *xInitial,
		    VEC *dir,
		    VEC *grad,
		    VEC *x_new,
		    double (*model)(VEC *,void *),
		    VEC *(*modeld)(VEC *,void *,VEC *),
		    void *stuff,
		    double *f_new,
		    double *alpha,
		    double fInitial,
		    double fPrimeInitial,
		    double a,
		    double b,
		    double f_a,
		    double fPrime_a,
		    double f_b,
		    double fPrime_b,
		    double rho,
		    double sigma,
		    double TolFun

		    )
{/* Fletcher, R. 1987 Practical methods of optimization. 2nd Edition. Page 35 */
  double eps = 1e-16;
  double tau2 = min(0.1,sigma);
  double tau3 = 0.5;
  double tol = TolFun/1000;

  int iter = 0;
  int iterMax = 100;
  while (iter < iterMax)
    {
      // printf("%d ",iter);
      iter = iter + 1;

      // Pick alpha in reduced bracket
      double brcktEndpntA = a + tau2*(b - a); 
      double brcktEndpntB = b - tau3*(b - a);

      // Find global minimizer in bracket [brcktEndpntA,brcktEndpntB] of 3rd-degree 
      // polynomial that interpolates f() and f'() at "a" and at "b".
      (*alpha) = pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,a,b,f_a,fPrime_a,f_b,fPrime_b);  

      // Evaluate f(alpha) and f'(alpha)
      VEC *dirtmp = v_get(xInitial->dim);
      sv_mlt(*alpha,dir,dirtmp);
      v_add(xInitial,dirtmp,x_new);

      double f_alpha = (*model)(x_new,stuff);

      *f_new = f_alpha;
      grad = (*modeld)(x_new,stuff,grad);

      double fPrime_alpha = in_prod(grad,dir);

      // Check if roundoff errors are stalling convergence.
      // Here the magnitude of a first order estimation on the
      // change in the objective is checked. The value of tol 
      // was chosen empirically through experimentation.
      if (fabs( ((*alpha) - a)*fPrime_a ) <= tol)
	return -2;  // No acceptable point could be found

      // Update bracket
      double aPrev = a; double bPrev = b; double f_aPrev = f_a; double f_bPrev = f_b; 
      double fPrime_aPrev = fPrime_a; double fPrime_bPrev = fPrime_b;
      if ((f_alpha > (fInitial + (*alpha)*rho*fPrimeInitial)) || (f_alpha >= f_a))
	{
	  a = aPrev; b = (*alpha);
	  f_a = f_aPrev; f_b = f_alpha;
	  fPrime_a = fPrime_aPrev; fPrime_b = fPrime_alpha;
	}
      else 
	{
  
	  if (fabs(fPrime_alpha) <= -sigma*fPrimeInitial)
	    return 0;	// Acceptable point found
	  a = (*alpha); f_a = f_alpha; fPrime_a = fPrime_alpha;
	  if ((b - a)*fPrime_alpha >= 0)
	    {
	      b = aPrev; f_b = f_aPrev; fPrime_b = fPrime_aPrev;
	    }
	  else
	    {
	      b = bPrev; f_b = f_bPrev; fPrime_b = fPrime_bPrev;
	    }
	}

      // Check if roundoff errors are stalling convergence
      if (fabs(b-a) < eps)
	return -2;	// No acceptable point could be found


    }

  // We reach this point if and only if maxFunEvals was reached
  return -1;
}

double pickAlphaWithinInterval(

			       double brcktEndpntA,
			       double brcktEndpntB,
			       double alpha1,
			       double alpha2,
			       double f1,
			       double fPrime1,
			       double f2,
			       double fPrime2

			       )
{
	// alpha = pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,alpha1,alpha2,f1,fPrime1,f2,fPrime2) finds 
	// a global minimizer alpha within the bracket [brcktEndpntA,brcktEndpntB] of the cubic polynomial 
	// that interpolates f() and f'() at alpha1 and alpha2. Here f(alpha1) = f1, f'(alpha1) = fPrime1, 
	// f(alpha2) = f2, f'(alpha2) = fPrime2.

  // Find interpolating Hermite polynomial in the z-space, 
  // where z = alpha1 + (alpha2 - alpha1)*z
  VEC *coeff = v_get(4);
  interpolatingCubic(alpha1,alpha2,f1,fPrime1,f2,fPrime2,coeff);

  // Convert bounds to the z-space
  double zlb = (brcktEndpntA - alpha1)/(alpha2 - alpha1);
  double zub = (brcktEndpntB - alpha1)/(alpha2 - alpha1);

  // Make sure zlb <= zub so that [zlb,zub] be an interval
  if (zlb > zub)	// swap zlb and zub
    {      
      double tmp = zlb;
      zlb = zub;
      zub = tmp;
    }

  // Minimize polynomial over interval [zlb,zub]
  double z = globalMinimizerOfPolyInInterval(zlb,zub,coeff);
  double alpha = alpha1 + z*(alpha2 - alpha1);
  return alpha;
}

void interpolatingCubic(

			double alpha1,
			double alpha2,
			double f1, 
			double fPrime1,
			double f2,
			double fPrime2,
			VEC *coeff

			)
{
	// coeff = interpolatingCubic(alpha1,alpha2,f1,fPrime1,f2,fPrime2) determines
	// the coefficients of the cubic polynomial that interpolates f and f' at alpha1 
	// and alpha2; that is, c(alpha1) = f1, c'(alpha1) = fPrime1, c(alpha2) = f2, 
	// c'(alpha2) = fPrime2.

  double deltaAlpha = alpha2 - alpha1;
  coeff->ve[3] = f1;
  coeff->ve[2] = deltaAlpha*fPrime1;
  coeff->ve[1] = 3.0*(f2 - f1) - (2.0*fPrime1 + fPrime2)*deltaAlpha;
  coeff->ve[0] = (fPrime1 + fPrime2)*deltaAlpha - 2.0*(f2 - f1);
}

double globalMinimizerOfPolyInInterval(

				       double lowerBound,
				       double upperBound,
				       VEC* coeff

				       )
{
	// alpha = globalMinimizerOfPolyInInterval(lowerBound,upperBound,coeff) finds a
	// global minimizer alpha in the interval lowerBound <= alpha <= upperBound of 
	// the cubic polynomial defined by the coefficients in the 4-vector coeff. 

	// Find stationary points 
  VEC *coeff_tmp = v_get(3);  
  coeff_tmp->ve[0] = 3*coeff->ve[0];
  coeff_tmp->ve[1] = 2*coeff->ve[1];
  coeff_tmp->ve[2] = coeff->ve[2];
  VEC *stationaryPoint = v_get(2);
  int flag = roots(coeff_tmp, stationaryPoint);

  // Which among the two endpoints has a lower polynomial value?
  double f1 = polyval(coeff, lowerBound);
  double f2 = polyval(coeff, upperBound);

  double alpha = lowerBound;
  double fmin = f1;
  if (f2 < f1)	
    {
      alpha = upperBound;
      fmin = f2;
    }

  // If any of the stationary points is feasible, update the current global minimizer.
  if (flag == 1)
    for (int i=0; i<2; i++)
      if ((lowerBound <= stationaryPoint->ve[i]) && (stationaryPoint->ve[i] <= upperBound))
	if (polyval(coeff,stationaryPoint->ve[i]) < fmin)
	  {
	    fmin = polyval(coeff,stationaryPoint->ve[i]);
	    alpha = stationaryPoint->ve[i];
	  }

  return alpha;
}

double polyval(

	       VEC *coeff, 
	       double x

	       )
{
  return coeff->ve[0]*pow(x, 3.0) + coeff->ve[1]*pow(x, 2.0) + coeff->ve[2]*x + coeff->ve[3];
}

int roots(

	  VEC *coef,
	  VEC *stationary_point

	  )
{
  if (coef->ve[0] == 0)
    return 1;
  double m11 = -coef->ve[1] / coef->ve[0];
  double m12 = -coef->ve[2] / coef->ve[0];
  double m21 = 1.0;
  double m22 = 0.0;

  double T = (m11 + m22);
  double D = m11 * m22 - m12 * m21;
  double dd = T*T/4.0 - D;
  if (dd < 0)
    return 0;	// complex values
  double d = sqrt( dd );
  stationary_point->ve[0] = T/2.0 + d;
  stationary_point->ve[1] = T/2.0 - d;
  return 1;
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
  double rhok =  1/in_prod(y,s);
  MAT *A1 = m_sub(m_ident(eye),sm_mlt(rhok,m_mlt(sMc,yMr,MNULL),MNULL),MNULL);
  MAT *A2 = m_sub(m_ident(eye),sm_mlt(rhok,m_mlt(yMc,sMr,MNULL),MNULL),MNULL);

  H =  m_add(m_mlt(A1,m_mlt(H,A2,MNULL),MNULL),sm_mlt(rhok,m_mlt(sMc,sMr,MNULL),MNULL),MNULL);

  return H;
}

VEC *get_obs(

	     VEC *dt,
	     VEC *xx,
	     int NdnF

	     )
{
  VEC *y = v_get(NdnF);
 
  double bw = get_bw(dt);

  VEC *zx = v_get(NdnF/2);

  double hw = linbin(dt,y,zx,bw);

  VEC *zetal = v_get(NdnF/2);
  zetal = kde(y,zetal,hw);

  VEC *obs = v_get(xx->dim);

  obs = regrid(zx,xx,zetal,obs);

  double qobs = Q(xx,obs);

  for (int i=0;i<obs->dim;i++)
    obs->ve[i]=obs->ve[i]/qobs;

  return obs;

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


double secondmodel(

		   VEC *x,
		   void *stuff

		   )
{

  struct SMS sms; 

  sms = * (struct SMS *) stuff;

  double rt = 0;
  for (int i=0;i<sms.length->dim;i++)
    rt += pow(sms.length->ve[i] - x->ve[1]*(1.0 - exp(-x->ve[0]*sms.age->ve[i])),2.0);
  return rt/sms.length->dim;

}

double thirdmodel(

		  VEC *x,
		  void *stuff

		  )
{

  struct TMS tms; 

  tms = * (struct TMS *) stuff;

  double rt = 0;

  for (int i=0;i<tms.xt->dim;i++)
    {
      double e = exp(-x->ve[0]*tms.tl->ve[i]);
      double e2 = exp(-2*x->ve[0]*tms.tl->ve[i]);
      double xt = tms.xt->ve[i];
      double xr = tms.xr->ve[i];
      double w = x->ve[1];
      double tmp1 = xt+e*xr+w*e2-w*e;
      double tmp2 = e2+1;
      rt += pow(e*tmp1/tmp2 - xr + w*(1-e),2.0) + pow(tmp1/tmp2 - xt,2.0);
    }

  return rt/tms.xt->dim;

}

VEC *thirdmodeld(

		 VEC *xp,
		 void *stuff,
		 VEC *grad

		 )
{

  struct TMS tms; 

  tms = * (struct TMS *) stuff;

  double rt1 = 0;
  for (int i=0;i<=tms.xt->dim;i++)
    {
      double t = tms.tl->ve[i];
      double x = tms.xt->ve[i];
      double y = tms.xr->ve[i];
      double k = xp->ve[0];
      double w = xp->ve[1];
      double tmp1 , tmp2 , tmp3 , tmp4 , tmp5 , tmp6 , tmp7 , tmp8 ; tmp1 = 1/exp(2*k*t) ; tmp2 = tmp1+1 ; tmp3 = 1/tmp2 ; tmp4 = 1/exp(k*t) ; tmp5 = tmp4*y+x-tmp4*w+tmp1*w ; tmp6 = 1/pow(tmp2,2) ; tmp7 = t*tmp4*w ; tmp8 = -t*tmp4*y+tmp7-2*t*tmp1*w ; rt1 += 2*(tmp3*tmp4*tmp5-y+(1-tmp4)*w)*(2*tmp5*tmp6*t/exp(3*k*t)+tmp3*tmp4*tmp8+tmp7-t*tmp3*tmp4*tmp5)+2*(tmp3*tmp5-x)*(tmp3*tmp8+2*t*tmp1*tmp6*tmp5) ;
    }

  double rt2 = 0;
  for (int i=0;i<=tms.xt->dim;i++)
    {
      double t = tms.tl->ve[i];
      double x = tms.xt->ve[i];
      double y = tms.xr->ve[i];
      double k = xp->ve[0];
      double w = xp->ve[1];

      double tmp1 , tmp2 , tmp3 , tmp4 , tmp5 , tmp6 ; tmp1 = 1/exp(2*k*t) ; tmp2 = 1/(tmp1+1) ; tmp3 = 1/exp(k*t) ; tmp4 = -tmp3 ; tmp5 = tmp4+tmp1 ; tmp6 = tmp3*y+x-tmp3*w+tmp1*w ; rt2 += 2*(tmp2*tmp5*tmp3+tmp4+1)*(tmp2*tmp3*tmp6-y+(tmp4+1)*w)+2*tmp2*tmp5*(tmp2*tmp6-x) ;
    }

  grad->ve[0] = rt1/tms.tl->dim;
  grad->ve[1] = rt2/tms.tl->dim;

  return grad;

}

VEC *secondmodeld(

		  VEC *x,
		  void *stuff,
		  VEC *grad

		  )
{

  struct SMS sms; 

  sms = * (struct SMS *) stuff;

  double rt1 = 0;
  for (int i=0;i<sms.length->dim;i++)
    rt1 = rt1 + -2*x->ve[1]*sms.age->ve[i]*(sms.length->ve[i] - x->ve[1]*(1.0 - exp(-x->ve[0]*sms.age->ve[i])))*exp(-x->ve[0]*sms.age->ve[i]);

  double rt2 = 0;
  for (int i=0;i<sms.length->dim;i++)
    rt2 = rt2 + 2*(sms.length->ve[i] - x->ve[1]*(1.0 - exp(-x->ve[0]*sms.age->ve[i])))*(exp(-x->ve[0]*sms.age->ve[i])-1.0);

  grad->ve[0] = rt1/sms.length->dim;
  grad->ve[1] = rt2/sms.length->dim;

  return grad;

}

double firstmodel(

		  VEC *x,
		  void *stuff

		  )
{

  struct FMS fms;

  fms = * (struct FMS *) stuff;

  VEC *v = v_get(fms.xx->dim);

  int Ndn = fms.xx->dim-1;
  
  VEC *w = v_get(Ndn+1);

  for (int i=0;i<=Ndn;i++)
    w->ve[i] = (x->ve[0] + x->ve[1]*exp(-pow(fms.xx->ve[i]-phi*iota1,2.)/(2*iota2*pow(phi,2.)))) / (fms.gp.kappa*(fms.gp.omega-fms.xx->ve[i]));

  for (int i=1;i<=Ndn;i++)
    v->ve[i]=exp(-pow(fms.xx->ve[i]-phi*iota1,2.)/(2*iota2*pow(phi,2.)))*(1./(fms.gp.omega-fms.xx->ve[i]))*exp(-Qn(fms.xx,w,i+1));
    
  v->ve[0] = exp(-pow(fms.xx->ve[0]-phi*iota1,2.)/(2*iota2*pow(phi,2.)))*(1./(fms.gp.omega));
  
  double qv = Q(fms.xx,v);

  for (int i=0;i<=Ndn;i++)
    v->ve[i] = v->ve[i]/qv;

  V_FREE(w);
    
  v = kullback(fms.obs,v);

  double rt = Q(fms.xx,v);

  V_FREE(v);

  return rt;

}

VEC *firstmodeld(

		 VEC *x,
		 void *stuff,
		 VEC *grad

		 )
{

  struct FMS fms;

  fms = * (struct FMS *) stuff;

  VEC *rt = v_get(2);

  VEC *v = v_get(fms.xx->dim);

  int Ndn = fms.xx->dim-1;
  
  VEC *w = v_get(Ndn+1);

  for (int i=0;i<=Ndn;i++)
    w->ve[i] = (x->ve[0] + x->ve[1]*exp(-pow(fms.xx->ve[i]-phi*iota1,2.)/(2*iota2*pow(phi,2.)))) / (fms.gp.kappa*(fms.gp.omega-fms.xx->ve[i]));

  for (int i=1;i<=Ndn;i++)
    v->ve[i]=exp(-pow(fms.xx->ve[i]-phi*iota1,2.)/(2*iota2*pow(phi,2.)))*(1./(fms.gp.omega-fms.xx->ve[i]))*exp(-Qn(fms.xx,w,i+1));
    
  v->ve[0] = exp(-pow(fms.xx->ve[0]-phi*iota1,2.)/(2*iota2*pow(phi,2.)))*(1./(fms.gp.omega));
  
  double B = Q(fms.xx,v);

  VEC *w2 = v_get(Ndn+1);

  for (int i=0;i<=Ndn;i++)
    w2->ve[i] = 1./(fms.gp.kappa*(fms.gp.omega-fms.xx->ve[i])); // used to have kappa

  VEC *w3 = v_get(Ndn+1);
  for (int i=0;i<=Ndn;i++)
    w3->ve[i] = exp(-pow(fms.xx->ve[i]-phi*iota1,2.)/(2*iota2*pow(phi,2.)))*(1./(fms.gp.omega-fms.xx->ve[i]))*exp(-Qn(fms.xx,w,i+1)) * Qn(fms.xx,w2,i+1);

  double C = Q(fms.xx,w3);

  for (int i=1;i<=Ndn;i++)
    v->ve[i] = (B*Qn(fms.xx,w2,i+1)-C ) / B;

  v->ve[0] = -C/B;
    
  for (int i=0;i<fms.obs->dim;i++)
    v->ve[i]=fms.obs->ve[i]*v->ve[i];

  rt->ve[0] = Q(fms.xx,v);

  Ndn = fms.xx->dim-1;
  
  for (int i=0;i<=Ndn;i++)
    w->ve[i] = (x->ve[0] + x->ve[1]*exp(-pow(fms.xx->ve[i]-phi*iota1,2.)/(2*iota2*pow(phi,2.)))) / (fms.gp.kappa*(fms.gp.omega-fms.xx->ve[i]));

  for (int i=1;i<=Ndn;i++)
    v->ve[i]=exp(-pow(fms.xx->ve[i]-phi*iota1,2.)/(2*iota2*pow(phi,2.)))*(1./(fms.gp.omega-fms.xx->ve[i]))*exp(-Qn(fms.xx,w,i+1));
    
  v->ve[0] = exp(-pow(fms.xx->ve[0]-phi*iota1,2.)/(2*iota2*pow(phi,2.)))*(1./(fms.gp.omega));
  
  B = Q(fms.xx,v);

  for (int i=0;i<=Ndn;i++)
    w2->ve[i] =  exp(-pow(fms.xx->ve[i]-phi*iota1,2.)/(2*iota2*pow(phi,2.)))/(fms.gp.kappa*(fms.gp.omega-fms.xx->ve[i]));

  for (int i=0;i<=Ndn;i++)
    w3->ve[i] = exp(-pow(fms.xx->ve[i]-phi*iota1,2.)/(2*iota2*pow(phi,2.)))*(1./(fms.gp.omega-fms.xx->ve[i]))*exp(-Qn(fms.xx,w,i+1)) * Qn(fms.xx,w2,i+1);

  C = Q(fms.xx,w3);

  for (int i=1;i<=Ndn;i++)
    v->ve[i] = (B*Qn(fms.xx,w2,i+1)-C ) / B;

  v->ve[0] = -C/B;

  V_FREE(w);
  V_FREE(w2);
  V_FREE(w3);
    
  for (int i=0;i<fms.obs->dim;i++)
    v->ve[i]=fms.obs->ve[i]*v->ve[i];

  rt->ve[1] = Q(fms.xx,v);

  V_FREE(v);

  return rt;

}

double linbin(

	      VEC *dt,
	      VEC *y,
	      VEC *zx,
	      double bw

	      )
{

  double from = dt->ve[0] - 3*bw;
  double to = dt->ve[7424] + 3*bw;
  double lo = from - 4*bw;
  double up = to + 4*bw;
  
  double step = (up-lo) / (double)y->dim;
  double ainc = 1./((double)dt->dim * step);
  double hw = bw / step;
  double fac1 = 32.0*pow(atan(1.0) * hw / (double)y->dim,2.0);

  double dlo1 = lo - step;

  for (int i=0;i<dt->dim;i++)
    {
      double wt = (dt->ve[i] - dlo1) / step;
      int jj = (int)wt;
      if (jj >= 1 && jj <= y->dim) 
	{
	  wt = wt - (double)jj;
	  double winc = wt * ainc;
	  int kk = jj+1;
	  if (jj==y->dim) kk=1;  // bug? Ndn ?
	  y->ve[jj-1]=y->ve[jj-1] + ainc - winc;
	  y->ve[kk-1]=y->ve[kk-1]+winc;
	}
    }

  for (int i=0;i<zx->dim;i++)
    zx->ve[i] = from + i*step*2;
  
  return hw;

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
  
  double fstQ = dt->ve[1855] + dt->ve[1856] / 2;
  double thrQ = dt->ve[5569] + dt->ve[5570] / 2;
  double IQR = thrQ - fstQ;
  
  IQR = IQR/1.34;

  double lo = (hi < IQR) ? hi : IQR;
  double bw = 0.9 * lo * pow((double)dt->dim,-0.2);

  PX_FREE(order);

  return bw;

}

VEC *kullback(

	      VEC *obs,
	      VEC *pred

	      )
{

  for (int i=0;i<obs->dim;i++)
    pred->ve[i]=obs->ve[i]*log(obs->ve[i]/(pred->ve[i]+1e-12) + 1e-12);

  return pred;
}


VEC *kde(

	 VEC *y,
	 VEC *zy,
	 double hw

	 )
{ /* 
    Kernel density estimation: 
    Rosenblatt (1956)
    Monro (1976)
    Silverman (1982) - main alg
    Jones and Lotwick (1984) - modification of alg
    Fan and Marron (1993) ? haven't read
  */

  double fac1 = 32.0*pow(atan(1.0) * hw / (double)y->dim,2.0);

  VEC *kim = v_get(y->dim);

  fft(y,kim);

  VEC *zim = v_get(zy->dim);

  for (int i=0;i<zy->dim;i++)
    {
      double rjfac = i*i*fac1;
      double bc = 1.0 - rjfac/(hw*hw*6.0);
      double fac = exp(-rjfac)/bc;
      zy->ve[i] = fac * y->ve[i];
      zim->ve[i] = fac * kim->ve[i];
    }

  ifft(zy,zim);

  for (int i=0;i<zy->dim;i++)
    zy->ve[i] = zy->ve[i]  - 0.0124531; // why do i need this?!

  V_FREE(kim);
  V_FREE(zim);
    
  return zy;

}


VEC *regrid(

	    VEC *zx,
	    VEC *xx,
	    VEC *zy,
	    VEC *obs

	    )
{ /* regrid from FFT grid to main grid */

  int i;
  for (i=0;i<xx->dim;i++)
    if (xx->ve[i] > zx->ve[0])
      break;

  int longidx=0;

  for (int j=i;j<xx->dim;j++)
    {

      if (zx->ve[longidx+1]<xx->ve[j])
        longidx = longidx + 1;

      if (longidx+1>=zy->dim)
	break;

      double m = (zy->ve[longidx+1] - zy->ve[longidx]) / (zx->ve[longidx+1] - zx->ve[longidx]);
      double c = zy->ve[longidx] - m*zx->ve[longidx];
      obs->ve[j] = m*xx->ve[j] + c;
      longidx = longidx + 1;

      if (longidx+1>=zy->dim)
	break;

    }

  return obs;
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
  
VEC *initial(

	     struct GP *gparams,
	     struct BP *bparams,
	     struct DP *dparams,
	     VEC *x,
	     VEC *u

	     )

{

  double inter = 27*bparams->alpha1*pow(gparams->kappa,2.)*gparams->omega+54*bparams->alpha2*pow(gparams->kappa,2.)*pow(gparams->omega,2.);
  double inter2 = pow(sqrt(pow(inter,2.) + 4*pow(-3*bparams->alpha1*gparams->kappa*gparams->omega - 3*pow(gparams->kappa,2.),3.)) + inter,1/3.);
  double Z = inter2 / (3*pow(2.,1/3.)) - pow(2.,1/3.) * (-3*bparams->alpha1*gparams->kappa*gparams->omega - 3*pow(gparams->kappa,2.)) / (3*inter2); 
  double ubar = (Z - dparams->beta - gparams->kappa) / dparams->gamma; 
  double vbar = (gparams->kappa*gparams->omega*ubar) / (dparams->beta+dparams->gamma*ubar+gparams->kappa);
  double wbar = (2*gparams->kappa*gparams->omega*vbar) / (dparams->beta+dparams->gamma*ubar+2*gparams->kappa);  

  for (int j=0;j<x->dim;j++) 
    u->ve[j] = (bparams->alpha1*vbar+bparams->alpha2*wbar)*pow(fmax(gparams->omega-x->ve[j],0.0),(dparams->beta+dparams->gamma*ubar)/gparams->kappa-1) / (gparams->kappa*pow(gparams->omega,(dparams->beta+dparams->gamma*ubar)/gparams->kappa));

  return u;
}

double g(

	 struct GP *params,
	 const double x

	 )
{ /* von-Bertalanffy growth */
  return params->kappa*(params->omega - x);
}

double b(

	 struct BP *params,
	 const double x

	 )
{ /* birth function */
  return params->alpha1*x + params->alpha2*pow(x,2.);
}
      
double themodeli(

		 void *stuff,
		 MAT *x,
		 MAT *u,
		 MAT *xhh,
		 MAT *xh,
		 MAT *xn,
		 MAT *uh,
		 VEC *Ui,
		 VEC *Uh,
		 VEC *Uhh,
		 IVEC *idxi

		 )
{
  struct TMI tmi;

  tmi = * (struct TMI *) stuff;
  VEC *xnt; VEC *un; 
  VEC *xht; VEC *uht;
  VEC *xhht; VEC *uhh;
  VEC *xt; VEC *ut;
  xnt = v_get(x->n+1);  un = v_get(x->n+1);
  xht = v_get(x->n+1);  uht = v_get(x->n+1);
  xhht = v_get(x->n+1); uhh = v_get(x->n+1);
  xt = v_get(x->n);      ut = v_get(x->n);
  VEC *C = v_get(x->m);

  for (int j=1;j<x->n;j++) 
    x->me[0][j] = h*j;
 
  set_row(u,0,initial((void *)&tmi.gp,(void *)&tmi.bp,(void *)&tmi.dp,get_row(x,0,xt),ut));
 
  Ui->ve[0] = Q(xt,ut);

  double gg = tmi.dp.gamma;
  double bb = tmi.dp.beta;
  double ii = tmi.dp.iota;
  double kk = tmi.gp.kappa;

  for (int i=1;i<x->m;i++)
    {

      double t = k*(i-tmi.I-1);
      double th = k*(i-tmi.I-.5);
      double thh = k*(i-tmi.I-.75);

      get_row(x,i-1,xt);
      get_row(u,i-1,ut);

      for (int j=1;j<=x->n;j++)
	{
	  xhht->ve[j] = xt->ve[j-1] + (k/4)*g(&tmi.gp,xt->ve[j-1]);
	  uhh->ve[j] = ut->ve[j-1]*exp(-(k/4)*zstar(bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1]));
	}

      Q2(&tmi.bp,&tmi.gp,xhht,uhh);
      Uhh->ve[i-1] = Q(xhht,uhh);
      set_row(xhh,i-1,xhht);

      for (int j=1;j<=x->n;j++)
	{
	  xht->ve[j] = xt->ve[j-1] + (k/2)*g(&tmi.gp,xhht->ve[j]);
	  uht->ve[j] = ut->ve[j-1]*exp(-(k/2)*zstar(bb,gg,kk,ii,thh,xhht->ve[j],Uhh->ve[i-1]));
	}

      Q2(&tmi.bp,&tmi.gp,xht,uht);
      Uh->ve[i-1] = Q(xht,uht);
      set_row(xh,i-1,xht);
      set_row(uh,i-1,uht);

      for (int j=1;j<=x->n;j++)
	{
	  xnt->ve[j] = xt->ve[j-1] + k*g(&tmi.gp,xht->ve[j]);
	  un->ve[j] = ut->ve[j-1]*exp(-k*zstar(bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1]));
	}

      Q2(&tmi.bp,&tmi.gp,xnt,un);
      Ui->ve[i] = Q(xnt,un);
      set_row(xn,i-1,xnt);
      
      VEC * ctm = v_get(x->n+1);
      for (int j=0;j<=x->n;j++)
	ctm->ve[j] = s(xnt->ve[j])*tmi.dp.iota*e(k*(i-tmi.I))*w(xnt->ve[j])*un->ve[j];
      C->ve[i] = Q(xnt,ctm)/1e3;
      V_FREE(ctm);
      
      int idx = idxselect(tmi.gp.omega,xnt);
      idxi->ive[i-1] = idx;

      set_row(x,i,idxremove(xnt,xt,idx));
      set_row(u,i,idxremove(un,ut,idx));

    }
  /*        
  printf("\n");
  for (int i=tmi.I;i<x->m;i++)
    printf("%f\n",C->ve[i]);
  printf("e\n");
  //exit(1);*/
  /*
  for (int i=tmi.I;i<x->m;i++)
    printf("%f\n",c(k*(i-tmi.I)));
  printf("e\n");
  
  exit(1);
  */

  VEC *chattmp = v_get(x->n);
  VEC *chat = v_get(x->m); 
  VEC *ttmp = v_get(x->m);

  for (int i=tmi.I;i<x->m;i++) {
    for (int j=0;j<x->n;j++)
      chattmp->ve[j] = s(x->me[i][j])*w(x->me[i][j])*tmi.dp.iota*e(k*(i-tmi.I))*u->me[i][j];

    chat->ve[i] = pow(Q(get_row(x,i,xt),chattmp)/1e3-c(k*(i-(tmi.I))),2.0);
    ttmp->ve[i] = k*(i-tmi.I);

  }

  double rt = Q(ttmp,chat);

  V_FREE(xt);
  V_FREE(xnt); 
  V_FREE(un); 
  V_FREE(xht); 
  V_FREE(uht);
  V_FREE(xhht); 
  V_FREE(uhh);
  V_FREE(C);
  V_FREE(xt); 
  V_FREE(ut);
  V_FREE(chattmp);
  V_FREE(chat);
  V_FREE(ttmp);

  return rt;

}

double themodelid(

		 void *stuff,
		 MAT *x,
		 MAT *u,
		 MAT *xhh,
		 MAT *xh,
		 MAT *xn,
		 MAT *uh,
		 VEC *Ui,
		 VEC *Uh,
		 VEC *Uhh,
		 IVEC *idxi

		 )
{
  struct TMI tmi;
  tmi = * (struct TMI *) stuff;

  MAT *p = m_get(x->m,x->n);

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
 
  double gg = tmi.dp.gamma;
  double bb = tmi.dp.beta;
  double ii = tmi.dp.iota;
  double kk = tmi.gp.kappa;

  VEC *Pi;
  Pi = v_get(x->m);
  
  Pi->ve[0] = Q(get_row(x,0,xt),get_row(p,0,pt));

  for (int i=1;i<x->m;i++)
    { 

      double t = k*(i-tmi.I-1);
      double th = k*(i-tmi.I-.5);
      double thh = k*(i-tmi.I-.75);

      get_row(x,i-1,xt);
      get_row(p,i-1,pt);
      get_row(xh,i-1,xht);
      get_row(xhh,i-1,xhht);
      get_row(uh,i-1,uht);
      get_row(xn,i-1,xnt);

      for (int j=1;j<=x->n;j++)
	ph->ve[j] = pt->ve[j-1]*exp(-(k/2)*zstar(bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1])) - exp(-(k/2)*zstar(bb,gg,kk,ii,t,xt->ve[j-1],Ui->ve[i-1]))*(k/2)*(s(xt->ve[j-1])*e(t)+gg*Pi->ve[i-1])*ut->ve[j-1];
	 
      Q2(&tmi.bp,&tmi.gp,xht,ph);
      double Ph = Q(xht,ph);

      for (int j=1;j<=x->n;j++)
	pn->ve[j] = pt->ve[j-1]*exp(-k*zstar(bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1])) - exp(-k*zstar(bb,gg,kk,ii,th,xht->ve[j],Uh->ve[i-1]))*k*(s(xht->ve[j])*e(th)+gg*Ph)*uht->ve[j]*exp((k/2)*zstar(bb,gg,kk,ii,thh,xhht->ve[j],Uhh->ve[i-1]));

      Q2(&tmi.bp,&tmi.gp,xnt,pn);
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

  VEC *chattmp = v_get(x->n);
  VEC *chattmp2 = v_get(x->n);
  VEC *chat = v_get(x->m); 
  VEC *ttmp = v_get(x->m);

  for (int i=tmi.I;i<x->m;i++) {
    for (int j=0;j<x->n;j++)
      {
	chattmp->ve[j] = s(x->me[i][j])*w(x->me[i][j])*tmi.dp.iota*e(k*(i-tmi.I))*u->me[i][j];
	chattmp2->ve[j] = s(x->me[i][j])*w(x->me[i][j])*e(k*(i-tmi.I))*u->me[i][j] + s(x->me[i][j])*w(x->me[i][j])*tmi.dp.iota*e(k*(i-tmi.I))*p->me[i][j];
      }

    double swieu = Q(get_row(x,i,xt),chattmp)/1e3;

    chat->ve[i] = 2*(swieu-c(k*(i-tmi.I)))*Q(get_row(x,i,xt),chattmp2)/1e3;
    ttmp->ve[i] = k*(i-tmi.I);

  }

  double rt = Q(ttmp,chat);

  M_FREE(p);
  V_FREE(xt); 
  V_FREE(xht);
  V_FREE(xnt); 
  V_FREE(ut); 
  V_FREE(uht);
  V_FREE(pt);
  V_FREE(chattmp);
  V_FREE(chattmp2);
  V_FREE(chat);
  V_FREE(ttmp);
  V_FREE(ph);
  V_FREE(pn);
  V_FREE(Pi);
  V_FREE(xhht);

  return rt;

}

/*
void phstep(

	    VEC *x,
	    VEC *xh,
	    VEC *u,
	    VEC *p,
	    VEC *ph,
	    struct BP *bpar,
	    struct GP *gpar,
	    struct DP *dpar,
	    double U,
	    double P,
	    double t

	    )

{

  for (int j=1;j<ph->dim;j++)
    ph->ve[j] = p->ve[j-1]*exp(-(k/2)*zstar(dpar->beta,dpar->gamma,gpar->kappa,dpar->iota,t,x->ve[j-1],U)) - exp(-(k/2)*zstar(dpar->beta,dpar->gamma,gpar->kappa,dpar->iota,t,x->ve[j-1],U))*(k/2)*(s(x->ve[j-1])*e(t)+dpar->gamma*P)*u->ve[j-1];//*exp(-(k/4)*zstar(dpar->beta,dpar->gamma,gpar->kappa,dpar->iota,t,x->ve[j-1],U));

  Q2(bpar,gpar,xh,ph);

}



void pstep(

	    VEC *x,
	    VEC *xh,
	    VEC *xn,
	    VEC *u,
	    VEC *uh,
	    VEC *p,
	    VEC *ph,
	    VEC *pn,
	    struct BP *bpar,
	    struct GP *gpar,
	    struct DP *dpar,
	    double U,
	    double Uh,
	    double Ph,
	    double th,
	    double thh,
	    double t

	    )

{

  VEC *xhh = v_get(xh->dim);
  VEC *uhh = v_get(xh->dim);

  double Uhh;

  for (int j=1;j<xh->dim;j++)
    {
      xhh->ve[j] = x->ve[j-1] + (k/4)*g(gpar,x->ve[j-1]);
      uhh->ve[j] = uhh->ve[j-1]*exp(-(k/4)*zstar(dpar->beta,dpar->gamma,gpar->kappa,dpar->iota,t,x->ve[j-1],U));
    }

  Q2(bpar,gpar,xhh,uhh);

  Uhh = Q(xhh,uhh);

  for (int j=1;j<pn->dim;j++)
    pn->ve[j] = p->ve[j-1]*exp(-k*zstar(dpar->beta,dpar->gamma,gpar->kappa,dpar->iota,th,xh->ve[j],Uh)) - exp(-k*zstar(dpar->beta,dpar->gamma,gpar->kappa,dpar->iota,th,xh->ve[j],Uh))*k*(s(xh->ve[j])*e(th)+dpar->gamma*Ph)*uh->ve[j]*exp(-(k/2)*zstar(dpar->beta,dpar->gamma,gpar->kappa,dpar->iota,thh,xhh->ve[j],Uhh));

  Q2(bpar,gpar,xn,pn);

  V_FREE(xhh);
  V_FREE(uhh);

}
*/

double zstar(

	     double beta,
	     double gamma,
	     double kappa,
	     double iota,
	     double t,
	     double x,
	     double U

	     )
{ /* death function: beta + gamma U + s(x)f(t) - kappa */

  return beta + gamma*U + s(x)*iota*e(t) - kappa;

}

double e(double t)
{

  if (t<0)
    return eff->ve[0];
  else
    {
      int idx = floor((t + (k/2) - 1e-12)/k); // better to use double epsilon?
      return eff->ve[idx];
    }
}

double c(double t)
{
  if (t<0)
    return cat->ve[0];
  else
    {
      int idx = floor((t + (k/2))/k);
      return cat->ve[idx];
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

	struct BP *bpar,
	struct GP *gpar,
	VEC *x,
	VEC *u

	)
{

  if (x->dim != u->dim)
    {
      printf("%d %d\n",x->dim,u->dim);
      error(E_SIZES,"Q2");
    }

  double rt = x->ve[1] * b(bpar,x->ve[1])*u->ve[1] - x->ve[0]*b(bpar,x->ve[1])*u->ve[1]; 

  for (int j=1;j<x->dim-1;j++) 
    rt = rt + (b(bpar,x->ve[j])*u->ve[j] + b(bpar,x->ve[j+1])*u->ve[j+1]) * (x->ve[j+1]-x->ve[j]);

  u->ve[0] = rt / (2*g(gpar,0.) + x->ve[0]*b(bpar,x->ve[0]) - x->ve[1]*b(bpar,x->ve[0])); 

}

void xhstep(

	    struct GP *gpar,
	    VEC *x,
	    VEC *xh

	    )
{

  //printf("%d %d\n",x->dim,xh->dim);

  if (x->dim != xh->dim-1)
    error(E_SIZES,"xhstep");

  xh->ve[0]=0;
  for (int i=1;i<xh->dim;i++)
    xh->ve[i]=x->ve[i-1]+(k/2)*g(gpar,x->ve[i-1]);
  
}

void uhstep(

	    VEC *x,
	    VEC *xh,
	    VEC *u,
	    VEC *uh,
	    struct BP *bpar,
	    struct GP *gpar,
	    struct DP *dpar,
	    double Qi,
	    double t

	    )

{ /* half step */

  for (int i=1;i<uh->dim;i++)
    uh->ve[i]=u->ve[i-1]*exp(-(k/2)*zstar(dpar->beta,dpar->gamma,gpar->kappa,dpar->iota,t,x->ve[i-1],Qi));

  Q2(bpar,gpar,xh,uh);

}

void xstep(

	   struct GP *gpar,
	   VEC *x,
	   VEC *xh,
	   VEC *xn

	   )
{ 

  if (x->dim != xh->dim-1)
    error(E_SIZES,"xstep x-xh");

  if (xh->dim != xn->dim)
    error(E_SIZES,"xstep xh-xn");

  xn->ve[0]=0;
  for (int i=1;i<xn->dim;i++)
    xn->ve[i]=x->ve[i-1]+k*g(gpar,xh->ve[i]);
  
}

void ustep(

	   VEC *xh,
	   VEC *u,
	   VEC *xn,
	   VEC *un,
	   struct BP *bpar,
	   struct GP *gpar,
	   struct DP *dpar,
	   double Qh,
	   double t
	   
	   )
{

  for (int i=1;i<un->dim;i++)
    un->ve[i]=u->ve[i-1]*exp(-k*zstar(dpar->beta,dpar->gamma,gpar->kappa,dpar->iota,t,xh->ve[i-1],Qh));

  Q2(bpar,gpar,xn,un);

}


  //x = m_get(I+1,J+1);
  //u = m_get(I+1,J+1);

  //pop(g,(void *)&gp,b,(void *)&bp,zstar,(void *)&dp,x,u);
  //printf("%f\n",u->me[I][J]);
  //> write.table(t(csskc$tl/10),file="mesch.dat",sep=" ",row.names=FALSE,col.names=FALSE)

  /*
  FILE *ifp;
  ifp = fopen("mesch.dat","r");
  int Ndt=7425;
  float data[Ndt];

  for (int i=0;i<Ndt;i++)
    fscanf(ifp,"%f", &data[i]);

  fclose(ifp);

  VEC *dt = v_get(Ndt);

  for (int i=0;i<Ndt;i++)
    dt->ve[i] = (double)data[i];

  VEC *xx;
  int Ndn=1000;

  xx = v_get(Ndn+1);

  double hh = 173./(float)Ndn;
  double f = .09; 

  for (int i=0;i<=Ndn;i++)
    xx->ve[i] = hh*i;
  xx->ve[Ndn] = xx->ve[Ndn]-1e-10;

  VEC *obs = get_obs(dt,xx,1024);

  struct FMS fms;
  fms.xx = v_get(Ndn+1);

  for (int i=0;i<=Ndn;i++)
    fms.xx->ve[i] = hh*i;
  fms.xx->ve[Ndn] = xx->ve[Ndn]-1e-10;

  fms.obs = get_obs(dt,xx,1024);

  fms.gp = gp;

  VEC *par = v_get(2);

  par->ve[0] = 0.5; //dp.beta;
  par->ve[1] = 0.5;

  double of = firstmodel(par,(void *)&fms); //objfun((void *)&gp,par,xx,obs);

  VEC *rt = v_get(2);
  rt = firstmodeld(par,(void *)&fms,rt);
  */

  //printf("%f\n",of);
  //printf("%.20f %f\n",rt->ve[0],rt->ve[1]);

  //bfgsrun(firstmodel,firstmodeld,par,(void *)&fms);

  //printf("\n");
  //for (int i=0;i<=Ndn;i++)
  //  printf("%f %f\n", xx->ve[i],v->ve[i]);
  //printf("e\n");

  //printf("\n");
  //for (int i=0;i<Ndn;i++)
  //  printf("%f %f\n",xx->ve[i],obs->ve[i]);
  //printf("e\n");  

  /*
  FILE *ifp2;
  ifp2 = fopen("meschla.dat","r");
  int Ndtla=6413;
  float datal[Ndtla];
  float dataa[Ndtla];

  for (int i=0;i<Ndtla;i++)
    fscanf(ifp2,"%f %f", &dataa[i], &datal[i]);

  fclose(ifp2);

  struct SMS sms;

  sms.length = v_get(Ndtla);
  sms.age = v_get(Ndtla);

  for (int i=0;i<Ndtla;i++)
    sms.length->ve[i] = datal[i];

  for (int i=0;i<Ndtla;i++)
    sms.age->ve[i] = dataa[i];
  */

  //double ts = secondmodel(z,(void *)&sms); 

  //VEC * gr2 = v_get(2);
  //gr2 = secondmodeld(z,(void *)&sms,gr2);

  /*
  printf("\n");
  for (int i=0;i<300;i++) 
    {
      double ag = (double)i / 10.0;
      double ln = gp.omega*(1-exp(-gp.kappa*ag));
      printf("%f %f\n",ag,ln);
    }
  printf("e");

  printf("\n");
  for (int i=0;i<sms.length->dim;i++) 
      printf("%f %f\n",sms.age->ve[i],sms.length->ve[i]);
  printf("e");
  */

  /*    
  printf("\n");
  for (int i=0;i<300;i++)
    {
      for (int j=0;j<300;j++)
	{
	  par->ve[0] = (double)i/1000.0;
	  par->ve[1] = 50.0 + (double)j;
	  of = secondmodel(par,(void *)&sms);
	  printf("%f ",of);
	}
      printf("\n");
    }
  printf("e\n");
  */

  //printf("%f %f\n",gr2->ve[0],gr2->ve[1]);
  
  /*
  FILE *ifp3;
  ifp3 = fopen("meschtr.dat","r");
  int Ndtr=5610;
  float dxt[Ndtr];
  float dxr[Ndtr];
  float dtl[Ndtr];

  for (int i=0;i<Ndtr;i++)
    fscanf(ifp3,"%f %f %f", &dxt[i], &dxr[i], &dtl[i]);

  fclose(ifp3);

  struct TMS tms;

  tms.xt = v_get(Ndtr);
  tms.xr = v_get(Ndtr);
  tms.tl = v_get(Ndtr);

  for (int i=0;i<Ndtr;i++)
    tms.xt->ve[i] = dxt[i];

  for (int i=0;i<Ndtr;i++)
    tms.xr->ve[i] = dxr[i];

  for (int i=0;i<Ndtr;i++)
    tms.tl->ve[i] = dtl[i];

  double bleh = thirdmodel(z,(void *)&tms); 
  */

  //printf("%f\n",bleh);  
  //bfgsrun(thirdmodel,thirdmodeld,z,(void *)&tms);
  /*
  printf("\n");

  for (int i=0;i<tms.xt->dim;i++)
    if (tms.tl->ve[i] > 1.95 && tms.tl->ve[i] < 2.05)
      printf("%f %f\n",tms.xt->ve[i],tms.xr->ve[i]);

  printf("e\n");
  */


/*void pop(

	 void *gparams,
	 void *bparams,
	 void *dparams,
	 MAT *x,
	 MAT *u

	 )
{

  VEC *xh; VEC *uh; VEC *xn; VEC *un;
  int J1 = x->n+1;
  xh = v_get(J1); uh = v_get(J1);
  xn = v_get(J1); un = v_get(J1);

  VEC *tmp1; VEC *tmp2;
  tmp1 = v_get(x->n);
  tmp2 = v_get(x->n);

  for (int j=1;j<x->n;j++) 
    x->me[0][j] = h*j;

  VEC *Qi;
  Qi = v_get(x->m);
 
  set_row(u,0,initial(gparams,bparams,dparams,get_row(x,0,tmp1),tmp2));
 
  Qi->ve[0] = Q(get_row(x,0,tmp1),get_row(u,0,tmp2));

  for (int i=1;i<x->m;i++)
    {

      xhstep(g,gparams,get_row(x,i-1,tmp1),xh);
      uhstep(zstar,dparams,g,gparams,b,bparams,get_row(x,i-1,tmp1),xh,get_row(u,i-1,tmp2),uh,Qi->ve[i-1],0.);

      double Qh = Q(xh,uh);

      xstep(g,gparams,get_row(x,i-1,tmp1),xh,xn);
      ustep(zstar,dparams,g,gparams,b,bparams,xh,get_row(u,i-1,tmp1),xn,Qh,un,0.);

      Qi->ve[i] = Q(xn,un);

      int idx = idxselect(gparams,xn);

      set_row(x,i,idxremove(xn,tmp1,idx));
      set_row(u,i,idxremove(un,tmp1,idx));

      //printf("%f\n",Qi->ve[i]);

    }

    }*/
