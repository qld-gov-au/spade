// SPADE
// Stock assessment using PArtial Differential Equations
// Alex Campbell 'ghostofsandy' 2015

#ifndef DEBUG
#  define DEBUG 0
#endif
#define _GNU_SOURCE
#define debug_print(...) do { if (DEBUG) fprintf(stderr, __VA_ARGS__); } while (0)

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

double themodeli(double, void *,MAT *,MAT *,MAT *,MAT *,VEC *,MAT *,IVEC *);
double themodelid(double,void *,MAT *,MAT *,MAT *,MAT *,VEC *,MAT *,IVEC *);
double zstar(double,double,double,double,double,double,double); 
double wstar(double,double,double,double,double);
double Q(VEC *,VEC *);
double Qn(VEC *,VEC *,int);
void Q2(struct BP *,struct GP *,VEC *,VEC *);
void xhstep(struct GP *,VEC *,VEC *);
void xstep(struct GP *,VEC *,VEC *,VEC *);
void uhstep(VEC *,VEC *,VEC *,VEC *,struct BP *,struct GP *,struct DP *,double,double);
void ustep(VEC *,VEC *,VEC *,VEC *, struct BP *,struct GP *,struct DP *,double,double);
void phstep(VEC *,VEC *,VEC *,VEC *,VEC *,struct BP *,struct GP *,struct DP *,double,double,double);
void pstep(VEC *,VEC *,VEC *,VEC *,VEC *,struct BP *,struct GP *,struct DP *,double,double,double);
VEC *initial(struct GP *,struct BP *,struct DP *,VEC *);
double idxselect(double,VEC *);
VEC *idxremove(VEC *,VEC *,int);
//void pop(double (*)(void *,double),void *,double (*)(void *,double),void *,double (*)(void *,double),void *,MAT *,MAT *);
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
  k = .05;

  //VEC *z = v_get(2);
  //z->ve[0] = gp.kappa;
  //z->ve[1] = gp.omega;

  cat = v_get(481);
  eff = v_get(481);

  for (int i=0;i<Nce;i++)
    {
      int idxr = floor((ti[i]+k/2)/k);
      //int idxa = floor(ti[i]/k);
      int idxr2 = floor((ti[i]+(k/4))/(k/2));
      cat->ve[idxr] += ct[i];
      eff->ve[idxr] += 1.0/(k/2);
    }

  //  cat = v_get(481);
  //  eff = v_get(481);

  /*
  printf("\n");
  for (int i=0;i<=480;i++)
    printf("%.1f\n",tmi.eff->ve[i]);
  printf("e");
  */

  double iota = 1e-4;

  qn1d(iota);

  //printf("%f %f\n",fa,fad);

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
  dp.gamma = 1e-6;
  dp.iota = iota;

  struct TMI tmi;
  tmi.gp = gp;
  tmi.dp = dp;
  tmi.bp = bp;

  tmi.I = (int)24/k;
  tmi.J = 400;

  MAT *x;
  MAT *u;

  int LI = 2*tmi.I + 1;
  int J = tmi.J+1;

  x = m_get(LI,J);
  u = m_get(LI,J);

  MAT *xh; MAT *uh; MAT *xn;

  xh = m_get(LI,J+1);
  uh = m_get(LI,J+1);
  xn = m_get(LI-1,J+1);

  VEC *Ui = v_get(LI);
  IVEC *idxi = iv_get(LI-1);

  double newiota = iota;
  double eps = 100;

  for (int i=0;i<20;i++)
    {

      iota = newiota;

      double fa = themodeli(iota,(void *)&tmi,x,u,xh,uh,Ui,xn,idxi);
      double fad = themodelid(iota,(void *)&tmi,x,u,xh,uh,Ui,xn,idxi);

      //printf("%f %f %f\n",fa,fad,eps);

      if (fabs(fad) < eps)
	break;
      else
	newiota = iota - 1e-19 * fad;

      printf("%d %f %f %f %f\n",i,newiota,iota,fa,fad);


    }

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
	     VEC *v

	     )

{

  double inter = 27*bparams->alpha1*pow(gparams->kappa,2.)*gparams->omega+54*bparams->alpha2*pow(gparams->kappa,2.)*pow(gparams->omega,2.);
  double inter2 = pow(sqrt(pow(inter,2.) + 4*pow(-3*bparams->alpha1*gparams->kappa*gparams->omega - 3*pow(gparams->kappa,2.),3.)) + inter,1/3.);
  double Z = inter2 / (3*pow(2.,1/3.)) - pow(2.,1/3.) * (-3*bparams->alpha1*gparams->kappa*gparams->omega - 3*pow(gparams->kappa,2.)) / (3*inter2); 
  double ubar = (Z - dparams->beta - gparams->kappa) / dparams->gamma; 
  double vbar = (gparams->kappa*gparams->omega*ubar) / (dparams->beta+dparams->gamma*ubar+gparams->kappa);
  double wbar = (2*gparams->kappa*gparams->omega*vbar) / (dparams->beta+dparams->gamma*ubar+2*gparams->kappa);  

  for (int j=0;j<v->dim;j++) 
    v->ve[j] = (bparams->alpha1*vbar+bparams->alpha2*wbar)*pow(gparams->omega-v->ve[j],(dparams->beta+dparams->gamma*ubar)/gparams->kappa-1) / (gparams->kappa*pow(gparams->omega,(dparams->beta+dparams->gamma*ubar)/gparams->kappa));

  return v;
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

		 double iota,
		 void *stuff,
		 MAT *x,
		 MAT *u,
		 MAT *xh,
		 MAT *uh,
		 VEC *Qi,
		 MAT *xn,
		 IVEC *idxi

		 )
{
  struct TMI tmi;

  tmi = * (struct TMI *) stuff;
  VEC *xn_tmp; VEC *un; VEC *xh_tmp; VEC *uh_tmp;
  int J1 = x->n+1;
  xn_tmp = v_get(J1); un = v_get(J1);
  xh_tmp = v_get(J1); uh_tmp = v_get(J1);

  VEC *tmp1; VEC *tmp2;
  tmp1 = v_get(x->n);
  tmp2 = v_get(x->n);

  for (int j=1;j<x->n;j++) 
    x->me[0][j] = h*j;
 
  set_row(u,0,initial((void *)&tmi.gp,(void *)&tmi.bp,(void *)&tmi.dp,get_row(x,0,tmp1)));
 
  Qi->ve[0] = Q(get_row(x,0,tmp1),get_row(u,0,tmp2));

  for (int i=1;i<x->m;i++)
    {

      xhstep(&tmi.gp,get_row(x,i-1,tmp1),xh_tmp);
      uhstep(get_row(x,i-1,tmp1),xh_tmp,get_row(u,i-1,tmp2),uh_tmp,&tmi.bp,&tmi.gp,&tmi.dp,Qi->ve[i-1],k*(i-(tmi.I+1)));

      double Qh = Q(xh_tmp,uh_tmp);

      set_row(xh,i,xh_tmp);
      set_row(uh,i,uh_tmp);

      xstep(&tmi.gp,get_row(x,i-1,tmp1),xh_tmp,xn_tmp);
      set_row(xn,i-1,xn_tmp);
      ustep(xh_tmp,get_row(u,i-1,tmp1),xn_tmp,un,&tmi.bp,&tmi.gp,&tmi.dp,Qi->ve[i-1],k*(i-(tmi.I+1)+.5));

      Qi->ve[i] = Q(xn_tmp,un);

      int idx = idxselect(tmi.gp.omega,xn_tmp);
      idxi->ive[i-1] = idx;

      set_row(x,i,idxremove(xn_tmp,tmp1,idx));
      set_row(u,i,idxremove(un,tmp1,idx));

    }

  VEC *chattmp = v_get(x->n);
  VEC *chat = v_get(x->m); 
  VEC *ttmp = v_get(x->m);

  for (int i=tmi.I+1;i<x->m;i++) {
    for (int j=0;j<x->n;j++)
      chattmp->ve[j] = exp(-pow(x->me[i][j]-phi*iota1,2.)/(2*iota2*pow(phi,2.)))*w(x->me[i][j])*iota*e(k*(i-tmi.I+1))*u->me[i][j];

    chat->ve[i] = pow(Q(get_row(x,i,VNULL),chattmp)-c(k*(i-(tmi.I+1))),2.0);
    ttmp->ve[i] = k*(i-(tmi.I+1));

  }

  double rt = Q(ttmp,chat);

  return rt;

}

double themodelid(

		 double iota,
		 void *stuff,
		 MAT *x,
		 MAT *u,
		 MAT *xh,
		 MAT *uh,
		 VEC *Ui,
		 MAT *xn,
		 IVEC *idxi

		 )
{
  struct TMI tmi;

  tmi = * (struct TMI *) stuff;

  tmi.dp.iota = iota;

  MAT *p;

  int LI = 2*tmi.I + 1;
  int J = tmi.J+1;

  p = m_get(LI,J);

  VEC *ph; VEC *pn;
  int J1 = x->n+1;
  ph = v_get(J1);
  pn = v_get(J1);

  VEC *tmp1; VEC *tmp2; VEC *tmp3; VEC *tmp4; VEC *tmp5; VEC *tmp6;
  tmp1 = v_get(x->n);
  tmp2 = v_get(J1);
  tmp3 = v_get(x->n);
  tmp4 = v_get(x->n);
  tmp5 = v_get(J1);
  tmp6 = v_get(J1);

  VEC *Pi;
  Pi = v_get(x->m);
  
  Pi->ve[0] = Q(get_row(x,0,tmp1),get_row(p,0,tmp3));

  for (int i=1;i<x->m;i++)
    {

      phstep(get_row(x,i-1,tmp1),get_row(xh,i-1,tmp2),get_row(p,i-1,tmp3),ph,get_row(u,i-1,tmp4),&tmi.bp,&tmi.gp,&tmi.dp,Ui->ve[i-1],Pi->ve[i-1],k*(i-(tmi.I+1)));

      double Uh = Q(get_row(xh,i,tmp2),get_row(uh,i,tmp5));
      double Ph = Q(get_row(xh,i,tmp2),ph);

      //printf("%f ",Ph);

      pstep(get_row(xh,i,tmp2),get_row(p,i-1,tmp1),get_row(xn,i-1,tmp6),pn,get_row(uh,i,tmp5),&tmi.bp,&tmi.gp,&tmi.dp,Uh,Ph,k*(i-(tmi.I+1)+.5));

      int idx = idxi->ive[i-1];
      set_row(p,i,idxremove(pn,tmp2,idx));

    }

  VEC *chattmp = v_get(x->n);
  VEC *chattmp2 = v_get(x->n);
  VEC *chat = v_get(x->m); 
  VEC *ttmp = v_get(x->m);

  for (int i=tmi.I+1;i<x->m;i++) {
    for (int j=0;j<x->n;j++)
      {
	chattmp->ve[j] = s(x->me[i][j])*w(x->me[i][j])*e(k*(i-tmi.I+1))*u->me[i][j];
	chattmp2->ve[j] = s(x->me[i][j])*w(x->me[i][j])*e(k*(i-tmi.I+1))*u->me[i][j] + s(x->me[i][j])*w(x->me[i][j])*iota*e(k*(i-tmi.I+1))*p->me[i][j];
      }

    chat->ve[i] = 2*(Q(get_row(x,i,VNULL),chattmp)-c(k*(i-(tmi.I+1))))*(Q(get_row(x,i,VNULL),chattmp2));
    ttmp->ve[i] = k*(i-(tmi.I+1));

  }

  double rt = Q(ttmp,chat);

  /*
  M_FREE(x);
  M_FREE(u);
  M_FREE(xh);
  M_FREE(uh);
  V_FREE(Ui);
  M_FREE(xn);
  IV_FREE(idxi);
  */

  return rt;

}


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
      int idx = floor((t + (k/4))/(k/2));
      return eff->ve[idx];
    }
}

double c(double t)
{
  if (t<0)
    return cat->ve[0];
  else
    {
      int idx = floor((t + (k/4))/(k/2));
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

void phstep(

	    VEC *x,
	    VEC *xh,
	    VEC *p,
	    VEC *ph,
	    VEC *u,
	    struct BP *bpar,
	    struct GP *gpar,
	    struct DP *dpar,
	    double Ui,
	    double Pi,
	    double t

	    )

{


  for (int i=1;i<ph->dim;i++)
    ph->ve[i]=p->ve[i-1]*exp(-(k/2)*zstar(dpar->beta,dpar->gamma,gpar->kappa,dpar->iota,t,x->ve[i-1],Ui)) - wstar(dpar->gamma,Pi,u->ve[i-1],x->ve[i-1],t);

  Q2(bpar,gpar,xh,ph);

}

double wstar(

	     double gamma,
	     double P,
	     double u,
	     double x,
	     double t

	     )
{
  return (s(x)*e(t) + gamma*P)*u;
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

void pstep(

	   VEC *xh,
	   VEC *p,
	   VEC *xn,
	   VEC *pn,
	   VEC *uh,
	   struct BP *bpar,
	   struct GP *gpar,
	   struct DP *dpar,
	   double Uh,
	   double Ph,
	   double t
	   
	   )
{

  for (int i=1;i<pn->dim;i++)
    pn->ve[i]=p->ve[i-1]*exp(-(k/2)*zstar(dpar->beta,dpar->gamma,gpar->kappa,dpar->iota,t,xh->ve[i-1],Uh)) - wstar(dpar->gamma,Ph,uh->ve[i-1],xh->ve[i-1],t);

  Q2(bpar,gpar,xn,pn);

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
 
  set_row(u,0,initial(gparams,bparams,dparams,get_row(x,0,tmp1)));
 
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
