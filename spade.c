#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <GL/freeglut.h>
#include "meschach/matrix.h"
#include "meschach/matrix2.h"

#include "trackball.h"

int N;  
int M;
int I;
int S;
int P;
double *z;
double *d1;
double *d2;
double *effort,*effort2;

int N2;
int S2;
int I2;
double *z2;

int Nparams1;
int Nparamsb;

double *x,*x2;
double *N2x;

VEC *lx,*lx2;

VEC *lb,*lb2;
MAT *A, *LU, *A2, *LU2;
PERM *pivot;

typedef struct {

  int nBlocks;
  double *effort;
  int nPoints;

} GSL_reg;

const gsl_multimin_fminimizer_type *T;

gsl_multimin_fminimizer *s1,*s2,*sb,*sb2;
gsl_vector *ss1, *ss2, *ssb, *q1, *q2, *qb,*qb2,*ssb2;
gsl_multimin_function f1,f2,fb,fb2;

GSL_reg reg;
    
size_t iter;
int status;
double size;

int ci;

double maxcurv;
int mcidx;

double thingy,thingy2;

int pt,new_pt;

int WindWidth, WindHeight;
int window;

double curvaturei(

		  double x

		 )

{
  
  return fabs(2*lx->ve[ci*4+2] + 6*lx->ve[ci*4+3] * x) / pow( 1 + pow(lx->ve[ci*4+1] + 2*lx->ve[ci*4+2]*x + 3*lx->ve[ci*4+3] * pow(x,2.0),2.0), 2.0/3.0 );

}

double ncurvaturei(

		  double x

		 )

{
  
  return -fabs(2*lx->ve[ci*4+2] + 6*lx->ve[ci*4+3] * x) / pow( 1 + pow(lx->ve[ci*4+1] + 2*lx->ve[ci*4+2]*x + 3*lx->ve[ci*4+3] * pow(x,2.0),2.0), 2.0/3.0 );

}


double curvaturei2(

		  double x

		 )

{
  
  return fabs(2*lx2->ve[ci*4+2] + 6*lx2->ve[ci*4+3] * x) / pow( 1 + pow(lx2->ve[ci*4+1] + 2*lx2->ve[ci*4+2]*x + 3*lx2->ve[ci*4+3] * pow(x,2.0),2.0), 2.0/3.0 );

}

double ncurvaturei2(

		  double x

		 )

{
  
  return -fabs(2*lx2->ve[ci*4+2] + 6*lx2->ve[ci*4+3] * x) / pow( 1 + pow(lx2->ve[ci*4+1] + 2*lx2->ve[ci*4+2]*x + 3*lx2->ve[ci*4+3] * pow(x,2.0),2.0), 2.0/3.0 );

}

double curvatureboth2(

		      gsl_vector *p,
		      void *param

		      )
{


  double d0 = gsl_vector_get(p,0);
  double d1 = gsl_vector_get(p,1);
  double dd1 = gsl_vector_get(p,2);

  lb2->ve[N2+S2-1+S2-1+S2-1] = d0;
  lb2->ve[N2+S2-1+S2-1+S2-1+2] = d1;
  lb2->ve[N2+S2-1+S2-1+S2-1+3] = dd1;

  for (int i=1;i<S2;i++)
    x2[i] = gsl_vector_get(p,i-1+3);

  A2 = m_get(4*S2,4*S2);

  double * fakeeffort = (double *) calloc(N2,sizeof(double));

  int j=0;
  for (int i=0;i<N2;i++)
    {

      if (x2[j] > N2x[i])
	{

          double tmpfakeeffort=0;
	  for (int l=0;l<100;l++)
	    {

	      double xx = l*(x2[j]-N2x[i])/100 + N2x[i];
	      int oldi = S;
	      while (x[oldi]>xx) { oldi--; }

	      tmpfakeeffort += ((x2[j]-N2x[i])/100) * (lx->ve[4*oldi] + lx->ve[4*oldi+1]*xx + lx->ve[4*oldi+2]*pow(xx,2.0) + lx->ve[4*oldi+3]*pow(xx,3.0));

	    }
	  fakeeffort[i] += tmpfakeeffort; ///(N2x[i]-x2[j]);
	}

      int counter=0;
      while ( x2[j+1] < N2x[i+1] ) { j++; counter++; }

      int k=0;
      
      for (k=0;k<counter;k++)
	{

          double tmpfakeeffort=0;
	  for (int l=0;l<100;l++)
	    {

	      double xx = l*(x2[j-counter+k+1]-x2[j-counter+k])/100 + x2[j-counter+k];
	      int oldi = S;
	      while (x[oldi]>xx) { oldi--; }

	      tmpfakeeffort += (x2[j-counter+k+1]-x2[j-counter+k])/100 * (lx->ve[4*oldi] + lx->ve[4*oldi+1]*xx + lx->ve[4*oldi+2]*pow(xx,2.0) + lx->ve[4*oldi+3]*pow(xx,3.0));

	    }
	  fakeeffort[i] += tmpfakeeffort; ///(x2[j-counter+k+1]-x2[j-counter+k]);
	}

      double tmpfakeeffort=0;
      for (int l=0;l<100;l++)
	{

	  double xx = l*(N2x[i+1]-x2[j-counter+k])/100 + x2[j-counter+k];
	  int oldi = S;
	  while (x[oldi]>xx) { oldi--; }

	  tmpfakeeffort += (N2x[i+1]-x2[j-counter+k])/100 * (lx->ve[4*oldi] + lx->ve[4*oldi+1]*xx + lx->ve[4*oldi+2]*pow(xx,2.0) + lx->ve[4*oldi+3]*pow(xx,3.0));

	}

      fakeeffort[i] += tmpfakeeffort; //x[i+1]-x2[j-counter+k]);

      j++; 
    }

  j=0;

  for (int i=0;i<N2;i++)
    {

      if (x2[j] > N2x[i])
	{
	  A2->me[i][(j-1)*4 + 0] += x2[j] - N2x[i];
	  A2->me[i][(j-1)*4 + 1] += .5*( pow(x2[j],2.0) - pow(N2x[i],2.0) );
	  A2->me[i][(j-1)*4 + 2] += (1.0/3) * ( pow(x2[j],3.0) - pow(N2x[i],3.0));
	  A2->me[i][(j-1)*4 + 3] += (1.0/4) * ( pow(x2[j],4.0) - pow(N2x[i],4.0));
	}
      
      int counter=0;
	        
      while ( x2[j+1] < N2x[i+1] ) { j++; counter++; }

      int k=0;
      
      for (k=0;k<counter;k++)
	{
	
	  A2->me[i][(j-1)*4 + 0] += x2[j-counter+k+1] - x2[j-counter+k];
	  A2->me[i][(j-1)*4 + 1] += .5*( pow(x2[j-counter+k+1],2.0) - pow(x2[j-counter+k],2.0) );
	  A2->me[i][(j-1)*4 + 2] += (1.0/3) * ( pow(x2[j-counter+k+1],3.0) - pow(x2[j-counter+k],3.0));
	  A2->me[i][(j-1)*4 + 3] += (1.0/4) * ( pow(x2[j-counter+k+1],4.0) - pow(x2[j-counter+k],4.0));

	}

      A2->me[i][j*4 + 0] += N2x[i+1]-x2[j-counter+k];
      A2->me[i][j*4 + 1] += .5*( pow(N2x[i+1],2.0) - pow(x2[j-counter+k],2.0) );
      A2->me[i][j*4 + 2] += (1.0/3) * ( pow(N2x[i+1],3.0) - pow(x2[j-counter+k],3.0));
      A2->me[i][j*4 + 3] += (1.0/4) * ( pow(N2x[i+1],4.0) - pow(x2[j-counter+k],4.0));

      j++; 
      	  
    }

  for (int j=0;j<S2-1;j++)
    {

      A2->me[N2+j][j*4 + 0] = 1;
      A2->me[N2+j][j*4 + 1] = x2[j+1];
      A2->me[N2+j][j*4 + 2] = pow(x2[j+1],2.0);
      A2->me[N2+j][j*4 + 3] = pow(x2[j+1],3.0);
      A2->me[N2+j][j*4 + 4] = -1;
      A2->me[N2+j][j*4 + 5] = -(x2[j+1]);
      A2->me[N2+j][j*4 + 6] = -pow(x2[j+1],2.0);
      A2->me[N2+j][j*4 + 7] = -pow(x2[j+1],3.0);
      
    }

  for (int j=0;j<S2-1;j++)
    {

      A2->me[N2+S2-1+j][j*4 + 1] = 1;
      A2->me[N2+S2-1+j][j*4 + 2] = 2*(x2[j+1]);
      A2->me[N2+S2-1+j][j*4 + 3] = 3*pow(x2[j+1],2.0);
      A2->me[N2+S2-1+j][j*4 + 5] = -1;
      A2->me[N2+S2-1+j][j*4 + 6] = -2*(x2[j+1]);
      A2->me[N2+S2-1+j][j*4 + 7] = -3*pow(x2[j+1],2.0);
      
    }

  for (int j=0;j<S2-1;j++)
    { 
      A2->me[N2+S2-1+S2-1+j][j*4 + 2] = 2;
      A2->me[N2+S2-1+S2-1+j][j*4 + 3] = 6*(x2[j+1]);
      A2->me[N2+S2-1+S2-1+j][j*4 + 6] = -2;
      A2->me[N2+S2-1+S2-1+j][j*4 + 7] = -6*(x2[j+1]);
    }
  
  A2->me[N2+S2-1+S2-1+S2-1][2] = 2; // second derivative at t=0

  A2->me[N2+S2-1+S2-1+S2-1+1][0] = 1;  // value at t=0;
  
  A2->me[N2+S2-1+S2-1+S2-1+2][N2+S2-1+S2-1+S2-1+3-2] = 1;  // first derivative at t=N2
  A2->me[N2+S2-1+S2-1+S2-1+2][N2+S2-1+S2-1+S2-1+3-1] = 2*N2x[N2]; // first derivative at t=N2
  A2->me[N2+S2-1+S2-1+S2-1+2][N2+S2-1+S2-1+S2-1+3] = 3*pow(N2x[N2],2.0);  // first derivative at t=N2

  A2->me[N2+S2-1+S2-1+S2-1+3][N2+S2-1+S2-1+S2-1+3-1] = 2;  // second derivative at t=N2
  A2->me[N2+S2-1+S2-1+S2-1+3][N2+S2-1+S2-1+S2-1+3] = 6*N2x[N2];   

  double pp = .9;
  for (int i=0;i<N2;i++)
    lb2->ve[i] = pp*fakeeffort[i] + (1-pp)*effort2[i];

  LU2 = m_get(4*S2,4*S2);
  LU2 = m_copy(A2,LU2);
  pivot = px_get(A2->m);
  LUfactor(LU2,pivot);

  lx2 = LUsolve(LU2,pivot,lb2,VNULL);

  double f_best = 0;

  for (ci=0;ci<S2;ci++)
    {
      
      double x_lower = x2[ci];
      double x_upper = x2[ci+1];
      double x_mid = (x2[ci+1] + x2[ci] )/2;
  
      double f_lower = curvaturei2(x_lower);
      double f_upper = curvaturei2(x_upper);
  
      double f_mid = curvaturei2(x_mid);

      int found = 0;
  
      do
	{
          
	  if (f_mid < f_lower && f_mid < f_upper)
	    {

	      if (f_lower > f_upper)
		f_mid = f_lower;
	      else
		f_mid = f_upper;

	      found = 1;
	    }
	  else
	    {

	      if (f_mid > f_lower && f_mid > f_upper)
		{

		  int status;
		  int iter=0, max_iter =100;
		  const gsl_min_fminimizer_type *T;
		  gsl_min_fminimizer *s;
	  
		  gsl_function F;

		  F.function = &ncurvaturei2;
		  F.params = 0;

		  T = gsl_min_fminimizer_brent;
		  s = gsl_min_fminimizer_alloc (T);
		  gsl_min_fminimizer_set (s, &F, x_mid, x_lower, x_upper);

		  printf ("using %s method\n",gsl_min_fminimizer_name (s));

		  printf ("%5s [%9s, %9s] %9s %9s\n","iter", "lower", "upper", "min", "err(est)");

		  printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lower, x_upper, x_mid, x_upper - x_lower);
	  
		  do
		    {
   
		      iter++;
		      status = gsl_min_fminimizer_iterate (s);

		      x_mid = gsl_min_fminimizer_x_minimum (s);
		      x_lower = gsl_min_fminimizer_x_lower (s);
		      x_upper = gsl_min_fminimizer_x_upper (s);

		      status = gsl_min_test_interval (x_lower, x_upper, 0.001, 0.0);

		      if (status == GSL_SUCCESS)
			printf ("Converged:\n");

		      printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lower, x_upper, x_mid, x_upper - x_lower);
	      
		    }
		  while (status == GSL_CONTINUE && iter < max_iter);

		  found = 1;

		  f_mid = curvaturei2(x_mid);
	  
		}
	      else
		{
	  
		  if ( f_mid > f_lower)
		    {

		      x_lower = x_mid;
		      f_lower = f_mid;
	      
		      x_mid = (x_lower + x_upper) / 2;
		      f_mid = curvaturei2(x_mid);

		      if ( (x_mid - x_lower) < 0.001)
			found = 1;
	      
		    }
		  else
		    {

		      x_upper = x_mid;
		      f_upper = f_mid;
	      
		      x_mid = (x_lower + x_upper) / 2;
		      f_mid = curvaturei2(x_mid);

		      if ( (x_mid - x_lower) < 0.001)
			found = 1;
		    }
		}
	    }
	}
      
      while (found ==0);

      if (f_mid > f_best)
	f_best = f_mid;

    }


  double f_pen = 0;
  for (int i=0;i<S2;i++)
    {
      for (int j=0;j<100;j++)
	{
	  double xx = x2[i]+j*(x2[i+1]-x2[i])/100.0;
	  double yy = lx2->ve[i*4] + lx2->ve[i*4+1]*xx + lx2->ve[i*4+2]*pow(xx,2.0) + lx2->ve[i*4+3] * pow(xx,3.0);

	  if (yy < 0)
	    f_pen += pow(10*yy,2.0);

	}
    }

  

  double  f_pen2 = 0;

  for (int i=1;i<=S2;i++)  
    f_pen2 += pow(1.0/(20*(x2[i]-x2[i-1])),4.0);

  //if (x[S-1] > N)
  //f_pen2 = 100*pow(N-x[S-1],2.0);

  printf("fpen: %lf fen2: %lf\n",f_pen,f_pen2);

  return f_best + f_pen + f_pen2;

}

double curvatureboth(

		 gsl_vector *p,
		 void *param

		 )
{

  double d0 = gsl_vector_get(p,0);
  double d1 = gsl_vector_get(p,1);
  double dd1 = gsl_vector_get(p,2);

  lb->ve[N+S-1+S-1+S-1] = d0;
  lb->ve[N+S-1+S-1+S-1+2] = d1;
  lb->ve[N+S-1+S-1+S-1+3] = dd1;

  for (int i=1;i<S;i++)
    x[i] = gsl_vector_get(p,i-1+3);

  MAT *nA = m_get(4*S,4*S);
  
  int j=0;

  for (int i=0;i<N;i++)
    {

      if (x[j] > i)
	{

	  if (x[j] > i+1)
	    {
	      nA->me[i][(j-1)*4 + 0] += (i+1) - i;
	      nA->me[i][(j-1)*4 + 1] += .5*( pow((i+1),2.0) - pow(i,2.0) );
	      nA->me[i][(j-1)*4 + 2] += (1.0/3) * ( pow((i+1),3.0) - pow(i,3.0));
	      nA->me[i][(j-1)*4 + 3] += (1.0/4) * ( pow((i+1),4.0) - pow(i,4.0));
	      
	    }
	  else
	    {
	    
	      nA->me[i][(j-1)*4 + 0] += x[j] - i;
	      nA->me[i][(j-1)*4 + 1] += .5*( pow(x[j],2.0) - pow(i,2.0) );
	      nA->me[i][(j-1)*4 + 2] += (1.0/3) * ( pow(x[j],3.0) - pow(i,3.0));
	      nA->me[i][(j-1)*4 + 3] += (1.0/4) * ( pow(x[j],4.0) - pow(i,4.0));
	    }
	}
      
      int counter=0;
	        
      while (j<S && x[j+1] < i+1 )
	{
	  j++;
	  counter++;
	}

      /*
      if (j==S)
	{
	  printf("problem?\n");
	  break;
	}
      */
      
      int k=0;
      
      for (k=0;k<counter;k++)
	{
	
	  nA->me[i][(j-1)*4 + 0] += x[j-counter+k+1] - x[j-counter+k];
	  nA->me[i][(j-1)*4 + 1] += .5*( pow(x[j-counter+k+1],2.0) - pow(x[j-counter+k],2.0) );
	  nA->me[i][(j-1)*4 + 2] += (1.0/3) * ( pow(x[j-counter+k+1],3.0) - pow(x[j-counter+k],3.0));
	  nA->me[i][(j-1)*4 + 3] += (1.0/4) * ( pow(x[j-counter+k+1],4.0) - pow(x[j-counter+k],4.0));
	}

      if (x[j] < i+1)
	{
	  nA->me[i][j*4 + 0] += i+1-x[j-counter+k];
	  nA->me[i][j*4 + 1] += .5*( pow(i+1,2.0) - pow(x[j-counter+k],2.0) );
	  nA->me[i][j*4 + 2] += (1.0/3) * ( pow(i+1,3.0) - pow(x[j-counter+k],3.0));
	  nA->me[i][j*4 + 3] += (1.0/4) * ( pow(i+1,4.0) - pow(x[j-counter+k],4.0));	  	  
	}

      j++;      
      	  
    }

  for (int j=0;j<S-1;j++)
    {

      nA->me[N+j][j*4 + 0] = 1;
      nA->me[N+j][j*4 + 1] = x[j+1];
      nA->me[N+j][j*4 + 2] = pow(x[j+1],2.0);
      nA->me[N+j][j*4 + 3] = pow(x[j+1],3.0);
      nA->me[N+j][j*4 + 4] = -1;
      nA->me[N+j][j*4 + 5] = -(x[j+1]);
      nA->me[N+j][j*4 + 6] = -pow(x[j+1],2.0);
      nA->me[N+j][j*4 + 7] = -pow(x[j+1],3.0);
      
    }

  for (int j=0;j<S-1;j++)
    {

      nA->me[N+S-1+j][j*4 + 1] = 1;
      nA->me[N+S-1+j][j*4 + 2] = 2*(x[j+1]);
      nA->me[N+S-1+j][j*4 + 3] = 3*pow(x[j+1],2.0);
      nA->me[N+S-1+j][j*4 + 5] = -1;
      nA->me[N+S-1+j][j*4 + 6] = -2*(x[j+1]);
      nA->me[N+S-1+j][j*4 + 7] = -3*pow(x[j+1],2.0);
      
    }

  for (int j=0;j<S-1;j++)
    { 
      nA->me[N+S-1+S-1+j][j*4 + 2] = 2;
      nA->me[N+S-1+S-1+j][j*4 + 3] = 6*(x[j+1]);
      nA->me[N+S-1+S-1+j][j*4 + 6] = -2;
      nA->me[N+S-1+S-1+j][j*4 + 7] = -6*(x[j+1]);
    }

  //nA->me[N+S-1+S-1+S-1+1][1] = 1;  // first derivative at t=0

  nA->me[N+S-1+S-1+S-1+1][0] = 1;  // value at t=0;
  
  nA->me[N+S-1+S-1+S-1][2] = 2; // second derivative at t=0
  
  nA->me[N+S-1+S-1+S-1+2][N+S-1+S-1+S-1+3-2] = 1;  // first derivative at t=N
  nA->me[N+S-1+S-1+S-1+2][N+S-1+S-1+S-1+3-1] = 2*N; // first derivative at t=N
  nA->me[N+S-1+S-1+S-1+2][N+S-1+S-1+S-1+3] = 3*pow(N,2.0);  // first derivative at t=N

  nA->me[N+S-1+S-1+S-1+3][N+S-1+S-1+S-1+3-1] = 2;  // second derivative at t=N
  nA->me[N+S-1+S-1+S-1+3][N+S-1+S-1+S-1+3] = 6*N;     

  LU = m_get(4*S,4*S);
  LU = m_copy(nA,LU);
  pivot = px_get(nA->m);
  LUfactor(LU,pivot);

  lx = LUsolve(LU,pivot,lb,VNULL);

  double f_best = 0;

  for (ci=0;ci<S;ci++)
    {
      
      double x_lower = x[ci];
      double x_upper = x[ci+1];
      double x_mid = (x[ci+1] + x[ci] )/2;
  
      double f_lower = curvaturei(x_lower);
      double f_upper = curvaturei(x_upper);
  
      double f_mid = curvaturei(x_mid);

      int found = 0;
  
      do
	{
          
	  if (f_mid < f_lower && f_mid < f_upper)
	    {

	      if (f_lower > f_upper)
		f_mid = f_lower;
	      else
		f_mid = f_upper;

	      found = 1;
	    }
	  else
	    {

	      if (f_mid > f_lower && f_mid > f_upper)
		{

		  int status;
		  int iter=0, max_iter =100;
		  const gsl_min_fminimizer_type *T;
		  gsl_min_fminimizer *s;
	  
		  gsl_function F;

		  F.function = &ncurvaturei;
		  F.params = 0;

		  T = gsl_min_fminimizer_brent;
		  s = gsl_min_fminimizer_alloc (T);
		  gsl_min_fminimizer_set (s, &F, x_mid, x_lower, x_upper);

		  //printf ("using %s method\n",gsl_min_fminimizer_name (s));

		  //printf ("%5s [%9s, %9s] %9s %9s\n","iter", "lower", "upper", "min", "err(est)");

		  //printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lower, x_upper, x_mid, x_upper - x_lower);
	  
		  do
		    {
   
		      iter++;
		      status = gsl_min_fminimizer_iterate (s);

		      x_mid = gsl_min_fminimizer_x_minimum (s);
		      x_lower = gsl_min_fminimizer_x_lower (s);
		      x_upper = gsl_min_fminimizer_x_upper (s);

		      status = gsl_min_test_interval (x_lower, x_upper, 0.001, 0.0);

		      //if (status == GSL_SUCCESS)
			//printf ("Converged:\n");

			//printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lower, x_upper, x_mid, x_upper - x_lower);
	      
		    }
		  while (status == GSL_CONTINUE && iter < max_iter);

		  found = 1;

		  f_mid = curvaturei(x_mid);
	  
		}
	      else
		{
	  
		  if ( f_mid > f_lower)
		    {

		      x_lower = x_mid;
		      f_lower = f_mid;
	      
		      x_mid = (x_lower + x_upper) / 2;
		      f_mid = curvaturei(x_mid);

		      if ( (x_mid - x_lower) < 0.001)
			found = 1;
	      
		    }
		  else
		    {

		      x_upper = x_mid;
		      f_upper = f_mid;
	      
		      x_mid = (x_lower + x_upper) / 2;
		      f_mid = curvaturei(x_mid);

		      if ( (x_mid - x_lower) < 0.001)
			found = 1;
		    }
		}
	    }
	}
      
      while (found ==0);

      if (f_mid > f_best)
	f_best = f_mid;

    }


  double f_pen = 0;
  for (int i=0;i<S;i++)
    {
      for (int j=0;j<100;j++)
	{
	  double xx = x[i]+j*(x[i+1]-x[i])/100.0;
	  double yy = lx->ve[i*4] + lx->ve[i*4+1]*xx + lx->ve[i*4+2]*pow(xx,2.0) + lx->ve[i*4+3] * pow(xx,3.0);

	  if (yy < 0)
	    f_pen += pow(10*yy,2.0);

	}
    }

  

  double  f_pen2 = 0;

  for (int i=1;i<=S;i++)  
    f_pen2 += pow(1.0/(20*(x[i]-x[i-1])),4.0);

  //if (x[S-1] > N)
  //f_pen2 = 100*pow(N-x[S-1],2.0);

  //printf("fpen: %lf fen2: %lf\n",f_pen,f_pen2);

  return f_best + f_pen + f_pen2;

}

double curvaturex(

		 gsl_vector *p,
		 void *param

		 )
{

  //x[1] = gsl_vector_get(p,0);
  //x[2] = gsl_vector_get(p,1);

  for (int i=1;i<S;i++)
    x[i] = gsl_vector_get(p,i-1);

  MAT *nA = m_get(4*S,4*S);
  //A = m_get(4*S,4*S);
  
  int j=0;

  for (int i=0;i<N;i++)
    {

      if (x[j] > i)
	{

	  if (x[j] > i+1)
	    {
	      nA->me[i][(j-1)*4 + 0] += (i+1) - i;
	      nA->me[i][(j-1)*4 + 1] += .5*( pow((i+1),2.0) - pow(i,2.0) );
	      nA->me[i][(j-1)*4 + 2] += (1.0/3) * ( pow((i+1),3.0) - pow(i,3.0));
	      nA->me[i][(j-1)*4 + 3] += (1.0/4) * ( pow((i+1),4.0) - pow(i,4.0));
	      
	    }
	  else
	    {
	    
	      nA->me[i][(j-1)*4 + 0] += x[j] - i;
	      nA->me[i][(j-1)*4 + 1] += .5*( pow(x[j],2.0) - pow(i,2.0) );
	      nA->me[i][(j-1)*4 + 2] += (1.0/3) * ( pow(x[j],3.0) - pow(i,3.0));
	      nA->me[i][(j-1)*4 + 3] += (1.0/4) * ( pow(x[j],4.0) - pow(i,4.0));
	    }
	}
      
      int counter=0;
	        
      while (j<S && x[j+1] < i+1 )
	{
	  j++;
	  counter++;
	}

      /*
      if (j==S)
	{
	  printf("problem?\n");
	  break;
	}
      */
      
      int k=0;
      
      for (k=0;k<counter;k++)
	{
	
	  nA->me[i][(j-1)*4 + 0] += x[j-counter+k+1] - x[j-counter+k];
	  nA->me[i][(j-1)*4 + 1] += .5*( pow(x[j-counter+k+1],2.0) - pow(x[j-counter+k],2.0) );
	  nA->me[i][(j-1)*4 + 2] += (1.0/3) * ( pow(x[j-counter+k+1],3.0) - pow(x[j-counter+k],3.0));
	  nA->me[i][(j-1)*4 + 3] += (1.0/4) * ( pow(x[j-counter+k+1],4.0) - pow(x[j-counter+k],4.0));
	}

      /*      if (x[j+1] > i+1)
	{

	  nA->me[i][j*4 + 0] += i+1-x[j-counter+k];
	  nA->me[i][j*4 + 1] += .5*( pow(i+1,2.0) - pow(x[j-counter+k],2.0) );
	  nA->me[i][j*4 + 2] += (1.0/3) * ( pow(i+1,3.0) - pow(x[j-counter+k],3.0));
	  nA->me[i][j*4 + 3] += (1.0/4) * ( pow(i+1,4.0) - pow(x[j-counter+k],4.0));
	  
	}
      else
      {*/

      if (x[j] < i+1)
	{
	  nA->me[i][j*4 + 0] += i+1-x[j-counter+k];
	  nA->me[i][j*4 + 1] += .5*( pow(i+1,2.0) - pow(x[j-counter+k],2.0) );
	  nA->me[i][j*4 + 2] += (1.0/3) * ( pow(i+1,3.0) - pow(x[j-counter+k],3.0));
	  nA->me[i][j*4 + 3] += (1.0/4) * ( pow(i+1,4.0) - pow(x[j-counter+k],4.0));	  	  
	}

      j++;      
      	  
    }

  for (int j=0;j<S-1;j++)
    {

      nA->me[N+j][j*4 + 0] = 1;
      nA->me[N+j][j*4 + 1] = x[j+1];
      nA->me[N+j][j*4 + 2] = pow(x[j+1],2.0);
      nA->me[N+j][j*4 + 3] = pow(x[j+1],3.0);
      nA->me[N+j][j*4 + 4] = -1;
      nA->me[N+j][j*4 + 5] = -(x[j+1]);
      nA->me[N+j][j*4 + 6] = -pow(x[j+1],2.0);
      nA->me[N+j][j*4 + 7] = -pow(x[j+1],3.0);
      
    }

  for (int j=0;j<S-1;j++)
    {

      nA->me[N+S-1+j][j*4 + 1] = 1;
      nA->me[N+S-1+j][j*4 + 2] = 2*(x[j+1]);
      nA->me[N+S-1+j][j*4 + 3] = 3*pow(x[j+1],2.0);
      nA->me[N+S-1+j][j*4 + 5] = -1;
      nA->me[N+S-1+j][j*4 + 6] = -2*(x[j+1]);
      nA->me[N+S-1+j][j*4 + 7] = -3*pow(x[j+1],2.0);
      
    }

  for (int j=0;j<S-1;j++)
    { 
      nA->me[N+S-1+S-1+j][j*4 + 2] = 2;
      nA->me[N+S-1+S-1+j][j*4 + 3] = 6*(x[j+1]);
      nA->me[N+S-1+S-1+j][j*4 + 6] = -2;
      nA->me[N+S-1+S-1+j][j*4 + 7] = -6*(x[j+1]);
    }

  //nA->me[N+S-1+S-1+S-1+1][1] = 1;  // first derivative at t=0

  nA->me[N+S-1+S-1+S-1+1][0] = 1;  // value at t=0;
  
  nA->me[N+S-1+S-1+S-1][2] = 2; // second derivative at t=0
  
  nA->me[N+S-1+S-1+S-1+2][N+S-1+S-1+S-1+3-2] = 1;  // first derivative at t=N
  nA->me[N+S-1+S-1+S-1+2][N+S-1+S-1+S-1+3-1] = 2*N; // first derivative at t=N
  nA->me[N+S-1+S-1+S-1+2][N+S-1+S-1+S-1+3] = 3*pow(N,2.0);  // first derivative at t=N

  nA->me[N+S-1+S-1+S-1+3][N+S-1+S-1+S-1+3-1] = 2;  // second derivative at t=N
  nA->me[N+S-1+S-1+S-1+3][N+S-1+S-1+S-1+3] = 6*N;  
    
  P = S-N;
  for (int i=1;i<P;i++)
    {
      nA->me[N+S-1+S-1+S-1+3+i][i*(S/P)*4] = 1;  // S bigger than 3
      nA->me[N+S-1+S-1+S-1+3+i][i*(S/P)*4+1] = N/(double)P;
      nA->me[N+S-1+S-1+S-1+3+i][i*(S/P)*4+2] = pow(N/(double)P,2.0);
      nA->me[N+S-1+S-1+S-1+3+i][i*(S/P)*4+3] = pow(N/(double)P,3.0);
    }


  LU = m_get(4*S,4*S);
  LU = m_copy(nA,LU);
  pivot = px_get(nA->m);
  LUfactor(LU,pivot);

  lx = LUsolve(LU,pivot,lb,VNULL);

  /*
  for (int i=0;i<S;i++)
    for (int j=0;j<100;j++)
      {
	double xx = x[i]+j*(x[i+1]-x[i])/100.0;	  
	printf("%lf ",lx->ve[i*4] + lx->ve[i*4+1]*xx + lx->ve[i*4+2]*pow(xx,2.0) + lx->ve[i*4+3] * pow(xx,3.0));
      }

  printf("\n");
  */
  
  double f_best = 0;

  for (ci=0;ci<S;ci++)
    {
      
      double x_lower = x[ci];
      double x_upper = x[ci+1];
      double x_mid = (x[ci+1] + x[ci] )/2;
  
      double f_lower = curvaturei(x_lower);
      double f_upper = curvaturei(x_upper);
  
      double f_mid = curvaturei(x_mid);

      int found = 0;
  
      do
	{
          
	  if (f_mid < f_lower && f_mid < f_upper)
	    {

	      if (f_lower > f_upper)
		f_mid = f_lower;
	      else
		f_mid = f_upper;

	      found = 1;
	    }
	  else
	    {

	      if (f_mid > f_lower && f_mid > f_upper)
		{

		  int status;
		  int iter=0, max_iter =100;
		  const gsl_min_fminimizer_type *T;
		  gsl_min_fminimizer *s;
	  
		  gsl_function F;

		  F.function = &ncurvaturei;
		  F.params = 0;

		  T = gsl_min_fminimizer_brent;
		  s = gsl_min_fminimizer_alloc (T);
		  gsl_min_fminimizer_set (s, &F, x_mid, x_lower, x_upper);

		  printf ("using %s method\n",gsl_min_fminimizer_name (s));

		  printf ("%5s [%9s, %9s] %9s %9s\n","iter", "lower", "upper", "min", "err(est)");

		  printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lower, x_upper, x_mid, x_upper - x_lower);
	  
		  do
		    {
   
		      iter++;
		      status = gsl_min_fminimizer_iterate (s);

		      x_mid = gsl_min_fminimizer_x_minimum (s);
		      x_lower = gsl_min_fminimizer_x_lower (s);
		      x_upper = gsl_min_fminimizer_x_upper (s);

		      status = gsl_min_test_interval (x_lower, x_upper, 0.001, 0.0);

		      if (status == GSL_SUCCESS)
			printf ("Converged:\n");

		      printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lower, x_upper, x_mid, x_upper - x_lower);
	      
		    }
		  while (status == GSL_CONTINUE && iter < max_iter);

		  found = 1;

		  f_mid = curvaturei(x_mid);
	  
		}
	      else
		{
	  
		  if ( f_mid > f_lower)
		    {

		      x_lower = x_mid;
		      f_lower = f_mid;
	      
		      x_mid = (x_lower + x_upper) / 2;
		      f_mid = curvaturei(x_mid);

		      if ( (x_mid - x_lower) < 0.001)
			found = 1;
	      
		    }
		  else
		    {

		      x_upper = x_mid;
		      f_upper = f_mid;
	      
		      x_mid = (x_lower + x_upper) / 2;
		      f_mid = curvaturei(x_mid);

		      if ( (x_mid - x_lower) < 0.001)
			found = 1;
		    }
		}
	    }
	}
      
      while (found ==0);

      if (f_mid > f_best)
	f_best = f_mid;

    }

  return f_best;

}

double curvature(

		 gsl_vector *p,
		 void *param

		 )
{

  double d0 = gsl_vector_get(p,0);
  double d1 = gsl_vector_get(p,1);
  double dd1 = gsl_vector_get(p,2);

  //printf("before: ");
  //for (int i=0;i<4*S;i++)
  //  printf(" %lf ",lb->ve[i]);
  //printf("\n");

  lb->ve[N+S-1+S-1+S-1+2] = d1;
  lb->ve[N+S-1+S-1+S-1+3] = dd1;

  for (int i=1;i<P;i++)
    lb->ve[N+S-1+S-1+S-1+3+i] = gsl_vector_get(p,1+i);

  //printf("after: ");
  //for (int i=0;i<4*S;i++)
  //  printf(" %lf ",lb->ve[i]);
  //printf("\n\n");

  lx = LUsolve(LU,pivot,lb,VNULL);

  double f_best = 0;

  for (ci=0;ci<S;ci++)
    {
      
      double x_lower = x[ci];
      double x_upper = x[ci+1];
      double x_mid = (x[ci+1] + x[ci] )/2;
  
      double f_lower = curvaturei(x_lower);
      double f_upper = curvaturei(x_upper);
  
      double f_mid = curvaturei(x_mid);

      int found = 0;
  
      do
	{
          
	  if (f_mid < f_lower && f_mid < f_upper)
	    {

	      if (f_lower > f_upper)
		f_mid = f_lower;
	      else
		f_mid = f_upper;

	      found = 1;
	    }
	  else
	    {

	      if (f_mid > f_lower && f_mid > f_upper)
		{

		  int status;
		  int iter=0, max_iter =100;
		  const gsl_min_fminimizer_type *T;
		  gsl_min_fminimizer *s;
	  
		  gsl_function F;

		  F.function = &ncurvaturei;
		  F.params = 0;

		  T = gsl_min_fminimizer_brent;
		  s = gsl_min_fminimizer_alloc (T);
		  gsl_min_fminimizer_set (s, &F, x_mid, x_lower, x_upper);

		  printf ("using %s method\n",gsl_min_fminimizer_name (s));

		  printf ("%5s [%9s, %9s] %9s %9s\n","iter", "lower", "upper", "min", "err(est)");

		  printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lower, x_upper, x_mid, x_upper - x_lower);
	  
		  do
		    {
   
		      iter++;
		      status = gsl_min_fminimizer_iterate (s);

		      x_mid = gsl_min_fminimizer_x_minimum (s);
		      x_lower = gsl_min_fminimizer_x_lower (s);
		      x_upper = gsl_min_fminimizer_x_upper (s);

		      status = gsl_min_test_interval (x_lower, x_upper, 0.001, 0.0);

		      if (status == GSL_SUCCESS)
			printf ("Converged:\n");

		      printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lower, x_upper, x_mid, x_upper - x_lower);
	      
		    }
		  while (status == GSL_CONTINUE && iter < max_iter);

		  found = 1;

		  f_mid = curvaturei(x_mid);
	  
		}
	      else
		{
	  
		  if ( f_mid > f_lower)
		    {

		      x_lower = x_mid;
		      f_lower = f_mid;
	      
		      x_mid = (x_lower + x_upper) / 2;
		      f_mid = curvaturei(x_mid);

		      if ( (x_mid - x_lower) < 0.001)
			found = 1;
	      
		    }
		  else
		    {

		      x_upper = x_mid;
		      f_upper = f_mid;
	      
		      x_mid = (x_lower + x_upper) / 2;
		      f_mid = curvaturei(x_mid);

		      if ( (x_mid - x_lower) < 0.001)
			found = 1;
		    }
		}
	    }
	}
      
      while (found ==0);

      if (f_mid > f_best)
	f_best = f_mid;

    }

  return f_best;

}
  


/*
  double block1 = 4 * sqrt(135+144 * pow(lx->ve[ci*4+1],2.0)) * pow(lx->ve[ci*4+3],2.0) - 96 * lx->ve[ci*4+1] * pow(lx->ve[ci*4+2],2.0) * lx->ve[ci*4+3]+16 * pow(lx->ve[ci*4+2],4.0) / (45 * pow(lx->ve[ci*4+3],2.0));

  double block2 = (12 * lx->ve[ci*4+1] * lx->ve[ci*4+3] - 4 * pow(lx->ve[ci*4 + 2],2.0)) / (45 * pow(lx->ve[ci*4+3],2.0));
    
  double x1 = -sqrt( -block1 - block2) / 2.0 - lx->ve[ci*4+2] / (3.0*lx->ve[ci*4+3]);

  double x2 = sqrt( -block1 - block2) / 2.0 - lx->ve[ci*4+2] / (3.0*lx->ve[ci*4+3]);

  double x3 = sqrt( block1 - block2) /2.0 - lx->ve[ci*4+2] / (3.0*lx->ve[ci*4+3]);
  */

double distance(

		gsl_vector *x,
		void *param

		)

{

  double d0 = gsl_vector_get(x,0);
  double d1 = gsl_vector_get(x,1);

  lb->ve[N+S-1+S-1+S-1] = d0;
  lb->ve[N+S-1+S-1+S-1+2] = d1;

  lx = LUsolve(LU,pivot,lb,VNULL);
  
  double dist=0;

  double nons = (double)N / (double) S;
  
  for (int i=0;i<S;i++)
    for (int j=0;j<100;j++)
      dist += pow(effort[i] - lx->ve[i*4] + lx->ve[i*4+1]*(i+j/100.0)*nons + lx->ve[i*4+2]*pow((i+j/100.0)*nons,2.0) + lx->ve[i*4+3] * pow((i+j/100.0)*nons,3.0) , 2.0);

  return dist;

}

void process_Normal_Keys(int key, int xlah, int ylah) 
{
     switch (key) 
    {    
       case 27 :      break;
       case 100 :

	 //status = gsl_multimin_fminimizer_iterate(s);
	 
	 printf("GLUT_KEY_LEFT %d\n",key);

	 //do
	 //{
	 //iter++;

	 status = gsl_multimin_fminimizer_iterate(s1);

      /*
      if (status) 
        break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 5e-2);

      if (status == GSL_SUCCESS)
        {
          printf ("converged to minimum at\n");
        }

      printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", 
              iter,
              gsl_vector_get (s->x, 0), 
              gsl_vector_get (s->x, 1), 
              //gsl_vector_get (s->x, 2), 
              s->fval, size);
    }
  while (status == GSL_CONTINUE && iter < 100);  
      */
	 
  //	 lb->ve[N+S-1+S-1+S-1] += .1;
	   
  //	 lx = LUsolve(LU,pivot,lb,VNULL);
	 
	 glutPostRedisplay();
	 
	 break;
	 
       case 102: printf("GLUT_KEY_RIGHT %d\n",key);

	 gsl_multimin_fminimizer_set (s1,&f1,q1,ss1);
	 
	 glutPostRedisplay();
	 
	 break;
	 
       case 101 : printf("GLUT_KEY_UP %d\n",key); 
 
	 qb = gsl_vector_alloc (Nparamsb);

	 for (int i=0;i<Nparams1;i++)
	   gsl_vector_set (qb, i, gsl_vector_get(s1->x,i));

	 for (int i=1;i<S;i++)
	   gsl_vector_set (qb, Nparams1 + i-1, x[i]);

	 gsl_multimin_fminimizer_set (sb,&fb,qb,ssb);

	 glutPostRedisplay();
	 
	 break;

    case 103 : printf("GLUT_KEY_DOWN %d\n",key);   
	 
      //	 lb->ve[N+S-1+S-1+S-1+2] -= .1;
	 	   	 
      // lx = LUsolve(LU,pivot,lb,VNULL);

      status = gsl_multimin_fminimizer_iterate(sb);

      printf("%lf %lf %lf %lf\n",x[1],x[2],x[3],sb->fval);

	 glutPostRedisplay();

    	 break;

    }

}

void display(void)
{

  double nons = (double)N / (double) S;
  
  glClear(GL_COLOR_BUFFER_BIT);

  glPointSize(5.0f);
  
  glBegin(GL_POINTS);
  
  glColor3f(0.0f,0.0f,1.0f);

  for (int i=0;i<I;i++)
    glVertex2f(i/((float)(I-.5)),z[i]);
  
  glEnd();

  glBegin(GL_POINTS);
  
  glColor3f(1.0f,0.0f,1.0f);

  for (int i=0;i<I2;i++)
    glVertex2f(i/((float)(I2-.5)),.7*z2[i]);
  
  glEnd();
    
  glBegin(GL_POINTS);
  
  glColor3f(0.0f,1.0f,0.0f);

  for (int i=0;i<S;i++)
    {
      for (int j=0;j<100;j++)
	{
	  double xx = x[i]+j*(x[i+1]-x[i])/100.0;
	  
	  glVertex2f(xx/N,lx->ve[i*4] + lx->ve[i*4+1]*xx + lx->ve[i*4+2]*pow(xx,2.0) + lx->ve[i*4+3] * pow(xx,3.0));
	}
    }
  
  glEnd();


  // new curve!
  glBegin(GL_POINTS);
  
  glColor3f(0.0f,1.0f,0.5f);

  for (int i=0;i<S2;i++)
    {
      for (int j=0;j<100;j++)
	{
	  double xx = x2[i]+j*(x2[i+1]-x2[i])/100.0;
	  
	  glVertex2f(xx/N,.7*(lx2->ve[i*4] + lx2->ve[i*4+1]*xx + lx2->ve[i*4+2]*pow(xx,2.0) + lx2->ve[i*4+3] * pow(xx,3.0)));
	}
    }
  
  glEnd();



  glBegin(GL_POINTS);
  
  glColor3f(1.0f,0.1f,0.1f);

  int i=0;
  int j=0;
  int counter=0;
  int k=0;

  for (int l=0;l<100;l++)
    {

      double xx = l*(N2x[i+1]-x2[j-counter+k])/100 + x2[j-counter+k];
      int oldi = S;
      while (x[oldi]>xx) { oldi--; }

      glVertex2f(xx/N,lx->ve[4*oldi] + lx->ve[4*oldi+1]*xx + lx->ve[4*oldi+2]*pow(xx,2.0) + lx->ve[4*oldi+3]*pow(xx,3.0));

    }

  glEnd();



  glBegin(GL_POINTS);
  
  glColor3f(0.2f,0.4f,0.1f);

  i=1;
  j=1;

  for (int l=0;l<100;l++)
    {

      double xx = l*(x2[j]-N2x[i])/100 + N2x[i];

      int oldi = S;
      while (x[oldi]>xx) { oldi--; }

      glVertex2f(xx/N,lx->ve[4*oldi] + lx->ve[4*oldi+1]*xx + lx->ve[4*oldi+2]*pow(xx,2.0) + lx->ve[4*oldi+3]*pow(xx,3.0));

    }

  glEnd();



  double * erk = (double *) calloc(100*S,sizeof(double));

  for (int i=0;i<S;i++)
    for (int j=0;j<100;j++)
      {
	double xx = x[i]+j*(x[i+1]-x[i])/100.0;
	erk[i*100+j] = (1/5.0)*fabs(2*lx->ve[i*4+2] + 6*lx->ve[i*4+3] * xx) / pow( 1 + pow(lx->ve[i*4+1] + 2*lx->ve[i*4+2]*xx + 3*lx->ve[i*4+3] * pow(xx,2.0),2.0), 2.0/3.0);


	//fabs(2*lx->ve[ci*4+2] + 6*lx->ve[ci*4+3] * x) / pow( 1 + pow(lx->ve[ci*4+1] + 2*lx->ve[ci*4+2]*x + 3*lx->ve[ci*4+3] * pow(x,2.0),2.0), 2.0/3.0 );

      }

  double maxc = 0;
  for (int i=0;i<100*S;i++)
    if (erk[i] > maxc)
      maxc = erk[i];

  double minc = 0;
  for (int i=0;i<100*S;i++)
    if (erk[i] > minc)
      minc = erk[i];
  
  glBegin(GL_POINTS);
  
  glColor3f(1.0f,0.0f,0.0f);

  for (int i=0;i<S;i++)
    for (int j=0;j<100;j++)
      {
	  double xx = x[i]+j*(x[i+1]-x[i])/100.0;	  	
	  glVertex2f(xx/N,erk[i*100+j]);
      }

  glColor3f(1.0f,0.0f,1.0f);
  
  double xx = 1.5;
  i=1;
  glVertex2f(xx/N, (1/5.0)*fabs(2*lx->ve[i*4+2] + 6*lx->ve[i*4+3] * xx) / pow( 1 + pow(lx->ve[i*4+1] + 2*lx->ve[i*4+2]*xx + 3*lx->ve[i*4+3] * pow(xx,2.0),2.0), 2.0/3.0));

  glEnd();
    
  glFlush ();
  glutSwapBuffers();
}

void idle(void)
{
  glutPostRedisplay();
}

void Reshape(int width, int height)
{
  
  glViewport(0, 0, width, height);
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  gluOrtho2D(-0.1,1.1,0,2.2);

  glMatrixMode(GL_MODELVIEW);
  
  WindWidth = width;
  WindHeight = height;
}

void init(void)
{
   glClearColor (0.0, 0.0, 0.0, 0.0);
   glShadeModel (GL_FLAT);

   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   
}

int main(int argc, char *argv[])
{
 
  WindWidth = 1800;
  WindHeight = 1800;

  glutInit(&argc, (char **)argv);
    
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
  glutInitWindowSize(1800, 1800);
  
  tbInit(GLUT_LEFT_BUTTON);
  tbAnimate(GL_TRUE);

  window = glutCreateWindow("smooth");
  glutSpecialFunc( process_Normal_Keys );

  glutDisplayFunc(display);
  glutReshapeFunc(Reshape);
  
  init();  

  FILE * fp = fopen("bc2aex.dat","r");

  int K;
  
  if (fscanf(fp,"%d",&K) <1)
    {
      printf("error reading I\n");
      exit(1);
    }

  double * raweffort = (double *) calloc(K,sizeof(double));

  for (int i=0;i<K;i++)
    if ( fscanf(fp, "%lf ",&raweffort[i]) < 1)
      {
	printf("error reading raweffort %d\n",i);
	exit(1);
      }

  double * begintime = (double *) calloc(K, sizeof(double));
  
  for (int i=0;i<K;i++)
    if ( fscanf(fp, "%lf ",&begintime[i]) < 1)
      {
	printf("error reading begintime %d\n",i);
	exit(1);
      }

  double * endtime = (double *) calloc(K, sizeof(double));
  
  for (int i=0;i<K;i++)
    if ( fscanf(fp, "%lf ",&endtime[i]) < 1)
      {
	printf("error reading endtime %d\n",i);
	exit(1);
      }
  
  fclose(fp);

  N = 3;  // no. data blocks

  double mintime=begintime[0];
  for (int i=1;i<K;i++)
    if (begintime[i] < mintime)
      mintime = begintime[i];

  double maxtime=endtime[0];
  for (int i=1;i<K;i++)
    if (endtime[i] > maxtime)
      maxtime = endtime[i];
  
  effort = (double *) calloc(N,sizeof(double));
  
  double tottime = maxtime - mintime + 1e-12;

  //printf("%lf %lf %lf\n",mintime,maxtime,tottime);
  
  for (int i=0;i<N;i++)
    for (int j=0;j<K;j++)
      {

	double periodstart = mintime + i*tottime/N;
	double periodend = mintime + (i+1)*tottime/N;
	
	if (begintime[j] >= periodstart && begintime[j] < periodend)
	  {

	    if (endtime[j] < periodend)
	      effort[i] += raweffort[j];
	    else	      
              effort[i] += raweffort[j] * (periodend - begintime[j]) / (endtime[j] - begintime[j]); 			      

	  }
      }		      
  
  M = 10;
  I = M*N;

  double toteffort = 0;
  for (int i=0;i<N;i++)
    toteffort += effort[i];

  double avgeffort = toteffort / N;
  
  for (int i=0;i<N;i++)
    effort[i] /= avgeffort;
  
  z = (double *) calloc(I,sizeof(double));  
  for ( int i=0;i<N;i++)
    for (int m=0;m<M;m++)
      z[i*M+m] = effort[i];

  S=4; // no. splines
  
  A = m_get(4*S,4*S);

  x = (double *) calloc(S+1,sizeof(double));

  for (int i=0;i<=S;i++)
    x[i] = i * ( (double)N / (double)S);

  int j=0;

  for (int i=0;i<N;i++)
    {

      if (x[j] > i)
	{
	  A->me[i][(j-1)*4 + 0] += x[j] - i;
	  A->me[i][(j-1)*4 + 1] += .5*( pow(x[j],2.0) - pow(i,2.0) );
	  A->me[i][(j-1)*4 + 2] += (1.0/3) * ( pow(x[j],3.0) - pow(i,3.0));
	  A->me[i][(j-1)*4 + 3] += (1.0/4) * ( pow(x[j],4.0) - pow(i,4.0));
	}
      
      int counter=0;
	        
      while ( x[j+1] < i+1 ) { j++; counter++; }

      int k=0;
      
      for (k=0;k<counter;k++)
	{	
	  A->me[i][(j-1)*4 + 0] += x[j-counter+k+1] - x[j-counter+k];
	  A->me[i][(j-1)*4 + 1] += .5*( pow(x[j-counter+k+1],2.0) - pow(x[j-counter+k],2.0) );
	  A->me[i][(j-1)*4 + 2] += (1.0/3) * ( pow(x[j-counter+k+1],3.0) - pow(x[j-counter+k],3.0));
	  A->me[i][(j-1)*4 + 3] += (1.0/4) * ( pow(x[j-counter+k+1],4.0) - pow(x[j-counter+k],4.0));

	}

      A->me[i][j*4 + 0] += i+1-x[j-counter+k];
      A->me[i][j*4 + 1] += .5*( pow(i+1,2.0) - pow(x[j-counter+k],2.0) );
      A->me[i][j*4 + 2] += (1.0/3) * ( pow(i+1,3.0) - pow(x[j-counter+k],3.0));
      A->me[i][j*4 + 3] += (1.0/4) * ( pow(i+1,4.0) - pow(x[j-counter+k],4.0));

      j++; 
      	  
    }

  for (int j=0;j<S-1;j++)
    {

      A->me[N+j][j*4 + 0] = 1;
      A->me[N+j][j*4 + 1] = x[j+1];
      A->me[N+j][j*4 + 2] = pow(x[j+1],2.0);
      A->me[N+j][j*4 + 3] = pow(x[j+1],3.0);
      A->me[N+j][j*4 + 4] = -1;
      A->me[N+j][j*4 + 5] = -(x[j+1]);
      A->me[N+j][j*4 + 6] = -pow(x[j+1],2.0);
      A->me[N+j][j*4 + 7] = -pow(x[j+1],3.0);
      
    }

  for (int j=0;j<S-1;j++)
    {

      A->me[N+S-1+j][j*4 + 1] = 1;
      A->me[N+S-1+j][j*4 + 2] = 2*(x[j+1]);
      A->me[N+S-1+j][j*4 + 3] = 3*pow(x[j+1],2.0);
      A->me[N+S-1+j][j*4 + 5] = -1;
      A->me[N+S-1+j][j*4 + 6] = -2*(x[j+1]);
      A->me[N+S-1+j][j*4 + 7] = -3*pow(x[j+1],2.0);
      
    }

  for (int j=0;j<S-1;j++)
    { 
      A->me[N+S-1+S-1+j][j*4 + 2] = 2;
      A->me[N+S-1+S-1+j][j*4 + 3] = 6*(x[j+1]);
      A->me[N+S-1+S-1+j][j*4 + 6] = -2;
      A->me[N+S-1+S-1+j][j*4 + 7] = -6*(x[j+1]);
    }
  
  A->me[N+S-1+S-1+S-1][2] = 2; // second derivative at t=0

  A->me[N+S-1+S-1+S-1+1][0] = 1;  // value at t=0;
  
  A->me[N+S-1+S-1+S-1+2][N+S-1+S-1+S-1+3-2] = 1;  // first derivative at t=N
  A->me[N+S-1+S-1+S-1+2][N+S-1+S-1+S-1+3-1] = 2*N; // first derivative at t=N
  A->me[N+S-1+S-1+S-1+2][N+S-1+S-1+S-1+3] = 3*pow(N,2.0);  // first derivative at t=N

  A->me[N+S-1+S-1+S-1+3][N+S-1+S-1+S-1+3-1] = 2;  // second derivative at t=N
  A->me[N+S-1+S-1+S-1+3][N+S-1+S-1+S-1+3] = 6*N;  

  lb = v_get(4*S);

  for (int i=0;i<N;i++)
    lb->ve[i] = effort[i];

  LU = m_get(4*S,4*S);
  LU = m_copy(A,LU);
  pivot = px_get(A->m);
  LUfactor(LU,pivot);

  lx = LUsolve(LU,pivot,lb,VNULL);

  //v_output(lx);
  
  T = gsl_multimin_fminimizer_nmsimplex2;

  s1 = NULL;
  sb = NULL;

  // Set initial step sizes to 0.1 

  Nparams1 = 3;
  Nparamsb = Nparams1 + S-1;

  ss1 = gsl_vector_alloc (Nparams1);
  gsl_vector_set_all (ss1, 0.1);

  f1.n = Nparams1;
  f1.f = &curvature; 
  f1.params = NULL;

  q1 = gsl_vector_alloc (Nparams1);
  qb = gsl_vector_alloc (Nparamsb);

  gsl_vector_set (q1, 0, 0);
  gsl_vector_set (q1, 1, 0);
  gsl_vector_set (q1, 2, 0);
  
  s1 = gsl_multimin_fminimizer_alloc(T,Nparams1);

  ssb = gsl_vector_alloc (Nparamsb);
  gsl_vector_set_all (ssb, 0.1);

  fb.n = Nparamsb;
  fb.f = &curvatureboth; 
  fb.params = NULL;

  sb = gsl_multimin_fminimizer_alloc(T,Nparamsb);

  qb = gsl_vector_alloc (Nparamsb);

  for (int i=0;i<Nparams1;i++)
    gsl_vector_set (qb, i, gsl_vector_get(q1,i));

  for (int i=1;i<S;i++)
    gsl_vector_set (qb, Nparams1 + i-1, x[i]);

  gsl_multimin_fminimizer_set (sb,&fb,qb,ssb);

  do {

    iter++;
    status = gsl_multimin_fminimizer_iterate(sb);
      
    if (status) 
      break;

    size = gsl_multimin_fminimizer_size (sb);
    status = gsl_multimin_test_size (size, 1e-3);

    if (status == GSL_SUCCESS)
      {
	//printf ("converged to minimum at\n");
      }

    //printf ("%5d %10.3e %10.3e %10.3e %10.3e f() = %7.3f size = %.3f\n", iter, gsl_vector_get (sb->x, 0), gsl_vector_get (sb->x, 1), gsl_vector_get (sb->x, 2), gsl_vector_get (sb->x, 3), sb->fval, size);
  }
  while (status == GSL_CONTINUE && iter < 100);

  /*  
  double * testeffort = (double *) calloc(N,sizeof(double));

  int id;
  for (int i=0;i<N;i++) 
    {
      for (int l=0;l<100;l++)
	{
	  double xx = l/100.0 + i;
	  if (xx<x[1])
	    id = 0;
	  if (xx>=x[1] && xx < x[2])
	    id = 1;
	  if (xx>=x[2] && xx < x[3])
	    id = 2;
	  if (xx>=x[3] && xx < x[4])
	    id = 3;

	  testeffort[i] += (1.0/100) * (lx->ve[4*id] + lx->ve[4*id+1]*xx + lx->ve[4*id+2]*pow(xx,2.0) + lx->ve[4*id+3]*pow(xx,3.0));
	}

    }
  */

  N2 = 4;  // no. data blocks
 
  effort2 = (double *) calloc(N2,sizeof(double));  
 
  N2x = (double *) calloc(N2+1,sizeof(double));

  for (int i=0;i<=N2;i++)
    N2x[i] = i * (double)N / (double)N2;

  
  for (int i=0;i<N2;i++)
    {

      double periodstart = mintime + i*tottime/N2;
      double periodend = mintime + (i+1)*tottime/N2;

    for (int j=0;j<K;j++)
      {

	
	if (begintime[j] >= periodstart && begintime[j] < periodend)
	  {

	    if (endtime[j] < periodend)
	      effort2[i] += raweffort[j];
	    else	      
              effort2[i] += raweffort[j] * (periodend - begintime[j]) / (endtime[j] - begintime[j]); 			      

	  }
      }	
    }
	      

  I2 = N2*M;
  
  toteffort = 0;
  for (int i=0;i<N2;i++)
    toteffort += effort2[i];

  avgeffort = toteffort / N2;
  
  for (int i=0;i<N2;i++)
    effort2[i] /= avgeffort;

  z2 = (double *) calloc(I2,sizeof(double));  
  for ( int i=0;i<N2;i++)
    for (int m=0;m<M;m++)
      z2[i*M+m] = effort2[i] * 1/.75;

  S2=5; // no. splines
  
  double gap=0;
  int idx=-1;
  for (int i=0;i<S;i++)
    if (x[i+1]-x[i] > gap)
      {
	gap = x[i+1] - x[i];
	idx = i;
      }

  x2 = (double *) calloc(S2+1,sizeof(double));  

  int i;
  for (i=0;i<=idx;i++)
    x2[i] = x[i];

  x2[i] = (x[i]+x[i-1])/2;

  for (i=idx+1;i<S2;i++)
    x2[i+1] = x[i];

  A2 = m_get(4*S2,4*S2);

  double * fakeeffort = (double *) calloc(N2,sizeof(double));

  j=0;
  for (int i=0;i<N2;i++)
    {

      if (x2[j] > N2x[i])
	{

          double tmpfakeeffort=0;
	  for (int l=0;l<100;l++)
	    {

	      double xx = l*(x2[j]-N2x[i])/100 + N2x[i];
	      int oldi = S;
	      while (x[oldi]>xx) { oldi--; }

	      tmpfakeeffort += ((x2[j]-N2x[i])/100) * (lx->ve[4*oldi] + lx->ve[4*oldi+1]*xx + lx->ve[4*oldi+2]*pow(xx,2.0) + lx->ve[4*oldi+3]*pow(xx,3.0));

	    }
	  fakeeffort[i] += tmpfakeeffort; ///(N2x[i]-x2[j]);
	}

      int counter=0;
      while ( x2[j+1] < N2x[i+1] ) { j++; counter++; }

      int k=0;
      
      for (k=0;k<counter;k++)
	{

          double tmpfakeeffort=0;
	  for (int l=0;l<100;l++)
	    {

	      double xx = l*(x2[j-counter+k+1]-x2[j-counter+k])/100 + x2[j-counter+k];
	      int oldi = S;
	      while (x[oldi]>xx) { oldi--; }

	      tmpfakeeffort += (x2[j-counter+k+1]-x2[j-counter+k])/100 * (lx->ve[4*oldi] + lx->ve[4*oldi+1]*xx + lx->ve[4*oldi+2]*pow(xx,2.0) + lx->ve[4*oldi+3]*pow(xx,3.0));

	    }
	  fakeeffort[i] += tmpfakeeffort; ///(x2[j-counter+k+1]-x2[j-counter+k]);
	}

      double tmpfakeeffort=0;
      for (int l=0;l<100;l++)
	{

	  double xx = l*(N2x[i+1]-x2[j-counter+k])/100 + x2[j-counter+k];
	  int oldi = S;
	  while (x[oldi]>xx) { oldi--; }

	  tmpfakeeffort += (N2x[i+1]-x2[j-counter+k])/100 * (lx->ve[4*oldi] + lx->ve[4*oldi+1]*xx + lx->ve[4*oldi+2]*pow(xx,2.0) + lx->ve[4*oldi+3]*pow(xx,3.0));

	}

      fakeeffort[i] += tmpfakeeffort; //x[i+1]-x2[j-counter+k]);

      j++; 
    }

  j=0;

  for (int i=0;i<N2;i++)
    {

      if (x2[j] > N2x[i])
	{
	  A2->me[i][(j-1)*4 + 0] += x2[j] - N2x[i];
	  A2->me[i][(j-1)*4 + 1] += .5*( pow(x2[j],2.0) - pow(N2x[i],2.0) );
	  A2->me[i][(j-1)*4 + 2] += (1.0/3) * ( pow(x2[j],3.0) - pow(N2x[i],3.0));
	  A2->me[i][(j-1)*4 + 3] += (1.0/4) * ( pow(x2[j],4.0) - pow(N2x[i],4.0));
	}
      
      int counter=0;
	        
      while ( x2[j+1] < N2x[i+1] ) { j++; counter++; }

      int k=0;
      
      for (k=0;k<counter;k++)
	{
	
	  A2->me[i][(j-1)*4 + 0] += x2[j-counter+k+1] - x2[j-counter+k];
	  A2->me[i][(j-1)*4 + 1] += .5*( pow(x2[j-counter+k+1],2.0) - pow(x2[j-counter+k],2.0) );
	  A2->me[i][(j-1)*4 + 2] += (1.0/3) * ( pow(x2[j-counter+k+1],3.0) - pow(x2[j-counter+k],3.0));
	  A2->me[i][(j-1)*4 + 3] += (1.0/4) * ( pow(x2[j-counter+k+1],4.0) - pow(x2[j-counter+k],4.0));

	}

      A2->me[i][j*4 + 0] += N2x[i+1]-x2[j-counter+k];
      A2->me[i][j*4 + 1] += .5*( pow(N2x[i+1],2.0) - pow(x2[j-counter+k],2.0) );
      A2->me[i][j*4 + 2] += (1.0/3) * ( pow(N2x[i+1],3.0) - pow(x2[j-counter+k],3.0));
      A2->me[i][j*4 + 3] += (1.0/4) * ( pow(N2x[i+1],4.0) - pow(x2[j-counter+k],4.0));

      j++; 
      	  
    }

  for (int j=0;j<S2-1;j++)
    {

      A2->me[N2+j][j*4 + 0] = 1;
      A2->me[N2+j][j*4 + 1] = x2[j+1];
      A2->me[N2+j][j*4 + 2] = pow(x2[j+1],2.0);
      A2->me[N2+j][j*4 + 3] = pow(x2[j+1],3.0);
      A2->me[N2+j][j*4 + 4] = -1;
      A2->me[N2+j][j*4 + 5] = -(x2[j+1]);
      A2->me[N2+j][j*4 + 6] = -pow(x2[j+1],2.0);
      A2->me[N2+j][j*4 + 7] = -pow(x2[j+1],3.0);
      
    }

  for (int j=0;j<S2-1;j++)
    {

      A2->me[N2+S2-1+j][j*4 + 1] = 1;
      A2->me[N2+S2-1+j][j*4 + 2] = 2*(x2[j+1]);
      A2->me[N2+S2-1+j][j*4 + 3] = 3*pow(x2[j+1],2.0);
      A2->me[N2+S2-1+j][j*4 + 5] = -1;
      A2->me[N2+S2-1+j][j*4 + 6] = -2*(x2[j+1]);
      A2->me[N2+S2-1+j][j*4 + 7] = -3*pow(x2[j+1],2.0);
      
    }

  for (int j=0;j<S2-1;j++)
    { 
      A2->me[N2+S2-1+S2-1+j][j*4 + 2] = 2;
      A2->me[N2+S2-1+S2-1+j][j*4 + 3] = 6*(x2[j+1]);
      A2->me[N2+S2-1+S2-1+j][j*4 + 6] = -2;
      A2->me[N2+S2-1+S2-1+j][j*4 + 7] = -6*(x2[j+1]);
    }
  
  A2->me[N2+S2-1+S2-1+S2-1][2] = 2; // second derivative at t=0

  A2->me[N2+S2-1+S2-1+S2-1+1][0] = 1;  // value at t=0;
  
  A2->me[N2+S2-1+S2-1+S2-1+2][N2+S2-1+S2-1+S2-1+3-2] = 1;  // first derivative at t=N2
  A2->me[N2+S2-1+S2-1+S2-1+2][N2+S2-1+S2-1+S2-1+3-1] = 2*N2x[N2]; // first derivative at t=N2
  A2->me[N2+S2-1+S2-1+S2-1+2][N2+S2-1+S2-1+S2-1+3] = 3*pow(N2x[N2],2.0);  // first derivative at t=N2

  A2->me[N2+S2-1+S2-1+S2-1+3][N2+S2-1+S2-1+S2-1+3-1] = 2;  // second derivative at t=N2
  A2->me[N2+S2-1+S2-1+S2-1+3][N2+S2-1+S2-1+S2-1+3] = 6*N2x[N2];  
 
  lb2 = v_get(4*S2);

  double p = 0;
  for (int i=0;i<N2;i++)
    lb2->ve[i] = effort2[i]; //p*fakeeffort[i] + (1-p)*effort2[i];

  lb2->ve[N2+S2-1+S2-1+S2-1] = gsl_vector_get(sb->x,0);
  lb2->ve[N2+S2-1+S2-1+S2-1+2] = gsl_vector_get(sb->x,1);
  lb2->ve[N2+S2-1+S2-1+S2-1+3] = gsl_vector_get(sb->x,2);

  //v_output(lb2);

  LU2 = m_get(4*S2,4*S2);
  LU2 = m_copy(A2,LU2);
  pivot = px_get(A2->m);
  LUfactor(LU2,pivot);

  lx2 = LUsolve(LU2,pivot,lb2,VNULL);

  //for (int i=0;i<=N2;i++)
  //  printf("%lf\n",N2x[i]);
  //exit(1);

  //printf("%lf %lf\n",effort2[1],effort2[1] * (N2x[2]-N2x[1]));
  //exit(1);
  /*
  printf("\n");

  for (int i=0;i<N2;i++)
    for(int j=0;j<100;j++)
      printf("%lf %lf\n",j*(N2x[i+1]-N2x[i])/100.0+N2x[i],effort2[i]* (1/ (N2x[i+1]-N2x[i])));

  printf("e\n");

  for (int i=0;i<S2;i++)
    for (int j=0;j<100;j++)
      {
	double xx = j*(x2[i+1]-x2[i])/100 + x2[i];
	printf("%lf %lf\n",xx,lx2->ve[4*i] + lx2->ve[4*i+1]*xx + lx2->ve[4*i+2]*pow(xx,2.0) + lx2->ve[4*i+3]*pow(xx,3.0));
      }

  exit(1);
  */

  /*

  i=0;
  j=0;
  int counter=0;
  int k=0;
  double tmpfakeeffort=0;
  id=0;
  for (int l=0;l<100;l++)
    {

      double xx = l*(N2x[i+1]-x2[j-counter+k])/100 + x2[j-counter+k];
      tmpfakeeffort += (N2x[i+1]-x2[j-counter+k])/100 * (lx2->ve[4*id] + lx2->ve[4*id+1]*xx + lx2->ve[4*id+2]*pow(xx,2.0) + lx2->ve[4*id+3]*pow(xx,3.0));

    }

  printf("should: %lf is: %lf\n",effort2[0],tmpfakeeffort);

  i=1;
  j=1;

 tmpfakeeffort=0;
 for (int l=0;l<100000;l++)
   {
     double xx = l*(x2[j]-N2x[i])/100000 + N2x[i];
     tmpfakeeffort += ((x2[j]-N2x[i])/100000) * (lx2->ve[0] + lx2->ve[0+1]*xx + lx2->ve[0+2]*pow(xx,2.0) + lx2->ve[0+3]*pow(xx,3.0));
   }

 for (int l=0;l<100000;l++)
   {
     double xx = l*(x2[j+1]-x[j])/100000 + x2[j];
     tmpfakeeffort += ((x2[j+1]-x2[j])/100000) * (lx2->ve[j*4] + lx2->ve[j*4+1]*xx + lx2->ve[j*4+2]*pow(xx,2.0) + lx2->ve[j*4+3]*pow(xx,3.0));
   }

 j=2;
 for (int l=0;l<100000;l++)
   {
     double xx = l*(N2x[i+1]-x2[j])/100000 + x2[j];
     tmpfakeeffort += ((N2x[i+1]-x2[j])/100000) * (lx2->ve[j*4] + lx2->ve[j*4+1]*xx + lx2->ve[j*4+2]*pow(xx,2.0) + lx2->ve[j*4+3]*pow(xx,3.0));
   }

  printf("should: %lf is: %lf\n",effort2[i],tmpfakeeffort);

  i=2;
  j=3;

 tmpfakeeffort=0;
 for (int l=0;l<100000;l++)
   {
     double xx = l*(x2[j]-N2x[i])/100000 + N2x[i];
     tmpfakeeffort += ((x2[j]-N2x[i])/100000) * (lx2->ve[(j-1)*4] + lx2->ve[(j-1)*4+1]*xx + lx2->ve[(j-1)*4+2]*pow(xx,2.0) + lx2->ve[(j-1)*4+3]*pow(xx,3.0));
   }

 for (int l=0;l<100000;l++)
   {
     double xx = l*(N2x[i+1]-x2[j])/100000 + x2[j];
     tmpfakeeffort += ((N2x[i+1]-x2[j])/100000) * (lx2->ve[j*4] + lx2->ve[j*4+1]*xx + lx2->ve[j*4+2]*pow(xx,2.0) + lx2->ve[j*4+3]*pow(xx,3.0));
   }

  printf("should: %lf is: %lf\n",effort2[i],tmpfakeeffort);




  i=3;
  j=4;

 tmpfakeeffort=0;
 for (int l=0;l<100000;l++)
   {
     double xx = l*(x2[j]-N2x[i])/100000 + N2x[i];
     tmpfakeeffort += ((x2[j]-N2x[i])/100000) * (lx2->ve[(j-1)*4] + lx2->ve[(j-1)*4+1]*xx + lx2->ve[(j-1)*4+2]*pow(xx,2.0) + lx2->ve[(j-1)*4+3]*pow(xx,3.0));
   }

 for (int l=0;l<100000;l++)
   {
     double xx = l*(N2x[i+1]-x2[j])/100000 + x2[j];
     tmpfakeeffort += ((N2x[i+1]-x2[j])/100000) * (lx2->ve[j*4] + lx2->ve[j*4+1]*xx + lx2->ve[j*4+2]*pow(xx,2.0) + lx2->ve[j*4+3]*pow(xx,3.0));
   }

  printf("should: %lf is: %lf\n",effort2[i],tmpfakeeffort);








  exit(1);
  */
  sb2 = NULL;

  // Set initial step sizes to 0.1 

  Nparams1 = 3;
  Nparamsb = Nparams1 + S2-1;

  ssb2 = gsl_vector_alloc (Nparamsb);

  gsl_vector_set_all (ssb2, 0.1);

  fb2.n = Nparamsb;
  fb2.f = &curvatureboth2; 
  fb2.params = NULL;

  sb2 = gsl_multimin_fminimizer_alloc(T,Nparamsb);
  qb2 = gsl_vector_alloc (Nparamsb);

  for (int i=0;i<Nparams1;i++)
    gsl_vector_set (qb2, i, gsl_vector_get(qb,i));

  for (int i=1;i<S2;i++)
    gsl_vector_set (qb2, Nparams1 + i-1, x2[i]);

  /*
  gsl_multimin_fminimizer_set (sb2,&fb2,qb2,ssb2);

  iter = 0;

  do {

    iter++;
    status = gsl_multimin_fminimizer_iterate(sb2);
      
    if (status) 
      break;

    size = gsl_multimin_fminimizer_size (sb2);
    status = gsl_multimin_test_size (size, 1e-3);

    if (status == GSL_SUCCESS)
      {
	printf ("converged to minimum at\n");
      }

    printf ("%5d %10.3e %10.3e %10.3e %10.3e f() = %7.3f size = %.3f\n", iter, gsl_vector_get (sb2->x, 0), gsl_vector_get (sb2->x, 1), gsl_vector_get (sb2->x, 2), gsl_vector_get (sb2->x, 3), sb2->fval, size);
  }
  while (status == GSL_CONTINUE && iter < 100);
  
  */


  glutMainLoop();

  return(0);
  
}

