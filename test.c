#include <math.h>
#include <stdlib.h>

#define A1 8.588e-5
#define A2 0.00144

double iota1=5.2;
double iota2=0.619;
double phi=17;
double eta1=1.703205e-5;
double eta2=2.9526;

double h = 1./160;

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

double e(

	  const double * ef,
	  double r,
	  double t,
    int Y

	  )
{

  if (t<0)
    {
      double cek = r/4;
      double ept5 = ef[(int)floor((.5 + (cek/2) - 1e-12)/cek)];
      double m = ept5/Y;

      return m*t+ ept5;
    }
  else
    {
      double cek = r/4;
      int idx = floor((t + (cek/2) - 1e-12)/cek);
      return ef[idx];
    }
}

double zstar(

	      const double *eff,
	      double b,
	      double g,
	      double kappa,
	      double i,
	      double t,
	      double x,
	      double U,
	      double r,
              int Y

	     )
{ // equation X in Y

  return b + g*U + s(x)*i*e(eff,r,t,Y) - kappa;

}

double g(

	 const double kappa,
	 const double w,
	 const double x

	 )
{ /* von-Bertalanffy growth */
  return kappa*(w - x);
}


double b(

	 const double a,
	 const double x

	 )
{ /* birth function */
  return a*(A1*x + A2*pow(x,2.));
}

double Q(
	 
	 const double * x,
	 const double * u,
	 int terminator

	 )
{ 

  double rt = 0.;

  for (int i=0;i<terminator;i++) 
    rt = rt + .5 * (x[i+1] - x[i]) * (u[i] + u[i+1]);
   
  return rt;

}

void Q2(

	const double a,
	const double k,
	const double w,
	const double *x,
	double *u,
	const int terminator
  
	)
{

  double rt = x[1] * b(a,x[1])*u[1] - x[0]*b(a,x[1])*u[1]; 

  for (int j=1;j<terminator;j++) 
    rt = rt + (b(a,x[j])*u[j] + b(a,x[j+1])*u[j+1]) * (x[j+1]-x[j]);

  u[0] = rt / (2*k*w + x[0]*b(a,x[0]) - x[1]*b(a,x[0]));

}

void initial(
	     
	     const double *parameters,
	     const double *x,
	     double *u,
	     const int terminator

	     )

{

  double a = parameters[0];
  double b = parameters[1];
  double g = parameters[2];
  double k = parameters[3];
  double w = parameters[4];

  double zeta = sqrt( 81*k*k*w*w*pow(a*A1+2*a*A2*w,2.) - 12*k*pow(a*A1*w+k,3.) );
  double eta = 9*a*A1*k*k*w + 18*a*A2*k*k*w*w + k*zeta;
  double Z = pow(eta,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*k*(a*A1*w+k) / pow(eta,1./3);

  double ubar = (Z - b - k) / g; 
  double vbar = (k*w*ubar) / (b+g*ubar+k);
  double wbar = (2*k*w*vbar) / (b+g*ubar+2*k);  

  for (int j=0;j<=terminator;j++) 
    u[j] = (a*A1*vbar+a*A2*wbar)*pow(w-(x[j]-1e-5),(b+g*ubar)/k-1) / (k*pow(w,(b+g*ubar)/k));

}

void loop(

	  double *parameters,
	  double *eff,
	  int J,
	  int Y,
          int SPY
	  
	  )
{

  int N = Y*SPY;
  int I = 2*N;
  int S = N;
  
  double k = 1./SPY;

  double **x; 
  double **u; 
  
  x = (double **) calloc((I+1),sizeof(double *));
  for (int i=0;i<=I;i++)
    x[i] = (double *) calloc( (J+i+1),sizeof(double));

  u = (double **) calloc((I+1),sizeof(double *));
  for (int i=0;i<=I;i++)
    u[i] = (double *) calloc( (J+i+1),sizeof(double));
  
  for (int j=0;j<=J;j++) 
    x[0][j] = h*j;

  int terminator = J;

  initial(parameters,x[0],u[0],terminator);

  double Ui = Q(x[0],u[0],terminator);
      
  double aa = parameters[0];
  double bb = parameters[1];
  double gg = parameters[2];
  double kk = parameters[3];
  double ww = parameters[4];
  double ii = parameters[5];

  double * xhht = (double *) calloc((I+1+J),sizeof(double));
  double * uhh = (double *) calloc((I+1+J),sizeof(double));
  double * xht = (double *) calloc((I+1+J),sizeof(double));
  double * uht = (double *) calloc((I+1+J),sizeof(double));
  double * xnt = (double *) calloc((I+1+J),sizeof(double));
  double * unt = (double *) calloc((I+1+J),sizeof(double));
  
  for (int i=1;i<=I;i++)
    {

      double t = k*(i-S-1);
      double th = k*(i-S-.5);
      double thh = k*(i-S-.75);

      int terminator = J + i-1;

      for (int j=0;j<=terminator;j++)
	{
	  xhht[j+1] = x[i-1][j] + (k/4)*g(kk,ww,x[i-1][j]);
	  uhh[j+1] = u[i-1][j]*exp(-(k/4)*zstar(eff,bb,gg,kk,ii,t,x[i-1][j],Ui,k,Y));
	}

      Q2(aa,kk,ww,xhht,uhh,terminator);
      double Uhh = Q(xhht,uhh,terminator);
                  
      for (int j=0;j<=terminator;j++)
	{
	  xht[j+1] = x[i-1][j] + (k/2)*g(kk,ww,xhht[j+1]);
	  uht[j+1] = u[i-1][j]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht[j+1],Uhh,k,Y));
	}

      Q2(aa,kk,ww,xht,uht,terminator);
      double Uh = Q(xht,uht,terminator);
         
      for (int j=0;j<=terminator;j++)
	{
	  xnt[j+1] = x[i-1][j] + k*g(kk,ww,xht[j+1]);
	  unt[j+1] = u[i-1][j]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht[j+1],Uh,k,Y));
        }
        
      Q2(aa,kk,ww,xnt,unt,terminator);

      for (int j=0;j<=terminator;j++)
	{
	  x[i][j] = xnt[j];
	  u[i][j] = unt[j];
        }
      
      Ui = Q(xnt,unt,terminator);      

    }

  free(xnt); 
  free(unt); 
  free(xht); 
  free(uht);
  free(xhht); 
  free(uhh);

}

int main(void)
{

  double *parameters = (double *) calloc(6,sizeof(double));

  parameters[0] = .1;
  parameters[1] = .1;
  parameters[2] = 1.1;
  parameters[3] = .1;
  parameters[4] = .1;
  parameters[5] = 160;
    
  int J = 400;
  int Y = 24;
  int SPY = 20;

  int N = Y*SPY;
  int I = 2*N;

  double *eff = (double *) calloc((2*I+1),sizeof(double));

  loop(parameters,eff,J,Y,SPY);
  
  return(0);
}
      
