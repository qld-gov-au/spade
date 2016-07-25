int a[256], b[256], c[256];


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

  x[J] -= 1e-5;

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
    u[j] = (a*A1*vbar+a*A2*wbar)*pow(w-x[j],(b+g*ubar)/k-1) / (k*pow(w,(b+g*ubar)/k));

}


int foo (int n, int x) {
   int i;

   /* feature: support for unknown loop bound  */
   /* feature: support for loop invariants  */
   for (i=0; i<n; i++){
      b[i] = x;
   }

   /* feature: general loop exit condition  */
   /* feature: support for bitwise operations  */
   while (n--){
      a[i] = b[i]&c[i]; i++;
   }

   return(0);
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

  double k = 1./SPY;

  double **x; 
  double **u; 
  
  x = (double **) malloc((I+1) * sizeof(double *));
  for (int i=0;i<=I;i++)
    x[i] = (double *) malloc( (J+i+1) * sizeof(double));

  u = (double **) malloc((I+1) * sizeof(double *));
  for (int i=0;i<=I;i++)
    u[i] = (double *) malloc( (J+i+1) * sizeof(double));
  
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

  double * xhht = (double *) malloc((I+1+J) * sizeof(double));
  double * uhh = (double *) malloc((I+1+J) * sizeof(double));
  double * xht = (double *) malloc((I+1+J) * sizeof(double));
  double * uht = (double *) malloc((I+1+J) * sizeof(double));
  
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
            
      set_row(xhh,i-1,xhht);
      
      for (int j=0;j<=terminator;j++)
	    {
	      xht->ve[j+1] = xt->ve[j] + (k/2)*g(kk,ww,xhht->ve[j+1]);
	      uht->ve[j+1] = ut->ve[j]*exp(-(k/2)*zstar(eff,bb,gg,kk,ii,thh,xhht->ve[j+1],Uhh->ve[i-1],k,Y));
	    }

      Q2(aa,kk,ww,xht,uht);
      Uh->ve[i-1] = Q(xht,uht);
         
      set_row(xh,i-1,xht);
      set_row(uh,i-1,uht);

      for (int j=0;j<=terminator;j++)
	    {
	      xnt->ve[j+1] = xt->ve[j] + k*g(kk,ww,xht->ve[j+1]);
	      unt->ve[j+1] = ut->ve[j]*exp(-k*zstar(eff,bb,gg,kk,ii,th,xht->ve[j+1],Uh->ve[i-1],k,Y));
        }
        
      Q2(aa,kk,ww,xnt,unt);

      Ui->ve[i] = Q(xnt,unt);      

      if(!SGNM) 
        {
          xnt = v_resize(xnt,xn->n);
          unt = v_resize(unt,un->n);

          set_row(x,i,xnt); // todo: test this.
          set_row(u,i,unt);
        }
      else
        {
          set_row(xn,i-1,xnt);
          set_row(un,i-1,unt);

          int idx = idxselect(ww,xnt);
          idxi->ive[i-1] = idx;
        
          xt = idxremove(xnt,xt,idx);
          ut = idxremove(unt,ut,idx);

          set_row(x,i,xt);
          set_row(u,i,ut);
        }
    }

  V_FREE(xt);
  V_FREE(xnt); 
  V_FREE(unt); 
  V_FREE(xht); 
  V_FREE(uht);
  V_FREE(xhht); 
  V_FREE(uhh);
  V_FREE(ut);


int main(void)
{

  foo(256,256);

  return(0);
}
      
