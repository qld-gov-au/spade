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
double *z;
double *d1;
double *d2;
double *effort;

double *x;

VEC *lx;

VEC *lb;
MAT *A, *LU;
PERM *pivot;

typedef struct {

  int nBlocks;
  double *effort;
  int nPoints;

} GSL_reg;

const gsl_multimin_fminimizer_type *T;

gsl_multimin_fminimizer *s;
gsl_vector *ss, *qq, *qq2;
gsl_multimin_function minex;

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

double curvature(

		 gsl_vector *p,
		 void *param

		 )
{

  double d0 = gsl_vector_get(p,0);
  double d1 = gsl_vector_get(p,1);

  lb->ve[N+S-1+S-1+S-1] = d0;
  lb->ve[N+S-1+S-1+S-1+2] = d1;

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

  double dd0 = gsl_vector_get(x,0);
  double dd1 = gsl_vector_get(x,1);

  lb->ve[N+S-1+S-1+S-1] = dd0;
  lb->ve[N+S-1+S-1+S-1+1] = dd1;

  lx = LUsolve(LU,pivot,lb,VNULL);
  
  double dist=0;

  double nons = (double)N / (double) S;
  
  for (int i=0;i<S;i++)
    for (int j=0;j<100;j++)
      dist += pow(effort[i] - lx->ve[i*4] + lx->ve[i*4+1]*(i+j/100.0)*nons + lx->ve[i*4+2]*pow((i+j/100.0)*nons,2.0) + lx->ve[i*4+3] * pow((i+j/100.0)*nons,3.0) , 2.0);

  return dist;

}

void process_Normal_Keys(int key, int x, int y) 
{
     switch (key) 
    {    
       case 27 :      break;
       case 100 :

	 //status = gsl_multimin_fminimizer_iterate(s);
	 
	 printf("GLUT_KEY_LEFT %d\n",key);
	 
	 thingy += .1;
	 
	 lb->ve[N+S-1+S-1+S-1+2] = thingy;
	   
	 lx = LUsolve(LU,pivot,lb,VNULL);

	 printf("%lf\n",curvature(qq,0));
	 
	 //         v_output(lx);	 
	 
	 glutPostRedisplay();
	 
	 break;
	 
       case 102: printf("GLUT_KEY_RIGHT %d\n",key);

	 thingy -= .1;
	 
	 lb->ve[N+S-1+S-1+S-1+2] = thingy;
	 	   	 
	 lx = LUsolve(LU,pivot,lb,VNULL);
	 printf("%lf\n",curvature(qq,0));

	 glutPostRedisplay();
	 
	 break;
	 
       case 101 : printf("GLUT_KEY_UP %d\n",key);  ; 

	 thingy2 += .1;
	 
	 lb->ve[N+S-1+S-1+S-1+3] = thingy2;
	 	   	 
	 lx = LUsolve(LU,pivot,lb,VNULL);
	 printf("%lf\n",curvature(qq,0));

	 glutPostRedisplay();
	 
	 break;

    case 103 : printf("GLUT_KEY_DOWN %d\n",key);  ; 

	 thingy2 -= .1;
	 
	 lb->ve[N+S-1+S-1+S-1+3] = thingy2;
	 	   	 
	 lx = LUsolve(LU,pivot,lb,VNULL);
	 printf("%lf\n",curvature(qq,0));

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
  
  glColor3f(0.0f,1.0f,0.0f);

  for (int i=0;i<S;i++)
    {
      for (int j=0;j<100;j++)
	{
	  glVertex2f(nons*(i*100+j)/(N*100.0),lx->ve[i*4] + lx->ve[i*4+1]*(i+j/100.0)*nons + lx->ve[i*4+2]*pow((i+j/100.0)*nons,2.0) + lx->ve[i*4+3] * pow((i+j/100.0)*nons,3.0));
	}
    }
  
  glEnd();

  double * erk = (double *) calloc(100*S,sizeof(double));

  for (int i=0;i<S;i++)
      for (int j=0;j<100;j++)		  
	erk[i*100+j] = (2*lx->ve[i*4+2] + 6*lx->ve[i*4+3] * (i+j/100.0)*nons) / pow( 1 + pow(lx->ve[i*4+1] + 2*lx->ve[i*4+2]*(i+j/100.0)*nons + 3*lx->ve[i*4+3] * pow((i+j/100.0)*nons,2.0),2.0), 2.0/3.0 );

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
      glVertex2f(nons*(i*100+j)/(N*100.0),erk[i*100+j]); 

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

  gluOrtho2D(-0.1,1.1,0,1.2);

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
  
  if (fscanf(fp,"%d",&N) <1)
    {
      printf("error reading N\n");
      exit(1);
    }

  effort = (double *) calloc(N,sizeof(double));

  for (int i=0;i<N;i++)
    if ( fscanf(fp, "%lf ",&effort[i]) < 1)
      {
	printf("error reading effort %d\n",i);
	exit(1);
      }
    
  fclose(fp);

  N = 2;  // no. data blocks
  M = 10;

  I = M*N;

  double maxeffort = 0;
  for (int i=0;i<N;i++)
    if (effort[i] > maxeffort)
      maxeffort = effort[i];

  for (int i=0;i<N;i++)
    effort[i] /= 1.5*maxeffort;
  
  z = (double *) calloc(I,sizeof(double));
  
  for ( int i=0;i<N;i++)
    for (int m=0;m<M;m++)
      z[i*M+m] = effort[i];

  S=3; // no. splines
  
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
	  /*
	  A->me[i][j*4 + 0] += x[j-counter+k+2] - x[j-counter+k+1];
	  A->me[i][j*4 + 1] += .5*( pow(x[j-counter+k+2],2.0) - pow(x[j-counter+k+1],2.0) );
	  A->me[i][j*4 + 2] += (1.0/3) * ( pow(x[j-counter+k+2],3.0) - pow(x[j-counter+k+1],3.0));
	  A->me[i][j*4 + 3] += (1.0/4) * ( pow(x[j-counter+k+2],4.0) - pow(x[j-counter+k+1],4.0));
	  */
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

  A->me[N+S-1+S-1+S-1][1] = 1;  // first derivative at t=0
  A->me[N+S-1+S-1+S-1+1][2] = 2; // second derivative at t=0

  A->me[N+S-1+S-1+S-1+2][N+S-1+S-1+S-1+3-2] = 1;  // first derivative at t=N
  A->me[N+S-1+S-1+S-1+2][N+S-1+S-1+S-1+3-1] = 2*N; // first derivative at t=N
  A->me[N+S-1+S-1+S-1+2][N+S-1+S-1+S-1+3] = 3*pow(N,2.0);  // first derivative at t=N

  A->me[N+S-1+S-1+S-1+3][N+S-1+S-1+S-1+3-1] = 2;  
  A->me[N+S-1+S-1+S-1+3][N+S-1+S-1+S-1+3] = 6*N;  
  
  //A->me[N+S-1+S-1+S-1+S][N+S-1+S-1+S-1+S-1] = 2; // second derivative at t=N
  //A->me[N+S-1+S-1+S-1+S][N+S-1+S-1+S-1+S] = 6*N; // second derivative at t=N
  
  // second derivative at the join points

  /*
  for (int j=0;j<S-1;j++)
    {

      A->me[N+S-1+S-1+S-1+j][j*4 + 2] = 2;
      A->me[N+S-1+S-1+S-1+j][j*4 + 3] = 6*(x[j+1]);

    }
  */
  

  //m_output(A);
  //exit(1);
  
  lb = v_get(4*S);

  for (int i=0;i<N;i++)
    lb->ve[i] = effort[i];

  //lb->ve[N+S-1+S-1+S-1] = 0;
  //lb->ve[N+S-1+S-1+S-1+1] = 0;
  
  LU = m_get(4*S,4*S);
  LU = m_copy(A,LU);
  pivot = px_get(A->m);
  LUfactor(LU,pivot);

  lx = LUsolve(LU,pivot,lb,VNULL);

  v_output(lx);
  
  thingy = 0; thingy2 = 0;

  
  T = gsl_multimin_fminimizer_nmsimplex2;

  s = NULL;

  // Set initial step sizes to 0.1 
  ss = gsl_vector_alloc (2);
  gsl_vector_set_all (ss, 0.1);

  minex.n = 2;
  minex.f = &curvature;
  minex.params = NULL;

  qq = gsl_vector_alloc (2);
  gsl_vector_set (qq, 0, 0);
  gsl_vector_set (qq, 1, 0);
  
  s = gsl_multimin_fminimizer_alloc(T,2);
  
  gsl_multimin_fminimizer_set (s,&minex,qq,ss);


  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      
      if (status) 
        break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-2);

      if (status == GSL_SUCCESS)
        {
          printf ("converged to minimum at\n");
        }

      printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", 
              iter,
              gsl_vector_get (s->x, 0), 
              gsl_vector_get (s->x, 1), 
              s->fval, size);
    }
  while (status == GSL_CONTINUE && iter < 100);
  
  /*
  
  z = (double *) calloc(I,sizeof(double));
  
  for ( int i=0;i<N;i++)
    for (int m=0;m<M;m++)
      z[i*M+m] = effort[i];

  printf("\n");
  for (int i=0;i<I;i++)
    printf("%lf\n",z[i]);
  
  d1 = calloc(I-1,sizeof(*d1));
  d2 = calloc(I-2,sizeof(*d2));

  qq = gsl_vector_alloc(3);  
  
  reg.nBlocks = N;
  reg.effort = calloc(N,sizeof(double));
  reg.nPoints = M;
  
  for (int i=0;i<N;i++)
    reg.effort[i] = effort[i];
  
  iter = 0;

  for (int i=0;i<I-1;i++)
    d1[i] = z[i+1]-z[i];
  
  for (int i=0;i<I-2;i++)
    d2[i] = fabs(d1[i+1]-d1[i]);
  
  maxcurv = 0;
  mcidx = -1;
  for (int i=0;i<I-2;i++)
    if (d2[i] > maxcurv)
      {
	maxcurv = d2[i];
	mcidx = i;
      }

  pt = mcidx + 1;

  new_pt = pt;

  //printf("%d\n",pt);

      if (pt > 1 && pt < 13)
	{

	  if (pt % 2 != 0)
	    {

	      //printf("odd\n");
	      gsl_vector_set (qq, 0, z[pt-3]);
	      gsl_vector_set (qq, 1, z[pt-1]);
	      gsl_vector_set (qq, 2, z[pt+1]);
	    }
	  else
	    {	   
	      gsl_vector_set (qq, 0, z[pt-2]);
	      gsl_vector_set (qq, 1, z[pt]);
	      gsl_vector_set (qq, 2, z[pt+2]);
	    }       	     	  
	}
      else
	{
	   printf("outside %d\n",pt);
	  exit(1);
	}

  
  
  //for (int iter=0;iter<10000;iter++)
  //  {  
	  
  // }             
  */
  
  
  glutMainLoop();

  return(0);
  
}


	 /*
      status = gsl_multimin_fminimizer_iterate(s);

      printf("%d\n",pt);
      
       if (pt % 2 != 0)
	 {
	   z[pt-3] = gsl_vector_get(s->x,0);
	   z[pt-2] = 2*effort[(pt-3)/2] - z[pt-3];
	   z[pt-1] = gsl_vector_get(s->x,1);
	   z[pt] = 2*effort[(pt-1)/2] - z[pt-1];
	   z[pt+1] = gsl_vector_get(s->x,2);
	   z[pt+2] = 2*effort[(pt+1)/2] - z[pt+1];       
	 }
       else
	 {
	   z[pt-2] = gsl_vector_get(s->x,0);
	   z[pt-1] = 2*effort[(pt-2)/2] - z[pt-2];
	   z[pt] = gsl_vector_get(s->x,1);
	   z[pt+1] = 2*effort[pt/2] - z[pt];
	   z[pt+2] = gsl_vector_get(s->x,2);
	   z[pt+3] = 2*effort[(pt+2)/2] - z[pt+2];       
	 }

      for (int i=0;i<I-1;i++)
	d1[i] = z[i+1]-z[i];
  
      for (int i=0;i<I-2;i++)
	d2[i] = fabs(d1[i+1]-d1[i]);
  
      maxcurv = 0;
      mcidx = -1;
      for (int i=0;i<I-2;i++)	
	if (d2[i] > maxcurv)
	  {
	    maxcurv = d2[i];
	    mcidx = i;
	  }

      new_pt = mcidx + 1;

      //printf("pt: %d\n",pt);

      if (new_pt != pt)
	{

	  //printf("in here\n");

	  pt = new_pt;

      if (pt ==1)
	pt += 1;

      if (pt == 12)
	pt -= 1;

	  
	  if (pt % 2 != 0)
	    {

	      //printf("odd\n");
	      gsl_vector_set (qq, 0, z[pt-3]);
	      gsl_vector_set (qq, 1, z[pt-1]);
	      gsl_vector_set (qq, 2, z[pt+1]);
	    }
	  else
	    {

	      //printf("even\n");	      
	    }
	 	  
	  gsl_vector_set_all (ss, 0.00001);
	  gsl_multimin_fminimizer_set (s,&minex,qq,ss);
	}
          
      if (pt ==1)
	pt += 1;

      if (pt == 12)
	pt -= 1;
	  
	  if (pt % 2 != 0)
	    {

	      //printf("odd\n");
	      gsl_vector_set (qq, 0, z[pt-3]);
	      gsl_vector_set (qq, 1, z[pt-1]);
	      gsl_vector_set (qq, 2, z[pt+1]);
	    }
	  else
	    {

	      //printf("even\n");	      
	      gsl_vector_set (qq, 0, z[pt-2]);
	      gsl_vector_set (qq, 1, z[pt]);
	      gsl_vector_set (qq, 2, z[pt+2]);
	    }

      
	  printf("\n\n");
	  for (int i=0;i<I;i++)
	    printf("%lf ",z[i]);
	  printf("\n");

	 */
