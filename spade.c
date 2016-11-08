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
double *z;
double *d1;
double *d2;
double *effort;

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
  
double maxcurv;
int mcidx;

double thingy,thingy2;

int pt,new_pt;

int WindWidth, WindHeight;
int window;

double curvature0(
		 double x,
		 void *param

		 )

{

  double * p = (double *) param;

  return fabs(2*p[2] + 6*p[3] * x) / pow( 1 + pow(p[1] + 2*p[2]*x + 3*p[3] * pow(x,2.0),2.0), 2.0/3.0 );

}
	       

void process_Normal_Keys(int key, int x, int y) 
{
     switch (key) 
    {    
       case 27 :      break;
       case 100 :

	 printf("GLUT_KEY_LEFT %d\n",key);

	 thingy += .1;
	 
	 lb->ve[N+N-1+N-1+N-1+2] = thingy;
	   
	 lx = LUsolve(LU,pivot,lb,VNULL);

	 //         v_output(lx);
	 
	 glutPostRedisplay();
	 
	 break;
	 
       case 102: printf("GLUT_KEY_RIGHT %d\n",key);

	 thingy -= .1;
	 
	 lb->ve[N+N-1+N-1+N-1+2] = thingy;
	 	   	 
	 lx = LUsolve(LU,pivot,lb,VNULL);

	 glutPostRedisplay();
	 
	 break;
	 
       case 101 : printf("GLUT_KEY_UP %d\n",key);  ; 

	 thingy2 += .1;
	 
	 lb->ve[N+N-1+N-1+N-1+1] = thingy2;
	 	   	 
	 lx = LUsolve(LU,pivot,lb,VNULL);

	 glutPostRedisplay();
	 
	 break;

    case 103 : printf("GLUT_KEY_DOWN %d\n",key);  ; 

	 thingy2 -= .1;
	 
	 lb->ve[N+N-1+N-1+N-1+1] = thingy2;
	 	   	 
	 lx = LUsolve(LU,pivot,lb,VNULL);

	 glutPostRedisplay();

    	 break;

    }

}

void display(void)
{
  
  glClear(GL_COLOR_BUFFER_BIT);

  glPointSize(5.0f);
  
  glBegin(GL_POINTS);
  
  glColor3f(0.0f,0.0f,1.0f);

  for (int i=0;i<I;i++)
    glVertex2f(i/((float)(I-.5)),z[i]);
  
  glEnd();
  
  glBegin(GL_POINTS);
  
  glColor3f(0.0f,1.0f,0.0f);


  for (int i=0;i<N;i++)
    {
      for (int j=0;j<100;j++)
	{
	  glVertex2f((i*100+j)/(N*100.0),lx->ve[i*4] + lx->ve[i*4+1]*(i+j/100.0) + lx->ve[i*4+2]*pow((i+j/100.0),2.0) + lx->ve[i*4+3] * pow((i+j/100.0),3.0));
	}
    }
  
  glEnd();

  double * erk = (double *) calloc(100*N,sizeof(double));

  for (int i=0;i<N;i++)
      for (int j=0;j<100;j++)		  
	erk[i*100+j] = fabs(2*lx->ve[i*4+2] + 6*lx->ve[i*4+3] * (i+j/100.0)) / pow( 1 + pow(lx->ve[i*4+1] + 2*lx->ve[i*4+2]*(i+j/100.0) + 3*lx->ve[i*4+3] * pow((i+j/100.0),2.0),2.0), 2.0/3.0 );

  double maxc = 0;
  for (int i=0;i<100*N;i++)
    if (erk[i] > maxc)
      maxc = erk[i];

  double minc = 0;
  for (int i=0;i<100*N;i++)
    if (erk[i] > minc)
      minc = erk[i];
  
  glBegin(GL_POINTS);
  
  glColor3f(1.0f,0.0f,0.0f);

  for (int i=0;i<N;i++)
    for (int j=0;j<100;j++)	
      glVertex2f((i*100+j)/200.0,erk[i*100+j]); 

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

  N = 3;
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

  A = m_get(4*N,4*N);
    
  for (int i=0;i<N;i++)
    {
    
      A->me[i][i*4 + 0] = 1;
      A->me[i][i*4 + 1] = .5*( pow(i+1,2.0) - pow(i,2.0) );
      A->me[i][i*4 + 2] = (1.0/3) * ( pow(i+1,3.0) - pow(i,3.0));
      A->me[i][i*4 + 3] = (1.0/4) * ( pow(i+1,4.0) - pow(i,4.0));
    
    }

  for (int i=0;i<N-1;i++)
    {

      A->me[N+i][i*4 + 0] = 1;
      A->me[N+i][i*4 + 1] = i+1;
      A->me[N+i][i*4 + 2] = pow(i+1,2.0);
      A->me[N+i][i*4 + 3] = pow(i+1,3.0);
      A->me[N+i][i*4 + 4] = -1;
      A->me[N+i][i*4 + 5] = -(i+1);
      A->me[N+i][i*4 + 6] = -pow(i+1,2.0);
      A->me[N+i][i*4 + 7] = -pow(i+1,3.0);
      
    }

  for (int i=0;i<N-1;i++)
    {

      A->me[N+N-1+i][i*4 + 1] = 1;
      A->me[N+N-1+i][i*4 + 2] = 2*(i+1);
      A->me[N+N-1+i][i*4 + 3] = 3*pow(i+1,2.0);
      A->me[N+N-1+i][i*4 + 5] = -1;
      A->me[N+N-1+i][i*4 + 6] = -2*(i+1);
      A->me[N+N-1+i][i*4 + 7] = -3*pow(i+1,2.0);
      
    }

  for (int i=0;i<N-1;i++)
    {

      A->me[N+N-1+N-1+i][i*4 + 2] = 2;
      A->me[N+N-1+N-1+i][i*4 + 3] = 6*(i+1);
      A->me[N+N-1+N-1+i][i*4 + 6] = -2;
      A->me[N+N-1+N-1+i][i*4 + 7] = -6*(i+1);
    }
    
  A->me[N+N-1+N-1+N-1][1] = 1;  // first derivative at t=0

  A->me[N+N-1+N-1+N-1+1][2] = 2; // second derivative at t=0

  A->me[N+N-1+N-1+N-1+2][N+N-1+N-1+N-1+1] = 2; // second derivative at t=N
  A->me[N+N-1+N-1+N-1+2][N+N-1+N-1+N-1+2] = 6*N; // second derivative at t=N
   
  lb = v_get(4*N);

  for (int i=0;i<N;i++)
    lb->ve[i] = effort[i];
  
  LU = m_get(4*N,4*N);
  LU = m_copy(A,LU);
  pivot = px_get(A->m);
  LUfactor(LU,pivot);

  
  lx = LUsolve(LU,pivot,lb,VNULL);

  v_output(lx);
  
  thingy = 0; thingy2 = 0;

  T = gsl_multimin_fminimizer_nmsimplex2;

  s = NULL;

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

  
  // Set initial step sizes to 0.1 
  ss = gsl_vector_alloc (3);
  gsl_vector_set_all (ss, 0.0001);

  minex.n = 3;
  minex.f = &regularise_new;
  minex.params = (void *) &reg;

  s = gsl_multimin_fminimizer_alloc(T,3);
  gsl_multimin_fminimizer_set (s,&minex,qq,ss);
  
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
	      gsl_vector_set (qq, 0, z[pt-2]);
	      gsl_vector_set (qq, 1, z[pt]);
	      gsl_vector_set (qq, 2, z[pt+2]);
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
