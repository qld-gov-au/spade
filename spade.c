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

double pfake;

int Nparams1;
int Nparamsb;

double *x,*x2;
double *N2x;

VEC *lx,*lx2;

VEC *lb,*lb2;
MAT *A, *LU, *A2, *LU2;
PERM *pivot, *pivot2;

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

int pt,new_pt;

int WindWidth, WindHeight;
int window;

void test0 (VEC *, double *, double *, double);
void test1 (VEC *, double *, double *, double);

void fakeeffort(VEC *, double *, double *, int, double *);
MAT *generateA(MAT *,double *,double *,int);
double curvature(double ,VEC *);

double maxCurvature(gsl_vector *,void *);

void processNormalKeys(int, int, int);

void process_Normal_Keys(int key, int xlah, int ylah) 
{
     switch (key) 
    {    
       case 27 :      break;
       case 100 :
	 
	 break;
	 
       case 102: printf("GLUT_KEY_RIGHT %d\n",key);

	 break;
	 
    case 101 : //printf("GLUT_KEY_UP %d\n",key); 
	 
	 break;

    case 103 : 

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

  for (int i=0;i<S;i++)
    {
      for (int j=0;j<100;j++)
	{
	  double xx = x[i]+j*(x[i+1]-x[i])/100.0;
	  
	  glVertex2f(xx/N,lx->ve[i*3] + lx->ve[i*3+1]*xx + lx->ve[i*3+2]*pow(xx,2.0));
	}
    }
  
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

  S=3; // no. splines
  
  A = m_get(3*S,3*S);

  x = (double *) calloc(S+1,sizeof(double));

  for (int i=0;i<=S;i++) {
    x[i] = i * ( (double)N / (double)S);
    printf("%f\n",x[i]);
  }
  
  int j=0;

  for (int i=0;i<N;i++)
    {
      A->me[i][i*3 + 0] += x[i+1] - x[i];
      A->me[i][i*3 + 1] += .5*( pow(x[i+1],2.0) - pow(x[i],2.0) );
      A->me[i][i*3 + 2] += (1.0/3) * ( pow(x[i+1],3.0) - pow(x[i],3.0));
    }
  
  for (int j=0;j<S-1;j++)
    {

      A->me[N+j][j*3 + 0] = 1;
      A->me[N+j][j*3 + 1] = x[j+1];
      A->me[N+j][j*3 + 2] = pow(x[j+1],2.0);
      A->me[N+j][j*3 + 3] = -1;
      A->me[N+j][j*3 + 4] = -(x[j+1]);
      A->me[N+j][j*3 + 5] = -pow(x[j+1],2.0);
      
    }

  for (int j=0;j<S-1;j++)
    {

      A->me[N+S-1+j][j*3 + 1] = 1;
      A->me[N+S-1+j][j*3 + 2] = 2*(x[j+1]);
      A->me[N+S-1+j][j*3 + 4] = -1;
      A->me[N+S-1+j][j*3 + 5] = -2*(x[j+1]);
      
    }

  for (int j=0;j<S-1;j++)
    { 
      A->me[N+S-1+S-1+j][j*3 + 2] = 2;
      A->me[N+S-1+S-1+j][j*3 + 5] = -2;
    }
  
  
  lb = v_get(3*S);

  for (int i=0;i<N;i++)
    lb->ve[i] = effort[i];

  LU = m_get(3*S,3*S);
  LU = m_copy(A,LU);
  pivot = px_get(A->m);
  LUfactor(LU,pivot);

  lx = LUsolve(LU,pivot,lb,VNULL);
  
  glutMainLoop();

  return(0);
  
}

