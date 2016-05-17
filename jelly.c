#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h>

#include "trackball.c"

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#endif

#define randmax 2146319488.000000
#define pi 3.141593

// mac: gcc -Wall jelly.c -framework GLUT -framework OpenGL -o jelly 
// nix: gcc -Wall jelly.c -lGL -lglut -lGLU -o jelly

// for debugging: ulimit -c unlimited 
// this makes sure a core is dumped
// then add -g to gcc to make sure debugging symbols are loaded
// then: gdb jelly core

//double hx[10000];
//double hy[10000];

int WindWidth, WindHeight;

#define N 200

double s[N];
double xs[N];
double ys[N];

typedef struct {
  double a;       // semi-major axis
  double b;       // semi-minor axis
  double theta;   // rotation angle
  double clx;      // closest x
  double cly;      // closest y
} ellipse;

ellipse ellipsi[N-1]; // = NULL;

double area = 0.05;
double ecce = 0.9;


double purpt_x[100];
double purpt_y[100];
int purpt_n=0;

double foc1_x;
double foc1_y;
double foc2_x;
double foc2_y;

void fMotion(int x, int y)
{
  tbMotion(x, y);
}

void fMouse(int button, int state, int x, int y)
{ 
  tbMouse(button, state, x, y);
}

void fdisplay() {

  int i;

  glColor3f(1,1,1);

  glBegin(GL_POINTS);
 
  for (i=1;i<N;i++) {
    glVertex2f(xs[i-1],xs[i]);
  }

  glEnd();

  glFlush();

  glColor3f(1,0,0);

  glBegin(GL_POINTS);

  for (i=1;i<N;i++) {
    glVertex2f(s[i-1],s[i]);
  }

  glEnd();

  glFlush();

  glColor3f(0,1,0);

  glBegin(GL_LINES);

  for (i=0;i<2;i++) { //(N-1);i++) {
    //glBegin(GL_LINES);
    glVertex2f(s[i],s[i+1]);
    glVertex2f(ellipsi[i].clx,ellipsi[i].cly);
    //printf("clx %f", ellipsi[i].clx);
    //printf("cly %f", ellipsi[i].cly);
    //glEnd();
  }

  glEnd();

  glFlush();
 
  glBegin(GL_POINTS);
  glColor3f(1,0,1);

  //glVertex2f(foc1_x,foc1_y);
  //glVertex2f(foc2_x,foc2_y);
  for (i=0;i<purpt_n;i++) {
    glVertex2f(purpt_x[i],purpt_y[i]);
  }

  glEnd();

  glFlush();

  glColor3f(0,0,1);

  glBegin(GL_POINTS);

  double epst = 2.f*pi/100.f;
  double tt=0.f;
  i=0;
  for (i=0;i<2;i++){ //(N-1);i++) {
    int j;
    for (j=0;j<100;j++) {
      double theta = ellipsi[i].theta;
      double aa = ellipsi[i].a;
      double bb = ellipsi[i].b;
      glVertex2f(s[i]+aa*cos(tt)*cos(theta)-bb*sin(tt)*sin(theta),s[i+1]+aa*cos(tt)*sin(theta)+bb*sin(tt)*cos(theta));  // thanks wikipedia!
      tt = tt+epst;
    }
  }

  glEnd();
  
}

void display(void) 
{
  glClear(GL_COLOR_BUFFER_BIT); // | GL_DEPTH_BUFFER_BIT);

#ifdef THREEDEE
  glPushMatrix();
  tbMatrix();
#endif
  //glutSolidTorus(1.0, 3.0, 3.0, 12.0);
  
  fdisplay();
  
#ifdef THREEDEE
  glPopMatrix();
#endif

   glFlush ();
   glutSwapBuffers();
}

void idle(void)
{
  glutPostRedisplay();
}

void Reshape(int width, int height)
{

#ifdef THREEDEE  
  tbReshape(width, height);
#endif

  glViewport(0, 0, width, height);
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
#ifdef THREEDEE
  gluPerspective(60.0, (GLfloat)height / (GLfloat)width, 1.0, 128.0);
#else
  gluOrtho2D(-1.2,1.2,-1.2,1.2);
#endif
  glMatrixMode(GL_MODELVIEW);
#ifdef THREEDEE
  glLoadIdentity();
  glTranslatef(0.0, 0.0, -10.0);
#endif

  WindWidth = width;
  WindHeight = height;
}

void init(void) 
{
   glClearColor (0.0, 0.0, 0.0, 0.0);
   glShadeModel (GL_FLAT);
}

void timf(int value)
{
  glutPostRedisplay();
  glutTimerFunc(1, timf, 0);
}

int main(int argc, char* argv[]) 
{
  //ellipsi = malloc((N-1) * sizeof( *ellipse ));
  //ellipse *elptr = (ellipse
  //if (NULL == ellipsi) {
  //  printf("ellipsi memory allocation failed\n");
  //  return 1;
  //}

  //void function(const double eccentricity, const double areaa, double &ta, double &tb) {

  //unsigned int iseed = (unsigned int)time(NULL);
  //srand (iseed);
  srand (5);

  int i;                                                                               
  for (i=0; i<N; i++) { ys[i] = rand() / randmax; }

  double a = 1.85;

  xs[0] = .2;
  s[0] = xs[0] + (ys[0]-.5)*.05;

  for (i=1; i<N; i++) {
    //xs[i] = a * xs[i-1] - a * pow(xs[i-1],2.f);
    xs[i] = 1.f - a * pow(xs[i-1],2.f);
    //printf("%f\n",xs[i]);
    s[i] = xs[i] + (ys[i]-.5)*.15;
  } 

  int e;
  for (e=0;e<(N-1);e++) {
    
    ellipsi[e].a = pow( pow(area/pi,2.f) / ( 1.f - pow(ecce,2.f) ) , 1.f/4.f );  // major axis from eccentricity and area
    ellipsi[e].b = area / (pi * ellipsi[e].a);
  
    // equation of the manifold: y = 1 - ax^2. equation of the distance (s1-x)^2 + (s2-f(x))^2
    // take derivative wrt x and root it
    // maxima output:
    /*
                        3   3        2       3    2        3       2
      x = expt(sqrt(16 a  s2  + (24 a  - 48 a ) s2  + (48 a  - 48 a  + 12 a) s2
      
            2   2       3       2                 3/2  3     s1   1      
      + 27 a  s1  - 16 a  + 24 a  - 12 a + 2)/(4 3    a ) + ----, -)
                                                               2  3
                                                            4 a
      
                                2               3   3        2       3    2        3       2
       - (a (2 s2 - 2) + 1)/(6 a  expt(sqrt(16 a  s2  + (24 a  - 48 a ) s2  + (48 a  - 48 a  + 12 a) s2
 
             2   2       3       2                 3/2  3     s1   1
       + 27 a  s1  - 16 a  + 24 a  - 12 a + 2)/(4 3    a ) + ----, -)) 
                                                                2  3
                                                             4 a
    */

    //prelims
    double a2 = pow(a,2.f);
    double s1 = s[e];
    double s2 = s[e+1];
    double p = (1.f/(2.f*a2)) + ((s2-1.f) / a);
    double q = s1 / (2.f*a2);

    double complex nthingy1a = -.5 - .5*sqrt(3.f)*I;
    double complex nthingy1b = -.5 + .5*sqrt(3.f)*I;
    double complex nprethingy2 = (27.f*pow(q,2.f) + 4.f*pow(p,3.f)) + 0.f*I;
    double complex exppo1 = .5 + 0.f*I;
    double complex nprprethingy2 = cpow(nprethingy2,exppo1);
    double complex quo1 = 2.f*pow(3.f,3.f/2.f) + 0.0f;
    double complex anopart = .5*q + 0.f*I;
    double complex pr1 = nprprethingy2 / quo1;
    double complex inner = pr1 + anopart;
    double complex exppo2 = 1.f/3.f + 0.f*I;
    double complex nthingy2 = cpow(inner,exppo2);    
    double complex prfirstrt = nthingy1a*nthingy2;
    double complex pcomplx = p + 0.f*I;
    double complex thrcomplx = 3.f + 0.f*I;
    double complex numrtr = pcomplx * nthingy1b;
    double complex demntr = thrcomplx * nthingy2;
    double complex subtr = numrtr / demntr;
    double complex firstrt = prfirstrt - subtr;
    double complex prfirstrt2 = nthingy1b*nthingy2;
    double complex numrtr2 = pcomplx * nthingy1a;
    double complex subtr2 = numrtr2 / demntr;
    double complex secndrt = prfirstrt2 - subtr2;
    double complex subtr3 = pcomplx / demntr;
    double complex thirdrt = nthingy2 - subtr3;

    printf("firstrt imag %f\n", cimag(firstrt));
    printf("secndrt imag %f\n", cimag(secndrt));
    printf("thirdrt imag %f\n", cimag(thirdrt));

    double numwewant;

    if (!(fabs(cimag(thirdrt))>0.f)){
      numwewant = creal(thirdrt);
      printf("third\n");
    }
    else if ((fabs(cimag(secndrt))<0.000001)){ 
      if (cimag(secndrt) >= 0.f) {numwewant = creal(secndrt); printf("second\n"); }
      else if ((fabs(cimag(firstrt)) < 0.000001) && (cimag(firstrt) >= 0.f)) {numwewant = creal(firstrt); printf("first (in second)\n");  // ok three real roots - choose closest
       double dist1 = pow((s1 - creal(firstrt)),2.f) + pow(s2 - (1.f - a*pow(creal(firstrt),2.f)),2.f); 
       double dist2 = pow((s1 - creal(secndrt)),2.f) + pow(s2 - (1.f - a*pow(creal(secndrt),2.f)),2.f); 
       double dist3 = pow((s1 - creal(thirdrt)),2.f) + pow(s2 - (1.f - a*pow(creal(thirdrt),2.f)),2.f);

       if  (dist1 < dist2) {
         if (dist1 < dist3) {
           numwewant = creal(firstrt);
         } else {
           numwewant = creal(thirdrt);
         }
       } else if (dist2 < dist3) {
         numwewant = creal(secndrt);
       } else {
         numwewant = creal(thirdrt);
       }
    }
      else printf("fuk2\n");
    } else if ((fabs(cimag(firstrt))<0.000001)){
      if (cimag(firstrt) >= 0.f) {numwewant = creal(firstrt); printf("first (in first)\n"); }
    } else
      printf("fuck\n");

    double clx =  numwewant;
    double cly = 1.f - a*pow(clx,2.f);

    printf("a2 %f\n",a2);
    printf("s1 %f\n",s1);
    printf("s2 %f\n",s2);

    ellipsi[e].clx = clx;
    ellipsi[e].cly = cly;

    ellipsi[e].theta = atan( (clx-s1) / (s2 - cly) );  // check sign

  }

  int hashset[N-1][N-1] = { { 0 } };

  // determinism code
  int k=0;
  int j;
  for (i=0;i<(N-1);i++) {
    double this_x = s[i];
    double this_y = s[i+1];
    double this_a = ellipsi[i].a;
    double this_theta = ellipsi[i].theta;

    // construct the focii:
    double hypot = this_a*ecce;
    foc1_x = cos(this_theta) * hypot;
    foc1_y = sin(this_theta) * hypot;
    foc2_x = -foc1_x;
    foc2_y = -foc1_y;

    foc1_x = foc1_x + this_x;
    foc1_y = foc1_y + this_y;
    foc2_x = foc2_x + this_x;
    foc2_y = foc2_y + this_y;

    for (j=0;j<(N-1);j++) {

      if (!(i==j)) {
      
        double cand_x = s[j];
        double cand_y = s[j+1];
     
        double dist2foc1 = sqrt( pow(cand_x - foc1_x,2.f) + pow(cand_y - foc1_y,2.f) );
        double dist2foc2 = sqrt( pow(cand_x - foc2_x,2.f) + pow(cand_y - foc2_y,2.f) );
      
        if ((dist2foc1 + dist2foc2) < 2*this_a) { // thanks dr math!

          hashset[i][j] = 1;

          if (i < 3) {
            purpt_x[k] = cand_x;
            purpt_y[k] = cand_y;
            //purpt_j[k] = j;
            purpt_n = purpt_n + 1;
            k = k+1;
            //printf("yes %d\n",j);
          }
            
        }
      }
    }
    //printf("end\n");
  }

  printf("ellipsi[0].clx %f\n",ellipsi[0].clx);
  
  printf("ellipsi[0].cly %f\n",ellipsi[0].cly);
  printf("s[0] %f\n", s[0]);
  printf("s[1] %f\n", s[1]);

  printf("ellipsi[0].theta %f\n",ellipsi[0].theta);

  printf("set intersection:\n");

  int thescore[N-1] = {0};

  //int j;
  for (i=1;i<N;i++) {
    //printf("hash1 %d ",hashset[0][i]);
    //printf("hash2 %d ",hashset[1][i]);
    for (j=0;j<(N-1);j++) {

      if ((hashset[j][i-1]==1) && (hashset[j+1][i]==1)) {
        thescore[j] = thescore[j]+1;
      }
    }
  }  
    
  for (j=0;j<(N-1);j++) {
    printf("%d ",thescore[j]);
  }


  WindWidth = 1000;
  WindHeight = 1000;

  glutInit(&argc, (char **)argv);
    
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB); //GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowSize(1000, 1000);

  glutCreateWindow("neo");
  glutDisplayFunc(display);
  glutReshapeFunc(Reshape);

#ifdef THREEDEE
  glutMouseFunc(fMouse);
  glutMotionFunc(fMotion);
  //glutIdleFunc(idle);

  init();

  tbInit(GLUT_LEFT_BUTTON);
  tbAnimate(GL_TRUE);

  glutTimerFunc(40, timf, 0); // Set up timer for 40ms, about 25 fps
#endif

  glutMainLoop();

  return(0);

}

    /*printf("area %f\n", area);
    printf("ecce %f\n", ecce);
    printf("pi %f\n", pi);
    printf("pow(area/pi,2.f) %f\n", pow(area/pi,2.f));
    printf("1 -e2 %f\n", 1.f-pow(ecce,2.f));
    printf("penu %f\n", pow(area/pi,2.f) / ( 1.f - pow(ecce,2.f) ) );
    printf("ellipsi a %f\n", ellipsi[e].a);
    printf("ellipsi b %f\n", ellipsi[e].b);*/

    //ecceareatoab(ecce,area,&ta,&tb);
    //double ta,tb;

  /*int l;
  double ha=0.9;
  double hb=-0.6;
  double hc=2.0;
  double hd=0.5;
  hx[0] = 0.0;
  hy[0] = 0.5;
  for (l=1;l<10000;l++) {
    hx[l] = pow(hx[l-1],2.0) - pow(hy[l-1],2.0) + ha*hx[l-1] + hb*hy[l-1];
    hy[l] = 2.0*hx[l-1]*hy[l-1] + hc*hx[l-1] + hd*hy[l-1];
    }*/

  //#endif
