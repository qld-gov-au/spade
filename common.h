#ifndef SPADE_COMMON_H
#define SPADE_COMMON_H

#include "meschach/matrix.h"

#define PTH 1

#define A1 8.588e-5
#define A2 0.00144
#define PARAMETER_COUNT 6

typedef struct {
  double stp;
  double ftol;
  double gtol;
  double xtol;
  double stpmin;
  double stpmax;
  int maxfev;
} OptimControl;
    
typedef struct {
  VEC *eff;
  VEC *cat;
  double **lf; 
  int n; 
  int *t_id; 
  int *t_sz;
  int I,J,S;
  double k;
} Data;

typedef struct {
  // The gradient function for a given parameter (e.g. grad_alpha)
  void (*grad) (void *args);

  // The initial value for this parameter
  Real value;

  // Whether this parameter should be predicted by the model (TRUE)
  // or if it should retain a fixed value (FALSE).
  int active;
  
  // The current value of the gradient 
  double gradient;
  
} Parameter;

typedef struct {
  Parameter alpha;
  Parameter beta;
  Parameter gamma;
  Parameter iota;
  Parameter kappa;
  Parameter omega;
  Parameter * parameter[PARAMETER_COUNT];
  int count;
} Parameters;

typedef struct {
  MAT *x,*u,*xhh,*xh,*xn,*uh,*un;
  VEC *Ui,*Uh,*Uhh;
  IVEC *idxi;
} Solve_Core_Args;

typedef struct {
  Data *d;
  VEC * eff;
  double k;
  int S;
  Solve_Core_Args *core_args;  
  Parameters *parameters;
} Grad_Args;

double iota1;
double iota2;
double phi;
double eta1;
double eta2;

double h;
int J;

#endif
