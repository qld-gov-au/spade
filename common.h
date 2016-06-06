#ifndef SPADE_COMMON_H
#define SPADE_COMMON_H

#include "meschach/matrix.h"
#include "parameters.h"

#define PTH 1

#define A1 8.588e-5
#define A2 0.00144

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
