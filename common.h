#ifndef SPADE_COMMON_H
#define SPADE_COMMON_H

#include "meschach/matrix.h"
#include "parameters.h"

#define PTH 1

#define A1 8.588e-5
#define A2 0.00144

typedef struct {
  Real stp;
  Real ftol;
  Real gtol;
  Real xtol;
  Real stpmin;
  Real stpmax;
  int maxfev;
} OptimControl;
    
typedef struct {
  VEC *eff;
  VEC *cat;
  Real **lf; 
  int n; 
  int *t_id; 
  int *t_sz;
  int I,J,S;
  Real k;
} Data;

typedef struct {
  MAT *x,*u,*xhh,*xh,*xn,*uh,*un;
  VEC *Ui,*Uh,*Uhh;
  IVEC *idxi;
} Solve_Core_Args;

typedef struct {
  Data *d;
  VEC * eff;
  Real k;
  int S;
  Solve_Core_Args *core_args;  
  Parameters *parameters;
} Grad_Args;

Real iota1;
Real iota2;
Real phi;
Real eta1;
Real eta2;

Real h;
int J;

void spade_v_output(VEC* vec);

void data_read_ce(char * data_file_name, Data * data, int * N, Real k);

void data_read_lf(char * data_file_name, Data * data, int N, Real k, int minfish);

void optim_control_read(char * optim_file_name, OptimControl * optim);
#endif
