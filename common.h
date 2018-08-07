// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#ifndef SPADE_COMMON_H
#define SPADE_COMMON_H

#include "meschach/matrix.h"
#include "parameters.h"

// Whether to use MESCHACH for the solvers
#define MESCHACH 0

// Whether to use a quarter step in the numerics
#define QUARTER 0

// Whether to use the new objective funtion or not
#define NEWOBJ 0

// SGNM: Selection Grid Nodes Method
// 
// Determines whether to create full-size matrices
// or to dynamically add and remove rows from the
// matrices.
//
// Set to 0 to use full-size matrices
// Set to 1 to dynamically add and remove rows
#define SGNM 0

// Whether to run grad_* methods concurrently
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
  int nK;          // number knots
  int B;
  VEC *knots_e;    // knots effort
  VEC *knots_c;    // knots catch
  VEC *splcoef_e;    // knots effort
  VEC *splcoef_c;    // knots catch
  int Nlf;
  VEC *ln;
  VEC *tl;
} NewData;

typedef struct {
  VEC *eff;
  VEC *cat;
  Real **lf; 
  int n; 
  int *t_id; 
  int *t_sz;
  int I,J,S;
  Real k;
  int Y; // Number of years of input data
  
} Data;

typedef struct {
  Data *d;
  VEC * eff;
  Real k;
  int S;
  Parameters *parameters;
} Grad_Args_No_MESCHACH;

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
int interactive_mode_requested;

void request_interactive_mode(int a);

void spade_v_output(VEC* vec);

void data_read_ce(char * data_file_name, Data * data, int * N, Real k);

void data_read_lf(char * data_file_name, Data * data, int N, Real k, int minfish);

void data_read_ce_new(char * data_file_name, NewData * newdata);

void data_read_lf_new(char * data_file_name, NewData * newdata);

void optim_control_read(char * optim_file_name, OptimControl * optim);
#endif
