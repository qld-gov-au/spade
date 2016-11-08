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
  VEC *knots_e;      // knots effort
  VEC *knots_c;      // knots catch
  VEC *splcoef_e;    // knots effort
  VEC *splcoef_c;    // knots catch
  int Nlf;
  VEC *ln;
  VEC *tl;
} NewData;

typedef struct {
  int I;
  int J;
  Real *cat;
  Real *eff;
  Real **p;
  Real *Qp;
} Da;

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
  Da *d;
  Real k;
  int N;
  Parameters *parameters;
} Grad_Args_2;

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


extern Real iota1;
extern Real iota2;
extern Real phi;
extern Real eta1;
extern Real eta2;

extern Real h;
extern int J;
extern int interactive_mode_requested;
extern Real k;

extern Da d;

extern int * idx;

extern Real A1; 
extern Real A2;

void request_interactive_mode(int a);

void spade_v_output(VEC* vec);

void data_read_ce(const char * data_file_name, Data * data, int * N, Real k);

void data_read_lf(const char * data_file_name, Data * data, int N, Real k, int minfish);

void data_read(const char * data_file_name);

void data_read_fast(const char * data_file_name);

void data_read_ce_new(const char * data_file_name, NewData * newdata);

void data_read_lf_new(const char * data_file_name, NewData * newdata);

void optim_control_read(const char * optim_file_name, OptimControl * optim);
#endif
