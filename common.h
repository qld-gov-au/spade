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
  Real stpMemax;
  int Memaxfev;
} OptimControl;

typedef struct {
  int nK;          // number knots
  int B;
  MeVEC *knots_e;    // knots effort
  MeVEC *knots_c;    // knots catch
  MeVEC *splcoef_e;    // knots effort
  MeVEC *splcoef_c;    // knots catch
  int Nlf;
  MeVEC *ln;
  MeVEC *tl;
} NewData;

typedef struct {
  MeVEC *eff;
  MeVEC *cat;
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
  MeVEC * eff;
  Real k;
  int S;
  Parameters *parameters;
} Grad_Args_No_MESCHACH;

typedef struct {
  MeMAT *x,*u,*xhh,*xh,*xn,*uh,*un;
  MeVEC *Ui,*Uh,*Uhh;
  IMeVEC *idxi;
} Solve_Core_Args;

typedef struct {
  Data *d;
  MeVEC * eff;
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

void spade_v_output(MeVEC* vec);

void data_read_ce(const char * data_file_name, Data * data, int * N, Real k);

void data_read_lf(const char * data_file_name, Data * data, int N, Real k, int minfish);

void data_read_ce_new(const char * data_file_name, NewData * newdata);

void data_read_lf_new(const char * data_file_name, NewData * newdata);

void optim_control_read(const char * optim_file_name, OptimControl * optim);
#endif
