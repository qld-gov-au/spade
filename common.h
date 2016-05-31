#ifndef SPADE_COMMON_H
#define SPADE_COMMON_H

#include "meschach/matrix.h"

#define PTH 1
#define PLOT 0
#define PLOTDERIV 0
#define PLOTSOLVP 0
#define GCHECK 0

#define A1 8.588e-5
#define A2 0.00144

struct DATA {
  VEC *eff;
  VEC *cat;
  double **lf; 
  int n; 
  int *t_id; 
  int *t_sz;
  int I,J,S;
  double k;
};

typedef struct {
  VEC * theta;
  struct DATA *dataptr;
  VEC * grad;

  MAT * p, *x, *u, *xhh, *xh, *xn, *uh, *un;
  VEC * Ui, *Uh, *Uhh;
  IVEC * idxi;
  VEC * eff;
  double k;
	int S;
} Solve_Args;

float kappa,omega;

double iota1;
double iota2;
double phi;
double eta1;
double eta2;

double h;
int J;

#endif
