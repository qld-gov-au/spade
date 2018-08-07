// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#ifndef SPADE_ARGS_H
#define SPADE_ARGS_H

#include "meschach/matrix.h"

#define PARAMETER_COUNT 6

typedef struct {
  // The gradient function for a given parameter (e.g. grad_alpha)
  void (*grad) (void *args);

  // The initial value for this parameter
  Real value;

  // Whether this parameter should be predicted by the model (TRUE)
  // or if it should retain a fixed value (FALSE).
  int active;
  
  // The current value of the gradient 
  Real gradient;

  // A string identifier for this parameter. This is used when
  // reading values from the command line.
  char * name;
  
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

int parameters_read(Parameters * parameters, int argc, char * argv[]);

VEC * parameters_to_vec(Parameters * parameters);

int parameter_read(Parameter * parameter, int argc, char * argv[]);

int parse_argument_value(int argc, char *argv[], int i, Real * value);
#endif
