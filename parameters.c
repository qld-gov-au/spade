// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#include <stdio.h>
#include "parameters.h"

// Attempts to read all parameter configurations from command line
// arguments. Returns TRUE if all parameter configurations have been specified,
// FALSE if any parameter configurations were unspecified.
int parameters_read(Parameters * parameters, int argc, char * argv[]) {
  for(int i = 0; i < parameters->count; i++) {
    if(!parameter_read(parameters->parameter[i], argc, argv)) {
      return 0;
    }
  }
  return 1;
}

// Attempts to read a specific parameter configuration from command line
// arguments. Returns TRUE if the parameter configuration was specified,
// FALSE if the parameter configuration was unspecified.
int parameter_read(Parameter * parameter, int argc, char * argv[]) {
  // Convert "alpha" to "-alpha"
  char arg_name[32];
  snprintf(arg_name, sizeof(arg_name), "%s%s", "-", parameter->name);

  // Convert "alpha" to "-alpha-disabled"
  char arg_name_disabled[32];
  snprintf(arg_name_disabled, sizeof(arg_name_disabled), "%s%s%s", "-", parameter->name, "-disabled");

  // Attempt to find an argument that matches either "-alpha"
  // or "-alpha-disabled" (whichever appears first)
  for(int i = 0; i < argc - 1; i++) {
    char *arg = argv[i];
    Real value;

    // If argument matches "-alpha"
    if(strcmp(arg, arg_name) == 0) {
      if(!parse_argument_value(argc, argv, i + 1, &value)) {
        return 0;
      }
      parameter->value = value;
      parameter->active = TRUE;
      return 1;
    }

    // If argument matches "-alpha-disabled"
    if(strcmp(arg, arg_name_disabled) == 0) {
      if(!parse_argument_value(argc, argv, i + 1, &value)) {
        return 0;
      }
      parameter->value = value;
      parameter->active = FALSE;
      return 1;
    }
  }

  return 0;
}

// Attempts to read a command line argument at the given index as a
// numeric (Real) value. Returns TRUE if the argument could be read,
// FALSE if the argument couldnot be read as a Real or if it was unspecified.
int parse_argument_value(int argc, char *argv[], int i, Real * value) {
  #if REAL == DOUBLE
    // lf
    if(sscanf(argv[i], "%lf", value) != 1) {
  #elif REAL == FLOAT
    // f
    if(sscanf(argv[i], "%f", value) != 1) {
  #elif REAL == LONGDOUBLE
    // Lf
    if(sscanf(argv[i], "%Lf", value) != 1) {
  #endif
    return 0;
  }

  return 1;
}

// Maps a Pararameters struct to a vector of all active parameter values
VEC * parameters_to_vec(Parameters * parameters) {
  // Determine the number of parameters which are active
  int activeParameterCount = 0;
  for(int i = 0; i < parameters->count; i++) {
    if(parameters->parameter[i]->active == TRUE) {
      activeParameterCount++;
    }
  }

  // Load all active parameter values into a theta vector
  VEC *theta = v_get(activeParameterCount);
  int iTheta = 0;
  for(int i = 0; i < parameters->count; i++) {
    if(parameters->parameter[i]->active == TRUE) {
      theta->ve[iTheta++] = parameters->parameter[i]->value;
    }
  }

  return theta;
}
