// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#include "arg.h"
#include "meschach/matrix.h"

// Attempts to read in an argument with the specified name into the specified value.
// Returns FALSE if the argument was missing or invalid. Returns TRUE if the argument was read.
int arg_read_int(char * name, int * value, int argc, char * argv[]) {
  // Convert "arg" to "-arg"
  char arg_name[32];
  snprintf(arg_name, sizeof(arg_name), "%s%s", "-", name);

  // Attempt to find an argument that matches "-arg"
  for(int i = 0; i < argc - 1; i++) {
    char *arg = argv[i];

    // If argument matches "-arg"
    if(strcmp(arg, arg_name) == 0) {
      if(sscanf(argv[i + 1], "%d", value) != 1) {
        return 0;
      }
      return 1;
    }
  }
  return 0;
}

// Attempts to read in an argument with the specified name into the specified value.
// Returns FALSE if the argument was missing or invalid. Returns TRUE if the argument was read.
int arg_read_real(char * name, Real * value, int argc, char * argv[]) {
  // Convert "arg" to "-arg"
  char arg_name[32];
  snprintf(arg_name, sizeof(arg_name), "%s%s", "-", name);

  // Attempt to find an argument that matches "-arg"
  for(int i = 0; i < argc - 1; i++) {
    char *arg = argv[i];

    // If argument matches "-arg"
    if(strcmp(arg, arg_name) == 0) {
      #if REAL == DOUBLE
        // lf
        if(sscanf(argv[i + 1], "%lf", value) != 1) {
      #elif REAL == FLOAT
        // f
        if(sscanf(argv[i + 1], "%f", value) != 1) {
      #elif REAL == LONGDOUBLE
        // Lf
        if(sscanf(argv[i + 1], "%Lf", value) != 1) {
      #endif
        return 0;
      }
      return 1;
    }
  }
  return 0;
}

// Attempts to read in an argument with the specified name into the specified value.
// Returns FALSE if the argument was missing or invalid. Returns TRUE if the argument was read.
int arg_read_string(char * name, char ** value, int argc, char * argv[]) {
  // Convert "arg" to "-arg"
  char arg_name[32];
  snprintf(arg_name, sizeof(arg_name), "%s%s", "-", name);

  // Attempt to find an argument that matches "-arg"
  for(int i = 0; i < argc - 1; i++) {
    char *arg = argv[i];

    // If argument matches "-arg"
    if(strcmp(arg, arg_name) == 0) {
      *value = argv[i + 1];
      return 1;
    }
  }
  return 0;
}
