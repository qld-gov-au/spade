#ifndef SPADE_ARG_H
#define SPADE_ARG_H

#include "meschach/matrix.h"

int arg_read_int(char * name, int * value, int argc, char * argv[]);

int arg_read_real(char * name, Real * value, int argc, char * argv[]);

int arg_read_string(char * name, char * value, int argc, char * argv[]);

#endif