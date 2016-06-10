#include "meschach/matrix.h"

Real iota1=5.2;
Real iota2=0.619;
Real phi=17;
Real eta1=1.703205e-5;
Real eta2=2.9526;

void spade_v_output(VEC* vec) {
  printf("\n\n");
  for(int i = 0; i < vec->dim; i++) {
    printf("%d:\t%Lf\n", i, vec->ve[i]);
  }
}