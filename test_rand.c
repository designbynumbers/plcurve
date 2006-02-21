#include "plCurve.h"
#include "stdio.h"

plcl_vector plcl_random_2();

#define POINTS 1000
int main () {
  int i;
  plcl_vector R;

  printf("VECT\n1000 2000 0\n");
  for (i = 0; i < 1000; i++) {
    printf("2 ");
  }
  printf("\n");
  for (i = 0; i < 1000; i++) {
    printf("0 ");
  }
  printf("\n");
  for (i = 0; i < 1000; i++) {
    R = plcl_random_2();
    printf("%lg %lg %lg\n",plcl_M_clist(0.9*R));
    printf("%lg %lg %lg\n",plcl_M_clist(R));
  }

  return 0;
}
