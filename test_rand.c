#include "plCurve.h"
#include "stdio.h"

#define POINTS 1000
int main () {
  int i;
  plcl_vector R;

  printf("VECT\n%d %d 0\n",POINTS,2*POINTS);
  for (i = 0; i < POINTS; i++) {
    printf("2 ");
  }
  printf("\n");
  for (i = 0; i < POINTS; i++) {
    printf("0 ");
  }
  printf("\n");
  for (i = 0; i < POINTS; i++) {
    R = plcl_random_vect();
    printf("%lg %lg %lg\n",plcl_M_clist(0.9*R));
    printf("%lg %lg %lg\n",plcl_M_clist(R));
  }

  return 0;
}
