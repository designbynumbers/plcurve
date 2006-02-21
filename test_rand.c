#include "plCurve.h"

plcl_vector plcl_random_2();

int main () {
  int i;
  plcl_vector R;

  for (i = 0; i < 1000000; i++) {
    R = plcl_random_vect();
    R = plcl_random_2();
  }

  return 0;
}
