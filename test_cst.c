#include <stdio.h>
#include <math.h>

#include "plCurve.h"

int main() {
  plCurve *L;
  int nv[2] = { 100, 2 };
  int open[2] = { FALSE, TRUE };
  int cc[2] = { 0, 0 };
  FILE *outfile;
  int i;
 
  L = plCurve_new(2,nv,open,cc,0,NULL);
  plcl_status_check();
  for (i = 0; i < 100; i++) {
    plCurve_M_set_vertex(L,0,i,sin(0.03*i),cos(0.03*i),0);
  }
  plCurve_M_set_vertex(L,1,0,0,0,0);
  plCurve_M_set_vertex(L,1,1,0,0,1);
  plCurve_set_constraint(L,0,0,100,PLCL_IN_PLANE,0,0,1,0,0,0);
  plcl_status_check();
  printf("Ncst: %d\n",L->ncst);
  outfile = fopen("test.out","w");
  plCurve_write(outfile,L);
  fclose(outfile);
  plcl_status_check();

  return 0;
}
