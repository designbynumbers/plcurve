#include <stdio.h>
#include <math.h>

#include "plCurve.h"

int main() {
  plCurve *L;
  int nv[2] = { 100, 2 };
  int open[2] = { FALSE, TRUE };
  int cc[2] = { 100, 1 };
  FILE *outfile;
  int i;
 
  L = plCurve_new(2,nv,open,cc,0,NULL);
  plcl_status_check();
  for (i = 0; i < 100; i++) {
    plCurve_M_set_vert(L,0,i,sin(0.03*i),cos(0.03*i),0);
    L->cp[0].clr[i].r = sin(0.03*i);
    L->cp[0].clr[i].g = 0.0;
    L->cp[0].clr[i].b = cos(0.03*i);
    L->cp[0].clr[i].alpha = 1.0;
  }
  L->cp[1].clr[0].r = 0.0;
  L->cp[1].clr[0].g = 1.0;
  L->cp[1].clr[0].b = 0.0;
  L->cp[1].clr[0].alpha = 1;
  plCurve_M_set_vert(L,1,0,0,0,0);
  plCurve_M_set_vert(L,1,1,0,0,1);
  plCurve_set_constraint(L,0,0,100,PLCL_IN_PLANE,0,0,1,0,0,0);
  plcl_status_check();
  plCurve_set_constraint(L,0,50,1,PLCL_ON_LINE,1,2,3,4,5,6);
  plcl_status_check();
  plCurve_set_constraint(L,1,0,1,PLCL_ON_LINE,3,2,3,4,5,6);
  plcl_status_check();
  plCurve_set_constraint(L,1,1,1,PLCL_ON_LINE,1,2,3,4,5,6);
  plcl_status_check();
  plCurve_set_unconstrained(L,0,75,5);
  plcl_status_check();
  plCurve_set_unconstrained(L,0,25,5);
  plcl_status_check();
  plCurve_set_constraint(L,0,32,1,PLCL_FIXED,3.2,0,3.2,0,0,3.2);
  plcl_status_check();
  plCurve_set_unconstrained(L,0,2,1);
  plcl_status_check();
  for (i = 0; i < 100; i++) {
    printf("Constraint %d: %d -- %lf\n",L->cp[0].vt[i].cst,
      i,plCurve_cst_fix(L,0,i));
    plcl_status_check();
  }
  outfile = fopen("test_cst.vect","w");
  plCurve_write(outfile,L);
  fclose(outfile);
  plcl_status_check();

  return 0;
}
