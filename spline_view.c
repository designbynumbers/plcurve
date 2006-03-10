#include <config.h>
#include <plCurve.h>

#ifdef HAVE_STDIO_H
  #include <stdio.h>
#endif
#ifdef HAVE_ASSERT_H
  #include <assert.h>
#endif
#ifdef HAVE_STDLIB_H
  #include <stdlib.h>
#endif

int main(void) {
  plCurve *L;
  plCurve_spline *spL;
  int nv[2] = { 3,4 };
  bool open[2] = { false, true };
  int cc[2] = { 1,4 };
  int cmp, vert;

  L = plCurve_new(2,nv,open,cc);
  assert(L != NULL);

  L->cp[0].vt[0] = plcl_build_vect(0,0,0);
  L->cp[0].vt[1] = plcl_build_vect(1,0,0);
  L->cp[0].vt[2] = plcl_build_vect(0,1,0);
  L->cp[1].vt[0] = plcl_build_vect(2,2,1);
  L->cp[1].vt[1] = plcl_build_vect(2,1,1);
  L->cp[1].vt[2] = plcl_build_vect(2,1,0);
  L->cp[1].vt[3] = plcl_build_vect(2,2,0);
  L->cp[0].clr[0] = plCurve_build_color(1.0,0.0,1.0,1.0);
  L->cp[1].clr[0] = plCurve_build_color(0.0,1.0,0.0,1.0);
  L->cp[1].clr[1] = plCurve_build_color(0.33,1.0,0.0,1.0);
  L->cp[1].clr[2] = plCurve_build_color(0.67,1.0,0.0,1.0);
  L->cp[1].clr[3] = plCurve_build_color(1.0,1.0,1.0,1.0);

  plCurve_fix_wrap(L);

  spL = plCurve_convert_to_spline(L,NULL);
  assert(spL != NULL);

  for (cmp = 0; cmp < L->nc; cmp++) {
    for (vert = -1; vert <= L->cp[cmp].nv; vert++) {
      fprintf(stderr,"check(plcl_vecteq(spL->cp[%d].vt2[%d],\n"
          "plcl_build_vect(%.35g,\n%.35g,\n%.35g)));\n",cmp,vert,
          plcl_M_clist(spL->cp[cmp].vt2[vert]));
    }
  }
  printf("(normalization World none)\n");
  printf("(geometry plCurve \n");
  plCurve_write(stdout,L);
  printf(")\n");
  fflush(stdout);

  plCurve_free(L);
  nv[0] = 3;
  L = plCurve_convert_from_spline(spL,nv);
  
  printf("(geometry spline \n");
  plCurve_write(stdout,L);
  printf(")\n");
  printf("(look-recenter spline) (look-encompass spline)\n");

  exit(EXIT_SUCCESS);
}
