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
  FILE *infile;
  plCurve *L;
  plCurve_spline *spL;
  int cmp, vert;
  plcl_vector v[2];
  plCurve_color clr = { 1.0, 0.5, 0.5, 1.0 };
  int err_num;
  char err_str[80];
  int nv[1] = { 3 };
  bool open[1] = { false };
  int cc[1] = { 3 };

  L = plCurve_new(1,nv,open,cc);
  assert(L != NULL);

  L->cp[0].vt[0] = plcl_build_vect(1,0,0);
  L->cp[0].vt[1] = plcl_build_vect(0,1,0);
  L->cp[0].vt[2] = plcl_build_vect(0,0,0);
  L->cp[0].clr[0] = plCurve_build_color(0.0,1.0,0.0,1.0);
  L->cp[0].clr[1] = plCurve_build_color(1.0,1.0,0.0,1.0);
  L->cp[0].clr[2] = plCurve_build_color(1.0,0.0,0.0,1.0);

  plCurve_fix_wrap(L);

  L->cp[0].cc = 1;
  L->cp[0].clr[0] = plCurve_build_color(0.0,0.0,1.0,1.0);

  spL = plCurve_convert_to_spline(L,NULL);
  assert(spL != NULL);

  L->cp[0].cc = 3;
  L->cp[0].clr[0] = plCurve_build_color(0.0,1.0,0.0,1.0);

  for (vert = 0; vert < 3; vert++) {
    v[0] = spL->cp[0].vt[vert];
    v[1] = plcl_vect_sum(spL->cp[0].vt[vert],spL->cp[0].vt2[vert]);
    plCurve_add_component(L,L->nc,2,true,1,v,&clr);
  }

  printf("(normalization World none)\n");
  printf("(geometry plCurve \n");
  plCurve_write(stdout,L);
  printf(")\n");

  L->cp[0].cc = 1;
  L->cp[0].clr[0] = plCurve_build_color(0.0,0.0,1.0,1.0);

  plCurve_free(L);
  nv[0] = 200;
  L = plCurve_convert_from_spline(spL,nv);
  
  printf("(geometry spline \n");
  plCurve_write(stdout,L);
  printf(")\n");
  printf("(look-recenter spline) (look-encompass spline)\n");

  exit(EXIT_SUCCESS);
}
