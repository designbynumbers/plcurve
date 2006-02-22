/*
 * $Id: test_cst.c,v 1.12 2006-02-22 22:54:11 ashted Exp $
 *
 * Test the constraint-handling code
 *
 */

#include <config.h>
#include <plCurve.h>

#ifdef HAVE_STDIO_H
  #include <stdio.h>
#endif
#ifdef HAVE_MATH_H
  #include <math.h>
#endif
#ifdef HAVE_STDBOOL_H
  #include <stdbool.h>
#endif

int main() {
  plCurve *L;
  int nv[2] = { 100, 2 };
  bool open[2] = { false, true };
  int cc[2] = { 100, 1 };
  FILE *outfile;
  int i;
 
  L = plCurve_new(2,nv,open,cc);
  plcl_status_check();
  for (i = 0; i < 100; i++) {
    L->cp[0].vt[i] = plcl_build_vect(sin(0.03*i),cos(0.03*i),0);
    L->cp[0].clr[i].r = sin(0.03*i);
    L->cp[0].clr[i].g = 0.0;
    L->cp[0].clr[i].b = cos(0.03*i);
    L->cp[0].clr[i].alpha = 1.0;
  }
  L->cp[1].clr[0].r = 0.0;
  L->cp[1].clr[0].g = 1.0;
  L->cp[1].clr[0].b = 0.0;
  L->cp[1].clr[0].alpha = 1;
  L->cp[1].vt[0] = plcl_build_vect(0,0,0);
  L->cp[1].vt[1] = plcl_build_vect(0,0,1);
  plCurve_constrain_to_plane(L,0,0,100,plcl_build_vect(0,0,1),0);
  plcl_status_check();
    if ((outfile = fopen("test_cst.vect","w")) != NULL) {
      plCurve_write(outfile,L);
      (void)fclose(outfile);
      plcl_status_check();
      printf("->"); (void)fgetc(stdin);
    }
  plCurve_constrain_to_line(L,0,50,1,plcl_build_vect(1,2,3),
                                     plcl_build_vect(4,5,6));
  plcl_status_check();
    if ((outfile = fopen("test_cst.vect","w")) != NULL) {
      plCurve_write(outfile,L);
      (void)fclose(outfile);
      plcl_status_check();
      printf("->"); (void)fgetc(stdin);
    }
  plCurve_constrain_to_line(L,1,0,1,plcl_build_vect(3,2,3),
                                    plcl_build_vect(4,5,6));
  plcl_status_check();
    if ((outfile = fopen("test_cst.vect","w")) != NULL) {
      plCurve_write(outfile,L);
      (void)fclose(outfile);
      plcl_status_check();
      printf("->"); (void)fgetc(stdin);
    }
  plCurve_constrain_to_line(L,0,50,1,plcl_build_vect(1,2,3),
                                     plcl_build_vect(4,5,6));
  plcl_status_check();
    if ((outfile = fopen("test_cst.vect","w")) != NULL) {
      plCurve_write(outfile,L);
      (void)fclose(outfile);
      plcl_status_check();
      printf("->"); (void)fgetc(stdin);
    }
  plCurve_unconstrain(L,0,75,5);
  plcl_status_check();
    if ((outfile = fopen("test_cst.vect","w")) != NULL) {
      plCurve_write(outfile,L);
      (void)fclose(outfile);
      plcl_status_check();
      printf("->"); (void)fgetc(stdin);
    }
  plCurve_unconstrain(L,0,25,5);
  plcl_status_check();
    if ((outfile = fopen("test_cst.vect","w")) != NULL) {
      plCurve_write(outfile,L);
      (void)fclose(outfile);
      plcl_status_check();
      printf("->"); (void)fgetc(stdin);
    }
  plCurve_set_fixed(L,0,32,1,plcl_build_vect(3.2,0,3.2));
  plcl_status_check();
    if ((outfile = fopen("test_cst.vect","w")) != NULL) {
      plCurve_write(outfile,L);
      (void)fclose(outfile);
      plcl_status_check();
      printf("->"); (void)fgetc(stdin);
    }
  plCurve_unconstrain(L,0,2,1);
  plcl_status_check();
    if ((outfile = fopen("test_cst.vect","w")) != NULL) {
      plCurve_write(outfile,L);
      (void)fclose(outfile);
      plcl_status_check();
      printf("->"); (void)fgetc(stdin);
    }
  plCurve_fix_cst(L);
  plcl_status_check();
    if ((outfile = fopen("test_cst.vect","w")) != NULL) {
      plCurve_write(outfile,L);
      (void)fclose(outfile);
      plcl_status_check();
      printf("->"); (void)fgetc(stdin);
    }
  plCurve_free(L);

  return 0;
}
