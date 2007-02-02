/*
 * $Id: run_tests.c,v 1.24 2007-02-02 17:51:12 ashted Exp $
 *
 * Test all of the library code.
 *
 */

#include <config.h>
#include <plCurve.h>

#ifdef HAVE_STDIO_H
  #include <stdio.h>
#endif
#ifdef HAVE_STDLIB_H
  #include <stdlib.h>
#endif
#ifdef HAVE_ASSERT_H
  #include <assert.h>
#endif
#ifdef HAVE_FLOAT_H
  #include <float.h>
#endif
#ifdef HAVE_MATH_H
  #include <math.h>
#endif
#ifdef HAVE_STRING_H
  #include <string.h>
#endif
/* To get the macro WEXITSTATUS */
#ifdef HAVE_WAIT_H
  #include <wait.h>
#endif
#ifdef HAVE_SYS_WAIT_H
  #include <sys/wait.h>
#endif
#ifdef HAVE_UNISTD_H
  #include <unistd.h>
#endif

#define require(C) \
  if (!(C)) { \
    fprintf(stderr,"\n--------------======================-------------\n" \
        "%s:%d: failed requirement `%s'\n\n",__FILE__,__LINE__,#C); \
    return false; \
  }

/* Compare two plCurves and make sure that they match */
static bool curves_match(const plCurve A, const plCurve B)
  /*@modifies nothing@*/ {
  int cmp, vert, clr;
  plc_constraint *cstA,*cstB;
  plc_vert_quant *qntA,*qntB;

  require(A.nc == B.nc);
  /* Check components */
  for (cmp = 0; cmp < A.nc; cmp++) {
    require(A.cp[cmp].nv == B.cp[cmp].nv);
    require(A.cp[cmp].open == B.cp[cmp].open);
    require(A.cp[cmp].cc == B.cp[cmp].cc);
    for (vert = -1; vert <= A.cp[cmp].nv; vert++) {
      if (!plc_M_vecteq(A.cp[cmp].vt[vert],B.cp[cmp].vt[vert])) {
        printf("Vertex difference %d:%d: %.20g,%.20g,%.20g\n",cmp,vert,
            plc_M_clist(plc_vect_diff(A.cp[cmp].vt[vert],
                                        B.cp[cmp].vt[vert])));
      }
      require(plc_M_vecteq(A.cp[cmp].vt[vert],B.cp[cmp].vt[vert]));
    }
    for (clr = 0; clr < A.cp[cmp].cc; clr++) {
      require(fabs(A.cp[cmp].clr[clr].r - B.cp[cmp].clr[clr].r) < DBL_EPSILON);
      require(fabs(A.cp[cmp].clr[clr].g - B.cp[cmp].clr[clr].g) < DBL_EPSILON);
      require(fabs(A.cp[cmp].clr[clr].b - B.cp[cmp].clr[clr].b) < DBL_EPSILON);
      require(fabs(A.cp[cmp].clr[clr].alpha -
                  B.cp[cmp].clr[clr].alpha) < DBL_EPSILON);
    }
  }
  /* Check constraints */
  cstA = A.cst;
  cstB = B.cst;
  while (cstA != NULL && cstB != NULL) {
    require(cstA->kind == cstB->kind);
    require(plc_M_vecteq(cstA->vect[0],cstB->vect[0]));
    require(plc_M_vecteq(cstA->vect[1],cstB->vect[1]));
    require(cstA->cmp == cstB->cmp);
    require(cstA->vert == cstB->vert);
    require(cstA->num_verts == cstB->num_verts);
    cstA = cstA->next;
    cstB = cstB->next;
  }
  require(cstA == NULL);
  require(cstB == NULL);
  /* Check quantifiers */
  qntA = A.quant;
  qntB = B.quant;
  while (qntA != NULL && qntB != NULL) {
    require(qntA->cmp == qntB->cmp);
    require(qntA->vert == qntB->vert);
    require(strcmp(qntA->tag,qntB->tag) == 0);
    require(fabs(qntA->quant - qntB->quant) < DBL_EPSILON);
    qntA = qntA->next;
    qntB = qntB->next;
  }
  require(qntA == NULL);
  require(qntB == NULL);
  return true;
}

#define list_csts(S) list_constraints(S,#S);
/*@unused@*/ static inline void list_constraints(plCurve *L, const char *str) {
  plc_constraint *cst;

  cst = L->cst;
  printf("Constraint list (%s):\n",str);
  while (cst != NULL) {
    printf("  %d %d %d %d %g %g %g %g %g %g\n",
        (int)cst->kind, cst->cmp, cst->vert, cst->num_verts,
        plc_M_clist(cst->vect[0]), plc_M_clist(cst->vect[1]));
    cst = cst->next;
  }
}

#define components 2
#define check(C) \
  assert(C); \
  printf("Passed %d\r",__LINE__); \
  (void)fflush(stdout);

int main(void) {
  plCurve S; /* The standard against which to measure */
  plCurve *L;
  int nv[components] = { 3, 4 };
  bool open[components] = { false, true };
  int cc[components] = { 1, 4 };
  char version[80];
  char revision[] = "$Revision: 1.24 $";
  plc_vert_quant *quant;
  int cmp, vert, ctr;
  double dist, temp_dbl;
  bool ok;
  plc_vector temp_vect, temp_vect2, zero_vect;
  plc_color temp_clr, temp_clr2;
  int sysret;
  FILE *filehandle;
  char *filename;
  int err_num;
  char err_str[80];
#define num_bad_vect_files 14
  int bad_read_results[num_bad_vect_files+1] = 
    { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 11, 11, 8 };
  double clens[2],diameter;
  plc_spline *spL;
  plc_constraint *temp_cst;

  zero_vect = plc_build_vect(0.0,0.0,0.0);

  /* Test plc_version */
  plc_version(NULL,0);
  version[0] = '\0';
  plc_version(version,sizeof(version));
  revision[strlen(revision)-2] = '\0';
  printf("run_tests %s (plCurve v. %s)\n", revision+11, version);
  check(strcmp(PACKAGE_VERSION,version) == 0);

  /* Check on plc_build_vect */
  temp_vect.c[0] = 10.0*rand()/RAND_MAX;
  temp_vect.c[1] = 10.0*rand()/RAND_MAX;
  temp_vect.c[2] = 10.0*rand()/RAND_MAX;
  temp_vect2 = plc_build_vect(temp_vect.c[0],temp_vect.c[1],temp_vect.c[2]);
  check(fabs(temp_vect.c[0] - temp_vect2.c[0]) < DBL_EPSILON);
  check(fabs(temp_vect.c[1] - temp_vect2.c[1]) < DBL_EPSILON);
  check(fabs(temp_vect.c[2] - temp_vect2.c[2]) < DBL_EPSILON);

  /* Test plc_build_color */
  temp_clr.r = 10.0*rand()/RAND_MAX;
  temp_clr.g = 10.0*rand()/RAND_MAX;
  temp_clr.b = 10.0*rand()/RAND_MAX;
  temp_clr.alpha = 10.0*rand()/RAND_MAX;
  temp_clr2 = plc_build_color(temp_clr.r,temp_clr.g,temp_clr.b,
                                  temp_clr.alpha);
  check(fabs(temp_clr.r - temp_clr2.r) < DBL_EPSILON);
  check(fabs(temp_clr.g - temp_clr2.g) < DBL_EPSILON);
  check(fabs(temp_clr.b - temp_clr2.b) < DBL_EPSILON);
  check(fabs(temp_clr.alpha - temp_clr2.alpha) < DBL_EPSILON);

  /* Check on vect_eq, plc_random_vect and plc_normalize */
  check(plc_vecteq(temp_vect,temp_vect2));
  temp_vect2 = plc_random_vect();
  check(!plc_vecteq(temp_vect,temp_vect2));
  ok = true;
  dist = plc_norm(plc_normalize_vect(temp_vect2,&ok)) - 1.0;
  check(ok);
  check(fabs(dist) < DBL_EPSILON);
  (void)plc_normalize_vect(zero_vect,&ok);
  check(!ok);

  /* Run the external tests */
  /* 1) Check to see if plc_normalize_vect fails on zero vectors with no bool
   *    passed in, as advertised. 
   */
  sysret = system("./exit_failure_tests 1");
  check(WEXITSTATUS(sysret) == EXIT_FAILURE);
  /* 2) Make sure the plc_component_div fails properly */
  sysret = system("./exit_failure_tests 2");
  check(WEXITSTATUS(sysret) == EXIT_FAILURE);
  /* 3) And make sure that the "failure program" isn't just failing */
  sysret = system("./exit_failure_tests 3");
  check(WEXITSTATUS(sysret) == EXIT_SUCCESS);

  /* Now to see if we can create plCurves.  Tests plc_new, plc_vmadd, and
   * plc_fix_wrap. */
  S.nc = components;
  /* Allocate two extra components which we will use when we test
   * _add_component and _drop component */
  S.cp = (plc_strand *)calloc((size_t)4,sizeof(plc_strand));
  check(S.cp != NULL);
  S.quant = NULL;
  S.cst = NULL;
  S.cp[0].nv = nv[0];
  S.cp[0].open = open[0];
  S.cp[0].cc = cc[0];
  S.cp[0].vt = (plc_vector *)calloc((size_t)5,sizeof(plc_vector));
  check(S.cp[0].vt != NULL);
  S.cp[0].vt++;
  S.cp[0].clr = (plc_color *)calloc((size_t)1,sizeof(plc_color));
  check(S.cp[0].clr != NULL);
  S.cp[0].vt[0] = zero_vect;
  S.cp[0].vt[1] = plc_build_vect(1.0,0.0,0.0);
  S.cp[0].vt[2] = plc_build_vect(0.0,1.0,0.0);
  S.cp[0].vt[-1] = S.cp[0].vt[2];
  S.cp[0].vt[3] = S.cp[0].vt[0];
  S.cp[0].clr[0] = plc_build_color(1.0,0.0,1.0,1.0);
  S.cp[1].nv = nv[1];
  S.cp[1].open = open[1];
  S.cp[1].cc = cc[1];
  S.cp[1].vt = (plc_vector *)calloc((size_t)6,sizeof(plc_vector));
  check(S.cp[1].vt != NULL);
  S.cp[1].vt++;
  S.cp[1].clr = (plc_color *)calloc((size_t)4,sizeof(plc_color));
  check(S.cp[1].clr != NULL);
  S.cp[1].vt[0] = plc_build_vect(1.0,1.0,2.0);
  S.cp[1].vt[1] = plc_build_vect(1.0,1.0,1.0);
  S.cp[1].vt[2] = plc_build_vect(1.0,1.0,0.0);
  S.cp[1].vt[3] = plc_build_vect(1.0,1.0,-1.0);
  S.cp[1].vt[4] = S.cp[1].vt[2];
  S.cp[1].vt[-1] = S.cp[1].vt[1];
  S.cp[1].clr[0] = plc_build_color(0.00,1.0,0.0,1.0);
  S.cp[1].clr[1] = plc_build_color(0.33,1.0,0.0,1.0);
  S.cp[1].clr[2] = plc_build_color(0.67,1.0,0.0,1.0);
  S.cp[1].clr[3] = plc_build_color(1.00,1.0,0.0,1.0);

  L = plc_new(components,nv,open,cc);
  L->cp[0].vt[0] = zero_vect;
  L->cp[0].vt[1] = plc_build_vect(1.0,0.0,0.0);
  L->cp[0].vt[2] = plc_build_vect(0.0,1.0,0.0);
  L->cp[1].vt[0] = plc_build_vect(1.0,1.0,2.0);
  L->cp[1].vt[1] = plc_build_vect(1.0,1.0,1.0);
  L->cp[1].vt[2] = plc_vmadd(L->cp[0].vt[1],1.0,L->cp[0].vt[2]);
  check(plc_M_vecteq(L->cp[1].vt[2],plc_build_vect(1.0,1.0,0.0)));
  L->cp[1].vt[3] = plc_build_vect(1.0,1.0,-1.0);
  plc_fix_wrap(L);
  L->cp[0].clr[0] = plc_build_color(1.0,0.0,1.0,1.0);
  L->cp[1].clr[0] = plc_build_color(0.00,1.0,0.0,1.0);
  L->cp[1].clr[1] = plc_build_color(0.33,1.0,0.0,1.0);
  L->cp[1].clr[2] = plc_build_color(0.67,1.0,0.0,1.0);
  L->cp[1].clr[3] = plc_build_color(1.00,1.0,0.0,1.0);
  check(curves_match(S,*L));

  /* Calculate the diameter of the vertex set */
  diameter = 0.0;
  diameter = plc_pointset_diameter(L);
  check(fabs(diameter - 3.0) < DBL_EPSILON);

  /* Check plc_cross_prod */
  temp_vect = plc_cross_prod(L->cp[0].vt[1],L->cp[0].vt[2]);
  check(fabs(temp_vect.c[0]) < DBL_EPSILON);
  check(fabs(temp_vect.c[1]) < DBL_EPSILON);
  check(fabs(temp_vect.c[2] - 1.0) < DBL_EPSILON);

  /* Check plc_component_div (3 ways, the 4th is checked externally above). */
  ok = true;
  temp_vect = plc_component_div(L->cp[1].vt[1],L->cp[1].vt[0],&ok);
  check(ok);
  check(plc_vecteq(temp_vect,plc_build_vect(1.0,1.0,0.5)));
  temp_vect = plc_component_div(L->cp[1].vt[1],L->cp[0].vt[0],&ok);
  check(!ok);
  temp_vect = plc_component_div(L->cp[1].vt[1],L->cp[1].vt[0],NULL);
  check(plc_vecteq(temp_vect,plc_build_vect(1.0,1.0,0.5)));

  /* Check plc_vweighted */
  temp_vect = plc_vweighted(0.7,L->cp[0].vt[1],L->cp[0].vt[2]);
  check(plc_M_vecteq(temp_vect,plc_build_vect(0.3,0.7,0.0)));

  /* Check plc_dot_prod */
  check(fabs(plc_dot_prod(L->cp[1].vt[0],L->cp[1].vt[3])) < DBL_EPSILON);

  /* Check plc_distance and plc_sq_dist */
  dist = plc_distance(L->cp[1].vt[3],L->cp[0].vt[2]) - sqrt(2.0);
  check(fabs(dist) < DBL_EPSILON);
  dist = plc_sq_dist(L->cp[0].vt[2],L->cp[1].vt[3]);
  check(fabs(dist - 2.0) < DBL_EPSILON);

  /* Check plc_num_edges */
  check(plc_num_edges(L) == 6);
  check(plc_num_verts(L) == 7);

  /* Check all the constraint-setting code!  Be very careful when editing 
   * the code between here and the "End of constraint-setting code checks"
   * as this code goes to great trouble to check all of the branches in that
   * code and that's a fairly fragile thing to do.
   */

  /* Start by doing nothing (it is one of the code paths :-) */
  plc_unconstrain(L,0,0,3);
  check(curves_match(S,*L));

  /* Allocate 10 constraint slots, we'll attach and detach them as needed */
  S.cst = (plc_constraint *)calloc((size_t)10,sizeof(plc_constraint));
  check(S.cst != NULL);
  check(S.cst->next == NULL);
  S.cst[0].kind = fixed;
  S.cst[0].vect[0] = L->cp[0].vt[0];
  S.cst[0].vect[1] = zero_vect;
  S.cst[0].cmp = 0;
  S.cst[0].vert = 0;
  S.cst[0].num_verts = 1;
  plc_set_fixed(L,0,0,L->cp[0].vt[0]);
  check(curves_match(S,*L));

  /*@-immediatetrans@*/
  S.cst[0].next = &(S.cst[1]);
  S.cst[1].kind = fixed;
  S.cst[1].vect[0] = L->cp[0].vt[0];
  S.cst[1].vect[1] = zero_vect;
  S.cst[1].cmp = 0;
  S.cst[1].vert = 1;
  S.cst[1].num_verts = 1;
  plc_set_fixed(L,0,1,L->cp[0].vt[0]);
  check(curves_match(S,*L));

  plc_set_fixed(L,0,1,L->cp[0].vt[0]);
  check(curves_match(S,*L));

  S.cst[0].next = NULL;
  S.cst[0].kind = line;
  S.cst[0].num_verts = 3;
  plc_constrain_to_line(L,0,0,3,L->cp[0].vt[0],zero_vect);
  check(curves_match(S,*L));

  S.cst[0].next = &(S.cst[1]);
  S.cst[0].num_verts = 2;
  S.cst[1].kind = line;
  S.cst[1].vect[0] = L->cp[0].vt[1];
  S.cst[1].vect[1] = zero_vect;
  S.cst[1].cmp = 0;
  S.cst[1].vert = 2;
  S.cst[1].num_verts = 1;
  S.cst[1].next = NULL;
  plc_constrain_to_line(L,0,2,1,L->cp[0].vt[1],zero_vect);
  check(curves_match(S,*L));

  S.cst[0].next = NULL;
  plc_unconstrain(L,0,2,1);
  check(curves_match(S,*L));

  S.cst[0].num_verts = 3;
  plc_constrain_to_line(L,0,2,1,L->cp[0].vt[0],zero_vect);
  check(curves_match(S,*L));

  S.cst[0].num_verts = 1;
  plc_unconstrain(L,0,1,2);
  check(curves_match(S,*L));

  S.cst[0].next = &(S.cst[1]);
  S.cst[1] = S.cst[0];
  S.cst[1].vert = 2;
  S.cst[1].next = NULL;
  plc_constrain_to_line(L,0,2,1,L->cp[0].vt[0],zero_vect);
  check(curves_match(S,*L));

  S.cst[0].next = NULL;
  S.cst[0].num_verts = 3;
  plc_constrain_to_line(L,0,1,1,L->cp[0].vt[0],zero_vect);
  check(curves_match(S,*L));

  S.cst[0].num_verts = 1;
  S.cst[2] = S.cst[0];
  S.cst[2].vert = 2;
  S.cst[0].next = &(S.cst[1]);
  S.cst[1].next = &(S.cst[2]);
  S.cst[1].kind = plane;
  S.cst[1].vect[0] = L->cp[0].vt[0];
  S.cst[1].vert = 1;
  plc_constrain_to_plane(L,0,1,1,L->cp[0].vt[0],0.0);
  /*@-compmempass@*/
  check(curves_match(S,*L));

  S.cst[0].next = NULL;
  S.cst[0].num_verts = 3;
  plc_constrain_to_line(L,0,1,1,L->cp[0].vt[0],zero_vect);
  check(curves_match(S,*L));

  S.cst[0].next = &(S.cst[2]);
  S.cst[2].kind = line;
  S.cst[2].vect[0] = plc_cross_prod(L->cp[0].vt[1],L->cp[0].vt[2]);
  check(plc_M_vecteq(S.cst[2].vect[0],plc_build_vect(0.0,0.0,1.0)));
  S.cst[2].vect[0] = L->cp[1].vt[2];
  S.cst[2].cmp = 1;
  S.cst[2].vert = 2;
  S.cst[2].num_verts = 2;
  plc_constrain_to_line(L,1,3,1,S.cst[2].vect[0],S.cst[2].vect[1]);
  plc_constrain_to_line(L,1,2,1,S.cst[2].vect[0],S.cst[2].vect[1]);
  check(curves_match(S,*L));

  /*@-nullret@*/
  S.cst[0] = S.cst[2];
  /*@=nullret@*/
  plc_unconstrain(L,0,0,3);
  check(curves_match(S,*L));

  S.cst[0].next = &(S.cst[2]);
  S.cst[0].kind = line;
  S.cst[0].vect[0] = S.cst[2].vect[0];
  S.cst[0].vect[1] = S.cst[0].vect[0];
  S.cst[0].cmp = 0;
  S.cst[0].vert = 2;
  S.cst[0].num_verts = 1;
  plc_constrain_to_line(L,0,2,1,S.cst[2].vect[0],S.cst[2].vect[0]);
  /*@-onlytrans@*/
  check(curves_match(S,*L));

  S.cst[1] = S.cst[0];
  /*@=onlytrans@*/
  S.cst[0].next = &(S.cst[1]);
  S.cst[0].vert = 0;
  plc_constrain_to_line(L,0,0,1,S.cst[2].vect[0],S.cst[2].vect[0]);
  check(curves_match(S,*L));

  S.cst[0].next = &(S.cst[2]);
  S.cst[0].num_verts = 3;
  S.cst[0].vect[1] = zero_vect;
  plc_constrain_to_line(L,0,0,3,S.cst[2].vect[0],zero_vect);
  check(curves_match(S,*L));

  S.cst[1] = S.cst[0];
  S.cst[0].next = &(S.cst[1]);
  S.cst[0].num_verts = 1;
  S.cst[1].vect[1] = S.cst[1].vect[0];
  S.cst[1].vert = 1;
  S.cst[1].num_verts = 2;
  plc_unconstrain(L,0,0,3);
  plc_constrain_to_line(L,0,0,1,S.cst[2].vect[0],zero_vect);
  plc_constrain_to_line(L,0,1,2,S.cst[0].vect[0],S.cst[0].vect[0]);
  check(curves_match(S,*L));

  S.cst[2].vert = 0;
  S.cst[2].num_verts = 4;
  plc_constrain_to_line(L,1,0,3,S.cst[2].vect[0],S.cst[2].vect[1]);
  check(curves_match(S,*L));

  S.cst[3] = S.cst[2];
  S.cst[2].next = &(S.cst[3]);
  S.cst[2].num_verts = 2;
  S.cst[3].vert = 3;
  S.cst[3].num_verts = 1;
  plc_unconstrain(L,1,2,1);
  check(curves_match(S,*L));

  /* Remove a specific type of constraint from all components (we check 
   * our code to remove all constraints just before cleanup below). */
  S.cst[0] = S.cst[1];
  S.cst[0].next = NULL;
  vert = plc_remove_constraint(L,line,S.cst[2].vect);
  check(vert == 4);
  check(curves_match(S,*L));

  /* A little more code needs covering in overrun_check */
  S.cst[1].kind = plane;
  S.cst[1].vect[0] = S.cp[0].vt[1];
  S.cst[1].vect[1] = zero_vect;
  S.cst[1].cmp = 1;
  S.cst[1].vert = 0;
  S.cst[1].num_verts = 3;
  S.cst[2].kind = plane;
  S.cst[2].vect[0] = zero_vect;
  S.cst[2].vect[1] = S.cp[0].vt[1];
  S.cst[2].cmp = 1;
  S.cst[2].vert = 3;
  S.cst[2].num_verts = 1;
  S.cst[0].next = &(S.cst[1]);
  S.cst[1].next = &(S.cst[2]);
  S.cst[2].next = NULL;
  plc_constrain_to_plane(L,1,0,1,S.cp[0].vt[1],0.0);
  plc_constrain_to_plane(L,1,2,2,zero_vect,1.0);
  plc_constrain_to_plane(L,1,1,2,S.cp[0].vt[1],0.0);
  check(curves_match(S,*L));
  S.cst[1].next = NULL;
  S.cst[1].num_verts = 4;
  plc_unconstrain(L,1,0,4);
  plc_constrain_to_plane(L,1,0,1,S.cp[0].vt[1],0.0);
  plc_constrain_to_plane(L,1,2,2,S.cp[0].vt[1],0.0);
  plc_constrain_to_plane(L,1,1,2,S.cp[0].vt[1],0.0);
  check(curves_match(S,*L));
  /* End of constraint-setting code checks */

  /* Now set somewhat reasonable constraints and get ready to test reading and
   * writing. */
  S.cst[0].kind = plane;
  plc_M_cross(S.cst[0].vect[0],S.cp[0].vt[1],S.cp[0].vt[2]);
  S.cst[0].vect[1] = S.cp[0].vt[0];
  S.cst[0].cmp = 0;
  S.cst[0].vert = 0;
  S.cst[0].num_verts = 1;
  S.cst[1] = S.cst[0];
  S.cst[1].vert = 2;
  S.cst[2].kind = fixed;
  S.cst[2].vect[0] = plc_build_vect(2.0,2.0,1.0);
  S.cst[2].vect[1] = zero_vect;
  S.cst[2].cmp = 1;
  S.cst[2].vert = 0;
  S.cst[2].num_verts = 1;
  S.cst[3].kind = line;
  S.cst[3].vect[0] = S.cst[0].vect[0];
  S.cst[3].vect[1] = plc_build_vect(2.0,1.0,3.7);
  S.cst[3].cmp = 1;
  S.cst[3].vert = 1;
  S.cst[3].num_verts = 2;
  S.cst[4] = S.cst[2];
  S.cst[4].vert = 3;
  S.cst[4].vect[0].c[2] = 0.0;
  /*@-usereleased*/
  S.cst[0].next = &(S.cst[1]);
  S.cst[1].next = &(S.cst[2]);
  S.cst[2].next = &(S.cst[3]);
  S.cst[3].next = &(S.cst[4]);
  S.cst[4].next = NULL;
  plc_constrain_to_plane(L,0,0,3,S.cst[0].vect[0],0.0);
  plc_unconstrain(L,0,1,1);
  plc_set_fixed(L,1,0,S.cst[2].vect[0]);
  plc_constrain_to_line(L,1,1,2,S.cst[3].vect[0],S.cst[3].vect[1]);
  plc_set_fixed(L,1,3,S.cst[4].vect[0]);
  check(curves_match(S,*L));
  /*@=usereleased*/

  /*@=immediatetrans@*/

  /* Now plc_check_cst, plc_component_mult and plc_vlincomb */
  plc_constrain_to_plane(L,0,0,1,S.cp[0].vt[1],0.0);
  plc_constrain_to_plane(L,0,1,1,S.cp[0].vt[2],0.0);
  dist = plc_check_cst(L) - sqrt(3.0);
  check(fabs(dist) < DBL_EPSILON);
  for (vert = 2; vert >= 0; vert--) {
    plc_constrain_to_plane(L,0,vert,3-vert,
        (vert == 2) ?
      plc_component_mult(S.cst[0].vect[0],plc_build_vect(-1.0,-1.0,-1.0)) :
      plc_vlincomb(-1.0,S.cp[0].vt[vert+1],1.0*vert,S.cp[0].vt[0]),0.0);
  }
  dist = plc_check_cst(L) - sqrt(3.0);
  check(fabs(dist) < DBL_EPSILON);

  /* Now put the constraint back the way it was */
  plc_constrain_to_plane(L,0,0,3,S.cst[0].vect[0],0.0);
  plc_unconstrain(L,0,1,1);
  check(curves_match(S,*L));

  /* And plc_fix_cst */
  plc_fix_cst(L);
  dist = plc_check_cst(L);
  check(fabs(dist) < DBL_EPSILON);
  S.cp[1].vt[0] = S.cst[2].vect[0];
  S.cp[1].vt[3] = S.cst[4].vect[0];
  S.cp[1].vt[1] = plc_vect_sum(S.cp[1].vt[1],S.cp[0].vt[1]);
  S.cp[1].vt[2] = plc_vect_sum(S.cp[1].vt[2],S.cp[0].vt[1]);
  S.cp[1].vt[-1] = S.cp[1].vt[1];
  S.cp[1].vt[4] = S.cp[1].vt[2];
  check(curves_match(S,*L));

  /* Now check _write and _read */ 
  filename = tmpnam(NULL);
  check(filename != NULL);
  filehandle = fopen(filename,"w");
  check(filehandle != NULL);
  plc_write(filehandle,L);
  sysret = fclose(filehandle);
  check(sysret == 0);
  plc_free(L);
  L = NULL;
  filehandle = fopen(filename,"r");
  check(filehandle != NULL);
  L = plc_read(filehandle,&err_num,err_str,sizeof(err_str));
  if (err_num != 0) {
    printf("Trouble reading file: %s\n",err_str);
  }
  check(err_num == 0);
  check(L != NULL);
  sysret = fclose(filehandle);
  filehandle = NULL;
  check(sysret == 0);
  sysret = unlink(filename);
  check(sysret == 0);
  check(curves_match(S,*L));

  /* Read invalid testfiles to exercise get_comment and PlCurve_read */
  filename = (char *)malloc((size_t)40);
  check(filename != NULL);
  for (ctr = 1; ctr <= num_bad_vect_files; ctr++) {
    (void)snprintf(filename,(size_t)40,"bad_vects/bad_%d.vect",ctr);
    filehandle = fopen(filename,"r");
    if (filehandle == NULL) {
      fprintf(stderr,"Unable to open bad_vects/bad_%d.vect\n",ctr);
    }
    check(filehandle != NULL);
    (void)plc_read(filehandle,&err_num,err_str,sizeof(err_str));
    check(err_num == bad_read_results[ctr]);
  }
  free(filename);
  filename = NULL;

  /* See if we can compute local minrad-curvature */
  dist = plc_MR_curvature(L,1,1) - 2.0;
  check(fabs(dist) < DBL_EPSILON);
  dist = plc_MR_curvature(L,1,0);
  check(fabs(dist) < DBL_EPSILON);

  /* Find a unit tangent vector */
  ok = true;
  temp_vect = plc_mean_tangent(L,0,0,&ok);
  check(ok);
  temp_dbl = sqrt(2.0)/2.0;
  temp_vect = plc_vect_diff(temp_vect,plc_build_vect(temp_dbl,-temp_dbl,0.0));
  check(plc_vecteq(temp_vect,zero_vect));
  temp_vect = plc_mean_tangent(L,1,0,&ok);
  check(ok);
  check(plc_vecteq(temp_vect,plc_build_vect(0.0,-1.0,0.0)));
  temp_vect = plc_mean_tangent(L,1,1,&ok);
  check(ok);
  check(plc_vecteq(temp_vect,plc_build_vect(0.0,-temp_dbl,-temp_dbl)));
  temp_vect = plc_mean_tangent(L,1,3,&ok);
  check(ok);
  check(plc_vecteq(temp_vect,plc_build_vect(0.0,1.0,0.0)));

  /* And the lengths of the components */
  dist = plc_arclength(L,clens);
  dist -= 5.0 + sqrt(2.0);
  clens[0] -= 2.0 + sqrt(2.0);
  clens[1] -= 3.0;
  check(fabs(dist) < DBL_EPSILON);
  check(fabs(clens[0]) < DBL_EPSILON);
  check(fabs(clens[1]) < DBL_EPSILON);

  /* Now subarc lengths */
  dist = plc_subarc_length(L,0,2,1);
  dist -= sqrt(2.0);
  check(fabs(dist) < DBL_EPSILON);
  dist = plc_subarc_length(L,0,0,2) - 1.0;
  check(fabs(dist) < DBL_EPSILON);
  dist = plc_subarc_length(L,1,0,3) - 3.0;
  check(fabs(dist) < DBL_EPSILON);

  /* Check out copying of plCurves */
  plc_free(L);
  L = plc_copy(&S);
  check(curves_match(S,*L));

  /* Now can we convert to a spline and back?  We remember that conversion
   * to splines removes constraints (and probably quantifiers) anyway. */
  ok = true;
  spL = plc_convert_to_spline(L,&ok);
  /* Overall data */
  check(ok);
  check(spL->nc == 2);
  check(spL->cp[0].open == false);
  check(spL->cp[1].open == true);
  check(spL->cp[0].ns == 3);
  check(spL->cp[1].ns == 4);
  check(spL->cp[0].cc == 1);
  check(spL->cp[1].cc == 4);
  /* Vertices */
  for (cmp = 0; cmp < L->nc; cmp++) {
    for (vert = -1; vert <= L->cp[cmp].nv; vert++) {
      check(plc_vecteq(spL->cp[cmp].vt[vert],L->cp[cmp].vt[vert]));
    }
  }
  /* Second derivatives */
  check(plc_vecteq(spL->cp[0].vt2[-1], plc_build_vect(0.0, 0.0, 0.0)));
  /*
  printf("%.16g %.16g %.16g\n",plc_M_clist(
        plc_vect_diff(spL->cp[0].vt2[0],
        plc_build_vect(2.4065244382070420, 1.3770871866841826, 0.0))));
  printf("%.16g %.16g %.16g\n",plc_M_clist(
        plc_vect_diff(spL->cp[0].vt2[1],
        plc_build_vect(-3.0556895635333689, 1.4884663137509202, 0.0))));
  printf("%.16g %.16g %.16g\n",plc_M_clist(
        plc_vect_diff(spL->cp[0].vt2[2],
        plc_build_vect(1.4884663137509202, -3.0556895635333694, 0.0))));
  printf("%.16g %.16g %.16g\n",plc_M_clist(
        plc_vect_diff(spL->cp[0].vt2[3],
        plc_build_vect(1.3770871866841824, 2.4065244382070418, 0.0))));
  */
  check(plc_vecteq(spL->cp[0].vt2[0],
        plc_build_vect(2.4065244382070420, 1.3770871866841826, 0.0)));
  check(plc_vecteq(spL->cp[0].vt2[1],
        plc_build_vect(-3.0556895635333689, 1.4884663137509202, 0.0)));
  check(plc_vecteq(spL->cp[0].vt2[2],
        plc_build_vect(1.4884663137509202, -3.0556895635333694, 0.0)));
  check(plc_vecteq(spL->cp[0].vt2[3],
        plc_build_vect(1.3770871866841824, 2.4065244382070418, 0.0)));
  check(plc_vecteq(spL->cp[1].vt2[-1], plc_build_vect(0.0, 0.0, 0.0)));
  check(plc_vecteq(spL->cp[1].vt2[0],
        plc_build_vect(0.0, -0.6666666666666666, 1.2)));
  check(plc_vecteq(spL->cp[1].vt2[1],
        plc_build_vect(0.0, 1.3333333333333333, -2.4)));
  check(plc_vecteq(spL->cp[1].vt2[2],
        plc_build_vect(0.0, 1.3333333333333332, 2.4000000000000004)));
  check(plc_vecteq(spL->cp[1].vt2[3],
        plc_build_vect(0.0, -0.6666666666666666, -1.2)));
  check(plc_vecteq(spL->cp[1].vt2[4], plc_build_vect(0.0, 0.0, 0.0)));
  check(fabs(spL->cp[0].svals[0] - 0.0) < DBL_EPSILON);
  dist = plc_distance(L->cp[0].vt[0],L->cp[0].vt[1]);
  check(fabs(spL->cp[0].svals[1] - dist) < DBL_EPSILON);
  check(fabs(spL->cp[1].svals[0] - 0.0) < DBL_EPSILON);
  dist = plc_distance(L->cp[1].vt[0],L->cp[1].vt[1]);
  check(fabs(spL->cp[1].svals[1] - dist) < DBL_EPSILON);
  dist = plc_distance(L->cp[1].vt[1],L->cp[1].vt[2]);
  check(fabs(spL->cp[1].svals[2] - spL->cp[1].svals[1] - dist) < DBL_EPSILON);
  dist = plc_distance(L->cp[1].vt[2],L->cp[1].vt[3]);
  check(fabs(spL->cp[1].svals[3] - spL->cp[1].svals[2] - dist) < DBL_EPSILON);
  check(spL->cp[0].cc == 1);
  check(fabs(spL->cp[0].clr[0].r - L->cp[0].clr[0].r) < DBL_EPSILON);
  check(fabs(spL->cp[0].clr[0].g - L->cp[0].clr[0].g) < DBL_EPSILON);
  check(fabs(spL->cp[0].clr[0].b - L->cp[0].clr[0].b) < DBL_EPSILON);
  check(fabs(spL->cp[0].clr[0].alpha-L->cp[0].clr[0].alpha) < DBL_EPSILON);
  check(plc_vecteq(spL->cp[1].vt[-1],L->cp[1].vt[-1]));
  check(plc_vecteq(spL->cp[1].vt[0],L->cp[1].vt[0]));
  check(plc_vecteq(spL->cp[1].vt[1],L->cp[1].vt[1]));
  check(plc_vecteq(spL->cp[1].vt[2],L->cp[1].vt[2]));
  check(plc_vecteq(spL->cp[1].vt[3],L->cp[1].vt[3]));
  check(plc_vecteq(spL->cp[1].vt[4],L->cp[1].vt[4]));
  check(spL->cp[1].cc == 4);
  check(fabs(spL->cp[1].clr[0].r - L->cp[1].clr[0].r) < DBL_EPSILON);
  check(fabs(spL->cp[1].clr[0].g - L->cp[1].clr[0].g) < DBL_EPSILON);
  check(fabs(spL->cp[1].clr[0].b - L->cp[1].clr[0].b) < DBL_EPSILON);
  check(fabs(spL->cp[1].clr[0].alpha - L->cp[1].clr[0].alpha) < DBL_EPSILON);
  check(fabs(spL->cp[1].clr[1].r - L->cp[1].clr[1].r) < DBL_EPSILON);
  check(fabs(spL->cp[1].clr[1].g - L->cp[1].clr[1].g) < DBL_EPSILON);
  check(fabs(spL->cp[1].clr[1].b - L->cp[1].clr[1].b) < DBL_EPSILON);
  check(fabs(spL->cp[1].clr[1].alpha - L->cp[1].clr[1].alpha) < DBL_EPSILON);
  check(fabs(spL->cp[1].clr[2].r - L->cp[1].clr[2].r) < DBL_EPSILON);
  check(fabs(spL->cp[1].clr[2].g - L->cp[1].clr[2].g) < DBL_EPSILON);
  check(fabs(spL->cp[1].clr[2].b - L->cp[1].clr[2].b) < DBL_EPSILON);
  check(fabs(spL->cp[1].clr[2].alpha - L->cp[1].clr[2].alpha) < DBL_EPSILON);
  check(fabs(spL->cp[1].clr[3].r - L->cp[1].clr[3].r) < DBL_EPSILON);
  check(fabs(spL->cp[1].clr[3].g - L->cp[1].clr[3].g) < DBL_EPSILON);
  check(fabs(spL->cp[1].clr[3].b - L->cp[1].clr[3].b) < DBL_EPSILON);
  check(fabs(spL->cp[1].clr[3].alpha - L->cp[1].clr[3].alpha) < DBL_EPSILON);

  /* Get ready for conversion back from spline */
  S.cp[0].vt[1] = plc_build_vect(1.0250988714690254, 0.11297231598867286, 0.0);
  S.cp[0].vt[2] = plc_build_vect(0.11297231598867283, 1.0250988714690255, 0.0);
  S.cp[0].vt[-1] = S.cp[0].vt[2];
  temp_cst = S.cst;
  S.cst = NULL;

  plc_free(L);
  L = NULL;
  L = plc_convert_from_spline(spL,nv);
  check(L != NULL);
  for (cmp = 0; cmp < L->nc; cmp++) {
    for (vert = 0; vert < L->cp[cmp].nv; vert++) {
    }
  }
  check(curves_match(S,*L));

  temp_vect = plc_sample_spline(spL,1,1.5);
  check(plc_vecteq(temp_vect,plc_build_vect(2.0,2.5/3.0,0.5)));

  plc_spline_free(spL);

  /* A little exercise in the spline code */
  spL = (plc_spline *)malloc(sizeof(plc_spline));
  check(spL != NULL);
  spL->nc = 1;
  spL->cp = NULL;
  /*@-nullstate@*/
  plc_spline_free(spL);
  /*@=nullstate@*/
  spL = NULL;
  plc_spline_free(spL);

  /* Now put things back where we expect them */
  S.cst = temp_cst;
  temp_cst = NULL;
  S.cp[0].vt[-1] = S.cp[0].vt[2] = plc_build_vect(1.0,0.0,0.0);
  S.cp[0].vt[1] = plc_build_vect(0.0,1.0,0.0);
  plc_free(L);
  L = plc_copy(&S);
  check(L != NULL);
  check(curves_match(S,*L));

  /* Now check to see if we can remove all constraints */
  /*@-kepttrans@*/
  free(S.cst);
  /*@=kepttrans@*/
  S.cst = NULL;
  plc_remove_all_constraints(L);
  check(curves_match(S,*L));

  /* And check on the copying of curves without constraints. */
  plc_free(L);
  L = plc_copy(&S);
  check(curves_match(S,*L));

  /* Can we add and remove components? */
  S.nc += 2;
  S.cp[2] = S.cp[1]; 
  S.cp[1] = S.cp[3] = S.cp[0]; 
  S.cp[0].cc = 0;
  plc_add_component(L, L->nc, L->cp[0].nv, L->cp[0].open, L->cp[0].cc,
      L->cp[0].vt, L->cp[0].clr);
  plc_add_component(L,0,L->cp[0].nv,L->cp[0].open,0,L->cp[0].vt,NULL);
  check(curves_match(S,*L));

  S.nc -= 2;
  S.cp[0].cc = S.cp[1].cc;
  S.cp[1] = S.cp[2];
  plc_drop_component(L,0);
  plc_drop_component(L,2);
  check(curves_match(S,*L));

  /* Eventually quantifier testing code goes here.  For now, put a blank 
   * quantifier in to make sure that plc_free cleans it up. */
  check(L->quant == NULL);
  L->quant = (plc_vert_quant *)calloc((size_t)1,sizeof(plc_vert_quant));
  check(L->quant != NULL);
  check(L->quant->next == NULL);

  /* Cleanup phase */
  plc_free(L);
  L = NULL;
  /* And just to make sure that it works */
  plc_free(L);

  while(S.quant != NULL) {
    quant = S.quant;
    S.quant = quant->next;
    free(quant);
    quant = NULL;
  }
  check(S.quant == NULL);
  check(S.cp[0].vt != NULL);
  S.cp[0].vt--;
  free(S.cp[0].vt);
  S.cp[0].vt = NULL;
  check(S.cp[0].vt == NULL); /* For the sake of splint */
  check(S.cp[0].clr != NULL);
  free(S.cp[0].clr);
  S.cp[0].clr = NULL;
  check(S.cp[0].clr == NULL); /* For the sake of splint */
  check(S.cp[1].vt != NULL);
  S.cp[1].vt--;
  free(S.cp[1].vt);
  S.cp[1].vt = NULL;
  check(S.cp[1].vt == NULL); /* For the sake of splint */
  check(S.cp[1].clr != NULL);
  free(S.cp[1].clr);
  S.cp[1].clr = NULL;
  check(S.cp[1].clr == NULL); /* For the sake of splint */
  free(S.cp);
  printf("Passed all tests.\n");
  return(EXIT_SUCCESS);

  /*@=compmempass@*/
}
