/*
 * $Id: run_tests.c,v 1.10 2006-03-02 05:32:57 ashted Exp $
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
    fprintf(stderr,"--------------======================-------------\n" \
        "%s:%d: failed requirement `%s'\n",__FILE__,__LINE__,#C); \
    return false; \
  }

/* Compare two plCurves and make sure that they match */
static bool curves_match(const plCurve A, const plCurve B)
  /*@modifies nothing@*/ {
  int cmp, vert, clr;
  plCurve_constraint *cstA,*cstB;
  plCurve_vert_quant *qntA,*qntB;

  require(A.nc == B.nc);
  /* Check components */
  for (cmp = 0; cmp < A.nc; cmp++) {
    require(A.cp[cmp].nv == B.cp[cmp].nv);
    require(A.cp[cmp].open == B.cp[cmp].open);
    require(A.cp[cmp].cc == B.cp[cmp].cc);
    for (vert = -1; vert <= A.cp[cmp].nv; vert++) {
      require(plcl_M_vecteq(A.cp[cmp].vt[vert],B.cp[cmp].vt[vert]));
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
    require(plcl_M_vecteq(cstA->vect[0],cstB->vect[0]));
    require(plcl_M_vecteq(cstA->vect[1],cstB->vect[1]));
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
  plCurve_constraint *cst;

  cst = L->cst;
  printf("Constraint list (%s):\n",str);
  while (cst != NULL) {
    printf("  %d %d %d %d %g %g %g %g %g %g\n",
        (int)cst->kind, cst->cmp, cst->vert, cst->num_verts,
        plcl_M_clist(cst->vect[0]), plcl_M_clist(cst->vect[1]));
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
  char revision[] = "$Revision: 1.10 $";
  plCurve_constraint *cst;
  plCurve_vert_quant *quant;
  int vert;
  double dist;
  bool ok;
  plcl_vector temp_vect, temp_vect2, zero_vect;
  plCurve_color temp_clr, temp_clr2;
  int sysret;
  FILE *outfile;
  char *outname;
  int err_num;
  char err_str[80];

  zero_vect = plcl_build_vect(0.0,0.0,0.0);

  /* Test plCurve_version */
  plCurve_version(NULL,0);
  version[0] = '\0';
  plCurve_version(version,sizeof(version));
  revision[strlen(revision)-2] = '\0';
  printf("run_tests %s (plCurve v. %s)\n", revision+11, version);
  check(strcmp(PACKAGE_VERSION,version) == 0);

  /* Check on plcl_build_vect */
  temp_vect.c[0] = 10.0*rand()/RAND_MAX;
  temp_vect.c[1] = 10.0*rand()/RAND_MAX;
  temp_vect.c[2] = 10.0*rand()/RAND_MAX;
  temp_vect2 = plcl_build_vect(temp_vect.c[0],temp_vect.c[1],temp_vect.c[2]);
  check(fabs(temp_vect.c[0] - temp_vect2.c[0]) < DBL_EPSILON);
  check(fabs(temp_vect.c[1] - temp_vect2.c[1]) < DBL_EPSILON);
  check(fabs(temp_vect.c[2] - temp_vect2.c[2]) < DBL_EPSILON);

  /* Test plCurve_build_color */
  temp_clr.r = 10.0*rand()/RAND_MAX;
  temp_clr.g = 10.0*rand()/RAND_MAX;
  temp_clr.b = 10.0*rand()/RAND_MAX;
  temp_clr.alpha = 10.0*rand()/RAND_MAX;
  temp_clr2 = plCurve_build_color(temp_clr.r,temp_clr.g,temp_clr.b,
                                  temp_clr.alpha);
  check(fabs(temp_clr.r - temp_clr2.r) < DBL_EPSILON);
  check(fabs(temp_clr.g - temp_clr2.g) < DBL_EPSILON);
  check(fabs(temp_clr.b - temp_clr2.b) < DBL_EPSILON);
  check(fabs(temp_clr.alpha - temp_clr2.alpha) < DBL_EPSILON);

  /* Check on vect_eq, plcl_random_vect and plcl_normalize */
  check(plcl_vecteq(temp_vect,temp_vect2));
  temp_vect2 = plcl_random_vect();
  check(!plcl_vecteq(temp_vect,temp_vect2));
  ok = true;
  dist = plcl_norm(plcl_normalize_vect(temp_vect2,&ok)) - 1.0;
  check(ok);
  check(fabs(dist) < DBL_EPSILON);
  (void)plcl_normalize_vect(zero_vect,&ok);
  check(!ok);

  /* Run the external tests */
  /* 1) Check to see if plcl_normalize_vect fails on zero vectors with no bool
   *    passed in, as advertised. 
   */
  sysret = system("./exit_failure_tests 1");
  check(WEXITSTATUS(sysret) == EXIT_FAILURE);
  /* 2) Make sure the plcl_component_div fails properly */
  sysret = system("./exit_failure_tests 2");
  check(WEXITSTATUS(sysret) == EXIT_FAILURE);
  /* 3) And make sure that the "failure program" isn't just failing */
  sysret = system("./exit_failure_tests 3");
  check(WEXITSTATUS(sysret) == EXIT_SUCCESS);

  /* Now to see if we can create plCurves.  Tests plCurve_new, plcl_vmadd, and
   * plCurve_fix_wrap. */
  S.nc = components;
  S.cp = (plCurve_strand *)calloc((size_t)2,sizeof(plCurve_strand));
  check(S.cp != NULL);
  S.quant = NULL;
  S.cst = NULL;
  S.cp[0].nv = nv[0];
  S.cp[0].open = open[0];
  S.cp[0].cc = cc[0];
  S.cp[0].vt = (plcl_vector *)calloc((size_t)5,sizeof(plcl_vector));
  check(S.cp[0].vt != NULL);
  S.cp[0].vt++;
  S.cp[0].clr = (plCurve_color *)calloc((size_t)1,sizeof(plCurve_color));
  check(S.cp[0].clr != NULL);
  S.cp[0].vt[0] = zero_vect;
  S.cp[0].vt[1] = plcl_build_vect(1.0,0.0,0.0);
  S.cp[0].vt[2] = plcl_build_vect(0.0,1.0,0.0);
  S.cp[0].vt[-1] = S.cp[0].vt[2];
  S.cp[0].vt[3] = S.cp[0].vt[0];
  S.cp[0].clr[0] = plCurve_build_color(1.0,0.0,1.0,1.0);
  S.cp[1].nv = nv[1];
  S.cp[1].open = open[1];
  S.cp[1].cc = cc[1];
  S.cp[1].vt = (plcl_vector *)calloc((size_t)6,sizeof(plcl_vector));
  check(S.cp[1].vt != NULL);
  S.cp[1].vt++;
  S.cp[1].clr = (plCurve_color *)calloc((size_t)4,sizeof(plCurve_color));
  check(S.cp[1].clr != NULL);
  S.cp[1].vt[0] = plcl_build_vect(1.0,1.0,2.0);
  S.cp[1].vt[1] = plcl_build_vect(1.0,1.0,1.0);
  S.cp[1].vt[2] = plcl_build_vect(1.0,1.0,0.0);
  S.cp[1].vt[3] = plcl_build_vect(1.0,1.0,-1.0);
  S.cp[1].vt[4] = S.cp[1].vt[2];
  S.cp[1].vt[-1] = S.cp[1].vt[1];
  S.cp[1].clr[0] = plCurve_build_color(0.00,1.0,0.0,1.0);
  S.cp[1].clr[1] = plCurve_build_color(0.33,1.0,0.0,1.0);
  S.cp[1].clr[0] = plCurve_build_color(0.67,1.0,0.0,1.0);
  S.cp[1].clr[1] = plCurve_build_color(1.00,1.0,0.0,1.0);

  L = plCurve_new(components,nv,open,cc);
  L->cp[0].vt[0] = zero_vect;
  L->cp[0].vt[1] = plcl_build_vect(1.0,0.0,0.0);
  L->cp[0].vt[2] = plcl_build_vect(0.0,1.0,0.0);
  L->cp[1].vt[0] = plcl_build_vect(1.0,1.0,2.0);
  L->cp[1].vt[1] = plcl_build_vect(1.0,1.0,1.0);
  L->cp[1].vt[2] = plcl_vmadd(L->cp[0].vt[1],1.0,L->cp[0].vt[2]);
  check(plcl_M_vecteq(L->cp[1].vt[2],plcl_build_vect(1.0,1.0,0.0)));
  L->cp[1].vt[3] = plcl_build_vect(1.0,1.0,-1.0);
  plCurve_fix_wrap(L);
  L->cp[0].clr[0] = plCurve_build_color(1.0,0.0,1.0,1.0);
  L->cp[1].clr[0] = plCurve_build_color(0.00,1.0,0.0,1.0);
  L->cp[1].clr[1] = plCurve_build_color(0.33,1.0,0.0,1.0);
  L->cp[1].clr[0] = plCurve_build_color(0.67,1.0,0.0,1.0);
  L->cp[1].clr[1] = plCurve_build_color(1.00,1.0,0.0,1.0);
  check(curves_match(S,*L));

  /* Check plcl_cross_prod */
  temp_vect = plcl_cross_prod(L->cp[0].vt[1],L->cp[0].vt[2]);
  check(fabs(temp_vect.c[0]) < DBL_EPSILON);
  check(fabs(temp_vect.c[1]) < DBL_EPSILON);
  check(fabs(temp_vect.c[2] - 1.0) < DBL_EPSILON);

  /* Check plcl_component_div (3 ways, the 4th is checked externally above). */
  ok = true;
  temp_vect = plcl_component_div(L->cp[1].vt[1],L->cp[1].vt[0],&ok);
  check(ok);
  check(plcl_vecteq(temp_vect,plcl_build_vect(1.0,1.0,0.5)));
  temp_vect = plcl_component_div(L->cp[1].vt[1],L->cp[0].vt[0],&ok);
  check(!ok);
  temp_vect = plcl_component_div(L->cp[1].vt[1],L->cp[1].vt[0],NULL);
  check(plcl_vecteq(temp_vect,plcl_build_vect(1.0,1.0,0.5)));

  /* Check plcl_vweighted */
  temp_vect = plcl_vweighted(0.7,L->cp[0].vt[1],L->cp[0].vt[2]);
  check(plcl_M_vecteq(temp_vect,plcl_build_vect(0.3,0.7,0.0)));

  /* Check plcl_dot_prod */
  check(fabs(plcl_dot_prod(L->cp[1].vt[0],L->cp[1].vt[3])) < DBL_EPSILON);

  /* Check plcl_distance and plcl_sq_dist */
  dist = plcl_distance(L->cp[1].vt[3],L->cp[0].vt[2]) - sqrt(2.0);
  check(fabs(dist) < DBL_EPSILON);
  dist = plcl_sq_dist(L->cp[0].vt[2],L->cp[1].vt[3]);
  check(fabs(dist - 2.0) < DBL_EPSILON);

  /* Check plCurve_num_edges */
  check(plCurve_num_edges(L) == 6);

  /* Check all the constraint-setting code!  Be very careful when editing 
   * the code between here and the "End of constraint-setting code checks"
   * as this code goes to great trouble to check all of the branches in that
   * code and that's a fairly fragile thing to do.
   */

  /* Start by doing nothing (it is one of the code paths :-) */
  plCurve_unconstrain(L,0,0,3);
  check(curves_match(S,*L));

  /* Allocate 10 constraint slots, we'll attach and detach them as needed */
  S.cst = (plCurve_constraint *)calloc((size_t)10,sizeof(plCurve_constraint));
  check(S.cst != NULL);
  check(S.cst->next == NULL);
  S.cst[0].kind = fixed;
  S.cst[0].vect[0] = L->cp[0].vt[0];
  S.cst[0].vect[1] = zero_vect;
  S.cst[0].cmp = 0;
  S.cst[0].vert = 0;
  S.cst[0].num_verts = 1;
  plCurve_set_fixed(L,0,0,L->cp[0].vt[0]);
  check(curves_match(S,*L));

  /*@-immediatetrans@*/
  S.cst[0].next = &(S.cst[1]);
  S.cst[1].kind = fixed;
  S.cst[1].vect[0] = L->cp[0].vt[0];
  S.cst[1].vect[1] = zero_vect;
  S.cst[1].cmp = 0;
  S.cst[1].vert = 1;
  S.cst[1].num_verts = 1;
  plCurve_set_fixed(L,0,1,L->cp[0].vt[0]);
  check(curves_match(S,*L));

  plCurve_set_fixed(L,0,1,L->cp[0].vt[0]);
  check(curves_match(S,*L));

  S.cst[0].next = NULL;
  S.cst[0].kind = line;
  S.cst[0].num_verts = 3;
  plCurve_constrain_to_line(L,0,0,3,L->cp[0].vt[0],zero_vect);
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
  plCurve_constrain_to_line(L,0,2,1,L->cp[0].vt[1],zero_vect);
  check(curves_match(S,*L));

  S.cst[0].next = NULL;
  plCurve_unconstrain(L,0,2,1);
  check(curves_match(S,*L));

  S.cst[0].num_verts = 3;
  plCurve_constrain_to_line(L,0,2,1,L->cp[0].vt[0],zero_vect);
  check(curves_match(S,*L));

  S.cst[0].num_verts = 1;
  plCurve_unconstrain(L,0,1,2);
  check(curves_match(S,*L));

  S.cst[0].next = &(S.cst[1]);
  S.cst[1] = S.cst[0];
  S.cst[1].vert = 2;
  S.cst[1].next = NULL;
  plCurve_constrain_to_line(L,0,2,1,L->cp[0].vt[0],zero_vect);
  check(curves_match(S,*L));

  S.cst[0].next = NULL;
  S.cst[0].num_verts = 3;
  plCurve_constrain_to_line(L,0,1,1,L->cp[0].vt[0],zero_vect);
  check(curves_match(S,*L));

  S.cst[0].num_verts = 1;
  S.cst[2] = S.cst[0];
  S.cst[2].vert = 2;
  S.cst[0].next = &(S.cst[1]);
  S.cst[1].next = &(S.cst[2]);
  S.cst[1].kind = plane;
  S.cst[1].vect[0] = L->cp[0].vt[0];
  S.cst[1].vert = 1;
  plCurve_constrain_to_plane(L,0,1,1,L->cp[0].vt[0],0.0);
  /*@-compmempass@*/
  check(curves_match(S,*L));

  S.cst[0].next = NULL;
  S.cst[0].num_verts = 3;
  plCurve_constrain_to_line(L,0,1,1,L->cp[0].vt[0],zero_vect);
  check(curves_match(S,*L));

  S.cst[0].next = &(S.cst[2]);
  S.cst[2].kind = line;
  S.cst[2].vect[0] = plcl_cross_prod(L->cp[0].vt[1],L->cp[0].vt[2]);
  check(plcl_M_vecteq(S.cst[2].vect[0],plcl_build_vect(0.0,0.0,1.0)));
  S.cst[2].vect[0] = L->cp[1].vt[2];
  S.cst[2].cmp = 1;
  S.cst[2].vert = 2;
  S.cst[2].num_verts = 2;
  plCurve_constrain_to_line(L,1,3,1,S.cst[2].vect[0],S.cst[2].vect[1]);
  plCurve_constrain_to_line(L,1,2,1,S.cst[2].vect[0],S.cst[2].vect[1]);
  check(curves_match(S,*L));

  /*@-nullret@*/
  S.cst[0] = S.cst[2];
  /*@=nullret@*/
  plCurve_unconstrain(L,0,0,3);
  check(curves_match(S,*L));

  S.cst[0].next = &(S.cst[2]);
  S.cst[0].kind = line;
  S.cst[0].vect[0] = S.cst[2].vect[0];
  S.cst[0].vect[1] = S.cst[0].vect[0];
  S.cst[0].cmp = 0;
  S.cst[0].vert = 2;
  S.cst[0].num_verts = 1;
  plCurve_constrain_to_line(L,0,2,1,S.cst[2].vect[0],S.cst[2].vect[0]);
  /*@-onlytrans@*/
  check(curves_match(S,*L));

  S.cst[1] = S.cst[0];
  /*@=onlytrans@*/
  S.cst[0].next = &(S.cst[1]);
  S.cst[0].vert = 0;
  plCurve_constrain_to_line(L,0,0,1,S.cst[2].vect[0],S.cst[2].vect[0]);
  check(curves_match(S,*L));

  S.cst[0].next = &(S.cst[2]);
  S.cst[0].num_verts = 3;
  S.cst[0].vect[1] = zero_vect;
  plCurve_constrain_to_line(L,0,0,3,S.cst[2].vect[0],zero_vect);
  check(curves_match(S,*L));

  S.cst[1] = S.cst[0];
  S.cst[0].next = &(S.cst[1]);
  S.cst[0].num_verts = 1;
  S.cst[1].vect[1] = S.cst[1].vect[0];
  S.cst[1].vert = 1;
  S.cst[1].num_verts = 2;
  plCurve_unconstrain(L,0,0,3);
  plCurve_constrain_to_line(L,0,0,1,S.cst[2].vect[0],zero_vect);
  plCurve_constrain_to_line(L,0,1,2,S.cst[0].vect[0],S.cst[0].vect[0]);
  check(curves_match(S,*L));

  S.cst[2].vert = 0;
  S.cst[2].num_verts = 4;
  plCurve_constrain_to_line(L,1,0,3,S.cst[2].vect[0],S.cst[2].vect[1]);
  check(curves_match(S,*L));

  S.cst[3] = S.cst[2];
  S.cst[2].next = &(S.cst[3]);
  S.cst[2].num_verts = 2;
  S.cst[3].vert = 3;
  S.cst[3].num_verts = 1;
  plCurve_unconstrain(L,1,2,1);
  check(curves_match(S,*L));

  /* End of constraint-setting code checks */

  /* Remove a specific type of constraint from all components (we check 
   * our code to remove all constraints just before cleanup below). */
  S.cst[0] = S.cst[1];
  S.cst[0].next = NULL;
  vert = plCurve_remove_constraint(L,line,S.cst[2].vect);
  check(vert == 4);
  check(curves_match(S,*L));

  /* Now set somewhat reasonable constraints and get ready to test reading and
   * writing. */
  S.cst[0].kind = plane;
  plcl_M_cross(S.cst[0].vect[0],S.cp[0].vt[1],S.cp[0].vt[2]);
  S.cst[0].vect[1] = S.cp[0].vt[0];
  S.cst[0].cmp = 0;
  S.cst[0].vert = 0;
  S.cst[0].num_verts = 1;
  S.cst[1] = S.cst[0];
  S.cst[1].vert = 2;
  S.cst[2].kind = fixed;
  S.cst[2].vect[0] = plcl_build_vect(2.0,2.0,1.0);
  S.cst[2].vect[1] = zero_vect;
  S.cst[2].cmp = 1;
  S.cst[2].vert = 0;
  S.cst[2].num_verts = 1;
  S.cst[3].kind = line;
  S.cst[3].vect[0] = S.cst[0].vect[0];
  S.cst[3].vect[1] = plcl_build_vect(2.0,1.0,3.7);
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
  plCurve_constrain_to_plane(L,0,0,3,S.cst[0].vect[0],0.0);
  plCurve_unconstrain(L,0,1,1);
  plCurve_set_fixed(L,1,0,S.cst[2].vect[0]);
  plCurve_constrain_to_line(L,1,1,2,S.cst[3].vect[0],S.cst[3].vect[1]);
  plCurve_set_fixed(L,1,3,S.cst[4].vect[0]);
  check(curves_match(S,*L));
  /*@=usereleased*/

  /*@=immediatetrans@*/

  /* Now plCurve_check_cst */
  dist = plCurve_check_cst(L) - sqrt(3.0);
  check(fabs(dist) < DBL_EPSILON);
  plCurve_fix_cst(L);
  dist = plCurve_check_cst(L);
  check(fabs(dist) < DBL_EPSILON);
  S.cp[1].vt[0] = S.cst[2].vect[0];
  S.cp[1].vt[3] = S.cst[4].vect[0];
  S.cp[1].vt[1] = plcl_vect_sum(S.cp[1].vt[1],S.cp[0].vt[1]);
  S.cp[1].vt[2] = plcl_vect_sum(S.cp[1].vt[2],S.cp[0].vt[1]);
  S.cp[1].vt[-1] = S.cp[1].vt[1];
  S.cp[1].vt[4] = S.cp[1].vt[2];
  check(curves_match(S,*L));

  /* Now check _write and _read */ 
  outname = tmpnam(NULL);
  check(outname != NULL);
  outfile = fopen(outname,"w");
  check(outfile != NULL);
  plCurve_write(outfile,L);
  sysret = fclose(outfile);
  check(sysret == 0);
  plCurve_free(L);
  L = NULL;
  outfile = fopen(outname,"r");
  check(outfile != NULL);
  L = plCurve_read(outfile,&err_num,err_str,sizeof(err_str));
  if (err_num != 0) {
    printf("Trouble reading file: %s\n",err_str);
  }
  check(err_num == 0);
  check(L != NULL);
  sysret = fclose(outfile);
  outfile = NULL;
  check(sysret == 0);
  sysret = unlink(outname);
  check(sysret == 0);
  plCurve_write(stdout,&S);
  plCurve_write(stdout,L);
  check(curves_match(S,*L));

  /*@=compmempass@*/

#if 0

  /* Try forcing the second component to close */
  S.cp[1].vt[0].c[2] -= 1.0/2.0;
  S.cp[1].vt[1].c[2] -= 1.0/6.0;
  S.cp[1].vt[2].c[2] += 1.0/6.0;
  S.cp[1].open = false;
  S.cp[1].nv--;
  S.cp[1].vt[-1] = S.cp[1].vt[2];
  S.cp[1].vt[3] = S.cp[1].vt[0];
  plCurve_force_closed(L);
  check(curves_match(S,*L));

  list_csts(L);
  cst = (plCurve_constraint *)calloc((size_t)1,sizeof(plCurve_constraint));
  check(cst != NULL);
  cst->kind = line;
  cst->vect[0] = plcl_build_vect(1.0,0.0,0.0);
  cst->vect[1] = cst->vect[0];
  cst->cmp = 0;
  cst->vert = 1;
  cst->num_verts = 1;
  cst->next = S.cst;
  S.cst = cst;
  plCurve_constrain_to_plane(L,0,0,3,plcl_build_vect(1.0,1.0,1.0),2.3);
  plCurve_constrain_to_line(L,0,1,1,S.cst->vect[0],S.cst->vect[1]);
  check(L->cst != NULL &&
         L->cst->kind == plane &&
         L->cst->cmp == 0);
  check(L->cst->next != NULL &&
         L->cst->next->kind == line &&
         L->cst->next->cmp == 0);
  check(L->cst->next->next != NULL &&
         L->cst->next->next->kind == plane &&
         L->cst->next->next->cmp == 0);
  check(L->cst->next->next->next != NULL &&
         L->cst->next->next->next->kind == fixed &&
         L->cst->next->next->next->cmp == 1);
  check(L->cst->next->next->next->next != NULL &&
         L->cst->next->next->next->next->kind == fixed &&
         L->cst->next->next->next->next->cmp == 1);
  check(L->cst->next->next->next->next->next == NULL);
  /* Should remove both first and third constraints */
  vert = plCurve_remove_constraint(L,plane,L->cst->vect);
  check(vert == 2);
  check(curves_match(S,*L));


#endif

  /* Now check to see if we can remove all constraints */
  /*@-kepttrans@*/
  free(S.cst);
  /*@=kepttrans@*/
  S.cst = NULL;
  plCurve_remove_all_constraints(L);
  check(curves_match(S,*L));

  /* Eventually quantifier testing code goes here */
  check(L->quant == NULL);
  L->quant = (plCurve_vert_quant *)calloc((size_t)1,sizeof(plCurve_vert_quant));
  check(L->quant != NULL);
  check(L->quant->next == NULL);

  /* Cleanup phase */
  plCurve_free(L);
  L = NULL;
  /* And just to make sure that it works */
  plCurve_free(L);

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
}
