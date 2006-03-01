/*
 * $Id: run_tests.c,v 1.8 2006-03-01 15:51:05 ashted Exp $
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
    fprintf(stderr,"%s:%d: failed requirement `%s'\n",__FILE__,__LINE__,#C); \
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
int main(void) {
  plCurve S; /* The standard against which to measure */
  plCurve *L;
  int nv[components] = { 3, 4 };
  bool open[components] = { false, true };
  int cc[components] = { 1, 2 };
  char version[80];
  char revision[] = "$Revision: 1.8 $";
  plCurve_constraint *cst;
  plCurve_vert_quant *quant;
  int vert;
  double dist;
  bool ok;
  plcl_vector temp_vect;
  int sysret;
  FILE *outfile;
  char *outname;
  int err_num;
  char err_str[80];

  plCurve_version(NULL,0);
  version[0] = '\0';
  plCurve_version(version,sizeof(version));
  revision[strlen(revision)-2] = '\0';
  printf("run_tests %s (plCurve v. %s)\n", revision+11, version);
  assert(strcmp(PACKAGE_VERSION,version) == 0);

  S.nc = 2;
  S.cp = (plCurve_strand *)calloc((size_t)2,sizeof(plCurve_strand));
  assert(S.cp != NULL);
  S.quant = NULL;
  S.cst = NULL;
  S.cp[0].nv = nv[0];
  S.cp[0].open = open[0];
  S.cp[0].cc = cc[0];
  S.cp[0].vt = (plcl_vector *)calloc((size_t)5,sizeof(plcl_vector));
  assert(S.cp[0].vt != NULL);
  S.cp[0].vt++;
  S.cp[0].clr = (plCurve_color *)calloc((size_t)1,sizeof(plCurve_color));
  assert(S.cp[0].clr != NULL);
  S.cp[0].vt[-1].c[0] = 0.0;
  S.cp[0].vt[-1].c[1] = 1.0;
  S.cp[0].vt[-1].c[2] = 0.0;
  S.cp[0].vt[0].c[0] = 0.0;
  S.cp[0].vt[0].c[1] = 0.0;
  S.cp[0].vt[0].c[2] = 0.0;
  S.cp[0].vt[1].c[0] = 1.0;
  S.cp[0].vt[1].c[1] = 0.0;
  S.cp[0].vt[1].c[2] = 0.0;
  S.cp[0].vt[2].c[0] = 0.0;
  S.cp[0].vt[2].c[1] = 1.0;
  S.cp[0].vt[2].c[2] = 0.0;
  S.cp[0].vt[3].c[0] = 0.0;
  S.cp[0].vt[3].c[1] = 0.0;
  S.cp[0].vt[3].c[2] = 0.0;
  S.cp[0].clr[0].r = 1.0;
  S.cp[0].clr[0].g = 0.0;
  S.cp[0].clr[0].b = 1.0;
  S.cp[0].clr[0].alpha = 1.0;
  S.cp[1].nv = nv[1];
  S.cp[1].open = open[1];
  S.cp[1].cc = cc[1];
  S.cp[1].vt = (plcl_vector *)calloc((size_t)6,sizeof(plcl_vector));
  assert(S.cp[1].vt != NULL);
  S.cp[1].vt++;
  S.cp[1].clr = (plCurve_color *)calloc((size_t)2,sizeof(plCurve_color));
  assert(S.cp[1].clr != NULL);
  S.cp[1].vt[-1].c[0] = 1.0;
  S.cp[1].vt[-1].c[1] = 1.0;
  S.cp[1].vt[-1].c[2] = 1.0;
  S.cp[1].vt[0].c[0] = 1.0;
  S.cp[1].vt[0].c[1] = 1.0;
  S.cp[1].vt[0].c[2] = 2.0;
  S.cp[1].vt[1].c[0] = 1.0;
  S.cp[1].vt[1].c[1] = 1.0;
  S.cp[1].vt[1].c[2] = 1.0;
  S.cp[1].vt[2].c[0] = 1.0;
  S.cp[1].vt[2].c[1] = 1.0;
  S.cp[1].vt[2].c[2] = 0.0;
  S.cp[1].vt[3].c[0] = 1.0;
  S.cp[1].vt[3].c[1] = 1.0;
  S.cp[1].vt[3].c[2] = -1.0;
  S.cp[1].vt[4].c[0] = 1.0;
  S.cp[1].vt[4].c[1] = 1.0;
  S.cp[1].vt[4].c[2] = 0.0;
  S.cp[1].clr[0].r = 0.0;
  S.cp[1].clr[0].g = 1.0;
  S.cp[1].clr[0].b = 0.0;
  S.cp[1].clr[0].alpha = 1.0;
  S.cp[1].clr[1].r = 1.0;
  S.cp[1].clr[1].g = 1.0;
  S.cp[1].clr[1].b = 0.0;
  S.cp[1].clr[1].alpha = 1.0;

  L = plCurve_new(components,nv,open,cc);
  L->cp[0].vt[0] = plcl_build_vect(0.0,0.0,0.0);
  L->cp[0].vt[1] = plcl_build_vect(1.0,0.0,0.0);
  L->cp[0].vt[2] = plcl_build_vect(0.0,1.0,0.0);
  L->cp[1].vt[0] = plcl_build_vect(1.0,1.0,2.0);
  L->cp[1].vt[1] = plcl_build_vect(1.0,1.0,1.0);
  L->cp[1].vt[2] = plcl_vmadd(L->cp[0].vt[1],1.0,L->cp[0].vt[2]);
  assert(plcl_M_vecteq(L->cp[1].vt[2],plcl_build_vect(1.0,1.0,0.0)));
  L->cp[1].vt[3] = plcl_build_vect(1.0,1.0,-1.0);
  plCurve_fix_wrap(L);
  L->cp[0].clr[0] = plCurve_build_color(1.0,0.0,1.0,1.0);
  L->cp[1].clr[0] = plCurve_build_color(0.0,1.0,0.0,1.0);
  L->cp[1].clr[1] = plCurve_build_color(1.0,1.0,0.0,1.0);
  assert(curves_match(S,*L));

  temp_vect = plcl_cross_prod(L->cp[0].vt[1],L->cp[0].vt[2]);
  assert(fabs(temp_vect.c[0]) < DBL_EPSILON);
  assert(fabs(temp_vect.c[1]) < DBL_EPSILON);
  assert(fabs(temp_vect.c[2] - 1.0) < DBL_EPSILON);

  ok = true;
  temp_vect = plcl_component_div(L->cp[1].vt[1],L->cp[1].vt[0],&ok);
  assert(ok);
  assert(plcl_vecteq(temp_vect,plcl_build_vect(1.0,1.0,0.5)));
  temp_vect = plcl_component_div(L->cp[1].vt[1],L->cp[0].vt[0],&ok);
  assert(!ok);
  temp_vect = plcl_component_div(L->cp[1].vt[1],L->cp[1].vt[0],NULL);
  assert(plcl_vecteq(temp_vect,plcl_build_vect(1.0,1.0,0.5)));

  temp_vect = plcl_vweighted(0.7,L->cp[0].vt[1],L->cp[0].vt[2]);
  assert(plcl_M_vecteq(temp_vect,plcl_build_vect(0.3,0.7,0.0)));

  /* Check to see if plcl_normalize_vect fails on zero vectors with no bool
   * passed in, as advertised. */
  sysret = system("./exit_failure_tests 1");
  assert(WEXITSTATUS(sysret) == EXIT_FAILURE);
  /* Make sure the plcl_component_div fails properly */
  sysret = system("./exit_failure_tests 2");
  assert(WEXITSTATUS(sysret) == EXIT_FAILURE);
  /* And make sure that the "failure program" isn't just failing */
  sysret = system("./exit_failure_tests 3");
  assert(WEXITSTATUS(sysret) == EXIT_SUCCESS);

  assert(fabs(plcl_dot_prod(L->cp[1].vt[0],L->cp[1].vt[3])) < DBL_EPSILON);

  dist = plcl_distance(L->cp[1].vt[3],L->cp[0].vt[2]) - sqrt(2.0);
  assert(fabs(dist) < DBL_EPSILON);

  dist = plcl_sq_dist(L->cp[0].vt[2],L->cp[1].vt[3]);
  assert(fabs(dist - 2.0) < DBL_EPSILON);

  S.cst = (plCurve_constraint *)calloc((size_t)1,sizeof(plCurve_constraint));
  assert(S.cst != NULL);
  assert(S.cst->next == NULL);
  S.cst->kind = on_line;
  S.cst->vect[0].c[0] = 0.0;
  S.cst->vect[0].c[1] = 0.0;
  S.cst->vect[0].c[2] = 1.0;
  S.cst->vect[1].c[0] = 1.0;
  S.cst->vect[1].c[1] = 2.0;
  S.cst->vect[1].c[2] = 0.0;
  S.cst->cmp = 1;
  S.cst->vert = 0;
  S.cst->num_verts = 4;
  plCurve_constrain_to_line(L,1,0,4,plcl_build_vect(0.0,0.0,1.0),
                                    plcl_build_vect(1.0,2.0,0.0));
  assert(curves_match(S,*L));

  cst = (plCurve_constraint *)calloc((size_t)1,sizeof(plCurve_constraint));
  assert(cst != NULL);
  assert(cst->next == NULL);
  cst->kind = in_plane;
  cst->vect[0].c[0] = 0.0;
  cst->vect[0].c[1] = 0.0;
  cst->vect[0].c[2] = 1.0;
  cst->vect[1].c[0] = 0.0;
  cst->vect[1].c[1] = 0.0;
  cst->vect[1].c[2] = 0.0;
  cst->cmp = 0;
  cst->vert = 0;
  cst->num_verts = 3;
  cst->next = S.cst;
  S.cst = cst;
  /* Build this one in two pieces */
  plCurve_constrain_to_plane(L,0,0,2,plcl_build_vect(0.0,0.0,1.0),0.0);
  plCurve_constrain_to_plane(L,0,1,2,plcl_build_vect(0.0,0.0,1.0),0.0);
  assert(curves_match(S,*L));

  /* Now knock a hole in the component 1 constraint and fill it up again */
  plCurve_constrain_to_line(L,1,2,1,plcl_build_vect(0.0,0.0,1.0),
                                    plcl_build_vect(1.0,1.0,2.0));
  plCurve_constrain_to_line(L,1,1,3,plcl_build_vect(0.0,0.0,1.0),
                                    plcl_build_vect(1.0,2.0,0.0));
  /* list_csts(&S) */
  /* list_csts(L) */

  assert(curves_match(S,*L));

  dist = plCurve_check_cst(L);
  assert(fabs(dist - 1.0) < DBL_EPSILON);
  plCurve_fix_cst(L);
  dist = plCurve_check_cst(L);
  assert(fabs(dist) < DBL_EPSILON);
  for (vert = 0; vert < S.cp[1].nv; vert++) {
    L->cp[1].vt[vert] =
      plcl_vect_sum(L->cp[1].vt[vert],plcl_build_vect(0.0,-1.0,0.0));
/*  printf("%g %g %g %g %g %g\n",
      plcl_M_clist(L->cp[1].vt[vert]),plcl_M_clist(S.cp[1].vt[vert])); */
  }
  assert(curves_match(S,*L));

  S.cst->next->kind = in_plane;
  S.cst->next->vect[0] = plcl_build_vect(0.0,2.0,0.0);
  S.cst->next->vect[1] = plcl_build_vect(2.0,0.0,0.0);
  plCurve_constrain_to_plane(L,1,0,4,plcl_build_vect(0.0,2.0,0.0),2.0);
  assert(curves_match(S,*L));

  dist = plCurve_check_cst(L);
  assert(fabs(dist) < DBL_EPSILON);

  S.cst->next->kind = in_plane;
  S.cst->next->vect[0] = plcl_build_vect(0.2,0.0,0.0);
  S.cst->next->vect[1] = plcl_build_vect(0.4,0.0,0.0);
  plCurve_constrain_to_plane(L,1,0,4,plcl_build_vect(0.2,0.0,0.0),0.4);
  assert(curves_match(S,*L));

  dist = plCurve_check_cst(L);
  assert(fabs(dist - 1.0) < DBL_EPSILON);
  plCurve_fix_cst(L);
  for (vert = 0; vert < S.cp[1].nv; vert++) {
    L->cp[1].vt[vert] =
      plcl_vect_diff(L->cp[1].vt[vert],plcl_build_vect(1.0,0.0,0.0));
  }
  assert(curves_match(S,*L));

  dist = plcl_norm(S.cp[1].vt[1]) - sqrt(3.0);
  assert(fabs(dist) < DBL_EPSILON);
  ok = true;
  dist = plcl_norm(plcl_normalize_vect(plcl_random_vect(),&ok)) - 1.0;
  assert(ok);
  assert(fabs(dist) < DBL_EPSILON);
  (void)plcl_normalize_vect(plcl_build_vect(0.0,0.0,0.0),&ok);
  assert(!ok);

  cst = (plCurve_constraint *)calloc((size_t)1,sizeof(plCurve_constraint));
  assert(cst != NULL);
  assert(cst->next == NULL);
  cst->kind = fixed;
  cst->vect[0] = plcl_vlincomb(3.0,L->cp[1].vt[1],-1.0,L->cp[1].vt[0]);
  /* vect[1] is already set to zeros by calloc */
  cst->cmp = 1;
  cst->vert = 0;
  cst->num_verts = 1;
  cst->next = S.cst->next;
  S.cst->next = cst;
  cst = S.cst->next->next;
  cst->vert++;
  cst->num_verts--;
  plCurve_set_fixed(L,1,0,plcl_build_vect(2.0,2.0,1.0));
  assert(curves_match(S,*L));

  cst = (plCurve_constraint *)calloc((size_t)1,sizeof(plCurve_constraint));
  assert(cst != NULL);
  assert(cst->next == NULL);
  cst->kind = fixed;
  cst->vect[0] = plcl_build_vect(2.0,2.0,0.0);
  /* vect[1] is already set to zeros by calloc */
  cst->cmp = 1;
  cst->vert = 3;
  cst->num_verts = 1;
  S.cst->next->next->next = cst;
  cst = S.cst->next->next;
  cst->num_verts--;
  plCurve_set_fixed(L,1,0,plcl_build_vect(2.0,2.0,1.0));
  plCurve_set_fixed(L,1,3, /* (2,2,1) times (1,1,0) componentwise */
      plcl_component_mult(S.cst->next->vect[0],L->cp[1].vt[2]));
  assert(curves_match(S,*L));

  /* Test read and write */
  plCurve_unconstrain(L,0,1,1);
  list_csts(L);
  outname = tmpnam(NULL);
  assert(outname != NULL);
  outfile = fopen(outname,"w");
  assert(outfile != NULL);
  plCurve_write(outfile,L);
  plCurve_write(stdout,L);
  sysret = fclose(outfile);
  assert(sysret == 0);
  plCurve_free(L);
  L = NULL;
  outfile = fopen(outname,"r");
  assert(outfile != NULL);
  L = plCurve_read(outfile,&err_num,err_str,sizeof(err_str));
  if (err_num != 0) {
    printf("Trouble reading file: %s\n",err_str);
  }
  assert(L != NULL);
  sysret = fclose(outfile);
  assert(sysret == 0);
  sysret = unlink(outname);
  assert(sysret == 0);
  assert(curves_match(S,*L));

  cst = S.cst->next->next;
  S.cst->next->next = S.cst->next->next->next;
  free(cst);
  cst = NULL;
  assert(cst == NULL);
  plCurve_unconstrain(L,1,1,2);
  assert(curves_match(S,*L));

  S.cp[1].vt[0] = S.cst->next->vect[0];
  S.cp[1].vt[3] = S.cst->next->next->vect[0];
  dist = plCurve_check_cst(L) - sqrt(3.0);
  assert(fabs(dist) < DBL_EPSILON);
  plCurve_fix_cst(L);
  assert(curves_match(S,*L));

  cst = S.cst;
  S.cst = S.cst->next;
  free(cst);
  cst = NULL;
  assert(cst == NULL);
  plCurve_unconstrain(L,0,0,3);
  assert(curves_match(S,*L));

  /* Try forcing the second component to close */
  S.cp[1].vt[0].c[2] -= 1.0/2.0;
  S.cp[1].vt[1].c[2] -= 1.0/6.0;
  S.cp[1].vt[2].c[2] += 1.0/6.0;
  S.cp[1].open = false;
  S.cp[1].nv--;
  S.cp[1].vt[-1] = S.cp[1].vt[2];
  S.cp[1].vt[3] = S.cp[1].vt[0];
  plCurve_force_closed(L);
  assert(curves_match(S,*L));

  list_csts(L);
  cst = (plCurve_constraint *)calloc((size_t)1,sizeof(plCurve_constraint));
  assert(cst != NULL);
  cst->kind = on_line;
  cst->vect[0] = plcl_build_vect(1.0,0.0,0.0);
  cst->vect[1] = cst->vect[0];
  cst->cmp = 0;
  cst->vert = 1;
  cst->num_verts = 1;
  cst->next = S.cst;
  S.cst = cst;
  plCurve_constrain_to_plane(L,0,0,3,plcl_build_vect(1.0,1.0,1.0),2.3);
  plCurve_constrain_to_line(L,0,1,1,S.cst->vect[0],S.cst->vect[1]);
  assert(L->cst != NULL &&
         L->cst->kind == in_plane &&
         L->cst->cmp == 0);
  assert(L->cst->next != NULL &&
         L->cst->next->kind == on_line &&
         L->cst->next->cmp == 0);
  assert(L->cst->next->next != NULL &&
         L->cst->next->next->kind == in_plane &&
         L->cst->next->next->cmp == 0);
  assert(L->cst->next->next->next != NULL &&
         L->cst->next->next->next->kind == fixed &&
         L->cst->next->next->next->cmp == 1);
  assert(L->cst->next->next->next->next != NULL &&
         L->cst->next->next->next->next->kind == fixed &&
         L->cst->next->next->next->next->cmp == 1);
  assert(L->cst->next->next->next->next->next == NULL);
  /* Should remove both first and third constraints */
  vert = plCurve_remove_constraint(L,in_plane,L->cst->vect);
  assert(vert == 2);
  assert(curves_match(S,*L));

  assert(plCurve_num_edges(L) == 6);

  while(S.cst != NULL) {
    cst = S.cst;
    S.cst = cst->next;
    free(cst);
    cst = NULL;
  }
  plCurve_remove_all_constraints(L);
  assert(curves_match(S,*L));

  plCurve_constrain_to_plane(L,0,0,2,plcl_build_vect(1.0,2.0,3.0),4.0);
  plCurve_constrain_to_plane(L,0,1,2,plcl_build_vect(1.0,2.0,3.0),4.0);
  plCurve_constrain_to_plane(L,0,1,1,plcl_build_vect(1.0,2.0,3.0),4.0);
  S.cst = (plCurve_constraint *)calloc((size_t)1,sizeof(plCurve_constraint));
  assert(S.cst != NULL);
  S.cst->kind = in_plane;
  S.cst->vect[0] = plcl_build_vect(1.0,2.0,3.0);
  S.cst->vect[1] = plcl_build_vect(4.0,0.0,0.0);
  S.cst->cmp = 0;
  S.cst->vert = 0;
  S.cst->num_verts = 3;
  assert(curves_match(S,*L));

  /* Eventually quantifier testing code goes here */
  assert(L->quant == NULL);
  L->quant = (plCurve_vert_quant *)calloc((size_t)1,sizeof(plCurve_vert_quant));
  assert(L->quant != NULL);
  assert(L->quant->next == NULL);

  /* Cleanup phase */
  plCurve_free(L);
  L = NULL;
  /* And just to make sure that it works */
  plCurve_free(L);

  while(S.cst != NULL) {
    cst = S.cst;
    S.cst = cst->next;
    free(cst);
    cst = NULL;
  }
  assert(S.cst == NULL);
  while(S.quant != NULL) {
    quant = S.quant;
    S.quant = quant->next;
    free(quant);
    quant = NULL;
  }
  assert(S.quant == NULL);
  assert(S.cp[0].vt != NULL);
  S.cp[0].vt--;
  free(S.cp[0].vt);
  S.cp[0].vt = NULL;
  assert(S.cp[0].vt == NULL); /* For the sake of splint */
  assert(S.cp[0].clr != NULL);
  free(S.cp[0].clr);
  S.cp[0].clr = NULL;
  assert(S.cp[0].clr == NULL); /* For the sake of splint */
  assert(S.cp[1].vt != NULL);
  S.cp[1].vt--;
  free(S.cp[1].vt);
  S.cp[1].vt = NULL;
  assert(S.cp[1].vt == NULL); /* For the sake of splint */
  assert(S.cp[1].clr != NULL);
  free(S.cp[1].clr);
  S.cp[1].clr = NULL;
  assert(S.cp[1].clr == NULL); /* For the sake of splint */
  free(S.cp);
  return(EXIT_SUCCESS);
}
