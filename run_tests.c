/*
 * $Id: run_tests.c,v 1.2 2006-02-26 02:33:21 ashted Exp $
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


int main(void) {
  plCurve S; /* The standard against which to measure */
  plCurve *L;
#define components 2
  int nv[components] = { 3, 4 };
  bool open[components] = { false, true };
  int cc[components] = { 1, 2 };
  char version[80];
  char revision[] = "$Revision: 1.2 $";

  plCurve_version(NULL,0);
  version[0] = '\0';
  plCurve_version(version,sizeof(version));
  revision[strlen(revision)-2] = '\0';
  printf("run_tests %s (plCurve v. %s)\n", revision+11, version);
  assert(strcmp(PACKAGE_VERSION,version) == 0);

  S.nc = 2;
  S.cp = (plCurve_strand *)calloc((size_t)2,sizeof(plCurve_strand));
  assert(S.cp != NULL);
  S.cst = NULL;
  S.quant = NULL;
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
  L->cp[1].vt[2] = plcl_build_vect(1.0,1.0,0.0);
  L->cp[1].vt[3] = plcl_build_vect(1.0,1.0,-1.0);
  plCurve_fix_wrap(L);
  L->cp[0].clr[0] = plCurve_build_color(1.0,0.0,1.0,1.0);
  L->cp[1].clr[0] = plCurve_build_color(0.0,1.0,0.0,1.0);
  L->cp[1].clr[1] = plCurve_build_color(1.0,1.0,0.0,1.0);
  
  assert(curves_match(S,*L));

  plCurve_free(L);

  /* Cleanup phase */
  assert(S.cp[0].vt != NULL);
  S.cp[0].vt--;
  free(S.cp[0].vt);
  S.cp[0].vt = NULL;
  assert(S.cp[0].clr != NULL);
  free(S.cp[0].clr);
  S.cp[0].clr = NULL;
  free(S.cp);
  return(EXIT_SUCCESS);
}
