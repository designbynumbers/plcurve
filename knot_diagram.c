#include <config.h>
#include <plCurve.h>

#ifdef HAVE_STDIO_H
  #include <stdio.h>
#endif
#ifdef HAVE_STDLIB_H
  #include <stdlib.h>
#endif
#ifdef HAVE_FLOAT_H
  #include <float.h>
#endif
#ifdef HAVE_MATH_H
  #include <math.h>
#endif
#ifdef HAVE_ASSERT_H
  #include <assert.h>
#endif
#ifdef HAVE_MALLOC_H
  #include <malloc.h>
#endif

/* Project the vertices onto the plane two which N is normal (returning them
 * as (x,y,0) vectors in a plCurve framework. */
plCurve *flatten(plCurve *L, plcl_vector N) {
  plCurve *F;
  int cmp;
  int vert;
  plcl_vector n;
  plcl_vector *vt;
  double a,ac,b,bc,r2;
  plcl_vector new_vect,cross;
  plcl_vector e3;
  double x,y,z;
  double cos_theta,sin_theta;

  F = plCurve_copy(L);
  n = plcl_normalize_vect(N,NULL); /* Bail out if we give it a bad vector */
  e3 = plcl_build_vect(0.0,0.0,1.0);

  if (n.c[1]*n.c[1]+n.c[2]*n.c[2] > DBL_EPSILON) {
    a = -n.c[1]/sqrt(n.c[1]*n.c[1]+n.c[2]*n.c[2]);
  } else {
    a = 1.0;
  }
  ac = sqrt(1.0 - a*a);
  r2 = a*n.c[1] - ac*n.c[2];
  /* Since n has length 1 (see above), this denominator isn't 0 */
  assert(n.c[0]*n.c[0] + r2*r2 > DBL_EPSILON);
  b = -n.c[0]/sqrt(n.c[0]*n.c[0]+r2*r2);
  bc = sqrt(1.0 - b*b);

  for (cmp = 0; cmp < F->nc; cmp++) {
    for (vert = 0; vert < F->cp[cmp].nv; vert++) {
      vt = &F->cp[cmp].vt[vert];
//    /* First, project onto the plane itself. */
//    *vt = plcl_cross_prod(n,plcl_cross_prod(*vt,n));
//    printf("%g %g %g -> %g %g %g\n",plcl_M_clist(L->cp[cmp].vt[vert]),
//        plcl_M_clist(F->cp[cmp].vt[vert]));
//    /* Then rotate that to the x,y plane */
//    new_vect.c[0] = b*(*vt).c[0] - a*bc*(*vt).c[1] + ac*bc*(*vt).c[2];
//    new_vect.c[1] = bc*(*vt).c[0] + a*b*(*vt).c[1] - b*ac*(*vt).c[2];
//    new_vect.c[2] = ac*(*vt).c[1] + a*(*vt).c[2];
//    *vt = new_vect;
      cross = plcl_cross_prod(n,e3);
      x = cross.c[0];
      y = cross.c[1];
      z = cross.c[2];
      cos_theta = plcl_dot_prod(n,e3);
      sin_theta = sqrt(1-cos_theta*cos_theta);
      /* Rotate space so that n is pointing along the z axis and just take 
         the first two components.  Formulas from
           http://en.wikipedia.org/wiki/Coordinate_rotation#Three_dimensions
       */
      new_vect.c[0] = 
        (cos_theta+(1-cos_theta)*x*x)*(*vt).c[0]+
        ((1-cos_theta)*x*y-sin_theta*z)*(*vt).c[1]+
        ((1-cos_theta)*x*z+sin_theta*y)*(*vt).c[2];
      new_vect.c[1] = 
        ((1-cos_theta)*y*z+sin_theta*z)*(*vt).c[0]+
        (cos_theta+(1-cos_theta)*y*y)*(*vt).c[1]+
        ((1-cos_theta)*y*z-sin_theta*x)*(*vt).c[2];
      new_vect.c[2] = 0.0;

      *vt = new_vect;

      printf("%g %g %g -> %g %g %g\n",plcl_M_clist(L->cp[cmp].vt[vert]),
          plcl_M_clist(F->cp[cmp].vt[vert]));
    }
  }

  return F;
}

/* Add a component which contains the convex hull of the other components
 * *which are presumed to be 2-d* */
void add_convex_hull(plCurve *L) {
  plCurve_strand *cp;
  plcl_vector *vt;
  int nv;
  int cmp,vert;
  plcl_vector first,next;

  for (nv = 0, cmp = 0; cmp < L->nc; cmp++) {
    nv += L->cp[cmp].nv;
  }
  vt = malloc(nv*sizeof(plcl_vector));
  nv = 0;

  for (cmp = 0; cmp < L->nc; cmp++) {
    for (vert = 0; vert < L->cp[cmp].nv; vert++) {
      /* Look for the lowest point */
    }
  }
  /* Area of the triangle made by three vertices.  Formula adapted from 
     http://www.cse.unsw.edu.au/~lambert/java/3d/triangle.html
   */
  // plCurve_add_component(L,nv,false,cc,vt,clr);
}

int main() {
  FILE *vectfile;
  plCurve *L, *F;
  char err_str[80];
  int err_num;

  vectfile = fopen("kl_3_1_I.vect","r");
  assert(vectfile != NULL);
  L = plCurve_read(vectfile,&err_num,err_str,80);
  fclose(vectfile);
  if (err_num != 0) {
    fprintf(stderr,"Error reading file %s:\n%s\n","kl_3_1_I.vect",err_str);
    exit(EXIT_FAILURE);
  }
  F = flatten(L,plcl_random_vect());
  assert(F != NULL);
  add_convex_hull(F);
  vectfile = fopen("kl_3_1_flat.vect","w");
  assert(vectfile != NULL);
  plCurve_write(vectfile,F);
  fclose(vectfile);
  exit(EXIT_SUCCESS);
}
