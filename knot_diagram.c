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

/* Project the vertices onto the plane to which N is normal (returning them
 * as (x,y,0) vectors in a plCurve framework. */
static plCurve *flatten(plCurve *L, plcl_vector N) {
  plCurve *F;
  int cmp;
  int vert;
  plcl_vector n;
  plcl_vector *vt;
  plcl_vector e3;
  double cos_theta,sin_theta;
  plcl_vector axle,u_para,u_perp,v_cross_u;
  bool ok;

  F = plCurve_copy(L);
  n = plcl_normalize_vect(N,NULL); /* Bail out if we give it a bad vector */
  e3 = plcl_build_vect(0.0,0.0,1.0);
  /* What should we rotate around? */
  axle = plcl_cross_prod(n,e3);
  ok = true;
  axle = plcl_normalize_vect(axle,&ok);
  if (ok) { 
    /* We actually do need to do some rotation */
    for (cmp = 0; cmp < F->nc; cmp++) {
      for (vert = 0; vert < F->cp[cmp].nv; vert++) {
        vt = &F->cp[cmp].vt[vert];
        /* Rotate space so that n is pointing along the z axis and just take 
           the first two components.  Formulas from
             http://en.wikipedia.org/wiki/Coordinate_rotation#Three_dimensions
         */
        cos_theta = plcl_dot_prod(n,e3);
        sin_theta = sqrt(1-cos_theta*cos_theta);
        u_para = plcl_scale_vect(plcl_dot_prod(*vt,axle),axle);
        u_perp = plcl_vect_diff(*vt,u_para);
        v_cross_u = plcl_cross_prod(*vt,axle);
        *vt = plcl_vect_sum(u_para,
            plcl_vlincomb(cos_theta,u_perp,sin_theta,v_cross_u));

        (*vt).c[2] = 0;

        printf("%g %g %g -> %g %g %g\n",plcl_M_clist(L->cp[cmp].vt[vert]),
            plcl_M_clist(F->cp[cmp].vt[vert]));
      }
    }
  }

  return F;
}

/* Add a component which contains the 2d convex hull of the other components
 * (looking only at the first 2 dimensions, of course).  We'll be using the
 * gift-wrapping algorithm. */
static void add_convex_hull(plCurve *L) {
  plcl_vector *vt;
  int nv,verts;
  int cmp,vert;
  plcl_vector first,next;
  plcl_vector x_vect,y_vect,cross;
  plCurve_color *clr;

  for (nv = 0, cmp = 0; cmp < L->nc; cmp++) {
    nv += L->cp[cmp].nv;
  }
  vt = malloc(nv*sizeof(plcl_vector));
  assert(vt != NULL);
  clr = malloc(nv*sizeof(plcl_vector));
  assert(clr != NULL);
  verts = 0;

  /* Look for the lowest point */
  first = L->cp[0].vt[0];
  for (cmp = 0; cmp < L->nc; cmp++) {
    for (vert = 0; vert < L->cp[cmp].nv; vert++) {
      if (first.c[1] - L->cp[cmp].vt[vert].c[1] > DBL_EPSILON) {
        first = L->cp[cmp].vt[vert];
      }
    }
  }
  vt[verts++] = first;

  /* Area of the triangle made by three vertices.  Formula adapted from 
       http://www.cse.unsw.edu.au/~lambert/java/3d/triangle.html
     which says that the area of the triangle defined by 
       p1 = (x1,y1), p2 = (x2,y2), p3 = (x3,y3) is

          1  | x1 x2 x3 |
         --- | y1 y2 y3 |
          2  |  1  1  1 |

       which looks to me a lot like

         1
         - (x1,x2,x3) x (y1,y2,y3)
         2 
   */

  next = plcl_build_vect(0.0,0.0,0.0);
  while (!plcl_vecteq(vt[0],next) && verts < nv) {
    /* There is one other vertex out there which, when it is the other end of a
     * line segment with "first" all of the other vectors are to the left of
     * the line.  Find it in as pedantic a way as possible, ye olde brute
     * search (and a particularly stupid one at that). */
    next = (plcl_vecteq(L->cp[0].vt[0],first)) ?
      L->cp[0].vt[1] : L->cp[0].vt[0];
    x_vect.c[0] = first.c[0];
    y_vect.c[0] = first.c[1];
    x_vect.c[1] = next.c[0];
    y_vect.c[1] = next.c[1];
    for (cmp = 0; cmp < L->nc; cmp++) {
      for (vert = 0; vert < L->cp[cmp].nv; vert++) {
        if (!plcl_vecteq(L->cp[cmp].vt[vert],first) &&
            !plcl_vecteq(L->cp[cmp].vt[vert],next)) {
          x_vect.c[2] = L->cp[cmp].vt[vert].c[0];
          y_vect.c[2] = L->cp[cmp].vt[vert].c[1];
          cross = plcl_cross_prod(x_vect,y_vect);
          if (- cross.c[0] - cross.c[1] - cross.c[2] > DBL_EPSILON) {
            /* cmp:vert is to the right of (or same as) next, take it instead */
            next = L->cp[cmp].vt[vert];
            x_vect.c[1] = x_vect.c[2];
            y_vect.c[1] = y_vect.c[2];
          }
        }
      }
    }
    if (plcl_vecteq(next,vt[verts-1])) {
      fprintf(stderr,"Looping...\n");
    }
    vt[verts++] = next;
    first = next;
  } /* while */
  verts--;
  for (vert = 0; vert < verts; vert++) {
    clr[vert].r = 1.0 - ((double)vert/(double)verts);
    clr[vert].g = 0.0;
    clr[vert].b = (double)vert/(double)verts;
    clr[vert].alpha = 1.0;
  }
  plCurve_add_component(L,verts,false,verts,vt,clr);
  free(vt);
  free(clr);
}

/* Rotate (in 2-dimensions) the first n-1 components based on the convex hull
 * in the nth component, so that the shortest "diameter" lies in the y
 * direction.  */
typedef struct vect_and_dist_type {
  double dist;
  plcl_vector dir;
} vect_and_dist;

static void rotate_2pic(plCurve *L) {
  plCurve_strand *cp;
  int cnt,cmp,vert;
  plcl_vector side_direction;
  vect_and_dist vals[1024];

  for (cnt=0; cnt < 1024; cnt++) {
    vals[cnt].dist = DBL_MAX;
  }
  cp = &(L->cp[L->nc-1]);
  for (vert = 0; vert < cp->nv; vert++) {
    side_direction =
      plcl_normalize_vect(plcl_vect_diff(cp->vt[vert+1],cp->vt[vert]),NULL);
  }
}

int main() {
  FILE *vectfile;
  plCurve *L, *F;
  char err_str[80];
  int err_num;

  vectfile = fopen("kl_3_1_I.vect","r");
  assert(vectfile != NULL);
  L = plCurve_read(vectfile,&err_num,err_str,80);
  assert(L != NULL);
  (void)fclose(vectfile);
  if (err_num != 0) {
    fprintf(stderr,"Error reading file %s:\n%s\n","kl_3_1_I.vect",err_str);
    exit(EXIT_FAILURE);
  }
  F = flatten(L,plcl_random_vect());
  assert(F != NULL);
  add_convex_hull(F);
  rotate_2pic(F);
  vectfile = fopen("kl_3_1_flat.vect","w");
  assert(vectfile != NULL);
  plCurve_write(vectfile,F);
  (void)fclose(vectfile);
  exit(EXIT_SUCCESS);
}
