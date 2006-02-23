/*
 *  Routines to create, destroy, read and write plCurves (and strands)
 * 
 *  $Id: plCurve.c,v 1.58 2006-02-23 12:33:06 ashted Exp $
 *
 */

/* Copyright 2004 The University of Georgia. */

/* This file is part of plCurve.
   
plCurve is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

plCurve is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with plCurve; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#include <config.h>
#include <plCurve.h>

#ifdef HAVE_STDLIB_H
  #include <stdlib.h>
#endif
#ifdef HAVE_MATH_H
  #include <math.h>
#endif
#ifdef HAVE_MALLOC_H
  #include <malloc.h>
#endif
#ifdef HAVE_STDARG_H
  #include <stdarg.h>
#endif
#ifdef HAVE_CTYPE_H
  #include <ctype.h>
#endif
#ifdef HAVE_STRING_H
  #include <string.h>
#endif
#ifdef HAVE_ASSERT_H
  #include <assert.h>
#endif
#ifdef HAVE_FLOAT_H
  #include <float.h>
#endif


/* Utility routines */
static inline int intmin(const int a, const int b) {
  return (a <= b) ? a : b;
}
static inline int intmax(const int a, const int b) {
  return (a >= b) ? a : b;
}

/*
 * Procedure allocates memory for a new plCurve. The number of components is
 * given by components. The number of vertices in each component shows up in
 * the buffer pointed to by nv.  The closed/open nature of each strand is given
 * in the array pointed to by open, and the number of colors per strand is 
 * found in the array cc.
 *
 * In each strand we allocate two extra vertices, at -1 and nv to make
 * "wrap-around" much simpler.
 *
 */
/*@only@*/ plCurve *plCurve_new(const int components, const int * const nv, 
                     const bool * const open, const int * const cc) {
  plCurve *L;
  int i;

  /* First, we check to see that the input values are reasonable. */

  assert(components > 0);
  assert(nv != NULL);
  assert(open != NULL);
  assert(cc != NULL);

  /* Now we attempt to allocate space for these components. */
  
  L = (plCurve *)malloc(sizeof(plCurve));
  assert(L != NULL);
  L->nc = components;
  L->cp = (plCurve_strand *)malloc(L->nc*sizeof(plCurve_strand));
  assert(L->cp != NULL);
  /*@+loopexec@*/
  for (i = 0; i < components; i++) {
    assert(nv[i] >= 1); /* Need to have at least one vertex */

    L->cp[i].open = open[i];
    L->cp[i].nv = nv[i];
    L->cp[i].vt = 
      (plcl_vector *)calloc((size_t)(nv[i]+2),sizeof(plcl_vector));
    assert(L->cp[i].vt != NULL);
    /*@-usedef@*/ /* Splint gets this one wrong */
    L->cp[i].vt++; /* so that L->cp[i].vt[-1] is a valid space */
    /*@=usedef@*/
    
    L->cp[i].cc = cc[i];
    L->cp[i].clr = 
      (plCurve_color *)calloc((size_t)cc[i],sizeof(plCurve_color));
    assert(L->cp[i].clr != NULL);
  }
  /*@=loopexec@*/

  /* Start off with no constraints */
  L->cst = NULL;

  /* Space for vertex quantifiers to be stored (someday) */
  L->quant = NULL;

  /*@-compdef@*/ /* Frustratingly, splint thinks that L->cp[i].nv, .open, .cc,
                    .vt and .clr don't get defined. */
  return L;
  /*@=compdef@*/
} /* plCurve_new */

/*
 * Find the closest point on the given line to the given point.
 *
 * coef should hold 6 doubles, which we will call a,b,c,d,e,f
 * the components of the vector point we will call x,y,z
 *
 * The line is given parametrically by (at + b, ct + d, et + f) and the
 * closest point then has
 *  
 *       a(x - b) + c(y - d) + e(z - f)
 *   t = ------------------------------
 *              a^2 + c^2 + e^2
 *
 */
static inline plcl_vector Closest_line_point(const plcl_vector point, 
                                             const plcl_vector tangent,
                                             const plcl_vector point_on_line) {
  plcl_vector ret_vect;
  double a = tangent.c[0];
  double b = point_on_line.c[0];
  double c = tangent.c[1];
  double d = point_on_line.c[1];
  double e = tangent.c[2];
  double f = point_on_line.c[2];
  double x = point.c[0];
  double y = point.c[1];
  double z = point.c[2];
  double t = (a*(x-b) + c*(y-d) + e*(z-f))/(a*a + c*c + e*e);

  /* Make sure that we aren't about to divide by zero */
  assert(a*a + c*c + e*e > DBL_EPSILON);

  ret_vect.c[0] = a*t + b;
  ret_vect.c[1] = c*t + d;
  ret_vect.c[2] = e*t + f;

  return ret_vect;
}

/*
 * Find the closest point on the given plane to the given point.
 *
 * coef should hold 4 doubles, which we will call a,b,c,d
 * the components of the vector point we will call x1,y1,z1
 *
 * The plane can be given by ax + by + cz = d or parameterized as
 *
 *          d - ax - by
 *   (x, y, -----------)
 *               c
 *
 * The closest point then has
 *
 *              a x1 + b y1 + c z1 - d
 *   x = x1 - a ----------------------
 *                 a^2 + b^2 + c^2
 *
 *              a x1 + b y1 + c z1 - d
 *   y = y1 - b ----------------------
 *                 a^2 + b^2 + c^2
 *
 */
static inline plcl_vector Closest_plane_point(const plcl_vector point, 
                                              const plcl_vector normal,
                                              const double d) {

  plcl_vector ret_vect;
  double a = normal.c[0];
  double b = normal.c[1];
  double c = normal.c[2];
  double x1 = point.c[0];
  double y1 = point.c[1];
  double z1 = point.c[2];
  int i0 = 0;
  int i1 = 1;
  int i2 = 2;
  double frac = (a*x1 + b*y1 + c*z1 - d)/(a*a + b*b + c*c);
  double x = x1 - a*frac;
  double y = y1 - b*frac;

  assert(a*a + b*b + c*c > DBL_EPSILON);

  /* 
   * We try to make sure that we aren't dividing by some tiny number later on.
   */
  if (c <= 1e-12 && c >= -1e-12) {
    if (b <= 1e-12 && b >= -1e-12) {
      i0 = 2; 
      i2 = 0;
      c = a;
      a = normal.c[2];
    } else {
      i1 = 2;
      i2 = 1;
      c = b;
      b = normal.c[2];
    }
  }

  /* And that we aren't dividing by zero now */

  ret_vect.c[i0] = x;
  ret_vect.c[i1] = y;
  ret_vect.c[i2] = (d - a*x - b*y)/c;

  return ret_vect;
}

/*
 * Check to see if the given constraint is satisfied and return value for 
 * how far off it is. 
 *
 */
double plCurve_check_cst(const plCurve * const L) {

  plcl_vector closest;
  plCurve_constraint *cst;
  double max_err = 0.0;
  int vert; 
  double sq_dist;  /* The squared distance */

  /* Sanity check */
  assert(L != NULL);

  cst = L->cst;
  while (cst != NULL) {
    assert(cst->kind == PLCL_FIXED ||
           cst->kind == PLCL_ON_LINE ||
           cst->kind == PLCL_IN_PLANE);
    if (cst->kind == PLCL_FIXED) {
      closest = cst->vect[0];
      /* PLCL_FIXED constraints only ever apply to one vertex */
      sq_dist = plcl_M_sq_dist(L->cp[cst->cmp].vt[cst->vert],closest);
      max_err = fmax(max_err,sq_dist);
    } else if (cst->kind == PLCL_ON_LINE) {
      for (vert = cst->vert; vert < cst->vert + cst->num_verts; vert++) {
        closest = Closest_line_point(L->cp[cst->cmp].vt[vert],
                                     cst->vect[0], cst->vect[1]);
        sq_dist = plcl_M_sq_dist(L->cp[cst->cmp].vt[vert],closest);
        max_err = fmax(max_err,sq_dist);
      } 
    } else if (cst->kind == PLCL_IN_PLANE) {
      for (vert = cst->vert; vert < cst->vert + cst->num_verts; vert++) {
        closest = Closest_plane_point(L->cp[cst->cmp].vt[vert],
                                      cst->vect[0], cst->vect[1].c[0]);
        sq_dist = plcl_M_sq_dist(L->cp[cst->cmp].vt[vert],closest);
        max_err = fmax(max_err,sq_dist);
      } 
    }
    cst = cst->next;
  }

  return sqrt(max_err);
} /* plCurve_check_cst */

/*
 * Fix any vertices which are out of compliance with their constraints.
 *
 */
void plCurve_fix_cst(plCurve * const L) {
  plcl_vector closest;
  plCurve_constraint *cst;
  int vert; 

  /* Sanity check */
  assert(L != NULL);

  cst = L->cst;
  while (cst != NULL) {
    assert(cst->kind == PLCL_FIXED ||
           cst->kind == PLCL_ON_LINE ||
           cst->kind == PLCL_IN_PLANE);
    if (cst->kind == PLCL_FIXED) {
      closest = cst->vect[0];
      /* PLCL_FIXED constraint is only ever allowed to be length 1 */
      L->cp[cst->cmp].vt[cst->vert] = closest;
    } else if (cst->kind == PLCL_ON_LINE) {
      for (vert = cst->vert; vert < cst->vert + cst->num_verts; vert++) {
        closest = Closest_line_point(L->cp[cst->cmp].vt[vert],
                                     cst->vect[0], cst->vect[1]);
        L->cp[cst->cmp].vt[vert] = closest;
      } 
    } else if (cst->kind == PLCL_IN_PLANE) {
      for (vert = cst->vert; vert < cst->vert + cst->num_verts; vert++) {
        closest = Closest_plane_point(L->cp[cst->cmp].vt[vert],
                                      cst->vect[0], cst->vect[1].c[0]);
        L->cp[cst->cmp].vt[vert] = closest;
      } 
    }
    cst = cst->next;
  }

  return;
} /* plCurve_fix_cst */

/*
 * Free the memory associated with a given plCurve.  We then set all the values
 * in the plCurve data structure to reflect the fact that the memory has been
 * freed.  We can call plCurve_free twice on the same plCurve pointer without
 * fear. 
 *
 */ 
void plCurve_free(/*@only@*/ /*@null@*/ plCurve *L) {
  int i;

  /* First, we check the input. */
  if (L == NULL) {
    return; /* Move along, nothing to see here */
  }

  /* Now we can get to work. */
  if (L->cp != NULL) {
    /*@+loopexec@*/
    for (i=0; i<L->nc; i++) {
      L->cp[i].nv = 0;
      if (L->cp[i].vt != NULL) {
        L->cp[i].vt--; /* undo our original vt++ (for wraparound) */
        free(L->cp[i].vt);
        L->cp[i].vt = NULL;
      }
 
      L->cp[i].cc = 0;
      if (L->cp[i].clr != NULL) {
        free(L->cp[i].clr);
      }
    }
    /*@=loopexec@*/
  
    /*@-compdestroy@*/ /* Splint thinks we aren't freeing .vt and .clr */
    free(L->cp);
    /*@=compdestroy@*/
    L->nc = 0;
  }

  if (L->cst != NULL) {
    free(L->cst);
    L->cst = NULL;
  }

  free(L);
  L = NULL;
} /* plCurve_free */

/*@only@*/ static inline 
plCurve_constraint *new_constraint(const int kind, const plcl_vector vect[], 
                                   const int cmp, const int vert, 
                                   const int num_verts, 
             /*@only@*/ /*@null@*/ plCurve_constraint *next) {
  plCurve_constraint *ncst;

  ncst = malloc(sizeof(plCurve_constraint));
  assert(ncst != NULL);
  ncst->kind = kind;
  ncst->vect[0] = vect[0];
  ncst->vect[1] = vect[1];
  ncst->cmp = cmp;
  ncst->vert = vert;
  ncst->num_verts = num_verts;
  ncst->next = next;

  return ncst;
}

/*
 * Remove any overlap between this (presumably new) constraint and the ones
 * which follow it. 
 *
 */
static inline void overrun_check(plCurve_constraint *cst) {

  plCurve_constraint *temp_cst;

  assert(cst != NULL);
  while (cst->next != NULL && 
         cst->next->cmp == cst->cmp &&
         cst->next->vert < cst->vert+cst->num_verts) {
    /* We overlap the next one, shrink or empty it */
    if (cst->next->vert + cst->next->num_verts <= 
        cst->vert       + cst->num_verts) {
      /* It has been completely subsumed and must be eliminated */
      temp_cst = cst->next;
      cst->next = cst->next->next; /* take it out of the list */
      free(temp_cst);
    } else {
      /* It just overlaps, either move the bottom up or join the two */
      assert(cst->next != NULL);
      if (cst->next->kind == cst->kind &&
          cst->kind != PLCL_FIXED && /* Never extend a fixed constraint */
          memcmp(cst->next->vect,cst->vect,sizeof(cst->vect)) == 0) {
        /* Same constraint so join the two */
        temp_cst = cst->next;
        cst->num_verts = 
          cst->next->vert+cst->next->num_verts - cst->vert;
        cst->next = cst->next->next;
        free(temp_cst);
      } else {
        /* Different constraint, just shorten the one ahead */
        cst->next->num_verts -= 
          cst->vert+cst->num_verts - cst->next->vert;
        cst->next->vert = cst->vert+cst->num_verts;
      }
    }
  }
} /* overrun_check */

/* Set a constraint on a vertex or run of vertices */
static inline void plCurve_set_constraint(plCurve * const L, const int cmp, 
                                          const int vert, const int num_verts, 
                                          const int kind, 
                                          const plcl_vector vect[]) {

  plCurve_constraint *cst,*pcst;  /* constraint, previous constraint */
  plCurve_constraint **pfn; /* Place for new constraint */

  /* First, we check the input. */
  assert(L != NULL);
  assert(cmp >= 0);
  assert(cmp < L->nc);
  assert(vert >= 0);
  assert(vert < L->cp[cmp].nv);
  assert(num_verts >= 1);
  assert(vert+num_verts <= L->cp[cmp].nv);
  assert(kind != PLCL_FIXED || num_verts == 1);

  /* Seek down the list 
   * Stop seeking when we either 
   *   1) run off the end of the list, or
   *   2) find one whose component is at least as large as ours and if it is
   *      the same as ours, has a range whose end extends to or beyond the
   *      beginning of ours. 
   */
  for (pcst = cst = L->cst; 
       cst != NULL && (cst->cmp < cmp || 
                       (cst->cmp == cmp && cst->vert+cst->num_verts <= vert));
       cst = cst->next) {
    pcst = cst;
  }
  if (pcst == cst) { /* Need to work at the head of the list */
    pfn = &L->cst;
  } else {
    assert(pcst != NULL);
    pfn = &pcst->next;
  }

  /* Now cst either points to a constraint which applies to the vertex in
   * question (or some vertex after it) or it points to NULL because there was
   * no such constraint found.  On the other hand, pcst either points to the
   * constraint just prior to the constraint which cst points to or else it is
   * NULL.  */
  if (pcst != NULL &&
      pcst->kind == kind &&
      kind != PLCL_FIXED && /* Never extend a fixed constraint */
      memcmp(pcst->vect,vect,sizeof(pcst->vect)) == 0 &&
      pcst->cmp == cmp &&
      pcst->vert+pcst->num_verts == vert) {
    /* We just passed one which has a range contiguous with the range we
       want to set and has the same attributes.  We'll just extend it. */
    pcst->num_verts += num_verts;
    overrun_check(pcst);
  } else if (cst != NULL &&
             cst->kind == kind &&
             kind != PLCL_FIXED && /* Never extend a fixed constraint */
             memcmp(pcst->vect,vect,sizeof(pcst->vect)) == 0 &&
             cst->cmp == cmp &&
             cst->vert <= vert+num_verts) {
    /* The one found can be extended to include our range. */
    cst->num_verts = intmax(vert+num_verts,cst->vert+cst->num_verts);
    cst->vert = intmin(vert,cst->vert);
    cst->num_verts -= cst->vert;
  } else if (cst == NULL || cst->cmp >cmp || cst->vert >= vert+num_verts) {
    /* Got to the end of the list without finding it or found one which deals
     * with vertices strictly beyond our range. */
    if (kind != PLCL_UNCST) {
      (*pfn) = new_constraint(kind, vect, cmp, vert, num_verts, cst);
    }
  } else {
    /* The one we found has a range that somehow overlaps to the one we have.
     * However, neither the one prior to it nor the one found are the same
     * constraint as we are trying to establish (otherwise we would have
     * extended them above).
     */
    if (cst->vert < vert) {
      /* Its range starts before ours, we'll have to leave part */
      if (cst->vert+cst->num_verts > vert+num_verts) {
        /* It also exteds past ours, we go in the middle */
        (*pfn) = new_constraint(cst->kind, cst->vect, cst->cmp, cst->vert,
                   vert - cst->vert, cst);
        assert((*pfn)->next != NULL);
        (*pfn)->next->num_verts -= vert + num_verts - (*pfn)->next->vert;
        (*pfn)->next->vert = vert + num_verts;
        if (kind != PLCL_UNCST) {
          (*pfn)->next = new_constraint(kind, vect, cmp, vert, num_verts, 
                                        (*pfn)->next);
        }
      } else {
        /* It's just before ours, but extends into ours, shorten it. */
        cst->num_verts = vert - cst->vert;
        /* And attach ours to the end */
        if (kind != PLCL_UNCST) {
          cst->next = 
            new_constraint(kind, vect, cmp, vert, num_verts, cst->next);
        }
        /* And check to see if we overran the next one */
        /*@-usereleased@*/
        cst = cst->next;
        /*@=usereleased@*/
        overrun_check(cst);
        /* Perhaps splint does have this right.  This code needs to be
         * re-thought. */
      }
    } else {
      /* It starts no earlier than ours, put ours in before it */
      (*pfn) = new_constraint(kind, vect, cmp, vert, num_verts, cst);
      overrun_check(*pfn);
    }
  }

  return;
} /* plCurve_set_constraint */

/* Now the four functions which call plCurve_set_constraint */
void plCurve_set_fixed(plCurve * const L, 
                       const int cmp, 
                       const int vert, 
                       const int num_verts, 
                       const plcl_vector point) {
  plcl_vector zv = {{0, 0, 0}};
  plcl_vector vect[2] = { point, zv };

  plCurve_set_constraint(L,cmp,vert,num_verts,PLCL_FIXED,vect);
}

void plCurve_constrain_to_line(plCurve * const L, 
                               const int cmp, 
                               const int vert, 
                               const int num_verts, 
                               const plcl_vector tangent, 
                               const plcl_vector point_on_line) {
  plcl_vector vect[2] = { tangent, point_on_line };

  plCurve_set_constraint(L,cmp,vert,num_verts,PLCL_ON_LINE,vect);
}

void plCurve_constrain_to_plane(plCurve * const L, 
                                const int cmp, 
                                const int vert, 
                                const int num_verts, 
                                const plcl_vector normal, 
                                const double dist_from_origin) {
  plcl_vector dz = {{dist_from_origin, 0, 0}};
  plcl_vector vect[2] = { normal, dz };

  plCurve_set_constraint(L,cmp,vert,num_verts,PLCL_IN_PLANE,vect);
}

void plCurve_unconstrain(plCurve * const L, const int cmp, 
                         const int vert, const int num_verts) {
        
  plcl_vector zv = {{0, 0, 0}};
  plcl_vector zero_vect[2] = {zv, zv};

  plCurve_set_constraint(L,cmp,vert,num_verts,PLCL_UNCST,zero_vect);
}

/*
 * Remove any constraint from the list which matches the data given.  Return
 * the number of vertices thus set unconstrained. 
 *
 */
int plCurve_remove_constraint(plCurve * const L, 
                              const int kind,
                              const plcl_vector vect[]) {
  plCurve_constraint **cst_pptr,*cst_ptr;
  int uncst = 0;

  assert(L != NULL);
  assert(kind == PLCL_FIXED || kind == PLCL_ON_LINE || kind == PLCL_IN_PLANE);

  /* This code gives splint fits (perhaps I just don't understand yet how to
   * explain to splint what is going on).  That's not unreasonable.  It's
   * tricky code.  Contact Ted if you think it has a bug. */
  cst_pptr = &L->cst;
  while ((*cst_pptr) != NULL) {
    if ((*cst_pptr)->kind == kind &&
        plcl_M_vecteq((*cst_pptr)->vect[0],vect[0]) &&
        plcl_M_vecteq((*cst_pptr)->vect[1],vect[1])) {
      cst_ptr = *cst_pptr; /* So we can free the space */
      *cst_pptr = (*cst_pptr)->next; /* Take it out of the list */
      uncst += cst_ptr->num_verts;
      /*@-compdestroy@*/
      free(cst_ptr);
      /*@=compdestroy@*/
      cst_ptr = NULL;
    } else {
      cst_pptr = &((*cst_pptr)->next);
    }
  }
  /*@-usereleased -compdef@*/
  return uncst;
  /*@=usereleased =compdef@*/
} /* plCurve_remove_constraint */

/* Remove all constraints */
void plCurve_remove_all_constraints(plCurve * const L) {
  plCurve_constraint *cst;

  /* Start with a sanity check */
  assert(L != NULL);

  cst = L->cst;
  while (cst != NULL) {
    L->cst = cst->next;
    free(cst);
    cst = L->cst;
  }
} /* plCurve_remove_all_constraints */

/*
 * Writes the plCurve to a file in Geomview VECT format.  The file format is:
 *
 * VECT                           # mandatory keyword
 * Ncomponents Nvertices Ncolors  # total number of components and vertices
 * Nv[0] ... Nv[NPolylines-1]     # number of vertices in each polyline 
 *                                # closed polylines use negative numbers
 * Nc[0] ... Nc[NPolylines-1]     # number of colors for each polyline 
 * # Constraint information
 * Vert[0] ... Vert[Nvertices-1]  # All the vertices, as triples of doubles
 *                                # with constraints listed in comments 
 * Color[0] ... Color[NColors]    # All the colors, in RGBA format
 *
 * Comments begin with #, and proceed to the end of the line. They are allowed
 * wherever a newline is allowed.
 *
 * We assume that file is open for writing.
 *
 */

void plCurve_write(FILE *file, plCurve * const L) {
  int i,j;                  /* Counters for the for loops */ 
  int nverts = 0;           /* Total number of vertices of all components */
  int colors = 0;           /* Total number of colors of all components */
  char outstr[80] = "";     /* So we can wrap the vertex lines */
  int num_cst = 0;          /* Tally of constraints */
  plCurve_constraint *cst;  /* Current constraint */
  plCurve_constraint *cst2;
  /* All the constraints (without duplicates) */
  plCurve_constraint *cst_list = NULL; 
  bool found;

  /* First, do a little sanity checking. */
  assert(file != NULL);
  assert(L != NULL);
  assert(L->nc >= 1);
  assert(L->cp != NULL);

  /* Now we begin work. */
  for(i=0;i<L->nc;i++) {
    nverts += L->cp[i].nv;
    colors += L->cp[i].cc;
  }

  /* We are ready to write the plCurve. */
  fprintf(file,"VECT \n");
  fprintf(file,"%d %d %d # Components Vertices Colors\n",L->nc,nverts,colors);
  
  for(i=0;i<L->nc;i++) {
    if (L->cp[i].open) {
      fprintf(file,"%d ",L->cp[i].nv); 
    } else {
      fprintf(file,"%d ",-L->cp[i].nv);
    }
  }
  fprintf(file,"# Vertices per Component\n");

  for(i=0;i<L->nc;i++) {
    fprintf(file,"%d ",L->cp[i].cc);
  }
  fprintf(file,"# Colors per Compoment\n");
  
  /* Slide the constraints, if any, in here */
  found = false;
  for (cst = L->cst; cst != NULL; cst = cst->next) {
    for (cst2 = cst_list; cst2 != NULL && !found; cst2 = cst2->next) { 
      if (cst->kind == cst2->kind &&
          plcl_M_vecteq(cst->vect[0],cst2->vect[0]) &&
          plcl_M_vecteq(cst->vect[1],cst2->vect[1])) {
        found = true;
      }
    }
  }
  fprintf(file,"# Constraints: %d\n",num_cst);
  cst = L->cst;
  num_cst = 0;
  while (cst != NULL) {
    num_cst++;
    assert(cst->kind == PLCL_FIXED ||
           cst->kind == PLCL_ON_LINE ||
           cst->kind == PLCL_IN_PLANE);
    if (cst->kind == PLCL_FIXED) {
      fprintf(file, "#  %d Fixed %lg %lg %lg\n",num_cst,
        plcl_M_clist(cst->vect[0]));
    } else if (cst->kind == PLCL_ON_LINE) {
      fprintf(file, "#  %d Line %lg %lg %lg %lg %lg %lg\n", num_cst,
        plcl_M_clist(cst->vect[0]), plcl_M_clist(cst->vect[1]));
    } else if (cst->kind == PLCL_IN_PLANE) {
      fprintf(file, "#  %d Plane %lg %lg %lg %lg\n", num_cst,
        plcl_M_clist(cst->vect[0]), cst->vect[1].c[0]);
    }
    cst = cst->next;
  }

  fprintf(file,"# Vertex coordinates\n");

  /* Now we write the vertex data . . . */
  cst = L->cst;
  num_cst = 1;
  for (i=0; i<L->nc; i++) {
    fprintf(file,"# Component %d\n",i);
    for (j=0; j<L->cp[i].nv; j++) {
      (void)snprintf(outstr,sizeof(outstr),"%.16g %.16g %.16g",
          plcl_M_clist(L->cp[i].vt[j]));
      while (cst != NULL &&
             (cst->cmp < i ||
              (cst->cmp == i && cst->vert+cst->num_verts <= j))) {
        cst = cst->next;
        num_cst++;
      }
      /* Now either we have nothing, or something yet to come, or something
       * currently applicable. */
      if (cst != NULL && cst->cmp == i && cst->vert <= j) {
        /* Something applicable */
        (void)snprintf(outstr,sizeof(outstr),"%s # Cst: %d",outstr,num_cst);
      }
      fprintf(file,"%s\n",outstr);
    }
  }

  fprintf(file,"# Colors (red green blue alpha)\n");

  /* . . . and the color data. */
  for (i=0; i < L->nc; i++) {
    fprintf(file,"# Component %d\n",i);
    for (j=0; j < L->cp[i].cc; j++) {
      fprintf(file,"%g %g %g %g\n", L->cp[i].clr[j].r, L->cp[i].clr[j].g,
        L->cp[i].clr[j].b, L->cp[i].clr[j].alpha);
    }
  }
      
  /* And we're done. */
  return;
}


/*
 * The next section of the library file includes some (private) procedures for
 * reading plCurve data reliably from Geomview VECT files.
 *
 */

/* 
 * get_comment acts just like skip_whitespace_and_comments but returns any
 * (possibly multi-line) comment in the given buffer.  
 *
 * If buflen > 0, get_comment makes sure that it contains a \0-terminated
 * string.
 * 
 */
#define end_str_and_return \
  if (i < (int)buflen) { buf[i] = '\0'; } \
  return;
#define add_char_to_buf(c) \
  if (i < (int)(buflen-1)) { buf[i++] = (char)c; } 

static inline void get_comment(FILE *infile, char *buf, const size_t buflen) {
  int i = 0;
  bool commentflag = false;
  bool skip_spaces = false; /* skip spaces after a hash, \n or other space */
  bool just_wrapped = false;
  int inchar;

  /* First, we check to make sure that infile looks legit. */
  assert(infile != NULL);
  
  while (true) {
    inchar = fgetc(infile);

    if (inchar == EOF) {
      end_str_and_return;
    } else if (inchar == (int)'#') { /* Started a comment. */
      if (commentflag) { /* We're already reading a comment */
        if (just_wrapped) { /* Looks like the comment is continued */
          skip_spaces = true;
          just_wrapped = false;
        } else {
          add_char_to_buf(inchar);
          skip_spaces = false;
        }
      } else {
        commentflag = skip_spaces = true;
      }
    } else if (inchar == (int)'\n') {
      if (commentflag) { 
        add_char_to_buf(' ') /* \n is replaced by space */
        skip_spaces = true; 
        just_wrapped = true;
      }
    } else if (isspace(inchar)) { /* blank */
      if (commentflag && !skip_spaces) {
        add_char_to_buf(' '); /* Any amount of any type of whitespace becomes
                                 one space. */
        skip_spaces = true;
      }
    } else {
      skip_spaces = false;
      if (commentflag && !just_wrapped) { /* Need another '#' after newline */
        add_char_to_buf(inchar);
      } else {
        (void)ungetc(inchar,infile);
        end_str_and_return;
      } 
    }
  }
}

/*
 * skip_whitespace_and_comments positions the file pointer on next
 * non-whitespace character, returning false if EOF happens first. We skip
 * anything between a # and a newline.
 * 
 */
static inline int skip_whitespace_and_comments(FILE *infile)
{
  int thischar;
  bool commentflag = false;

  /* First, we check to make sure that infile looks legit. */
  assert(infile != NULL);
  
  /* Now we start to work. */
  for(;;) {
    thischar = fgetc(infile);

    if (thischar == EOF) {  
      /* Reached end of file before a non-space, non-comment */
      return 0;
    } else if (thischar == (int)'#') { /* Started a comment. */
      commentflag = true;
    } else if (thischar == (int)'\n' && commentflag) { /* End a comment. */
      commentflag = false;
    } else if (!isspace(thischar) && !commentflag) { /* Found a hit! */
      (void)ungetc(thischar,infile);
      return 1;
    } /* It must have been a space or a non-space in a comment. */
  }
}

/* Procedure scans for nfloats floating point (or double) numbers, ignoring  *
 * whitespace and comments between them. We expect the variable length       *
 * arguments to contain a collection of pointers to doubles. If not, there's *
 * trouble.                                                                  */
static inline int scandoubles(FILE *infile,int ndoubles, ...)
{
  int nconverted = 0,i;
  va_list ap;
  double *thisdouble;

  /* First, we check for overall sanity. */
  assert(infile != NULL);
  assert(ndoubles >= 1);

  va_start(ap,ndoubles);

  /* Now we're ready to work. */

  for (i=0;i<ndoubles;i++) {    /* We expect to exit from the loop by */
                                /* returning, but this is a safety.   */
    if (skip_whitespace_and_comments(infile) == 0) { /* Failed */
      va_end(ap);              
      return nconverted;
    }

    thisdouble = va_arg(ap,double *);
    if (fscanf(infile,"%lf",thisdouble) != 1) { /* We couldn't scan. */
      va_end(ap);              
      return nconverted;        /* So give up here */
    } else {                    /* Else record our victory */
      nconverted++;
    }
  }
  va_end(ap);

  return nconverted;
}

/* Procedure scans for nints integers, ignoring whitespace and     *
 * comments between them. We expect the variable length arguments  *
 * to contain a collection of pointers to ints. If not,            *
 * there's trouble.                                                */
static inline int scanints(FILE *infile,int nints, ...)
{
  int nconverted = 0,i;
  va_list ap;
  int *thisint;

  /* First, we check for overall sanity. */
  assert(infile != NULL);
  assert(nints >= 1);

  va_start(ap,nints);

  /* Now we're ready to work. */
  for (i=0;i<nints;i++) {       /* We expect to exit from the loop by */
                                /* returning, but this is a safety.   */
    if (skip_whitespace_and_comments(infile) == 0) { /* Failed */
      va_end(ap);              
      return nconverted;
    }
    thisint = va_arg(ap,int *);

    if (fscanf(infile,"%d",thisint) != 1) {     /* We couldn't scan. */
      va_end(ap);              
      return nconverted;        /* So give up here */
    } else {                    /* Else record our victory */
      nconverted++;
    }
  }
  va_end(ap);
  return nconverted;
}

/* 
 * Touchup the "extra" vertices at each end of the component strands which are
 * used to implement "wraparound".
 *
 */
void plCurve_fix_wrap(plCurve * const L) {
  int i,nv;

  for (i = 0; i < L->nc; i++) {
    nv = L->cp[i].nv;
    if (L->cp[i].open) {
      /* fold it back on itself: v_{-1} = v_1 and v_{nv} = v_{nv-2} */
      L->cp[i].vt[-1] = L->cp[i].vt[1];
      L->cp[i].vt[nv] = L->cp[i].vt[nv-2];
    } else {
      /* wrap it around: v_{-1} = v_{nv-1} and v_{nv} = v_0 */
      L->cp[i].vt[-1] = L->cp[i].vt[nv-1];
      L->cp[i].vt[nv] = L->cp[i].vt[0];
    }
  }
}

/*
 * Read a Geomview VECT file and create a plCurve.  Color information is not
 * preserved.  File is assumed to be open for reading. Returns either a pointer
 * to a newly allocated plCurve structure (don't forget to FREE it!) or NULL on
 * failure. 
 *
 */
/*@only@*/ /*@null@*/ plCurve *plCurve_read(FILE *file, 
                                            int *error_num, 
                                            char error_str[], 
                                            size_t error_str_len) {
  plCurve *L;
  int nverts = 0, ncomp = 0, ncolors = 0;
  int *nvarray, *ccarray;
  bool *open; 
  int i, j;
  int nv;
  char comment[256] = ""; /* Space to store per-vertex comments */
  int ds; /* Doubles scanned */
  
  *error_num = 0;

  /* First, we check for the 'VECT' keyword. */
  if (fscanf(file," VECT ") == EOF) {
    *error_num = PLCL_E_NO_VECT;
    (void)snprintf(error_str,error_str_len,
        "plCurve_read: Couldn't find VECT keyword.\n");
    return NULL;
  }

  /* Now we read the three integers giving vertices, components, and colors. */

  if (scanints(file,3,&ncomp,&nverts,&ncolors) != 3) {
    *error_num = PLCL_E_BAD_CVC_LINE;
    (void)snprintf(error_str,error_str_len,
      "plCurve_read: Couldn't parse <ncomp> <nverts> <ncolors> line.\n");
    return NULL;
  }

  if (ncomp <= 0) {
    *error_num = PLCL_E_BAD_CMP_NUM;
    (void)snprintf(error_str,error_str_len,
        "plCurve_read: VECT file defines plCurve with %d compoments.\n",ncomp);
    return NULL;
  } 

  /* We now try to read the array of numbers of vertices. */

  nvarray = (int *)malloc(ncomp*sizeof(int));
  assert(nvarray != NULL);
  nvarray[0] = 0;
  open    = (bool *)malloc(ncomp*sizeof(bool));
  assert(open != NULL);
  ccarray = (int *)malloc(ncomp*sizeof(int));
  assert(ccarray != NULL);
  ccarray[0] = 0;

  /*@+loopexec@*/
  for(i=0; i<ncomp; i++) {
    if (scanints(file,1,&(nvarray[i])) != 1) {
      *error_num = PLCL_E_BAD_CVRT_LINE;
      (void)snprintf(error_str,error_str_len,
          "plCurve_read: Couldn't parse number of vertices in component %d.\n",
          i);    
      free(nvarray);
      free(open);
      free(ccarray);
      return NULL;
    }
    /* A negative number of vertices indicates a CLOSED component. */
    open[i] = (nvarray[i] >= 0);
    nvarray[i] = abs(nvarray[i]);
  }
  /*@=loopexec@*/

  /* We have set nvarray and open and are ready to read the color data.  */

  /*@+loopexec@*/
  for(i=0; i<ncomp; i++) {
    if (scanints(file,1,&(ccarray[i])) != 1) {
      *error_num = PLCL_E_BAD_CLR_LINE;
      (void)snprintf(error_str,error_str_len,
          "plCurve_read: Couldn't parse <ncolors> for component %d.\n", i);   
      free(nvarray);
      free(open);
      free(ccarray);
      return NULL;
    }
  }
  /*@=loopexec@*/

  /* We now allocate the plCurve data structure. */

  L = plCurve_new(ncomp,nvarray,open,ccarray);

  /* done with temorary arrays */
  free(nvarray);
  free(open);
  free(ccarray);

  assert(L != NULL);

  /* And get ready to read the actual data. */

  for(i = 0; i < ncomp; i++) {
    nv = L->cp[i].nv;
    for(j = 0; j < nv; j++) {
      if (scandoubles(file,3,plcl_M_clist(&L->cp[i].vt[j])) != 3) {
        plCurve_free(L);
        *error_num = PLCL_E_BAD_VERT_LINE;
        (void)snprintf(error_str,error_str_len,
          "plCurve_read: Couldn't parse <x> <y> <z> data for vertex %d of "
          "component %d.\n",j,i);
        return NULL;
      }
      get_comment(file,comment,sizeof(comment));
    }
  }
  /* Now set the "wrap-around" vertices */
  plCurve_fix_wrap(L);

  /*
   * And next the colors. Unfortunately, to really comply with 
   *   the Geomview standard here we have to be kind of careful. 
   */
  for (i=0; i < ncomp; i++) {
    for (j=0; j < L->cp[i].cc; j++) {
      if ((ds = scandoubles(file,4, &L->cp[i].clr[j].r, &L->cp[i].clr[j].g,
           &L->cp[i].clr[j].b, &L->cp[i].clr[j].alpha)) != 4) {
        plCurve_free(L);
        *error_num = PLCL_E_BAD_COLOR;
        (void)snprintf(error_str,error_str_len,
          "plCurve_read: Couldn't parse color %d in component %d of "
          "plCurve (%d).\n",j,i,ds);
        return NULL;
      }
    }
  }

  return L;
}

#define strand_edges(P) (((P).open) ? (P).nv-1 : (P).nv)
/* 
 *   Return the total number of edges in plCurve. 
 */
int plCurve_num_edges(plCurve * const L) 
{
  int i, edges = 0;

  for (i=0; i<L->nc; i++) {
    edges += strand_edges(L->cp[i]);
  }
  return edges;
}

/* Compute the curvature of L at vertex vt of component cp */

double plCurve_curvature(plCurve * const L, const int comp, const int vert) {
  double      kappa;
  plcl_vector in,out;
  double      normin, normout;
  double      dot_prod,cross_prod_norm;

  /* We start with some initializations. */
  
  plCurve_fix_wrap(L);
  
  /* Now we work. */

  in = plcl_vect_diff(L->cp[comp].vt[vert],L->cp[comp].vt[vert-1]);
  normin = plcl_M_norm(in);
  out = plcl_vect_diff(L->cp[comp].vt[vert+1],L->cp[comp].vt[vert]);
  normout = plcl_M_norm(out);
      
  dot_prod = plcl_M_dot(in,out);
  cross_prod_norm = plcl_M_norm(plcl_cross_prod(in,out));

  /* Check for infinite curvature condition */
  assert(normin*normout + dot_prod > DBL_EPSILON);  /* PERHAPS RETURN ERROR? */
      
  kappa = (2*cross_prod_norm)/(normin*normout + dot_prod);
  kappa /= (normin - normout < DBL_EPSILON) ? normin : normout;

  return kappa;
}

/* 
 * Duplicate a plCurve and return the duplicate.
 *
 */
plCurve *plCurve_copy(plCurve * const L) {
  plCurve *nL;
  int *nv, *ccarray;
  bool *open;
  int cnt,cnt2;
  plCurve_constraint *cst,*cst2;

  assert(L->nc >= 1);

  nv = (int *)malloc((L->nc)*sizeof(int));
  assert(nv != NULL);
  open = (bool *)malloc((L->nc)*sizeof(bool));
  assert(open != NULL);
  ccarray = (int *)malloc((L->nc)*sizeof(int));
  assert(ccarray != NULL);

  /*@+loopexec@*/
  for (cnt = 0; cnt < L->nc; cnt++) {
    nv[cnt] = L->cp[cnt].nv;
    open[cnt] = L->cp[cnt].open;
    ccarray[cnt] = L->cp[cnt].cc;
  }
  /*@=loopexec@*/
  nL = plCurve_new(L->nc,nv,open,ccarray);
  assert(nL->cst == NULL);

  free(ccarray);
  free(open);
  free(nv);

  /* Constraints */
  cst = L->cst;
  if (cst != NULL) {
    nL->cst = new_constraint(cst->kind, cst->vect, cst->cmp, cst->vert, 
                             cst->num_verts, NULL);
    cst2 = nL->cst;
    while (cst->next != NULL) {
      cst = cst->next;
      assert(cst2->next == NULL);
      cst2->next = new_constraint(cst->kind, cst->vect, cst->cmp, cst->vert,
                                  cst->num_verts, NULL);
      cst2 = cst2->next;
    }
  } else {
    nL->cst = NULL;
  }

  for (cnt = 0; cnt < L->nc; cnt++) {
    /*
     * Copy the vertices (including the "hidden" ones, so we don't have to call
     * plCurve_fix_wrap).
     */
    for (cnt2 = -1; cnt2 <= L->cp[cnt].nv; cnt2++) { 
      nL->cp[cnt].vt[cnt2] = L->cp[cnt].vt[cnt2];
    }
    for (cnt2 = 0; cnt2 < L->cp[cnt].cc; cnt2++) {
      nL->cp[cnt].clr[cnt2] = L->cp[cnt].clr[cnt2];
    }
  }

  return nL;
}

/*
 * Compute a (unit) tangent vector to L at vertex vert of component cmp. 
 *
 */
plcl_vector plCurve_tangent_vector(plCurve * const L, int cmp, int vert,
                                   bool *ok) {
  plcl_vector in, out, tan;

  if (L->cp[cmp].open) {
    if (vert == 0) {
      tan = plcl_vect_diff(L->cp[cmp].vt[1], L->cp[cmp].vt[0]);

      return plcl_normalize_vect(tan,ok);
    } else if (vert == L->cp[cmp].nv-1) {
      tan = plcl_vect_diff(L->cp[cmp].vt[L->cp[cmp].nv-1],
                           L->cp[cmp].vt[L->cp[cmp].nv-2]);

      return plcl_normalize_vect(tan,ok);
    }
  }

  /* We now know that either we are on a closed 
     component, or we are not at an endpoint.   */
  
   in = plcl_normalize_vect(
     plcl_vect_diff(L->cp[cmp].vt[vert+1],L->cp[cmp].vt[vert]),ok
   );

   out = plcl_normalize_vect(
     plcl_vect_diff(L->cp[cmp].vt[vert],L->cp[cmp].vt[vert-1]),ok
   );

   plcl_M_vweighted(tan,0.5,in,out);
   return plcl_normalize_vect(tan,ok);
} /* plCurve_tangent_vector */

/* Procedure computes the length of each component of the plCurve,
   and fills in the array of doubles "component_lengths", which 
   must be as long as L->nc. It returns the total length. We assume
   that fix_wrap has been called. */
double plCurve_arclength(const plCurve * const L,
                         double *component_lengths)
{
  double tot_length;
  int cmp; 
  int nv = 0;
  int vert;
  plCurve_strand *cp;

  tot_length = 0;
  for (cmp = 0; cmp < L->nc; cmp++) {

    component_lengths[cmp] = 0;
    cp = &L->cp[cmp];
    nv = (cp->open) ? cp->nv-1 : cp->nv;

    for (vert = 0; vert < nv; vert++) {
      component_lengths[cmp] += plcl_M_distance(cp->vt[vert+1],cp->vt[vert]);
    }

    tot_length += component_lengths[cmp];
  }

  return tot_length;
} /* plCurve_arclength */

/* Procedure reports the arclength distance from the given vertex */
/* to the 0th vertex of the given component of L. */
double plCurve_parameter(const plCurve * const L,const int cmp,const int vert)
{
  double tot_length;
  int v,nv;
  plCurve_strand *cp;
  plcl_vector temp_vect;

  assert(L != NULL);
  assert(cmp >= 0);
  assert(cmp < L->nc);
  assert(vert >= 0);
  assert(vert < L->cp[cmp].nv);

  tot_length = 0;
  cp = &L->cp[cmp];
  nv = (cp->open) ? cp->nv-1 : cp->nv;

  for (v = 0; v < vert; v++) {
    temp_vect = cp->vt[v+1];
    plcl_M_sub_vect(temp_vect,cp->vt[v]);
    tot_length += plcl_M_norm(temp_vect);
  }

  return tot_length;
}

/*
 * This procedure closes all open components of plCurve by distributing a small
 * change of all vertices of each such component. It also changes the "open"
 * flag and calls fix_wrap. We remove one vertex in this process. 
 */
void plCurve_force_closed(plCurve * const L)
{
  int i, cmp;
  plcl_vector diff;

  for (cmp=0;cmp < L->nc;cmp++) {
    if (L->cp[cmp].open == true) {  /* Isolate the open components. */

      /* Compute the error in closure */
      diff = L->cp[cmp].vt[L->cp[cmp].nv-1];   
      plcl_M_sub_vect(diff,L->cp[cmp].vt[0]);

      for (i=0;i<L->cp[cmp].nv;i++) {
        plcl_M_vlincomb(L->cp[cmp].vt[i],
          1.0,L->cp[cmp].vt[i],
         -1.0*i/(L->cp[cmp].nv-1),diff);
      }

      /* We claim to have moved the last vertex on top of the first. */
      diff = plcl_vect_diff(L->cp[cmp].vt[0],
                            L->cp[cmp].vt[L->cp[cmp].nv-1]);
      assert(plcl_M_dot(diff,diff) < DBL_EPSILON);

      /* Thus we eliminate the last vertex. */
      L->cp[cmp].nv--;
      L->cp[cmp].open = false;
    }
  }

  plCurve_fix_wrap(L);
}

/* 
 * Either return or display the library version number.
 *
 */
void plCurve_version(char *version) {
  if (version == NULL) {
    printf("plCurve Version: %s\n",PACKAGE_VERSION);
  } else {
    (void)snprintf(version,sizeof(version),PACKAGE_VERSION);
  }
}
