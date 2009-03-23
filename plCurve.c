/*
 *  Routines to create, destroy, read and write plCurves (and strands)
 *
 *  $Id: plCurve.c,v 1.90 2007-09-21 21:08:43 cantarel Exp $
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
/*@only@*/ plCurve *plc_new(const int          components,
                             const int  * const nv,
                             const bool * const open,
                             const int  * const cc) {
  plCurve *L;
  int i;

  /* First, we check to see that the input values are reasonable. */

  assert(components > 0);
  assert(nv != NULL);
  assert(open != NULL);
  assert(cc != NULL);
  
  for(i=0;i<components;i++) { assert(cc[i] >= 0); assert(cc[i] <= nv[i]); }

  /* Now we attempt to allocate space for these components. */

  L = (plCurve *)malloc(sizeof(plCurve));
  assert(L != NULL);
  L->nc = components;
  L->cp = (plc_strand *)malloc(L->nc*sizeof(plc_strand));
  assert(L->cp != NULL);
  /*@+loopexec@*/
  for (i = 0; i < components; i++) {
    assert(nv[i] >= 1); /* Need to have at least one vertex */

    L->cp[i].open = open[i];
    L->cp[i].nv = nv[i];
    L->cp[i].vt =
      (plc_vector *)calloc((size_t)(nv[i]+2),sizeof(plc_vector));
    assert(L->cp[i].vt != NULL);
    /*@-usedef@*/ /* Splint gets this one wrong */
    L->cp[i].vt++; /* so that L->cp[i].vt[-1] is a valid space */
    /*@=usedef@*/

    L->cp[i].cc = cc[i];
    L->cp[i].clr =
      (plc_color *)calloc((size_t)cc[i],sizeof(plc_color));
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
} /* plc_new */

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
static inline plc_vector Closest_line_point(const plc_vector point,
                                             const plc_vector tangent,
                                             const plc_vector point_on_line) {
  plc_vector ret_vect;
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
static inline plc_vector Closest_plane_point(const plc_vector point,
                                              const plc_vector normal,
                                              const double d) {

  plc_vector ret_vect;
  double a = normal.c[0];
  double b = normal.c[1];
  double c = normal.c[2];
  int i0 = 0;
  int i1 = 1;
  int i2 = 2;
  double frac;
  double x;
  double y;

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

  frac = (a*point.c[i0] + b*point.c[i1] + c*point.c[i2] - d)/(a*a + b*b + c*c);
  x = point.c[i0] - a*frac;
  y = point.c[i1] - b*frac;
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
double plc_check_cst(const plCurve * const L) {

  plc_vector closest;
  plc_constraint *cst;
  double max_err = 0.0;
  int vert;
  double sq_dist;  /* The squared distance */

  /* Sanity check */
  assert(L != NULL);

  cst = L->cst;
  while (cst != NULL) {
    assert(cst->kind == fixed || cst->kind == line || cst->kind == plane);
    if (cst->kind == fixed) {
      closest = cst->vect[0];
      /* PLC_FIXED constraints only ever apply to one vertex */
      sq_dist = plc_M_sq_dist(L->cp[cst->cmp].vt[cst->vert],closest);
      max_err = fmax(max_err,sq_dist);
    } else if (cst->kind == line) {
      for (vert = cst->vert; vert < cst->vert + cst->num_verts; vert++) {
        closest = Closest_line_point(L->cp[cst->cmp].vt[vert],
                                     cst->vect[0], cst->vect[1]);
        sq_dist = plc_M_sq_dist(L->cp[cst->cmp].vt[vert],closest);
        max_err = fmax(max_err,sq_dist);
      }
    } else if (cst->kind == plane) {
      for (vert = cst->vert; vert < cst->vert + cst->num_verts; vert++) {
        closest = Closest_plane_point(L->cp[cst->cmp].vt[vert],
                                      cst->vect[0], cst->vect[1].c[0]);
        sq_dist = plc_M_sq_dist(L->cp[cst->cmp].vt[vert],closest);
        max_err = fmax(max_err,sq_dist);
      }
    }
    cst = cst->next;
  }

  return sqrt(max_err);
} /* plc_check_cst */

/*
 * Fix any vertices which are out of compliance with their constraints.
 *
 */
void plc_fix_cst(plCurve * const L) {
  plc_vector closest;
  plc_constraint *cst;
  int vert;

  /* Sanity check */
  assert(L != NULL);

  cst = L->cst;
  while (cst != NULL) {
    assert(cst->kind == fixed || cst->kind == line || cst->kind == plane);
    if (cst->kind == fixed) {
      closest = cst->vect[0];
      /* PLC_FIXED constraint is only ever allowed to be length 1 */
      L->cp[cst->cmp].vt[cst->vert] = closest;
    } else if (cst->kind == line) {
      for (vert = cst->vert; vert < cst->vert + cst->num_verts; vert++) {
        closest = Closest_line_point(L->cp[cst->cmp].vt[vert],
                                     cst->vect[0], cst->vect[1]);
        L->cp[cst->cmp].vt[vert] = closest;
      }
    } else if (cst->kind == plane) {
      for (vert = cst->vert; vert < cst->vert + cst->num_verts; vert++) {
        closest = Closest_plane_point(L->cp[cst->cmp].vt[vert],
                                      cst->vect[0], cst->vect[1].c[0]);
        L->cp[cst->cmp].vt[vert] = closest;
      }
    }
    cst = cst->next;
  }

  plc_fix_wrap(L);

  return;
} /* plc_fix_cst */

/*
 * Free the memory associated with a given plCurve.  We then set all the values
 * in the plCurve data structure to reflect the fact that the memory has been
 * freed.  We can call plc_free twice on the same plCurve pointer without
 * fear.
 *
 */
void plc_free(/*@only@*/ /*@null@*/ plCurve *L) {
  int i;
  plc_constraint *cst;
  plc_vert_quant *qnt;

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
        L->cp[i].clr = NULL;
      }
    }
    /*@=loopexec@*/

    /*@-compdestroy@*/ /* Splint isn't sure we are freeing .vt and .clr */
    free(L->cp);
    /*@=compdestroy@*/
    L->cp = NULL;
    L->nc = 0;
  }

  while (L->cst != NULL) {
    cst = L->cst;
    L->cst = cst->next;
    free(cst);
    cst = NULL;
  }

  while (L->quant != NULL) {
    qnt = L->quant;
    L->quant = qnt->next;
    free(qnt);
    qnt = NULL;
  }

  free(L);
} /* plc_free */

/*@only@*/ static inline
plc_constraint *new_constraint(const plc_cst_kind kind,
                                   const plc_vector vect[],
                                   const int cmp,
                                   const int vert,
                                   const int num_verts,
             /*@only@*/ /*@null@*/ plc_constraint *next) {
  plc_constraint *ncst;

  ncst = malloc(sizeof(plc_constraint));
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
static inline void overrun_check(plc_constraint *cst) {

  plc_constraint *temp_cst;

  while (cst->next != NULL &&
         cst->next->cmp == cst->cmp &&
         cst->next->vert < cst->vert+cst->num_verts) {
    /* We overlap (or almost overlap) the next one, shrink or empty it */
    if (cst->next->vert + cst->next->num_verts <=
        cst->vert       + cst->num_verts) {
      /* It has been completely subsumed and must be eliminated */
      temp_cst = cst->next;
      cst->next = cst->next->next; /* take it out of the list */
      free(temp_cst);
      temp_cst = NULL;
    } else {
      /* It just overlaps, either move the bottom up or join the two */
      assert(cst->next != NULL);
      if (cst->next->kind == cst->kind &&
          cst->kind != fixed && /* Never extend a fixed constraint */
          memcmp(cst->next->vect,cst->vect,sizeof(cst->vect)) == 0) {
        /* Same constraint so join the two */
        cst->num_verts = cst->next->vert+cst->next->num_verts - cst->vert;
        temp_cst = cst->next;
        cst->next = cst->next->next;
        free(temp_cst);
        temp_cst = NULL;
      } else {
        /* Different constraint, just shorten the one ahead */
        cst->next->num_verts -=
          cst->vert+cst->num_verts - cst->next->vert;
        cst->next->vert = cst->vert+cst->num_verts;
      }
    }
  }
  /* Now combine adjacent constraints if possible */
  while (cst->next != NULL &&
         cst->next->vert == cst->vert+cst->num_verts &&
         cst->next->kind == cst->kind &&
         cst->kind != fixed &&
         memcmp(cst->next->vect,cst->vect,sizeof(cst->vect)) == 0) {
    /* Same constraint so join the two */
    cst->num_verts = cst->next->vert+cst->next->num_verts - cst->vert;
    temp_cst = cst->next->next;
    free(cst->next);
    cst->next = temp_cst;
  }
} /* overrun_check */

/* Set a constraint on a vertex or run of vertices */
static inline void plc_set_constraint(plCurve * const L,
                                          const int          cmp,
                                          const int          vert,
                                          const int          num_verts,
                                          const plc_cst_kind kind,
                                          const plc_vector vect[]) {

  plc_constraint *cst,*pcst; /* constraint, previous constraint */
  plc_constraint **pfn; /* Place for new constraint */
  plc_constraint *temp_cst;

  /* First, we check the input. */
  assert(L != NULL);
  assert(cmp >= 0);
  assert(cmp < L->nc);
  assert(vert >= 0);
  assert(vert < L->cp[cmp].nv);
  assert(num_verts >= 1);
  assert(vert+num_verts <= L->cp[cmp].nv);
  assert(kind != fixed || num_verts == 1);

  /* Seek down the list
   * Stop seeking when we either
   *   1) run off the end of the list, or
   *   2) find one whose component is at least as large as ours and if it is
   *      the same as ours, has a range whose end extends to or beyond the
   *      beginning of ours.
   */
  pcst = NULL;
  cst = L->cst;
  for (; cst != NULL && (cst->cmp < cmp ||
                        (cst->cmp == cmp && cst->vert+cst->num_verts <= vert));
       cst = cst->next) {
    pcst = cst;
  }
  if (pcst == NULL) { /* Need to work at the head of the list */
    pfn = &L->cst;
  } else {
    pfn = &pcst->next;
  }

  /* Now cst either points to a constraint which applies to the vertex in
   * question (or some vertex after it) or it points to NULL because there was
   * no such constraint found.  On the other hand, pcst either points to the
   * constraint just prior to the constraint which cst points to or else it is
   * NULL.  */
#define disp_cond(C) \
  if (C) { printf(#C "\n"); }
  if (pcst != NULL &&
      kind != fixed && /* Never extend a fixed constraint */
      pcst->kind == kind &&
      pcst->cmp == cmp &&
      pcst->vert+pcst->num_verts >= vert &&
      memcmp(pcst->vect,vect,sizeof(pcst->vect)) == 0) {
    /* We just passed one which has a range contiguous with the range we
       want to set and has the same attributes.  We'll just extend it. */
    pcst->num_verts = intmax(vert+num_verts,pcst->vert+pcst->num_verts);
    pcst->vert = intmin(vert,pcst->vert);
    pcst->num_verts -= pcst->vert;
    overrun_check(pcst);
  } else if (cst != NULL &&
             kind != fixed && /* Never extend a fixed constraint */
             cst->kind == kind &&
             cst->cmp == cmp &&
             cst->vert <= vert+num_verts &&
             memcmp(cst->vect,vect,sizeof(cst->vect)) == 0) {
    /* The one found can be extended to include our range. */
    cst->num_verts = intmax(vert+num_verts,cst->vert+cst->num_verts);
    cst->vert = intmin(vert,cst->vert);
    cst->num_verts -= cst->vert;
    overrun_check(cst);
  } else if (cst == NULL || cst->cmp > cmp || cst->vert >= vert+num_verts) {
    /* Got to the end of the list without finding it or found one which deals
     * with vertices strictly beyond our range. */
    if (kind != unconstrained) {
      (*pfn) = new_constraint(kind, vect, cmp, vert, num_verts, cst);
    }
  } else {
    /* The one we found has a range that somehow overlaps to the one we have.
     * However, neither the one prior to it nor the one found are the same
     * constraint as we are trying to establish (otherwise we would have
     * extended them above).
     */
    assert(cst != NULL);
    if (cst->vert < vert) {
      /* Its range starts before ours, we'll have to leave part */
      if (cst->vert+cst->num_verts > vert+num_verts) {
        /* It also extends past ours, we go in the middle */
        (*pfn) = new_constraint(cst->kind, cst->vect, cst->cmp, cst->vert,
                   vert - cst->vert, cst);
        cst = (*pfn)->next;
        assert(cst != NULL);
        cst->num_verts -= vert + num_verts - cst->vert;
        cst->vert = vert + num_verts;
        if (kind != unconstrained) {
          (*pfn)->next = new_constraint(kind, vect, cmp, vert, num_verts, cst);
        }
      } else {
        /* It's just before ours, but extends into ours, shorten it. */
        cst->num_verts = vert - cst->vert;
        /* And attach ours to the end */
        temp_cst = cst->next;
        cst->next = new_constraint(kind, vect, cmp, vert, num_verts, temp_cst);
        /* And check to see if we overran the next one */
        overrun_check(cst->next);
        if (kind == unconstrained) {
          /* Now remove the "unconstraint" */
          temp_cst = cst->next;
          cst->next = cst->next->next;
          free(temp_cst);
          temp_cst = NULL;
        }
      }
    } else {
      /* It starts no earlier than ours, put ours in before it */
      temp_cst = new_constraint(kind, vect, cmp, vert, num_verts, cst);
      (*pfn) = temp_cst;
      overrun_check(temp_cst);
      if (kind == unconstrained) {
        /* Now remove the "unconstraint" */
        temp_cst = (*pfn);  /* This line here for splint's benefit */
        (*pfn) = temp_cst->next;
        free(temp_cst);
      }
    }
  /*@-branchstate@*/
  }
  /*@=branchstate@*/

  /* Splint believes that L->cst was freed and is yet accessible :-( */
  /*@-usereleased -compdef@*/
  return;
  /*@=usereleased =compdef@*/
} /* plc_set_constraint */

/* Now the four functions which call plc_set_constraint */
void plc_set_fixed(plCurve * const L,
                       const int          cmp,
                       const int          vert,
                       const plc_vector point) {
  plc_vector zv = {{ 0.0, 0.0, 0.0 }};
  plc_vector vect[2] = { point, zv };

  plc_set_constraint(L,cmp,vert,1,fixed,vect);
}

void plc_constrain_to_line(plCurve * const L,
                               const int          cmp,
                               const int          vert,
                               const int          num_verts,
                               const plc_vector tangent,
                               const plc_vector point_on_line) {
  plc_vector vect[2] = { tangent, point_on_line };

  plc_set_constraint(L,cmp,vert,num_verts,line,vect);
}

void plc_constrain_to_plane(plCurve * const L,
                                const int          cmp,
                                const int          vert,
                                const int          num_verts,
                                const plc_vector normal,
                                const double dist_from_origin) {
  plc_vector dz = {{ dist_from_origin, 0.0, 0.0 }};
  plc_vector vect[2] = { normal, dz };

  plc_set_constraint(L,cmp,vert,num_verts,plane,vect);
}

void plc_unconstrain(plCurve * const L, const int cmp,
                         const int vert, const int num_verts) {

  plc_vector zv = {{ 0.0, 0.0, 0.0 }};
  plc_vector zero_vect[2] = {zv, zv};

  plc_set_constraint(L,cmp,vert,num_verts,unconstrained,zero_vect);
}

/*
 * Remove any constraint from the list which matches the data given.  Return
 * the number of vertices thus set unconstrained.
 *
 */
int plc_remove_constraint(plCurve * const L,
                              const plc_cst_kind kind,
                              const plc_vector vect[]) {
  plc_constraint **cst_pptr,*cst_ptr;
  int uncst = 0;

  assert(L != NULL);
  assert(kind == fixed || kind == line || kind == plane);

  /* This code gives splint fits (perhaps I just don't understand yet how to
   * explain to splint what is going on).  That's not unreasonable.  It's
   * tricky code.  Contact Ted if you think it has a bug. */
  cst_pptr = &L->cst;
  while ((*cst_pptr) != NULL) {
    if ((*cst_pptr)->kind == kind &&
        plc_M_vecteq((*cst_pptr)->vect[0],vect[0]) &&
        plc_M_vecteq((*cst_pptr)->vect[1],vect[1])) {
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
} /* plc_remove_constraint */

/* Remove all constraints */
void plc_remove_all_constraints(plCurve * const L) {
  plc_constraint *cst;

  /* Start with a sanity check */
  assert(L != NULL);

  cst = L->cst;
  while (cst != NULL) {
    L->cst = cst->next;
    free(cst);
    cst = L->cst;
  }
} /* plc_remove_all_constraints */

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

void plc_write(FILE *outfile, plCurve * const L) {
  int i,j;                  /* Counters for the for loops */
  int nverts = 0;           /* Total number of vertices of all components */
  int colors = 0;           /* Total number of colors of all components */
  char outstr[80] = "";     /* So we can wrap the vertex lines */
  char outstr2[80] = "";    /* Again for wrapping */
  int cst_num;              /* This constraint */
  int **cst_nums = NULL;    /* All constraint numbers */
  plc_constraint *cst;  /* Current constraint */
  plc_constraint *cst2;
  /* All the constraints (without duplicates) */
  plc_constraint *cst_list = NULL;
  bool found;

  /* First, do a little sanity checking. */
  assert(outfile != NULL);
  assert(L != NULL);
  assert(L->nc >= 1);
  assert(L->cp != NULL);

  cst_nums = (int **)calloc((size_t)L->nc,sizeof(int *));
  assert(cst_nums != NULL);

  /* Now we begin work. */
  for (i=0; i<L->nc; i++) {
    nverts += L->cp[i].nv;
    colors += L->cp[i].cc;
    /* Make space for the constraint numbers and set them all to zero */
    cst_nums[i] = (int *)calloc((size_t)L->cp[i].nv,sizeof(int));
    assert(cst_nums[i] != NULL);
  }

  /* We are ready to write the plCurve. */
  fprintf(outfile,"VECT \n");
  fprintf(outfile,
    "%d %d %d # Components Vertices Colors\n",L->nc,nverts,colors);

  for(i=0;i<L->nc;i++) {
    if (L->cp[i].open) {
      fprintf(outfile,"%d ",L->cp[i].nv);
    } else {
      fprintf(outfile,"%d ",-L->cp[i].nv);
    }
  }
  fprintf(outfile,"# Vertices per Component\n");

  for(i=0;i<L->nc;i++) {
    fprintf(outfile,"%d ",L->cp[i].cc);
  }
  fprintf(outfile,"# Colors per Compoment\n");

  /* Slide the constraints, if any, in here */
  for (cst = L->cst; cst != NULL; cst = cst->next) {
    assert(cst->kind == fixed || cst->kind == line || cst->kind == plane);
    found = false;
    cst_num = 0;
    for (cst2 = cst_list; cst2 != NULL && !found; cst2 = cst2->next) {
      cst_num++;
      if (cst->kind == cst2->kind &&
          plc_M_vecteq(cst->vect[0],cst2->vect[0]) &&
          plc_M_vecteq(cst->vect[1],cst2->vect[1])) {
        found = true;
      }
    }
    if (!found) {
      /* It's a new one, add it to the list */
      cst_list = new_constraint(cst->kind, cst->vect, cst->cmp, cst->vert,
                                cst->num_verts, cst_list);
      cst_num++;
      if (cst->kind == fixed) {
        fprintf(outfile, "# Constraint %d Fixed %g %g %g\n", cst_num,
          plc_M_clist(cst->vect[0]));
      } else if (cst->kind == line) {
        fprintf(outfile, "# Constraint %d Line %g %g %g %g %g %g\n", cst_num,
          plc_M_clist(cst->vect[0]), plc_M_clist(cst->vect[1]));
      } else if (cst->kind == plane) {
      fprintf(outfile, "# Constraint %d Plane %g %g %g %g\n", cst_num,
          plc_M_clist(cst->vect[0]), cst->vect[1].c[0]);
      }
    }
    /* Now mark it down in the list */
    for (i = cst->vert; i < cst->vert + cst->num_verts; i++) {
      cst_nums[cst->cmp][i] = cst_num;
    }
  }

  fprintf(outfile,"# Vertex coordinates\n");

  /* Now we write the vertex data . . . */
  for (i=0; i<L->nc; i++) {
    fprintf(outfile,"# Component %d\n",i);
    for (j=0; j<L->cp[i].nv; j++) {
      (void)snprintf(outstr,sizeof(outstr),"%.16g %.16g %.16g",
          plc_M_clist(L->cp[i].vt[j]));
      if (cst_nums[i][j] != 0) {
        strcpy(outstr2,outstr);
        (void)snprintf(outstr,sizeof(outstr),
                       "%s # Cst: %d",outstr2,cst_nums[i][j]);
      }
      /* Here is where we will eventually write out the quantifiers, wrapping
       * as necessary */
      /* Now complete the string. */
      fprintf(outfile,"%s\n",outstr);
    }
  }

  fprintf(outfile,"# Colors (red green blue alpha)\n");

  /* . . . and the color data. */
  for (i=0; i < L->nc; i++) {
    fprintf(outfile,"# Component %d\n",i);
    for (j=0; j < L->cp[i].cc; j++) {
      fprintf(outfile,"%.16g %.16g %.16g %.16g\n", L->cp[i].clr[j].r,
          L->cp[i].clr[j].g, L->cp[i].clr[j].b, L->cp[i].clr[j].alpha);
    }
  }

  /* Clean up temporary storage */
  for (i=0; i < L->nc; i++) {
    free(cst_nums[i]);
    cst_nums[i] = NULL;
  }
  free(cst_nums);
  cst_nums = NULL;


  /* Don't forget to free the list of constraints */

  cst = cst_list;
  while (cst != NULL) {
    cst_list = cst->next;
    free(cst);
    cst = cst_list;
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
 * get_comment acts like skip_whitespace_and_comments but returns any (possibly
 * multi-line) comment in a buffer. 
 *
 */
#define end_str_and_return \
  buf[i] = '\0'; \
  return buf;

#define add_char_to_buf(c) \
  if (i == buflen-1) { \
    buflen += 256; \
    buf = realloc(buf,buflen); \
    assert(buf != NULL); \
  } \
  buf[i++] = (char)c;

/*@only@*/ static inline char *get_comment(FILE *infile) {
  bool comment_flag = false;
  bool skip_spaces = false; /* skip spaces after a hash, \n or other space */
  bool just_wrapped = false;
  size_t i;
  int inchar;
  char *buf;
  size_t buflen;

  /* First, we check to make sure that infile looks legit. */
  assert(infile != NULL);

  buflen = (size_t)256;
  buf = (char *)malloc(buflen);
  assert(buf != NULL);
  buf[0] = '\0';

  i = 0;
  while (true) {
    inchar = fgetc(infile);

    if (inchar == EOF) {
      end_str_and_return;
    } else if (inchar == (int)'#') { /* Started a comment. */
      if (comment_flag) { /* We're already reading a comment */
        if (just_wrapped) { /* Looks like the comment is continued */
          skip_spaces = true;
          just_wrapped = false;
        } else {
          add_char_to_buf(inchar);
          skip_spaces = false;
        }
      } else {
        comment_flag = skip_spaces = true;
      }
    } else if (inchar == (int)'\n') {
      if (comment_flag) {
        add_char_to_buf(' ') /* \n is replaced by space */
        skip_spaces = true;
        just_wrapped = true;
      }
    } else if (isspace(inchar)) { /* blank */
      if (comment_flag && !skip_spaces) {
        add_char_to_buf(' '); /* Any amount of any type of whitespace becomes
                                 one space. */
        skip_spaces = true;
      }
    } else {
      skip_spaces = false;
      if (comment_flag && !just_wrapped) { /* Need another '#' after newline */
        add_char_to_buf(inchar);
      } else {
        (void)ungetc(inchar,infile);
        end_str_and_return;
      }
    }
  }
  end_str_and_return;
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
 * used to implement "wraparound".  Also, if the number of vertices is less
 * than 3, require the component to be considered open.
 *
 */
void plc_fix_wrap(plCurve * const L) {
  int          i,nv;

  for (i = 0; i < L->nc; i++) {
    nv = L->cp[i].nv;
    if (nv < 3) {
      L->cp[i].open = true;
    }
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
/*@only@*/ /*@null@*/ plCurve *plc_read(FILE *file,
                                  /*@out@*/ int *error_num,
                                  /*@out@*/ char error_str[],
                                            size_t error_str_len) {
  plCurve *L;
  int nverts = 0, ncomp = 0, ncolors = 0;
  int *nvarray, *ccarray;
  bool *open;
  int i, j;
  int nv;
  char *comment; 
  int ds; /* Doubles scanned */
  /*@temp@*/ char *place;
  char *space = " ";
  size_t parsed;
  plc_constraint *cst_list;
  size_t cst_list_len = (size_t)2;
  int cst_num, prev_cst_num, prev_cst_vert;

  assert(file != NULL);
  assert(error_str_len > 0);
  error_str[0] = '\0';
  *error_num = 0;

  /* First, we check for the 'VECT' keyword. */
  if (fscanf(file," VECT ") == EOF) {
    *error_num = PLC_E_NO_VECT;
    (void)snprintf(error_str,error_str_len,
        "plc_read: Couldn't find VECT keyword.\n");
    return NULL;
  }

  /* Now we read the three integers giving vertices, components, and colors. */

  if (scanints(file,3,&ncomp,&nverts,&ncolors) != 3) {
    *error_num = PLC_E_BAD_CVC_LINE;
    (void)snprintf(error_str,error_str_len,
      "plc_read: Couldn't parse <ncomp> <nverts> <ncolors> line.\n");
    return NULL;
  }

  if (ncomp <= 0) {
    *error_num = PLC_E_BAD_CMP_NUM;
    (void)snprintf(error_str,error_str_len,
        "plc_read: VECT file defines plCurve with %d compoments.\n",ncomp);
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
      *error_num = PLC_E_BAD_CVRT_LINE;
      (void)snprintf(error_str,error_str_len,
          "plc_read: Couldn't parse number of vertices in component %d.\n",
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
      *error_num = PLC_E_BAD_CLR_LINE;
      (void)snprintf(error_str,error_str_len,
          "plc_read: Couldn't parse <ncolors> for component %d.\n", i);
      free(nvarray);
      free(open);
      free(ccarray);
      return NULL;
    }
  }
  /*@=loopexec@*/

  /* We now allocate the plCurve data structure. */

  L = plc_new(ncomp,nvarray,open,ccarray);

  /* done with temorary arrays */
  free(nvarray); free(open); free(ccarray);
  nvarray = NULL; open = NULL; ccarray = NULL;

  assert(L != NULL);

  /* Now get the constraint definitions, if any */

  cst_list = (plc_constraint *)
    calloc(cst_list_len,sizeof(plc_constraint));
  assert(cst_list != NULL);
  comment = get_comment(file);
/* If we don't know about strcasestr, use strstr instead, losing case
 * insensitivity.  Because some systems don't have strcasestr, we need to 
 * write our checks with the case we wrote out to the file. */
  #define strcasestr strstr

  while ((place = strcasestr(comment,"Constraint ")) != NULL) {
    /* Move the word "Constraint to the beginning of the string */
    memmove(comment,place,strlen(place)+1);
    /* Now start tokenizing */
    place = strtok(comment, space); /* "Constraint" */
#define place_not_null \
    if (place == NULL) { \
      plc_free(L); \
      assert(cst_list->next == NULL); \
      free(cst_list); \
      free(comment); \
      *error_num = PLC_E_BAD_CST_LINE; \
      (void)snprintf(error_str,error_str_len, \
        "plc_read: Couldn't parse constraint line.\n"); \
      return NULL; \
    }
    place_not_null;
    parsed = strlen(place)+1;
    place = strtok(NULL, space); /* The constraint number */
    place_not_null;
    parsed += strlen(place)+1;
    if (sscanf(place,"%d",&cst_num) != 1 ||
        cst_num < 1) {
      plc_free(L);
      assert(cst_list->next == NULL);
      free(cst_list);
      free(comment);
      *error_num = PLC_E_BAD_CST_NUM;
      (void)snprintf(error_str,error_str_len,
        "plc_read: Couldn't parse constraint number from '%s'.\n",place);
      return NULL;
    }
    if (cst_num+1 > (int)cst_list_len) {
      cst_list_len = (size_t)(10*((cst_num+11) % 10)); /* Allocate by 10s */
      assert(cst_list->next == NULL);
      cst_list = realloc(cst_list,cst_list_len*sizeof(plc_constraint));
      assert(cst_list != NULL);
    }
    place = strtok(NULL, space); /* The constraint type */
    place_not_null;
    parsed += strlen(place)+1;
    if (strcasestr(place,"Fixed") != NULL) {
      cst_list[cst_num].kind = fixed;
      if (sscanf(&comment[parsed],"%lf %lf %lf",
            &cst_list[cst_num].vect[0].c[0],
            &cst_list[cst_num].vect[0].c[1],
            &cst_list[cst_num].vect[0].c[2]) != 3) {
        plc_free(L);
        assert(cst_list->next == NULL);
        free(cst_list);
        free(comment);
        *error_num = PLC_E_BAD_CST_NUMS;
        (void)snprintf(error_str,error_str_len,
          "plc_read: Bad numbers for constraint %d.\n",cst_num);
        return NULL;
      }
      cst_list[cst_num].vect[1] = plc_build_vect(0.0, 0.0, 0.0);
      for (i=0; i < 3; i++) {
        place = strtok(NULL, space); /* Step over the numbers */
        place_not_null;
        parsed += strlen(place)+1;
      }
      memmove(comment,&comment[parsed],strlen(&comment[parsed])+1);
    } else if (strcasestr(place,"Line") != NULL) {
      cst_list[cst_num].kind = line;
      if (sscanf(&comment[parsed],"%lf %lf %lf %lf %lf %lf",
            &cst_list[cst_num].vect[0].c[0],
            &cst_list[cst_num].vect[0].c[1],
            &cst_list[cst_num].vect[0].c[2],
            &cst_list[cst_num].vect[1].c[0],
            &cst_list[cst_num].vect[1].c[1],
            &cst_list[cst_num].vect[1].c[2]) != 6) {
        plc_free(L);
        assert(cst_list->next == NULL);
        free(cst_list);
        free(comment);
        *error_num = PLC_E_BAD_CST_NUMS;
        (void)snprintf(error_str,error_str_len,
            "plc_read: Bad numbers for constraint %d.\n",cst_num);
        return NULL;
      }
      for (i=0; i < 6; i++) {
        place = strtok(NULL, space); /* Step over the numbers */
        place_not_null;
        parsed += strlen(place)+1;
      }
      memmove(comment,&comment[parsed],strlen(&comment[parsed])+1);
    } else if (strcasestr(place,"Plane") != NULL) {
      cst_list[cst_num].kind = plane;
      if (sscanf(&comment[parsed],"%lf %lf %lf %lf",
            &cst_list[cst_num].vect[0].c[0],
            &cst_list[cst_num].vect[0].c[1],
            &cst_list[cst_num].vect[0].c[2],
            &cst_list[cst_num].vect[1].c[0]) != 4) {
        plc_free(L);
        assert(cst_list->next == NULL);
        free(cst_list);
        free(comment);
        *error_num = PLC_E_BAD_CST_NUMS;
        (void)snprintf(error_str,error_str_len,
          "plc_read: Bad numbers for constraint %d.\n",cst_num);
        return NULL;
      }
      cst_list[cst_num].vect[1].c[1] = 0.0;
      cst_list[cst_num].vect[1].c[2] = 0.0;
      for (i=0; i < 4; i++) {
        place = strtok(NULL, space); /* Step over the numbers */
        place_not_null;
        parsed += strlen(place)+1;
      }
      memmove(comment,&comment[parsed],strlen(&comment[parsed])+1);
    } else {
      plc_free(L);
      assert(cst_list->next == NULL);
      free(cst_list);
      free(comment);
      *error_num = PLC_E_BAD_CST_KIND;
      (void)snprintf(error_str,error_str_len,
        "plc_read: Unrecognized constraint kind '%s'.\n",place);
      return NULL;
    }
  }
  free(comment);

  /* And get ready to read the actual data. */

  for (i = 0; i < ncomp; i++) {
    prev_cst_num = prev_cst_vert = cst_num = 0;
    nv = L->cp[i].nv;
    for (j = 0; j < nv; j++) {
      if (scandoubles(file,3,plc_M_clist(&L->cp[i].vt[j])) != 3) {
        plc_free(L);
        assert(cst_list->next == NULL);
        free(cst_list);
        *error_num = PLC_E_BAD_VERT_LINE;
        (void)snprintf(error_str,error_str_len,
          "plc_read: Couldn't parse <x> <y> <z> data for vertex %d of "
          "component %d.\n",j,i);
        return NULL;
      }
      comment = get_comment(file);
      cst_num = 0;
      if ((place = strcasestr(comment,"Cst: ")) != NULL) {
        memmove(comment,place,strlen(place)+1);
        place = strtok(comment,space); /* Cst: */
        place_not_null;
        if (sscanf(&comment[5],"%d",&cst_num) != 1) {
          plc_free(L);
          assert(cst_list->next == NULL);
          free(cst_list);
          free(comment);
          *error_num = PLC_E_BAD_CST_NUM;
          (void)snprintf(error_str,error_str_len,
            "plc_read: Bad numbers for constraint %d.\n",cst_num);
          return NULL;
        }
      }
      free(comment);
      if (cst_num != prev_cst_num) {
        if (prev_cst_num != 0) {
          plc_set_constraint(L, i, prev_cst_vert, j - prev_cst_vert,
              cst_list[prev_cst_num].kind, cst_list[prev_cst_num].vect);
        }
        prev_cst_num = cst_num;
        prev_cst_vert = j;
      }
    }
    if (prev_cst_num != 0) {
      plc_set_constraint(L, i, prev_cst_vert, j - prev_cst_vert,
          cst_list[cst_num].kind, cst_list[cst_num].vect);
    }
  }
  /* Now set the "wrap-around" vertices */
  plc_fix_wrap(L);

  /*
   * And next the colors. Unfortunately, to really comply with
   *   the Geomview standard here we have to be kind of careful.
   */
  for (i=0; i < ncomp; i++) {
    for (j=0; j < L->cp[i].cc; j++) {
      if ((ds = scandoubles(file,4, &L->cp[i].clr[j].r, &L->cp[i].clr[j].g,
           &L->cp[i].clr[j].b, &L->cp[i].clr[j].alpha)) != 4) {
        plc_free(L);
        assert(cst_list->next == NULL);
        free(cst_list);
        *error_num = PLC_E_BAD_COLOR;
        (void)snprintf(error_str,error_str_len,
          "plc_read: Couldn't parse color %d in component %d of "
          "plCurve (%d).\n",j,i,ds);
        return NULL;
      }
    }
  }

  assert(cst_list->next == NULL);
  free(cst_list);
  return L;
}

#define strand_edges(P) (((P).open) ? (P).nv-1 : (P).nv)
/*
 *   Return the total number of edges in plCurve.
 */
int plc_num_edges(const plCurve * const L) /*@modifies nothing@*/
{
  int i, edges = 0;

  for (i=0; i<L->nc; i++) {
    edges += strand_edges(L->cp[i]);
  }
  return edges;
}

/* Return total number of edges in plCurve, store #edges for each 
   strand in component_edges. */

int plc_edges(const plCurve * const L,
	      /*@null@*/ /*@out@*/ int *component_edges)

{
 int i, edges = 0;

 for (i=0; i<L->nc; i++) {

   edges += strand_edges(L->cp[i]);
   
   if (component_edges != NULL) { component_edges[i] = strand_edges(L->cp[i]); }

 }
 return edges;
} 
     
/* 
 * Return the total number of vertices in plCurve.
 */
int plc_num_verts(const plCurve * const L) /*@modifies nothing@*/
{
  int cmp, verts = 0;
  
  for (cmp = 0; cmp < L->nc; cmp++) {
    verts += L->cp[cmp].nv;
  }
  return verts;
}

/* Compute an index between 0 and plc_num_verts(L) - 1 for a (cp,vt) pair in a plCurve, 
   using full wraparound addressing for closed components and repeating the last or first vertex 
   for open ones. */
int plc_vertex_num(const plCurve * const L, int cp, int vt)
{
  int wrapVt=0,i;

  if (L->cp[cp].open) {

    if (vt < 0) { wrapVt = 0; }
    else if (vt > L->cp[cp].nv-1) { wrapVt = L->cp[cp].nv-1; }
    else { wrapVt = vt; }

  } else {

    for(wrapVt = vt;wrapVt < 0;wrapVt += L->cp[cp].nv);  // Make sure we start with a non-negative index (for modular arithmetic)
    wrapVt = wrapVt % L->cp[cp].nv;
  
  }

  for(i=0;i<cp;i++) { wrapVt += L->cp[i].nv; }

  return wrapVt;
}

int plc_cp_num (const plCurve * const L, int wrapVt)

{
  int i;
  assert(wrapVt < plc_num_verts(L));

  for(i=0;i<L->nc;i++) {

    wrapVt -= L->cp[i].nv;

    if (wrapVt < 0) { return i; }

  }

  fprintf(stderr,"plCurve: Library error in plc_cp_num.\n");
  exit(1); // We should never reach this point.

}

int plc_vt_num (const plCurve * const L, int wrapVt)

{
  int i;
  assert(wrapVt < plc_num_verts(L));

  for(i=0;i<L->nc;i++) {

    if (wrapVt < L->cp[i].nv) { return wrapVt; }
    wrapVt -= L->cp[i].nv;

  }

  fprintf(stderr,"plCurve: Library error in plc_vt_num.\n");
  exit(1);

}
  
/* Compute the MinRad-based curvature of L at vertex vt of component cp */
double plc_MR_curvature(plCurve * const L, const int cmp, const int vert) {

  double      kappa;
  plc_vector in,out;
  double      normin, normout;
  double      dot_prod,cross_prod_norm;

  if (L->cp[cmp].open && 
    /* There is no curvature at the endpoints of an open strand */
      (vert == 0 || vert == L->cp[cmp].nv-1)) {
    return 0.0;
  }

  in = plc_vect_diff(L->cp[cmp].vt[vert],L->cp[cmp].vt[vert-1]);
  normin = plc_M_norm(in);
  out = plc_vect_diff(L->cp[cmp].vt[vert+1],L->cp[cmp].vt[vert]);
  normout = plc_M_norm(out);

  dot_prod = plc_M_dot(in,out);
  cross_prod_norm = plc_M_norm(plc_cross_prod(in,out));

  /* Check for infinite curvature condition */
  assert(normin*normout + dot_prod > DBL_EPSILON); 

  kappa = (2*cross_prod_norm)/(normin*normout + dot_prod);
  kappa /= (normin - normout < DBL_EPSILON) ? normin : normout;

  return kappa;
}

/*
 * Duplicate a plCurve and return the duplicate.
 *
 */
plCurve *plc_copy(const plCurve * const L) {
  plCurve *nL;
  int *nv, *ccarray;
  bool *open;
  int cnt,cnt2;
  plc_constraint *cst,*cst2;

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
  nL = plc_new(L->nc,nv,open,ccarray);
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
     * plc_fix_wrap).
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
 * Compute a (unit) tangent vector to L at vertex vert of component cmp by
 * taking the incoming tangent and outgoing tangent and averaging them *with
 * their lengths taken into account*.  That is, the average of (0.0,0.0,6.0)
 * and (8.0,0.0,0.0) is (4.0,0.0,3.0) which normalizes to (0.8,0.0,0.6)
 *
 */
plc_vector plc_mean_tangent(const plCurve * const L, const int cmp,
                                 const int vert, bool *ok) {

  if (L->cp[cmp].open) {
    /* For open strands, take the only tangent we have. */
    if (vert == 0) {
      return plc_normalize_vect(
          plc_vect_diff(L->cp[cmp].vt[vert+1],L->cp[cmp].vt[vert]),
          ok);
    } else if (vert == L->cp[cmp].nv-1) {
      return plc_normalize_vect(
          plc_vect_diff(L->cp[cmp].vt[vert],L->cp[cmp].vt[vert-1]),
          ok);
    }
  }

  /* We now know that either we are on a closed
     component, or we are not at an endpoint.   */

  return plc_normalize_vect(
      plc_vect_diff(L->cp[cmp].vt[vert+1],L->cp[cmp].vt[vert-1]),
      ok);
} /* plc_tangent_vector */

/*
 * Procedure computes the length of each component of the plCurve, and fills in
 * the array of doubles "component_lengths" (if it exists), which must be as
 * long as L->nc. It returns the total length. We assume that fix_wrap has been
 * called. 
 *
 */
double plc_arclength(const plCurve * const L,
    /*@null@*/ /*@out@*/ double *component_lengths)
{
  double tot_length,c_length;
  int cmp;
  int nv = 0;
  int vert;
  plc_strand *cp;

  assert(L != NULL);

  tot_length = 0.0;
  /*@+loopexec@*/
  for (cmp = 0; cmp < L->nc; cmp++) {

    c_length = 0.0;
    cp = &L->cp[cmp];
    nv = (cp->open) ? cp->nv-1 : cp->nv;

    for (vert = 0; vert < nv; vert++) {
      c_length += plc_M_distance(cp->vt[vert+1],cp->vt[vert]);
    }

    tot_length += c_length;
    if (component_lengths != NULL) {
      component_lengths[cmp] = c_length;
    }
  }
  /*@=loopexec@*/

  return tot_length;
} /* plc_arclength */

/*
 * Report the arclength distance from one vertex to another.  On closed
 * strands, give the shortest of the two options. 
 *
 */
double plc_subarc_length(const plCurve * const L, const int cmp, 
                             const int vert1, const int vert2)
{
  double length1, length2;
  int v, v1, v2, start, end;
  plc_strand *cp;

  assert(L != NULL);
  assert(cmp >= 0);
  assert(cmp < L->nc);
  assert(vert1 >= 0);
  assert(vert1 < L->cp[cmp].nv);
  assert(vert2 >= 0);
  assert(vert2 < L->cp[cmp].nv);

  cp = &L->cp[cmp];

  if (vert1 <= vert2) {
    v1 = vert1;
    v2 = vert2;
  } else {
    v2 = vert1;
    v1 = vert2;
  }

  if (cp->open) {
    start = v1;
    end = v2;
  } else {
    start = 0;
    end = cp->nv;
  }

  length1 = length2 = 0.0;
  for (v = start; v < end; v++) {
    if (v >= v1 && v < v2) {
      length1 += plc_M_norm(plc_vect_diff(cp->vt[v+1],cp->vt[v]));
    } else {
      length2 += plc_M_norm(plc_vect_diff(cp->vt[v+1],cp->vt[v]));
    }
  }

  return ((cp->open) ? length1 : ((length1 <= length2) ? length1 : length2));
}

/*
 * Either return or display the library version number.
 *
 */
void plc_version(/*@null@*/ char *version, size_t strlen) {
  if (version == NULL) {
    printf("plCurve Version: %s\n",PACKAGE_VERSION);
  } else {
    (void)snprintf(version,strlen,PACKAGE_VERSION);
  }
}

/* Put 4 doubles together into a color */
plc_color plc_build_color(const double r,
                                  const double g,
                                  const double b,
                                  const double alpha) {
  plc_color C = { r, g, b, alpha };
  return C;
}

/* Add a component with nv vertices read from the array vt and cc colors read
 * from the array clr.  The parameter add_as tells what the new component
 * number will be (any components above it will be moved up by 1. */
void plc_add_component(plCurve *L, const int add_as, const int nv, 
                      const bool open, const int cc,
                      const plc_vector * const vt,
           /*@null@*/ const plc_color * const clr) {
  plc_strand *cp;

  assert(L != NULL);
  assert(nv >= 1);
  assert(vt != NULL);
  assert(cc >= 0);
  assert(cc == 0 || clr != NULL);
  assert(add_as >= 0);
  assert(add_as <= L->nc);

  L->nc++;
  /*@-compdestroy@*/
  L->cp = realloc(L->cp,L->nc*sizeof(plc_strand));
  /*@=compdestroy@*/
  assert(L->cp != NULL);
  if (L->nc - add_as > 1) {
    memmove(&(L->cp[add_as+1]),&(L->cp[add_as]),
        (L->nc-add_as-1)*sizeof(plc_strand));
  }
  cp = &L->cp[add_as];
  cp->open = open;

  cp->nv = nv;
  cp->vt = (plc_vector *)calloc((size_t)(nv+2),sizeof(plc_vector));
  assert(cp->vt != NULL);
  /*@-usedef@*/ /* Splint gets this one wrong */
  cp->vt++; /* so that cp->vt[-1] is a valid space */
  /*@=usedef@*/

  memmove(cp->vt,vt,nv*sizeof(plc_vector));

  cp->cc = cc;
  if (clr != NULL) {
    cp->clr = (plc_color *)calloc((size_t)cc,sizeof(plc_color));
    assert(cp->clr != NULL);
    memmove(cp->clr,clr,cc*sizeof(plc_color));
  } else {
    cp->clr = NULL;
  }
/*@-nullstate@*/
  plc_fix_wrap(L);
}
/*@=nullstate@*/

void plc_drop_component(plCurve *L, const int cmp) {
  int cnt;

  assert(L != NULL);
  assert(cmp >= 0);
  assert(cmp < L->nc);
  assert(L->nc > 1);

  L->nc--;
  L->cp[cmp].vt--;
  free(L->cp[cmp].vt);
  free(L->cp[cmp].clr);
  for (cnt = cmp; cnt < L->nc; cnt++) {
    L->cp[cnt] = L->cp[cnt+1];
  }
  /* Leave the memory sitting as the cost of the extra few bytes is small 
     compared to the cost of a realloc.  It will get cleaned up at the end and
     may get reused before then. */
  /*@-compdef -usereleased@*/
}
/*@=compdef =usereleased@*/

void plc_resize_colorbuf(plCurve *L,const int cp, const int cc)
/* Change the size of the color buffer for a plCurve, preserving existing data if it exists */
{
  int i;
  
  if (cc == 0) { /* We are supposed to eliminate the color buffer, if it exists. */

    if (L->cp[cp].cc == 0 && L->cp[cp].clr == NULL) { /* Nothing to do */

      return;

    } else {
    
      L->cp[cp].cc = 0;
      free(L->cp[cp].clr); 
      L->cp[cp].clr = NULL;

      return;

    }

  } else {

    if (L->cp[cp].cc == 0) { /* There is no existing data, and we are adding some. */

      L->cp[cp].cc = cc;
      L->cp[cp].clr = calloc(sizeof(plc_color),cc);
      assert(L->cp[cp].clr != NULL);

      return;

    } else if (L->cp[cp].cc == 1) { /* There is a single color for the component already, copy it to all locations. */

      L->cp[cp].cc = cc;
      L->cp[cp].clr = realloc(L->cp[cp].clr,sizeof(plc_color)*cc);
      assert(L->cp[cp].clr != NULL);

      for(i=1;i<cc;i++) {

	L->cp[cp].clr[i] = L->cp[cp].clr[0];

      }

      return;

    } else { /* There is more than one vertex worth of information
		present. We keep what we get under realloc, as we have
		no meaningful way to downsample color information. */

      L->cp[cp].cc = cc;
      L->cp[cp].clr = realloc(L->cp[cp].clr,sizeof(plc_color)*cc);
      assert(L->cp[cp].clr != NULL);

      return;

    }

  }

  printf("plc_resize_colorbuf: Library bug. Please report to cantarel@math.uga.edu.\n");

}

    

void mpsverts(plCurve *L,int cmp,int vt,plc_vector mpsverts[2])

     /* Computes the pair of minrad preserving spline verts which 
	should replace the vertex vt on component cmp of L. We assume
	that the vertex before and after vt are legal. */

{
  /* Take care of the trivial cases first. */

  if (L->cp[cmp].open) {

    if (vt == 0) {

      mpsverts[0] = L->cp[cmp].vt[vt];
      mpsverts[1] = plc_vweighted(0.25,L->cp[cmp].vt[vt],L->cp[cmp].vt[vt+1]);

      return;

    } else if (vt == L->cp[cmp].nv-1) {

      mpsverts[0] = plc_vweighted(0.75,L->cp[cmp].vt[vt-1],L->cp[cmp].vt[vt]);
      mpsverts[1] = L->cp[cmp].vt[vt];

      return;

    }

  }

  /* Now we work on the real stuff. */

  plc_vector v0,v1,v2,in,out;
  double mr,s0,s1,normin,normout,rad,curv,mc,dot_prod,cross_prod_norm;
  double sin2a,cos2a,tana,dist;

  v0 = L->cp[cmp].vt[vt-1]; v1 = L->cp[cmp].vt[vt]; v2 = L->cp[cmp].vt[vt+1];
  in = plc_vect_diff(v1,v0); out = plc_vect_diff(v2,v1);
  normin = plc_norm(in); normout = plc_norm(out);

  if (fabs(normin) < 1e-12 || fabs(normout) < 1e-12) {

    /* If there are doubled vertices, we probably won't make the situation worse by adding more. */

    mpsverts[0] = v1;
    mpsverts[1] = v1;

    return;
    
  }

  /* We now use the minrad formula. */

  dot_prod = plc_M_dot(in,out);
  cross_prod_norm = plc_norm(plc_cross_prod(in,out));

  curv = 2*cross_prod_norm/(normin*normout + dot_prod);      
  rad = (normin*normout + dot_prod)/(2*cross_prod_norm);

  /* Minrad is the least of normin * rad and normout * rad. */
  /* We keep track of the reciprocal of minrad as well so that */
  /* we can deal effectively with straight segments. */

  mc = (normin < normout) ? (1/normin)*curv : (1/normout)*curv;
  mr = (normin < normout) ? normin * rad : normout * rad;

  if (mc < 1e-3) {

    mr = 1e3;  /* If minrad is too high, we simply compute as if it was 1000 */

    /* We could have tested mr directly, instead of computing mc, but if the segments */
    /* v0v1 and v1v2 are in a straight line, you get mr = nan and it doesn't work. */

  }

  /* a is the angle between the contact point of the circle and the mpsvects */

  if (normin < normout) {

    sin2a = (normin/2.0) / sqrt( pow(normin/2.0,2.0) + pow(mr,2.0) );
    cos2a = (mr) / sqrt( pow(normin/2.0,2.0) + pow(mr,2.0) );
    tana =  sin2a/(1 + cos2a);

    dist = (normin/2.0 - mr*tana); /* Actual distance from v1 to mpsverts. */

  } else {

    sin2a = (normout/2.0) / sqrt( pow(normout/2.0,2.0) + pow(mr,2.0) );
    cos2a = (mr) / sqrt( pow(normout/2.0,2.0) + pow(mr,2.0) );
    tana =  sin2a/(1 + cos2a);

    dist = (normout/2.0 - mr*tana); /* Actual distance from v1 to mpsverts. */

  }

  s0 = dist/normin; s1 = dist/normout;
    
  /* Now we generate the verts. */

  mpsverts[0] = plc_vweighted(s0,v1,v0);
  mpsverts[1] = plc_vweighted(s1,v1,v2);

}

plCurve *plc_double_verts(plCurve *L)
{

  plCurve *newL;
  int cmp,vt;
  int *newVerts,*newCc;
  bool *newOpen;

  //  plc_fix_cst(L);
  
  newVerts = calloc(L->nc,sizeof(int)); assert(newVerts != NULL); 
  newCc = calloc(L->nc,sizeof(int)); assert(newCc != NULL);
  newOpen = calloc(L->nc,sizeof(bool)); assert(newOpen != NULL);

  for(cmp=0;cmp<L->nc;cmp++) {

    newVerts[cmp] = 2*(L->cp[cmp].nv);

    /* Colors are a bit tricky, since there are either 
       0 (no color info)
       1 (all verts the same color)
       
       or

       L->cp[cmp].nv (colors change along the curve)
    */

    if (L->cp[cmp].cc == 0 || L->cp[cmp].cc == 1) { 

      newCc[cmp] = L->cp[cmp].cc; /* Just preserve info */

    } else {

      assert(L->cp[cmp].cc == L->cp[cmp].nv); /* We will spline colors, too. */
      newCc[cmp] = newVerts[cmp];

    }

    newOpen[cmp] = L->cp[cmp].open;

  }

  newL = plc_new(L->nc,newVerts,newOpen,newCc);
  assert(newL != NULL);

  /* We now attempt to translate the constraints from L to newL */

  plc_constraint *thisCst;

  for(thisCst = L->cst;thisCst != NULL;thisCst = thisCst->next) {

    if (thisCst->kind == fixed) {

      plc_set_fixed(newL,thisCst->cmp,thisCst->vert*2,thisCst->vect[0]);

    } else if (thisCst->kind == line) {

      plc_constrain_to_line(newL,thisCst->cmp,thisCst->vert*2,
			    2*thisCst->num_verts-1,thisCst->vect[0],thisCst->vect[1]);
    } else if (thisCst->kind == plane) {

      if (thisCst->vert == L->cp[thisCst->cmp].nv-1) { // Last has a special meaning 

	plc_constrain_to_plane(newL,thisCst->cmp,thisCst->vert*2+1,2*thisCst->num_verts-1,
			       thisCst->vect[0],thisCst->vect[1].c[0]);

      } else {

	plc_constrain_to_plane(newL,thisCst->cmp,thisCst->vert*2,2*thisCst->num_verts-1,
			       thisCst->vect[0],thisCst->vect[1].c[0]);
	
      }

    } 

  }

  /* Now we are ready to interpolate vertices */

  plc_fix_wrap(L); /* Just paranoia. This should do nothing. */

  for(cmp=0;cmp<L->nc;cmp++) {

    for(vt=0;vt<L->cp[cmp].nv;vt++) {

      plc_vector mps[2];
      mpsverts(L,cmp,vt,mps);

      newL->cp[cmp].vt[2*vt]   = mps[0]; 
      newL->cp[cmp].vt[2*vt+1] = mps[1];
     
    } 

    if (newL->cp[cmp].cc == 1) {

      newL->cp[cmp].clr[0] = L->cp[cmp].clr[0];

    } else {

      for(vt=0;vt<L->cp[cmp].nv-1;vt++) {

	newL->cp[cmp].clr[2*vt+1] = newL->cp[cmp].clr[2*vt] = L->cp[cmp].clr[vt];
	
      }

    }

  }

  plc_fix_wrap(newL);
  //  plc_fix_cst(newL);

  free(newVerts); 
  free(newCc);
  free(newOpen);

  return newL;

}
    

double plc_pointset_diameter(const plCurve * const L) {
  int cmp, vert, cmp2, vert2;
  double diameter, dist;

  assert(L != NULL);

  diameter = 0.0;
  for (cmp = 0; cmp < L->nc; cmp++) {
    for (vert = 0; vert < L->cp[cmp].nv; vert++) {
      for (cmp2 = cmp; cmp2 < L->nc; cmp2++) {
        for (vert2 = (cmp == cmp2) ? vert : 0; 
             vert2 < L->cp[cmp2].nv; vert2++) {
          dist = plc_M_distance(L->cp[cmp].vt[vert],
                                L->cp[cmp2].vt[vert2]);
          if (dist > diameter) {
            diameter = dist;
          }
        }
      }
    }
  }
  return diameter;
}

void plc_scale( plCurve *link, const double alpha)

     /* Scale link by alpha, including constraints, if any! */

{
  int cp, vt;
  plc_constraint *tc; 

  /* First, scale the vertices. */

  for(cp=0;cp<link->nc;cp++) {
    for(vt=0;vt<link->cp[cp].nv;vt++) {
      plc_M_scale_vect(alpha,link->cp[cp].vt[vt]);
    }
  }

  plc_fix_wrap(link);
  
  /* Now deal with the constraints. */

  for(tc=link->cst;tc!=NULL;tc=tc->next) {  /* Follow the linked list of csts. */

    if (tc->kind == unconstrained) { 
    } else if (tc->kind == fixed) {
      plc_M_scale_vect(alpha,tc->vect[0]);
    } else if (tc->kind == line) {
      plc_M_scale_vect(alpha,tc->vect[1]);
    } else if (tc->kind == plane) {
      plc_M_scale_vect(alpha,tc->vect[1]);
    }

  }
   
}   

