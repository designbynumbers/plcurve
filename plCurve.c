/*
 *  Routines to create, destroy, read and write plCurves (and strands)
 * 
 *  $Id: plCurve.c,v 1.49 2006-02-16 20:28:18 ashted Exp $
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif
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

#include <plCurve.h>

/* Utility routines */
static inline double doublemin(const double a, const double b) {
  return (a <= b) ? a : b;
}
static inline double doublemax(const double a, const double b) {
  return (a >= b) ? a : b;
}
static inline int intmin(const int a, const int b) {
  return (a <= b) ? a : b;
}
static inline int intmax(const int a, const int b) {
  return (a >= b) ? a : b;
}

/*
 * Set up a new strand.  Pl should point to an *ALREADY ALLOCATED* strand (but
 * with an unallocated space for vertices).  The number of vertices is given in
 * nv, the number of colors in cc, and open is set to TRUE or FALSE depending
 * on whether the strand is open or closed. 
 *
 * We allocate two extra vertices, at -1 and nv to make "wrap-around" much 
 * simpler.
 */
static inline void strand_new(plCurve_strand *Pl,int nv, int open, int cc) {
 
  if (nv < 1) {
    plcl_error_num = PLCL_E_TOO_FEW_VERTS;
    sprintf(plcl_error_str,
      "strand_new: Can't create a strand with %d vertices.\n",nv);
    return;
  }

  Pl->open = open;
  Pl->nv = nv;
  if ((Pl->vt = 
       (plcl_vector *)calloc((nv+2),sizeof(plcl_vector))) == NULL) {
    plcl_error_num = PLCL_E_CANT_ALLOC;
    sprintf(plcl_error_str,
      "strand_new: Can't allocate space for %d vertices in strand_new.\n",nv);
    return;
  }
  Pl->vt++; /* so that Pl->vt[-1] is a valid space */
  
  Pl->cc = cc;
  if ((Pl->clr = (plCurve_color *)calloc(cc,sizeof(plCurve_color))) == NULL) {
    plcl_error_num = PLCL_E_CANT_ALLOC;
    sprintf(plcl_error_str,
      "strand_new: Can't allocate space for %d colors in strand_new.\n",cc);
    return;
  }
}

/*
 * Procedure allocates memory for a new plCurve. The number of components is
 * given by components. The number of vertices in each component shows up in
 * the buffer pointed to by nv.  The closed/open nature of each strand is given
 * in the array pointed to by open, and the number of colors per strand is 
 * found in the array cc.
 *
 */
plCurve *plCurve_new(const int components, const int * const nv, 
                     const int * const open, const int * const cc) {
  plCurve *L;
  int i;

  /* First, we check to see that the input values are reasonable. */

  plcl_error_num = plcl_error_str[0] = 0;
  if (components < 0) {
    plcl_error_num = PLCL_E_TOO_FEW_COMPS;
    sprintf(plcl_error_str,
      "plCurve_new: Can't create a plCurve with %d components.\n",components);
    return NULL;
  }

  if (components > 0) {
    if (nv == NULL) {
      plcl_error_num = PLCL_E_NULL_PTR;
#ifdef HAVE_STRLCPY
      strlcpy(plcl_error_str,"plCurve_new: nv is NULL.\n",
        sizeof(plcl_error_str));
#else
      strlncpy(plcl_error_str,"plCurve_new: nv is NULL.\n",
        sizeof(plcl_error_str)-1);
      plcl_error_str[sizeof(plcl_error_str)-1] = '\0';
#endif
      return NULL;
    }
    if (open == NULL) {
      plcl_error_num = PLCL_E_NULL_PTR;
#ifdef HAVE_STRLCPY
      strlcpy(plcl_error_str,"plCurve_new: open is NULL.\n",
        sizeof(plcl_error_str));
#else
      strlncpy(plcl_error_str,"plCurve_new: open is NULL.\n",
        sizeof(plcl_error_str)-1);
      plcl_error_str[sizeof(plcl_error_str)-1] = '\0';
#endif
      return NULL;
    }
    if (cc == NULL) {
      plcl_error_num = PLCL_E_NULL_PTR;
#ifdef HAVE_STRLCPY
      strlcpy(plcl_error_str,"plCurve_new: cc is NULL.\n",
        sizeof(plcl_error_str));
#else
      strlncpy(plcl_error_str,"plCurve_new: cc is NULL.\n",
        sizeof(plcl_error_str)-1);
      plcl_error_str[sizeof(plcl_error_str)-1] = '\0';
#endif
      return NULL;
    }
  }

  /* Now we attempt to allocate space for these components. */
  
  if ((L = (plCurve *)malloc(sizeof(plCurve))) == NULL) {
    plcl_error_num = PLCL_E_CANT_ALLOC;
    sprintf(plcl_error_str,
      "plCurve_new: Could not allocate space for plCurve.\n");
    return NULL;
  }
  L->nc = components;
  if (L->nc > 0) {
    if ((L->cp = (plCurve_strand *)
      malloc(L->nc*sizeof(plCurve_strand))) == NULL) {
      plcl_error_num = PLCL_E_CANT_ALLOC;
      sprintf(plcl_error_str,
        "plCurve_new: Can't allocate array of strand ptrs.\n");
      return NULL;
    }
  }

  for (i = 0; i < L->nc; i++) {
    strand_new(&L->cp[i],nv[i],open[i],cc[i]);
  }

  /* Start off with no constraints */
  L->cst = NULL;

  return L;
} /* plCurve_new */

/*
 * Find the closest point on the given line to the given point.
 *
 * coef should hold 6 doulbes, which we will call a,b,c,d,e,f
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
                                             const double coef[6]) {
  plcl_vector ret_vect;
  double a = coef[0];
  double b = coef[1];
  double c = coef[2];
  double d = coef[3];
  double e = coef[4];
  double f = coef[5];
  double x = point.c[0];
  double y = point.c[1];
  double z = point.c[2];
  double denom = (a*a + c*c + e*e);
  double t = (a*(x-b) + c*(y-d) + e*(z-f))/(a*a + c*c + e*e);

  /* Make sure that we aren't about to divide by zero */
  if (denom == 0.0) {
    plcl_error_num = PLCL_E_BAD_CST;
#ifdef HAVE_STRLCPY
    strlcpy(plcl_error_str, "Closest_line_point: Constraint corrupted. "
      "a, c, and e are all 0.\n", sizeof(plcl_error_str));
#else
    strncpy(plcl_error_str, "Closest_line_point: Constraint corrupted. "
      "a, c, and e are all 0.\n", sizeof(plcl_error_str)-1);
    plcl_error_str[sizeof(plcl_error_str)-1] = '\0';
#endif
    plcl_M_set_vect(ret_vect,0,0,0);
    return ret_vect;
  }

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
                                              const double coef[6]) {

  plcl_vector ret_vect;
  double a = coef[0];
  double b = coef[1];
  double c = coef[2];
  double d = coef[3];
  double x1 = point.c[0];
  double y1 = point.c[1];
  double z1 = point.c[2];
  int i0 = 0;
  int i1 = 1;
  int i2 = 2;
  double denom = (a*a + b*b + c*c);
  double frac = (a*x1 + b*y1 + c*z1 - d)/denom;
  double x = x1 - a*frac;
  double y = y1 - b*frac;

  /* 
   * We need to make sure that we aren't dividing by some tiny number later on.
   */
  if (c <= 1e-12 && c >= -1e-12) {
    if (b <= 1e-12 && b >= -1e-12) {
      i0 = 2; 
      i2 = 0;
      c = a;
      a = coef[2];
    } else {
      i1 = 2;
      i2 = 1;
      c = b;
      b = coef[2];
    }
  }

  /* And that we aren't dividing by zero now */
  if (denom == 0) {
    plcl_error_num = PLCL_E_BAD_CST;
#ifdef HAVE_STRLCPY
    strlcpy(plcl_error_str, "Closest_plane_point: Constraint corrupted. "
      "First 3 coefficients are 0.\n", sizeof(plcl_error_str));
#else
    strncpy(plcl_error_str, "Closest_plane_point: Constraint corrupted. "
      "First 3 coefficients are 0.\n", sizeof(plcl_error_str)-1);
    plcl_error_str[sizeof(plcl_error_str)-1] = '\0';
#endif
    plcl_M_set_vect(ret_vect,0,0,0);
    return ret_vect;
  }

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

  plcl_error_num = plcl_error_str[0] = 0;

  /* Sanity check */
  if (L == NULL) {
    plcl_error_num = PLCL_E_NULL_PTR;
    sprintf(plcl_error_str, 
      "plCurve_cst_check: Called with NULL pointer.\n");
    return -1;
  }
  cst = L->cst;
  while (cst != NULL) {
    if (cst->kind == PLCL_FIXED) {
      closest.c[0] = cst->coef[0];
      closest.c[1] = cst->coef[1];
      closest.c[2] = cst->coef[2];
      /* PLCL_FIXED constraints only ever apply to one vertex */
      sq_dist = plcl_M_sq_dist(L->cp[cst->cmp].vt[cst->vert],closest);
      max_err = doublemax(max_err,sq_dist);
    } else if (cst->kind == PLCL_ON_LINE) {
      for (vert = cst->vert; vert < cst->vert + cst->num_verts; vert++) {
        closest = Closest_line_point(L->cp[cst->cmp].vt[vert],cst->coef);
        if (plcl_error_num != 0) {
          return -1;
        }
        sq_dist = plcl_M_sq_dist(L->cp[cst->cmp].vt[vert],closest);
        max_err = doublemax(max_err,sq_dist);
      } 
    } else if (cst->kind == PLCL_IN_PLANE) {
      for (vert = cst->vert; vert < cst->vert + cst->num_verts; vert++) {
        closest = Closest_plane_point(L->cp[cst->cmp].vt[vert],cst->coef);
        if (plcl_error_num != 0) {
          return -1;
        }
        sq_dist = plcl_M_sq_dist(L->cp[cst->cmp].vt[vert],closest);
        max_err = doublemax(max_err,sq_dist);
      } 
    } else {
      plcl_error_num = PLCL_E_BAD_CST_KIND;
      sprintf(plcl_error_str,
        "plCurve_cst_fix: Unknown constraint kind: %d.\n",cst->kind);
      return -1;
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

  plcl_error_num = plcl_error_str[0] = 0;

  /* Sanity check */
  if (L == NULL) {
    plcl_error_num = PLCL_E_NULL_PTR;
    sprintf(plcl_error_str, 
      "plCurve_cst_check: Called with NULL pointer.\n");
    return;
  }
  cst = L->cst;
  while (cst != NULL) {
    if (cst->kind == PLCL_FIXED) {
      closest.c[0] = cst->coef[0];
      closest.c[1] = cst->coef[1];
      closest.c[2] = cst->coef[2];
      /* PLCL_FIXED constraint is only ever allowed to be length 1 */
      L->cp[cst->cmp].vt[cst->vert] = closest;
    } else if (cst->kind == PLCL_ON_LINE) {
      for (vert = cst->vert; vert < cst->vert + cst->num_verts; vert++) {
        closest = Closest_line_point(L->cp[cst->cmp].vt[vert],cst->coef);
        if (plcl_error_num != 0) {
          return;
        }
        L->cp[cst->cmp].vt[vert] = closest;
      } 
    } else if (cst->kind == PLCL_IN_PLANE) {
      for (vert = cst->vert; vert < cst->vert + cst->num_verts; vert++) {
        closest = Closest_plane_point(L->cp[cst->cmp].vt[vert],cst->coef);
        if (plcl_error_num != 0) {
          return;
        }
        L->cp[cst->cmp].vt[vert] = closest;
      } 
    } else {
      plcl_error_num = PLCL_E_BAD_CST_KIND;
      sprintf(plcl_error_str,
        "plCurve_cst_fix: Unknown constraint kind: %d.\n",cst->kind);
      return;
    }
    cst = cst->next;
  }

  return;
} /* plCurve_fix_cst */

/*
 * Free the memory used to hold vertices in a given strand (not the memory of
 * the strand itself).  We then set all the values in the strand data structure
 * to reflect the fact that the memory has been freed.  We can call strand_free
 * twice on the same strand pointer without fear. 
 *
 */ 
static inline void strand_free(plCurve_strand *Pl) {
  
  if (Pl == NULL) {
    return;
  }

  Pl->nv = 0;
  if (Pl->vt != NULL) {
    Pl->vt--; /* undo our original vt++ (for wraparound) */
    free(Pl->vt);
    Pl->vt = NULL;
  }
 
  Pl->cc = 0;
  if (Pl->clr != NULL) {
    free(Pl->clr);
  }
} /* strand_free */

/*
 * Free the memory associated with a given plCurve.  We then set all the values
 * in the plCurve data structure to reflect the fact that the memory has been
 * freed.  We can call plCurve_free twice on the same plCurve pointer without
 * fear. 
 *
 */ 
void plCurve_free(plCurve *L) {
  int i;

  plcl_error_num = plcl_error_str[0] = 0;

  /* First, we check the input. */
  if (L == NULL) {
    return; /* Move along, nothing to see here */
  }

  /* Now we can get to work. */
  if (L->cp != NULL) {
    for (i=0; i<L->nc; i++) {
      strand_free(&L->cp[i]); /* strand_free is ok, even if L->cp[i] is NULL */
    }
  
    free(L->cp);
    L->nc = 0;
  }

  if (L->cst != NULL) {
    free(L->cst);
    L->cst = NULL;
  }

  free(L);
  L = NULL;
} /* plCurve_free */

static inline 
plCurve_constraint *new_constraint(const int kind, const double coef[6], 
                                   const int cmp, const int vert, 
                                   const int num_verts, 
                                   plCurve_constraint *next) {
  plCurve_constraint *ncst;

  if ((ncst = malloc(sizeof(plCurve_constraint))) == NULL) {
    plcl_error_num = PLCL_E_CANT_ALLOC;
    sprintf(plcl_error_str, 
      "new_constraint: Couldn't allocate space for new constraint.\n");
    return NULL;
  }
  ncst->kind = kind;
  memcpy(ncst->coef,coef,sizeof(ncst->coef));
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
      if (cst->next->kind == cst->kind &&
 //       cst->kind != PLCL_FIXED && /* Never extend a fixed constraint */
          memcmp(cst->next->coef,cst->coef,sizeof(cst->coef)) == 0) {
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
void plCurve_set_constraint(plCurve * const L, const int cmp, 
                            const int vert, const int num_verts, 
                            const int kind, const double coef0,
                            const double coef1, const double coef2, 
                            const double coef3, const double coef4, 
                            const double coef5) {
  double coef[6];
  plCurve_constraint *cst,*pcst;  /* constraint, previous constraint */
  plCurve_constraint **pfn; /* Place for new constraint */

  plcl_error_num = plcl_error_str[0] = 0;

  /* First, we check the input. */
  if (L == NULL) {
    plcl_error_num = PLCL_E_NULL_PTR;
    sprintf(plcl_error_str, 
      "plCurve_set_constraint: Called with NULL pointer.\n");
    return;
  }
  if (cmp < 0 || cmp >= L->nc) {
    plcl_error_num = PLCL_E_BAD_COMPONENT;
    sprintf(plcl_error_str,
      "plCurve_set_constraint: Component value out of range (0..%d): %d.\n",
      L->nc-1, cmp);
    return;
  }
  if (vert < 0 || vert >= L->cp[cmp].nv) {
    plcl_error_num = PLCL_E_BAD_VERTEX;
    sprintf(plcl_error_str,
      "plCurve_set_constraint: Vertex value out of range (0..%d): %d.\n",
      L->cp[cmp].nv-1, vert);
    return;
  }
  if (num_verts < 1) {
    plcl_error_num = PLCL_E_TOO_FEW_VERTS;
    sprintf(plcl_error_str,
      "plCurve_set_constraint: Can't set constraint on %d vertices.\n",
      num_verts);
    return;
  }
  if (L->cp[cmp].nv < vert+num_verts) {
    plcl_error_num = PLCL_E_TOO_MANY_VRTS;
    sprintf(plcl_error_str, "plCurve_set_constraint: Component only has %d "
      "verts, can't set %d--%d.\n", L->cp[cmp].nv, vert, vert+num_verts-1);
    return;
  }
  if (kind == PLCL_FIXED && num_verts > 1) {
    plcl_error_num = PLCL_E_TOO_MANY_VRTS;
    sprintf(plcl_error_str, "plCurve_set_constraint: Cannot fix %d "
      "consecutive vertices to the same spot.\n", num_verts);
    return;
  }

  /* Put these in one spot for easier passage */
  coef[0] = coef0;
  coef[1] = coef1;
  coef[2] = coef2;
  coef[3] = coef3;
  coef[4] = coef4;
  coef[5] = coef5;

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
    pfn = &pcst->next;
  }

  /* Now cst either points to a constraint which applies to the vertex in
   * question (or some vertex after it) or it points to NULL because there was
   * no such constraint found.  On the other hand, pcst either points to the
   * constraint just prior to the constraint which cst points to or else it is
   * NULL.  */
  if (pcst != NULL &&
      pcst->kind == kind &&
//    kind != PLCL_FIXED && /* Never extend a fixed constraint */
      memcmp(pcst->coef,coef,sizeof(coef)) == 0 &&
      pcst->cmp == cmp &&
      pcst->vert+pcst->num_verts == vert) {
    /* We just passed one which has a range contiguous with the range we
       want to set and has the same attributes.  We'll just extend it. */
    pcst->num_verts += num_verts;
    overrun_check(pcst);
  } else if (cst != NULL &&
             cst->kind == kind &&
//           kind != PLCL_FIXED && /* Never extend a fixed constraint */
             memcmp(pcst->coef,coef,sizeof(coef)) == 0 &&
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
      (*pfn) = new_constraint(kind, coef, cmp, vert, num_verts, cst);
      if (plcl_error_num != 0) { return; }
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
        (*pfn) = new_constraint(cst->kind, cst->coef, cst->cmp, cst->vert,
                   vert - cst->vert, cst);
        if (plcl_error_num != 0) { return; }
        if (kind != PLCL_UNCST) {
          (*pfn)->next = new_constraint(kind, coef, cmp, vert, num_verts, cst);
          if (plcl_error_num != 0) { return; }
        }
        cst->num_verts -= vert + num_verts - cst->vert;
        cst->vert = vert + num_verts;
      } else {
        /* It's just before ours, but extends into ours, shorten it. */
        cst->num_verts = vert - cst->vert;
        /* And attach ours to the end */
        if (kind != PLCL_UNCST) {
          cst->next = 
            new_constraint(kind, coef, cmp, vert, num_verts, cst->next);
          if (plcl_error_num != 0) { return; }
        }
        /* And check to see if we overran the next one */
        cst = cst->next;
        overrun_check(cst);
      }
    } else {
      /* It starts no earlier than ours, put ours in before it */
      (*pfn) = new_constraint(kind, coef, cmp, vert, num_verts, cst);
      if (plcl_error_num != 0) { return; }
      overrun_check(*pfn);
    }
  }

  return;
} /* plCurve_set_constraint */

/*
 * Remove a constraint from the list of constraints returning the number of
 * vertices thus set unconstrained. 
 *
 */
int plCurve_remove_constraint(plCurve * const L, plCurve_constraint *cst) {
  plCurve_constraint *cst_ptr;
  int uncst = 0;

  plcl_error_num = plcl_error_str[0] = 0;

  if (L == NULL || cst == NULL) {
    plcl_error_num = PLCL_E_NULL_PTR;
    sprintf(plcl_error_str, 
      "plCurve_remove_constraint: Called with NULL pointer.\n");
    return -1;
  }

  if (L->cst == cst) {
    /* Top one on the list */
    L->cst = cst->next;
    uncst = cst->num_verts;
    free(cst);
    return uncst;
  } else {
    cst_ptr = L->cst;
    while (cst_ptr->next != NULL &&
           cst_ptr->next != cst) {
      cst_ptr = cst_ptr->next;
    }
    if (cst_ptr->next == cst) {
      cst_ptr->next = cst->next;
      uncst = cst->num_verts;
      free(cst);
      return uncst;
    } else {
      plcl_error_num = PLCL_E_BAD_CST;
      sprintf(plcl_error_str, "plCurve_remove_constraint: Constraint not in "
        "list for plCurve given.\n");
      return -1;
    }
  }
} /* plCurve_remove_constraint */

/* Remove all constraints */
void plCurve_remove_all_constraints(plCurve * const L) {
  plCurve_constraint *cst;

  plcl_error_num = plcl_error_str[0] = 0;

  if (L == NULL) {
    plcl_error_num = PLCL_E_NULL_PTR;
    sprintf(plcl_error_str, 
      "plCurve_remove_all_constraints: Called with NULL pointer.\n");
    return;
  }

  cst = L->cst;
  while (cst != NULL) {
    L->cst = cst->next;
    free(cst);
    cst = L->cst;
  }
} /* plCurve_remove_all_constraints */

/* Set a vertex to the desired triple.  */
inline void plCurve_set_vert(plCurve * const L, const int cmp, const int vert,
                             const double x, const double y, const double z)
{
  plcl_error_num = plcl_error_str[0] = 0;

  /* First, we check the input. */
  if (L == NULL) {
    plcl_error_num = PLCL_E_NULL_PTR;
    sprintf(plcl_error_str, "plCurve_set_vert: Called with NULL pointer.\n");
    return;
  }
  if (cmp < 0 || cmp >= L->nc) {
    plcl_error_num = PLCL_E_BAD_COMPONENT;
    sprintf(plcl_error_str,
      "plCurve_set_vert: Component value out of range (0..%d): %d.\n",
      L->nc-1, cmp);
    return;
  }
  if (vert < 0 || vert >= L->cp[cmp].nv) {
    plcl_error_num = PLCL_E_BAD_VERTEX;
    sprintf(plcl_error_str,
      "plCurve_set_vert: Vertex value out of range (0..%d): %d.\n",
      L->cp[cmp].nv-1, vert);
    return;
  }
  L->cp[cmp].vt[vert].c[0] = x;
  L->cp[cmp].vt[vert].c[1] = y;
  L->cp[cmp].vt[vert].c[2] = z;

  return;
}

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
 * Returns TRUE if successful write, FALSE otherwise.
 *
 */

int plCurve_write(FILE *file, plCurve * const L) {
  int i,j;                  /* Counters for the for loops */ 
  int nverts = 0;           /* Total number of vertices of all components */
  int colors = 0;           /* Total number of colors of all components */
  char outstr[80] = "";     /* So we can wrap the vertex lines */
  plCurve_constraint *cst;  /* Current constraint */
  int cst_cnt;              /* Tally of constraints */

  plcl_error_num = plcl_error_str[0] = 0;

  /* First, do a little sanity checking. */
  if (file == NULL) {
    plcl_error_num = PLCL_E_NULL_PTR;
    sprintf(plcl_error_str,"plCurve_write: Passed NULL pointer as file.\n");
    return -1;
  }

  if (L == NULL) {
    plcl_error_num = PLCL_E_NULL_PTR;
    sprintf(plcl_error_str,
      "plCurve_write: Passed NULL pointer as plCurve.\n");
    return -1;
  }

  if (L->nc < 0) {
    plcl_error_num = PLCL_E_TOO_FEW_COMPS;
    sprintf(plcl_error_str,
      "plCurve_write: plCurve corrupted. L.nc = %d.\n",L->nc);
    return -1;
  }

  if (L->nc > 0 && L->cp == NULL) {
    plcl_error_num = PLCL_E_NULL_PTR;
#ifdef HAVE_STRLCPY
    strlcpy(plcl_error_str,
      "plCurve_write: plCurve corrupted. L.cp is NULL.\n",
      sizeof(plcl_error_str));
#else
    strncpy(plcl_error_str,
      "plCurve_write: plCurve corrupted. L.cp is NULL.\n",
      sizeof(plcl_error_str)-1);
    plcl_error_str[sizeof(plcl_error_str)-1] = '\0';
#endif
    return -1;
  }

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
  fprintf(file,"# Constraints\n");
  cst = L->cst;
  cst_cnt = 0;
  while (cst != NULL) {
    cst_cnt++;
    if (cst->kind == PLCL_FIXED) {
      fprintf(file, "#  %d Fixed %lg %lg %lg\n",cst_cnt, cst->coef[0],
        cst->coef[1], cst->coef[2]);
    } else if (cst->kind == PLCL_ON_LINE) {
      fprintf(file, "#  %d Line %lg %lg %lg %lg %lg %lg\n", cst_cnt,
        cst->coef[0], cst->coef[1], cst->coef[2],
        cst->coef[3], cst->coef[4], cst->coef[5]);
    } else if (cst->kind == PLCL_IN_PLANE) {
      fprintf(file, "#  %d Plane %lg %lg %lg %lg\n", cst_cnt,
        cst->coef[0], cst->coef[1], cst->coef[2], cst->coef[3]);
    } else {
      plcl_error_num = PLCL_E_BAD_CST_KIND;
      sprintf(plcl_error_str,
        "plCurve_write: Unknown constraint kind: %d.\n",cst->kind);
      return -1;
    }
    cst = cst->next;
  }

  fprintf(file,"# Vertex coordinates\n");

  /* Now we write the vertex data . . . */
  cst = L->cst;
  cst_cnt = 1;
  for (i=0; i<L->nc; i++) {
    fprintf(file,"# Component %d\n",i);
    for (j=0; j<L->cp[i].nv; j++) {
      sprintf(outstr,"%.16g %.16g %.16g", plcl_M_clist(L->cp[i].vt[j]));
      while (cst != NULL &&
             (cst->cmp < i ||
              (cst->cmp == i && cst->vert+cst->num_verts <= j))) {
        cst = cst->next;
        cst_cnt++;
      }
      /* Now either we have nothing, or something yet to come, or something
       * currently applicable. */
      if (cst != NULL && cst->cmp == i && cst->vert <= j) {
        /* Something applicable */
        sprintf(outstr,"%s # Cst: %d",outstr,cst_cnt);
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
  return 0;
}


/*
 * The next section of the library file includes some (private) procedures for
 * reading plCurve data reliably from Geomview VECT files.
 *
 */

/*
 * skip_whitespace_and_comments positions the file pointer on next
 * non-whitespace character, returning FALSE if EOF happens first. We skip
 * anything between a # and a newline.
 * 
 */
static inline int skip_whitespace_and_comments(FILE *infile)
{
  int thischar,commentflag = {FALSE};

  /* First, we check to make sure that infile looks legit. */
  if (infile == NULL) {
    plcl_error_num = PLCL_E_NULL_PTR;
    sprintf(plcl_error_str,
      "skip_whitespace_and_comments: infile is a null pointer.\n");
    return -1;
  }
  
  /* Now we start to work. */
  for(;;) {
    thischar = fgetc(infile);

    if (thischar == EOF) {  
      /* Reached end of file before a non-space, non-comment */
      return 0;
    } else if (thischar == '#') { /* Started a comment. */
      commentflag = TRUE;
    } else if (thischar == '\n' && commentflag) { /* End a comment. */
      commentflag = FALSE;
    } else if (!isspace(thischar) && !commentflag) { /* Found a hit! */
      ungetc(thischar,infile);
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

  if (infile == NULL) {
    plcl_error_num = PLCL_E_NULL_PTR;
    sprintf(plcl_error_str,"scandoubles: infile is a null pointer.\n");
    return -1;
  }

  if (ndoubles < 1) {
    plcl_error_num = PLCL_E_TOO_FEW_DBLS;
    sprintf(plcl_error_str,
      "scandoubles: ndoubles (%d) is less than one.\n", ndoubles);
    return -1;
  }

  va_start(ap,ndoubles);

  /* Now we're ready to work. */

  for (i=0;i<ndoubles;i++) {    /* We expect to exit from the loop by */
                                /* returning, but this is a safety.   */
    if (skip_whitespace_and_comments(infile) == 0) { /* Failed */
      return nconverted;
    }

    thisdouble = va_arg(ap,double *);
    if (fscanf(infile,"%lf",thisdouble) != 1) { /* We couldn't scan. */
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

  if (infile == NULL) {
    plcl_error_num = PLCL_E_NULL_PTR;
    sprintf(plcl_error_str,"scanints: infile is a null pointer.\n");
    return -1;
  }

  if (nints < 1) {
    plcl_error_num = PLCL_E_TOO_FEW_INTS;
    sprintf(plcl_error_str,"scanints: nints (%d) is less than one.\n",nints);
    return -1;
  }

  va_start(ap,nints);

  /* Now we're ready to work. */
  for (i=0;i<nints;i++) {       /* We expect to exit from the loop by */
                                /* returning, but this is a safety.   */
    if (skip_whitespace_and_comments(infile) == 0) { /* Failed */
      return nconverted;
    }
    thisint = va_arg(ap,int *);

    if (fscanf(infile,"%d",thisint) != 1) {     /* We couldn't scan. */
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

  plcl_error_num = plcl_error_str[0] = 0;

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
plCurve *plCurve_read(FILE *file) 
{
  plCurve *L;
  int nverts, ncomp, ncolors;
  int *nvarray, *open, *ccarray;
  int i, j;
  int nv;
  
  plcl_error_num = plcl_error_str[0] = 0;

  /* First, we check for the 'VECT' keyword. */
  if (fscanf(file," VECT ") == EOF) {
    plcl_error_num = PLCL_E_NO_VECT;
    sprintf(plcl_error_str,"plCurve_read: Couldn't find VECT keyword.\n");
    return NULL;
  }

  /* Now we read the three integers giving vertices, components, and colors. */

  if (scanints(file,3,&ncomp,&nverts,&ncolors) != 3) {
    plcl_error_num = PLCL_E_BAD_CVC_LINE;
    sprintf(plcl_error_str,
      "plCurve_read: Couldn't parse <ncomp> <nverts> <ncolors> line.\n");
    return NULL;
  }

  /* We now try to read the array of numbers of vertices. */

  nvarray = (int *)calloc(ncomp,sizeof(int));
  open    = (int *)calloc(ncomp,sizeof(int));
  ccarray = (int *)calloc(ncomp,sizeof(int));

  for(i=0;i<ncomp;i++) {
    if (scanints(file,1,&(nvarray[i])) != 1) {
      plcl_error_num = PLCL_E_BAD_CVRT_LINE;
      sprintf(plcl_error_str,"plCurve_read: Couldn't parse number"
              "of vertices in component %d.\n",i);    
      return NULL;
    }
    /* A negative number of vertices indicates a CLOSED component. */
    open[i] = (nvarray[i] >= 0);
    nvarray[i] = abs(nvarray[i]);
  }

  /* We have set nvarray and open and are ready to read the color data.  */

  for(i=0;i<ncomp;i++) {
    if (scanints(file,1,&(ccarray[i])) != 1) {
      plcl_error_num = PLCL_E_BAD_CLR_LINE;
      sprintf(plcl_error_str,"plCurve_read: Couldn't parse <ncolors>"
        "for component %d.\n", i);   
      return NULL;
    }
  }

  /* We now allocate the plCurve data structure. */

  L = plCurve_new(ncomp,nvarray,open,ccarray);

  /* done with temorary arrays */
  free(nvarray);
  free(open);
  free(ccarray);

  if (L == NULL) {   /* If we don't have this much memory, then return NULL. */
    sprintf(plcl_error_str,"plCurve_read: Error in plCurve_new: %d.\n",
      plcl_error_num);
    plcl_error_num = PLCL_E_NEW_FAILED;
    return NULL;
  }

  /* And get ready to read the actual data. */

  for(i = 0; i < ncomp; i++) {
    nv = L->cp[i].nv;
    for(j = 0; j < nv; j++) {
      if (scandoubles(file,3,plcl_M_clist(&L->cp[i].vt[j])) != 3) {
        plCurve_free(L);
        plcl_error_num = PLCL_E_BAD_VERT_LINE;
        sprintf(plcl_error_str,"plCurve_read: Couldn't parse "
          " <x> <y> <z> data for vertex %d of component %d.\n",j,i);
        return NULL;
      }
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
      if (scandoubles(file,4, &L->cp[i].clr[j].r, &L->cp[i].clr[j].g,
           &L->cp[i].clr[j].b, &L->cp[i].clr[j].alpha) != 4) {
        plcl_error_num = PLCL_E_BAD_COLOR;
        sprintf(plcl_error_str,"plCurve_read: Couldn't parse color %d "
          "in component %d of plCurve.\n",j,i);
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

  plcl_error_num = plcl_error_str[0] = 0;

  for (i=0;i<L->nc;i++) {
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
  
  plcl_error_num = plcl_error_str[0] = 0;
  plCurve_fix_wrap(L);
  
  /* Now we work. */

  in = plcl_vect_diff(L->cp[comp].vt[vert],L->cp[comp].vt[vert-1]);
  normin = plcl_M_norm(in);
  out = plcl_vect_diff(L->cp[comp].vt[vert+1],L->cp[comp].vt[vert]);
  normout = plcl_M_norm(out);
      
  dot_prod = plcl_M_dot(in,out);
  cross_prod_norm = plcl_M_norm(plcl_cross_prod(in,out));

  if (normin*normout + dot_prod < 1e-12) {
    plcl_error_num = PLCL_E_INF_KAPPA;
    sprintf(plcl_error_str,"plCurve_curvature: kappa not finite "
      "at vertex %d of component %d.\n",comp,vert);
    return -1.0;
  }
      
  kappa = (2*cross_prod_norm)/(normin*normout + dot_prod);
  kappa /= (normin < normout) ? normin : normout;

  return kappa;
}

/* 
 * Duplicate a plCurve and return the duplicate.
 *
 */
plCurve *plCurve_copy(plCurve * const L) {
  plCurve *nL;
  int *nv,*open,*ccarray;
  int cnt,cnt2;

  plcl_error_num = plcl_error_str[0] = 0;

  if ((nv = (int *)malloc((L->nc)*sizeof(int))) == NULL ||
      (open = (int *)malloc((L->nc)*sizeof(int))) == NULL ||
      (ccarray = (int *)malloc((L->nc)*sizeof(int))) == NULL) {
    plcl_error_num = PLCL_E_CANT_ALLOC;
    sprintf(plcl_error_str,
      "plCurve_copy: Unable to malloc space for new plCurve.\n");
    return NULL;
  }
  for (cnt = 0; cnt < L->nc; cnt++) {
    nv[cnt] = L->cp[cnt].nv;
    open[cnt] = L->cp[cnt].open;
    ccarray[cnt] = L->cp[cnt].cc;
  }
  nL = plCurve_new(L->nc,nv,open,ccarray);

  fprintf(stderr,"NEED TO HANDLE CONSTRAINTS ON COPY\n");
  exit(1);

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

  free(ccarray);
  free(open);
  free(nv);

  return nL;
}

/*
 * Compute a (unit) tangent vector to L at vertex vert of component cmp. 
 *
 */
plcl_vector plCurve_tangent_vector(plCurve * const L,int cmp, int vert) {
  plcl_vector in, out, tan;

  plcl_error_num = plcl_error_str[0] = 0;

  if (L->cp[cmp].open) {
    if (vert == 0) {
      tan = plcl_vect_diff(L->cp[cmp].vt[1], L->cp[cmp].vt[0]);

      return plcl_normalize_vect(tan);
    } else if (vert == L->cp[cmp].nv-1) {
       tan = plcl_vect_diff(L->cp[cmp].vt[L->cp[cmp].nv-1],
                            L->cp[cmp].vt[L->cp[cmp].nv-2]);

       return plcl_normalize_vect(tan);
    }
  }

  /* We now know that either we are on a closed 
     component, or we are not at an endpoint.   */
  
   in = plcl_normalize_vect(
     plcl_vect_diff(L->cp[cmp].vt[vert+1],L->cp[cmp].vt[vert])
   );

   out = plcl_normalize_vect(
     plcl_vect_diff(L->cp[cmp].vt[vert],L->cp[cmp].vt[vert-1])
   );

   plcl_M_vweighted(tan,0.5,in,out);
   return plcl_normalize_vect(tan);
} /* plCurve_tangent_vector */

/* Procedure computes the length of each component of the plCurve,
   and fills in the array of doubles "component_lengths", which 
   must be as long as L->nc. It returns the total length. We assume
   that fix_wrap has been called. */
double plCurve_arclength(const plCurve * const L,double *component_lengths)
{
  double tot_length;
  int cmp, nv, vert;
  plCurve_strand *cp;

  plcl_error_num = plcl_error_str[0] = 0;

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

  plcl_error_num = plcl_error_str[0] = 0;

  if (L == NULL) {
    plcl_error_num = PLCL_E_NULL_PTR;
    sprintf(plcl_error_str,
      "plCurve_parameter: Passed NULL pointer as plCurve.\n");
    return -1;
  }

  if (cmp < 0 || cmp >= L->nc) {
    plcl_error_num = PLCL_E_BAD_COMPONENT;
    sprintf(plcl_error_str,
      "plCurve_parameter: Component value out of range (0..%d): %d.\n",
      L->nc-1, cmp);
    return -1;
  }

  if (vert < 0 || vert >= L->cp[cmp].nv) {
    plcl_error_num = PLCL_E_BAD_VERTEX;
    sprintf(plcl_error_str,
      "plCurve_parameter: Vertex value out of range (0..%d): %d.\n",
      L->cp[cmp].nv-1, vert);
    return -1;
  }

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

  plcl_error_num = plcl_error_str[0] = 0;

  for (cmp=0;cmp < L->nc;cmp++) {
    if (L->cp[cmp].open == TRUE) {  /* Isolate the open components. */

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
      assert(plcl_M_norm(diff) < 1e-10);

      /* Thus we eliminate the last vertex. */
      L->cp[cmp].nv--;
      L->cp[cmp].open = FALSE;
    }
  }

  plCurve_fix_wrap(L);
}

/*
 * Check plcl_error_num, report if there is an error or warning and exit if it
 * is an error.
 *
 */
inline void plcl_status_check() {
  if (plcl_error_num != 0) {
    fprintf(stderr,"%s",plcl_error_str);
    if (plcl_error_num > 0) {
      exit(plcl_error_num);
    }
  }
}

/* 
 * Either return or display the library version number.
 *
 */
inline void plCurve_version(char *version) {
  if (version == NULL) {
    printf("plCurve Version: %s\n",PACKAGE_VERSION);
  } else {
    sprintf(version,PACKAGE_VERSION);
  }
}
