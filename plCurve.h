/* plCurve.h.  Generated by configure.  */
/*
 *
 * Data structures and prototypes for the plCurve library
 *
 *  $Id: plCurve.h,v 1.42 2006-03-01 23:12:46 ashted Exp $
 *
 */

/* Copyright 2004 The University of Georgia */

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

/*@-exportlocal@*/
#ifndef PLCURVE_H
#define PLCURVE_H

#if (__cplusplus || c_plusplus)
extern "C" {
#endif

/*
 * Take our chances that stdio.h and stdbool.h exist.  We need stdio to define
 * the FILE type and stdbool to define bool.  If stdbool doesn't exist, one
 * option is to just
 *   typedef int bool;
 * and
 *   #define true 1
 *   #define false 0
 */
#include <stdio.h>
#include <stdbool.h>

/* Define 3-space vectors */
typedef struct plcl_vector_type {
  double c[3];
} plcl_vector;

typedef struct plCurve_color_type {
  double r;
  double g;
  double b;
  double alpha;
} plCurve_color;

typedef struct plCurve_strand_type {
  int            nv;     /* Number of vertices */
  bool           open;   /* This is an "open" strand (with distinct ends) */
  int            cc;     /* Color count (number of colors) */
  plcl_vector   *vt;     /* Actual vertices */
  plCurve_color *clr;    /* Colors */
} plCurve_strand;

/* Curve constraint kind */
typedef enum plCurve_constraint_kind {
  unconstrained = 0,
  fixed,
  line,
  plane,
} plCurve_cst_kind;

typedef struct plCurve_constraint_type {
  plCurve_cst_kind kind;    /* What kind of constraint */
  plcl_vector      vect[2]; /* Vectors to define plane, line or fixed point */
  int              cmp;     /* Component */
  int              vert;    /* Starting vertex */
  int              num_verts; /* Length of run */
  /*@only@*/ /*@null@*/ struct plCurve_constraint_type *next;
} plCurve_constraint;

typedef struct plCurve_vert_quant_type { /* Vertex quantifiers */
  int              cmp;    /* Component */
  int              vert;   /* Vertex */
  char             tag[4]; /* 3-character tag */
  double           quant;  /* Quantifier */
  /*@only@*/ /*@null@*/ struct plCurve_vert_quant_type *next;
} plCurve_vert_quant;

typedef struct plCurve_type {
  int             nc;                              /* Number of components */
  plCurve_strand *cp;                              /* Components */
  /*@only@*/ /*@null@*/ plCurve_constraint *cst;   /* Constraints */
  /*@only@*/ /*@null@*/ plCurve_vert_quant *quant; /* per-vertex quantifiers */
} plCurve;

/* PlCurve_spline types */
typedef struct plCurve_spline_strand_type {
  bool            open;     /* This is an "open" strand (with distinct ends) */
  int             ns;       /* Number of samples used to build spline. */
  double         *svals;    /* s values at samples */
  plcl_vector    *vt;       /* positions at these s values */
  plcl_vector    *vt2;      /* second derivatives at these s values */
  int             cc;
  plCurve_color  *clr;      /* color values at samples */
  plCurve_color  *clr2;     /* second derivatives at these s values */
} plCurve_spline_strand;

typedef struct plCurve_spline_type {
  int                    nc;     /* Number of components */
  plCurve_spline_strand *cp;     /* Components */
} plCurve_spline;

/*
 * Prototypes for vector routines.
 *
 */

plcl_vector plcl_vect_sum(plcl_vector A,plcl_vector B);   /* A + B */
plcl_vector plcl_vect_diff(plcl_vector A,plcl_vector B);  /* A - B */
plcl_vector plcl_cross_prod(plcl_vector A,plcl_vector B); /* A x B */
plcl_vector plcl_scale_vect(double s,plcl_vector A);      /* sA */
plcl_vector plcl_normalize_vect(const plcl_vector V,
                                /*@null@*/ bool *ok);/*V/|V|*/
plcl_vector plcl_random_vect(void);

/* Translate three doubles into a vector */
plcl_vector plcl_build_vect(const double x, const double y, const double z);

/* Multiply or divide two ordered triples componetwise */
plcl_vector plcl_component_mult(plcl_vector A, plcl_vector B);
plcl_vector plcl_component_div(plcl_vector A, plcl_vector B,
                               /*@null@*/ bool *ok);

/* Return a linear combination: a*A + b*B */
plcl_vector plcl_vlincomb(double a,plcl_vector A, double b,plcl_vector B);
plcl_vector plcl_vmadd(plcl_vector A, double s, plcl_vector B); /* A + sB */
plcl_vector plcl_vweighted(double s, plcl_vector A, plcl_vector B);
  /* (1-s)A + sB */

/* Different vector measurements */
double plcl_dot_prod(plcl_vector A,plcl_vector B);
double plcl_norm(plcl_vector A);
double plcl_distance(plcl_vector A, plcl_vector B);
double plcl_sq_dist(plcl_vector A, plcl_vector B);

/* Do two vectors match ? */
bool plcl_vecteq(plcl_vector A, plcl_vector B);

/*
 * Macros for vector work (require careful programming)
 *
 */

#define plcl_M_dot(A,B)      \
  ((A).c[0]*(B).c[0] + (A).c[1]*(B).c[1] + (A).c[2]*(B).c[2])

#define plcl_M_norm(A)       \
  sqrt(plcl_M_dot((A),(A)))

#define plcl_M_add_vect(A,B)    \
  (A).c[0] += (B).c[0]; (A).c[1] += (B).c[1]; (A).c[2] += (B).c[2];

#define plcl_M_sub_vect(A,B)    \
  (A).c[0] -= (B).c[0]; (A).c[1] -= (B).c[1]; (A).c[2] -= (B).c[2];

#define plcl_M_scale_vect(s,V)  \
  (V).c[0] *= s; (V).c[1] *= s; (V).c[2] *= s;

/* Add a multiple of B to A */
#define plcl_M_vmadd(A,s,B) \
  (A).c[0] += (s)*(B).c[0]; (A).c[1] += (s)*(B).c[1]; (A).c[2] += (s)*(B).c[2];

/* A becomes a linear combination of B and C */
#define plcl_M_vlincomb(A,s,B,t,C) \
  (A).c[0] = s*(B).c[0] + t*(C).c[0]; \
  (A).c[1] = s*(B).c[1] + t*(C).c[1]; \
  (A).c[2] = s*(B).c[2] + t*(C).c[2];

/* A = B + s(C-B)                               *
 * equivalent to                                *
 *   A = C; vsub(A,B); vsmult(s,A); vsadd(A,B); */
#define plcl_M_vweighted(A,s,B,C)  \
  plcl_M_vlincomb(A,(1.0-s),B,s,C)

#define plcl_M_component_mult(A,B) \
  (A).c[0] *= (B).c[0]; \
  (A).c[1] *= (B).c[1]; \
  (A).c[2] *= (B).c[2];

#define plcl_M_component_div(A,B) \
  (A).c[0] /= (B).c[0]; \
  (A).c[1] /= (B).c[1]; \
  (A).c[2] /= (B).c[2];

#define plcl_M_distance(A,B) \
  plcl_norm(plcl_vect_diff((A),(B)));

/* The squared distance from A to B */
#define plcl_M_sq_dist(A,B) \
  ((A).c[0]-(B).c[0])*((A).c[0]-(B).c[0])+ \
  ((A).c[1]-(B).c[1])*((A).c[1]-(B).c[1])+ \
  ((A).c[2]-(B).c[2])*((A).c[2]-(B).c[2]);

/* The coordinates of a vector, as a list */
#define plcl_M_clist(A) \
  A.c[0], A.c[1], A.c[2]

/* Are two vectors equal? */
#define plcl_M_vecteq(A,B) \
  ((A).c[0] - (B).c[0] < DBL_EPSILON && -((A).c[0]-(B).c[0]) < DBL_EPSILON && \
   (A).c[1] - (B).c[1] < DBL_EPSILON && -((A).c[1]-(B).c[1]) < DBL_EPSILON && \
   (A).c[2] - (B).c[2] < DBL_EPSILON && -((A).c[2]-(B).c[2]) < DBL_EPSILON)

/*
 * Prototypes for routines to deal with plCurves.
 *
 */

/* Build a new plCurve (with associated strands) */
/*@only@*/ plCurve *plCurve_new(const int components,
                                const int * const nv,
                                const bool * const open,
                                const int * const cc);

/* Free the plCurve (and strands) */
void plCurve_free(/*@only@*/ /*@null@*/ plCurve *L);

/* Set a constraint on a vertex or run of vertices */
void plCurve_set_fixed(plCurve * const L,
                       const int          cmp,
                       const int          vert,
                       const plcl_vector point);

void plCurve_constrain_to_line(plCurve * const L,
                               const int          cmp,
                               const int          vert,
                               const int          num_verts,
                               const plcl_vector tangent,
                               const plcl_vector point_on_line);

void plCurve_constrain_to_plane(plCurve * const L,
                                const int          cmp,
                                const int          vert,
                                const int          num_verts,
                                const plcl_vector normal,
                                const double dist_from_origin);

void plCurve_unconstrain(plCurve * const L, const int cmp,
                         const int vert, const int num_verts);

/* Remove a constraint from the list of constraints returning the number of
 * vertices thus set unconstrained.  */
int plCurve_remove_constraint(plCurve * const L,
                              const plCurve_cst_kind kind,
                              const plcl_vector vect[]);

/* Remove all constraints */
void plCurve_remove_all_constraints(plCurve * const L);

/* Read plCurve data from a file */
/*@only@*/ /*@null@*/ plCurve *plCurve_read(FILE *file,
                                  /*@out@*/ int *error_num,
                                  /*@out@*/ char error_str[],
                                            size_t error_str_len);

/* Write plCurve data to a file */
void plCurve_write(FILE *outfile, plCurve * const L);

/* Fix the "hidden vertices" for easy handling of closed components */
void plCurve_fix_wrap(plCurve * const L);

/* Count the edges in a plCurve (correctly handling open/closed) */
int plCurve_num_edges(plCurve * const L);

/* Compute the (minrad-based) curvature of L at vertex vt of component cp */
double plCurve_curvature(plCurve * const L, const int cmp, const int vert);

/* Copy a plCurve */
plCurve *plCurve_copy(plCurve * const L);

/* Compute tangent vector */
plcl_vector plCurve_tangent_vector(const plCurve * const L,
                                   const int cmp,
                                   const int vert,
                                   bool *ok);

/* Find the arclength of a plCurve. */
double plCurve_arclength(const plCurve * const L,double *component_lengths);

/* Find the arclength position of a point on a plCurve. */
double plCurve_parameter(const plCurve * const L,
                         const int cmp,
                         const int vert);

/* Force a plCurve to close as gently as possible. */
void plCurve_force_closed(plCurve * const L);

/* Return how far a constraint is from being satisfied (sup norm). */
double plCurve_check_cst(const plCurve * const L);

/* Fix all the vertices which are out of compliance with their constraints. */
void plCurve_fix_cst(plCurve * const L);

/* Either return (if given a char *) or print out the library version number */
void plCurve_version(/*@null@*/ char *version, size_t strlen);

/* Put 4 doubles together into a color */
plCurve_color plCurve_build_color(const double r,
                                  const double g,
                                  const double b,
                                  const double alpha);

/* Allocate new spline. */
plCurve_spline *plCurve_spline_new(const int          components,
                                   const int  * const ns,
                                   const bool * const open,
                                   const int  * const cc);

/* Free memory for spline. */
void plCurve_spline_free(/*@only@*/ /*@null@*/ plCurve_spline *L);

/* Convert plCurve to spline representation. */
plCurve_spline *plCurve_convert_to_spline(plCurve * const L, bool *ok);

/* Convert splined curve to plCurve (with resampling). */
plCurve *plCurve_convert_from_spline(const plCurve_spline * const spL,
                                     const int          * const nv);

/* Samples a spline at a particular s value. */
plcl_vector plCurve_sample_spline(const plCurve_spline * const spL,
                                  const int          cmp,
                                  double s);
/* Define the error codes */
#define PLCL_E_NO_VECT       1
#define PLCL_E_BAD_CVC_LINE  2
#define PLCL_E_BAD_CMP_NUM   3
#define PLCL_E_BAD_CVRT_LINE 4
#define PLCL_E_BAD_CLR_LINE  5
#define PLCL_E_BAD_COLOR     6
#define PLCL_E_BAD_VERT_LINE 7
#define PLCL_E_BAD_CST_NUM   8
#define PLCL_E_BAD_CST_LINE  9
#define PLCL_E_BAD_CST_KIND  10
#define PLCL_E_BAD_CST_NUMS  11
#if (__cplusplus || c_plusplus)
};
#endif
#endif /* PLCURVE_H */
/*@=exportlocal@*/
