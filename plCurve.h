/* plCurve.h.  Generated by configure.  */
/*
 *
 * Data structures and prototypes for the plCurve library
 *
 *  $Id: plCurve.h,v 1.57 2007-07-12 15:27:52 cantarel Exp $
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
typedef struct plc_vector_type {
  double c[3];
} plc_vector;

typedef struct plc_color_type {
  double r;
  double g;
  double b;
  double alpha;
} plc_color;

typedef struct plc_strand_type {
  int            nv;     /* Number of vertices */
  bool           open;   /* This is an "open" strand (with distinct ends) */
  int            cc;     /* Color count (number of colors) */
  plc_vector   *vt;     /* Actual vertices */
  plc_color    *clr;    /* Colors */
} plc_strand;

/* Curve constraint kind */
typedef enum plc_cst_kind_type {
  unconstrained = 0,
  fixed,
  line,
  plane,
} plc_cst_kind;

typedef struct plc_constraint_type {
  plc_cst_kind kind;    /* What kind of constraint */
  plc_vector   vect[2]; /* Vectors to define plane, line or fixed point */
  int           cmp;     /* Component */
  int           vert;    /* Starting vertex */
  int           num_verts; /* Length of run */
  /*@only@*/ /*@null@*/ struct plc_constraint_type *next;
} plc_constraint;

  /* The vectors in vect change meaning depending on the type of constraint:

  fixed : vect[0] is the fixed point, vect[1] is ignored 
  line:   vect[0] is the tangent vector, vect[1] is a point on the line
  plane:  vect[0] is the normal vector, vect[1].x is the distance from origin.

  */

typedef struct plc_vert_quant_type { /* Vertex quantifiers */
  int              cmp;    /* Component */
  int              vert;   /* Vertex */
  char             tag[4]; /* 3-character tag */
  double           quant;  /* Quantifier */
  /*@only@*/ /*@null@*/ struct plc_vert_quant_type *next;
} plc_vert_quant;

typedef struct plc_type {
  int             nc;                           /* Number of components */
  plc_strand *cp;                           /* Components */
  /*@only@*/ /*@null@*/ plc_constraint *cst;   /* Constraints */
  /*@only@*/ /*@null@*/ plc_vert_quant *quant; /* per-vertex quantifiers */
} plCurve;

/* PlCurve_spline types */
typedef struct plc_spline_strand_type {
  bool         open;     /* This is an "open" strand (with distinct ends) */
  int          ns;       /* Number of samples used to build spline. */
  double      *svals;    /* s values at samples */
  plc_vector *vt;       /* positions at these s values */
  plc_vector *vt2;      /* _second_ derivatives at these s values */
  int          cc;       /* Number of colors */
  plc_color  *clr;      /* color values at samples */
} plc_spline_strand;

typedef struct plc_spline_type {
  int                 nc;     /* Number of components */
  plc_spline_strand *cp;     /* Components */
} plc_spline;

/*
 * Prototypes for vector routines.
 *
 */

plc_vector plc_vect_sum(plc_vector A,plc_vector B);   /* A + B */
plc_vector plc_vect_diff(plc_vector A,plc_vector B);  /* A - B */
plc_vector plc_cross_prod(plc_vector A,plc_vector B); /* A x B */
plc_vector plc_scale_vect(double s,plc_vector A);      /* sA */
plc_vector plc_normalize_vect(const plc_vector V,
                                /*@null@*/ bool *ok);     /* V/|V| */

plc_vector plc_random_vect(void);

/* Translate three doubles into a vector */
plc_vector plc_build_vect(const double x, const double y, const double z);

/* Multiply or divide two ordered triples componetwise */
plc_vector plc_component_mult(plc_vector A, plc_vector B);
plc_vector plc_component_div(plc_vector A, plc_vector B,
                               /*@null@*/ bool *ok);

/* Return a linear combination: a*A + b*B */
plc_vector plc_vlincomb(double a,plc_vector A, double b,plc_vector B);
plc_vector plc_vmadd(plc_vector A, double s, plc_vector B); /* A + sB */
plc_vector plc_vweighted(double s, plc_vector A, plc_vector B); /* (1-s)A + sB */

/* More sophisticated geometric operations */
plc_vector plc_circumcenter(plc_vector A, plc_vector B, plc_vector C,
			    double *circumradius,bool *ok);
  /* Returns the center and radius of the circle through three points. 
     If circumradius is NULL, then it won't be written to. */
plc_vector plc_normal(plc_vector A, plc_vector B, plc_vector C, bool *ok);
  /* Returns the normal vector to plane defined by A,B,C. */

plc_vector plc_3plane_intersection(plc_vector N0, plc_vector p0,
				   plc_vector N1, plc_vector p1,
				   plc_vector N2, plc_vector p2,
				   bool *ok);
/* Returns the intersection point of 3 planes. */ 					

double plc_tetrahedron_inradius(plc_vector A,plc_vector B,plc_vector C,plc_vector D);
/* Returns the inradius of a tetrahedron. */	
					  
/* Different vector measurements */
double plc_dot_prod(plc_vector A,plc_vector B);
double plc_norm(plc_vector A);
double plc_distance(plc_vector A, plc_vector B);
double plc_sq_dist(plc_vector A, plc_vector B);
double plc_angle(plc_vector A, plc_vector B, bool *ok);

/* Do two vectors match ? */
bool plc_vecteq(plc_vector A, plc_vector B);

/*
 * Macros for vector work (require careful programming)
 *
 */

#define plc_vect_copy(A,B) \
  (A) = (B)

#define plc_M_dot(A,B)      \
  ((A).c[0]*(B).c[0] + (A).c[1]*(B).c[1] + (A).c[2]*(B).c[2])

#define plc_M_cross(A,B,C) \
  (A).c[0] = (B).c[1] * (C).c[2] - (B).c[2] * (C).c[1]; \
  (A).c[1] = (B).c[2] * (C).c[0] - (B).c[0] * (C).c[2]; \
  (A).c[2] = (B).c[0] * (C).c[1] - (B).c[1] * (C).c[0];

#define plc_M_norm(A)       \
  sqrt(plc_M_dot((A),(A)))

#define plc_M_add_vect(A,B)    \
  (A).c[0] += (B).c[0]; (A).c[1] += (B).c[1]; (A).c[2] += (B).c[2];

#define plc_M_sub_vect(A,B)    \
  (A).c[0] -= (B).c[0]; (A).c[1] -= (B).c[1]; (A).c[2] -= (B).c[2];

#define plc_M_scale_vect(s,V)  \
  (V).c[0] *= s; (V).c[1] *= s; (V).c[2] *= s;

/* Add a multiple of B to A */
#define plc_M_vmadd(A,s,B) \
  (A).c[0] += (s)*(B).c[0]; (A).c[1] += (s)*(B).c[1]; (A).c[2] += (s)*(B).c[2];

/* A becomes a linear combination of B and C */
#define plc_M_vlincomb(A,s,B,t,C) \
  (A).c[0] = s*(B).c[0] + t*(C).c[0]; \
  (A).c[1] = s*(B).c[1] + t*(C).c[1]; \
  (A).c[2] = s*(B).c[2] + t*(C).c[2];

/* A = B + s(C-B)                               *
 * equivalent to                                *
 *   A = C; vsub(A,B); vsmult(s,A); vsadd(A,B); */
#define plc_M_vweighted(A,s,B,C)  \
  plc_M_vlincomb(A,(1.0-s),B,s,C)

#define plc_M_component_mult(A,B) \
  (A).c[0] *= (B).c[0]; \
  (A).c[1] *= (B).c[1]; \
  (A).c[2] *= (B).c[2];

#define plc_M_component_div(A,B) \
  (A).c[0] /= (B).c[0]; \
  (A).c[1] /= (B).c[1]; \
  (A).c[2] /= (B).c[2];

#define plc_M_distance(A,B) \
  plc_norm(plc_vect_diff((A),(B)));

/* The squared distance from A to B */
#define plc_M_sq_dist(A,B) \
  ((A).c[0]-(B).c[0])*((A).c[0]-(B).c[0])+ \
  ((A).c[1]-(B).c[1])*((A).c[1]-(B).c[1])+ \
  ((A).c[2]-(B).c[2])*((A).c[2]-(B).c[2]);

/* The coordinates of a vector, as a list */
#define plc_M_clist(A) \
  A.c[0], A.c[1], A.c[2]

/* Are two vectors equal? */
#define plc_M_vecteq(A,B) \
  (  (A).c[0] - (B).c[0]  <= 2*DBL_EPSILON && \
   -((A).c[0] - (B).c[0]) <= 2*DBL_EPSILON && \
     (A).c[1] - (B).c[1]  <= 2*DBL_EPSILON && \
   -((A).c[1] - (B).c[1]) <= 2*DBL_EPSILON && \
     (A).c[2] - (B).c[2]  <= 2*DBL_EPSILON && \
   -((A).c[2] - (B).c[2]) <= 2*DBL_EPSILON)

/*
 * Prototypes for routines to deal with plCurves.
 *
 */

/* Build a new plCurve (with associated strands) */
/*@only@*/ plCurve *plc_new(const int components,
                            const int * const nv,
                            const bool * const open,
                            const int * const cc);

/* Free the plCurve (and strands) */
void plc_free(/*@only@*/ /*@null@*/ plCurve *L);

/* Add a component to the curve which will become component number add_as. */
void plc_add_component(plCurve *L, const int add_as, const int nv, 
                       const bool open, const int cc,
                       const plc_vector * const vt,
            /*@null@*/ const plc_color  * const clr);

/* And remove one */
void plc_drop_component(plCurve *L, const int cmp);

/* Change the size of the color buffer for a plCurve, preserving existing data if it exists */
void plc_resize_colorbuf(plCurve *L, const int cp, const int cc);

/* Set a constraint on a vertex or run of vertices */
void plc_set_fixed(plCurve * const L,
                   const int cmp,
                   const int vert,
                   const plc_vector point);

void plc_constrain_to_line(plCurve * const L,
                           const int cmp,
                           const int vert,
                           const int num_verts,
                           const plc_vector tangent,
                           const plc_vector point_on_line);

void plc_constrain_to_plane(plCurve * const L,
                            const int cmp,
                            const int vert,
                            const int num_verts,
                            const plc_vector normal,
                            const double dist_from_origin);

void plc_unconstrain(plCurve * const L, const int cmp,
                     const int vert, const int num_verts);

/* Remove a constraint from the list of constraints returning the number of
 * vertices thus set unconstrained.  */
int plc_remove_constraint(plCurve * const L,
                          const plc_cst_kind kind,
                          const plc_vector vect[]);

/* Remove all constraints */
void plc_remove_all_constraints(plCurve * const L);

/* Read plCurve data from a file */
/*@only@*/ /*@null@*/ plCurve *plc_read(FILE *file,
                              /*@out@*/ int *error_num,
                              /*@out@*/ char error_str[],
                                        size_t error_str_len);

/* Write plCurve data to a file */
void plc_write(FILE *outfile, plCurve * const L);

/* Fix the "hidden vertices" for easy handling of closed components */
void plc_fix_wrap(plCurve * const L);

/* Count the edges in a plCurve (correctly handling open/closed) */
/* Deprecated in versions > 1.3 in favor of plc_edges call below. */
int plc_num_edges(const plCurve * const L);

/* Count edges in plCurve, returning total and storing #edges for */
/* each component in component_edges if this is non-NULL. */
int plc_edges(const plCurve * const L,
/*@null@*/ /*@out@*/ int *component_edges);

/* Count the vertices in a plCurve */
int plc_num_verts(const plCurve * const L);

/* Compute an index between 0 and plc_num_verts(L) - 1 for a (cp,vt) pair in a plCurve, 
   using full wraparound addressing for closed components and repeating the last or first vertex 
   for open ones. */
int plc_vertex_num(const plCurve * const L, int cp, int vt);

/* Compute the MinRad-based curvature of L at vertex vt of component cp */
double plc_MR_curvature(plCurve * const L, const int cmp, const int vert);

/* Copy a plCurve */
plCurve *plc_copy(const plCurve * const L);

/* Compute average of inward and outward tangents (and normalize) */
plc_vector plc_mean_tangent(const plCurve * const L, const int cmp,
                            const int vert, bool *ok);

/* Find the arclength of a plCurve. Total arclength is returned, arclength */
/* of individual strands stored in component_lengths if this pointer is non-NULL. */
double plc_arclength(const plCurve * const L,
/*@null@*/ /*@out@*/ double *component_lengths);

/* Find the arclength distance from one vertex to another.  On closed
 * strands, give the shortest of the two options.  */
double plc_subarc_length(const plCurve * const L, const int cmp,
                         const int vert1, const int vert2);

/* Return how far a constraint is from being satisfied (sup norm). */
double plc_check_cst(const plCurve * const L);

/* Fix all the vertices which are out of compliance with their constraints. */
void plc_fix_cst(plCurve * const L);

/* Either return (if given a char *) or print out the library version number */
void plc_version(/*@null@*/ char *version, size_t strlen);

/* Put 4 doubles together into a color */
plc_color plc_build_color(const double r, const double g,
                          const double b, const double alpha);

/* Allocate new spline. */
plc_spline *plc_spline_new(const int          components,
                           const int  * const ns,
                           const bool * const open,
                           const int  * const cc);

/* Free memory for spline. */
void plc_spline_free(/*@only@*/ /*@null@*/ plc_spline *L);

/* Convert plCurve to spline representation. */
plc_spline *plc_convert_to_spline(plCurve * const L, bool *ok);

/* Convert splined curve to plCurve (with resampling). */
plCurve *plc_convert_from_spline(const plc_spline * const spL,
                                 const int * const nv);

/* Samples a spline at a particular s value. */
plc_vector plc_sample_spline(const plc_spline * const spL,
                             const int cmp,
                             double s);

/* Doubles the number of vertices of L by inserting new vertices at midpoints
   of edges. Attempts to preserve constraints. */
plCurve *plc_double_verts(plCurve * L);


/* Calculate the diameter of the plCurve, thinking of the vertices as 
   a set of points in R^3 */
double plc_pointset_diameter(const plCurve * const L);

/* Scale a plCurve (and its' constraints!) by a factor. */
void plc_scale( plCurve *link, const double alpha);  

/* Define the error codes */
#define PLC_E_NO_VECT       1
#define PLC_E_BAD_CVC_LINE  2
#define PLC_E_BAD_CMP_NUM   3
#define PLC_E_BAD_CVRT_LINE 4
#define PLC_E_BAD_CLR_LINE  5
#define PLC_E_BAD_COLOR     6
#define PLC_E_BAD_VERT_LINE 7
#define PLC_E_BAD_CST_NUM   8
#define PLC_E_BAD_CST_LINE  9
#define PLC_E_BAD_CST_KIND  10
#define PLC_E_BAD_CST_NUMS  11
#define PLC_E_BAD_CC        12

#if (__cplusplus || c_plusplus)
};
#endif
#endif /* PLCURVE_H */
/*@=exportlocal@*/
