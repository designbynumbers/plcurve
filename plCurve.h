/*
 *
 * Data structures and prototypes for the plCurve library
 *
 *  $Id: plCurve.h,v 1.17 2006-02-12 19:37:01 ashted Exp $
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

#ifndef __PLCURVE_H
#define __PLCURVE_H

#if (__cplusplus || c_plusplus)
extern "C" {
#endif

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* We need to define FILE */
#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif 

/* Define TRUE and FALSE */
#ifndef FALSE
#define FALSE (1 == 0)
#endif /* FALSE */
#ifndef TRUE
#define TRUE (1 == 1)
#endif /* TRUE */

/* Variables for reporting errors */
int  plcl_error_num;
char plcl_error_str[80];

/* Define 3-space vectors */
typedef struct plcl_vector_type {  
  double c[3];
  int cst; /* Which constraint, if any, this vertex satisfies */
} plcl_vector;

typedef struct plCurve_color_type {
  double r;
  double g;
  double b;
  double alpha;
} plCurve_color;

typedef struct plCurve_pline_type {
  int            open;   /* This is an "open" pline (with distinct ends) */
  int            nv;     /* Number of vertices */
  int            cc;     /* Color count (number of colors) */
  plcl_vector   *vt;     /* Actual vertices */
  plCurve_color *clr;    /* Colors */
} plCurve_pline;

/* Curve constraint kind */
#define PLCL_FIXED     1  /* Vertex is not allowed to move */
#define PLCL_ON_LINE   2  /* Vertex must lie on the given line */
#define PLCL_IN_PLANE  3  /* Vertex must lie in the given plane */

typedef struct plCurve_constraint_type {
  int    kind;    /* What kind of constraint */
  double coef[6]; /* Coefficients for defining plane or line */
} plCurve_constraint;
  
typedef struct plCurve_type {   
  int nc;                       /* Number of components */
  plCurve_pline *cp;            /* Components */
  int ncst;                     /* Number of constraints */
  plCurve_constraint *cst;      /* The constraints themselves */
} plCurve;

/* PlCurve_spline types */
typedef struct linklib_spline_pline_type {
  int             open;     /* This is an "open" pline (with distinct ends) */
  int             ns;       /* Number of samples used to build spline. */
  double         *svals;    /* s values at samples */
  plcl_vector    *vt;       /* positions at these s values */
  plcl_vector    *vt2;      /* second derivatives at these s values */
  int             cc;
  plCurve_color  *clr;      /* color values at samples */
  plCurve_color  *clr2;     /* second derivatives at these s values */
} linklib_spline_pline;

typedef struct plCurve_spline_type {    
  int nc;                       /* Number of components */
  linklib_spline_pline *cp;     /* Components */
} plCurve_spline;

/*
 * Prototypes for vector routines. 
 *
 */

inline plcl_vector plcl_vect_sum(plcl_vector A,plcl_vector B);   /* A + B */
inline plcl_vector plcl_vect_diff(plcl_vector A,plcl_vector B);  /* A - B */
inline plcl_vector plcl_cross_prod(plcl_vector A,plcl_vector B); /* A x B */
inline plcl_vector plcl_scale_vect(double s,plcl_vector A);      /* sA */
inline plcl_vector plcl_normalize_vect(const plcl_vector V);     /* V / |V| */
inline plcl_vector plcl_random_vect();                         

/* Multiply or divide two ordered triples componetwise */
inline plcl_vector plcl_component_mult(plcl_vector A,plcl_vector B);
inline plcl_vector plcl_component_div(plcl_vector A,plcl_vector B);

/* Return a linear combination: a*A + b*B */
inline plcl_vector plcl_vlincomb(double a,plcl_vector A,
                                 double b,plcl_vector B);

inline double plcl_dot_prod(plcl_vector A,plcl_vector B);
inline double plcl_norm(plcl_vector A);

/*
 * Macros for vector work (requires careful programming)
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
#define linklib_vmadd(A,s,B) \
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

#define plCurve_M_set_vertex(L,cmp,vertex,x,y,z) \
  L->cp[cmp].vt[vertex].c[0] = x; \
  L->cp[cmp].vt[vertex].c[1] = y; \
  L->cp[cmp].vt[vertex].c[2] = z;

/* 
 * Prototypes for routines to deal with plCurves.
 *
 */

/* Build a new link (with associated plines) */
plCurve *plCurve_new(int components, const int *nv, 
                     const int *open, const int *cc, 
                     const int ncst, const plCurve_constraint *cst);

/* Free the link (and plines) */
void plCurve_free(plCurve *L);

/* Set a vertex */
inline void plCurve_set_vertex(plCurve *L, const int cmp, const int vertex,
                               const double x, const double y, const double z);

/* Set a constraint on a vertex or run of vertices */
inline int plCurve_set_constraint(plCurve *L, const int cmp, 
                                  const int vertex, const int num_verts, 
                                  const int kind, const int coef1,
                                  const int coef2, const int coef3,
                                  const int coef4, const int coef5,
                                  const int coef6);

/* Set a vertex or run of vertices to unconstrained */
inline void plCurve_set_unconstrained(const plCurve *L, const int cmp,
                                      const int vertex, const int num_verts);

/* Read link data from a file */
plCurve *plCurve_read(FILE *infile);

/* Write link data to a file */
int plCurve_write(FILE *outfile, const plCurve *L);

/* Fix the "hidden vertices" for easy handling of closed components */
void plCurve_fix_wrap(const plCurve *L);

/* Count the edges in a link (correctly handling open/closed) */
int plCurve_num_edges(const plCurve *L);

/* Compute the (minrad-based) curvature of L at vertex vt of component cp */
double plCurve_curvature(const plCurve *L, const int cp, const int vt);

/* Copy a link */
plCurve *plCurve_copy(const plCurve *L);

/* Compute tangent vector */
plcl_vector plCurve_tangent_vector(plCurve *link,int cp, int vt);

/* Find the arclength of a link. */
double plCurve_arclength(plCurve *L,double *component_lengths);

/* Find the arclength position of a point on a link. */
double plCurve_parameter(plCurve *L,int cmp,int vertex);

/* Force a plCurve to close as gently as possible */
void plCurve_force_closed(plCurve *link);

/* Return a value for how far a constraint is from being satisfied */
double plCurve_cst_check(const plCurve *L, const int cmp, const int vertex);

/* 
 * Force a constraint to be satisfied and return a value for how far from being
 * satisfied it was 
 */
double plCurve_cst_fix(const plCurve *L, const int cmp, const int vertex);

/* Check plcl_error_num, report on nonzero, terminate on positive */
inline void plcl_status_check();

/* Either return (if given a char *) or print out the library version number */
inline void plCurve_version(char *version);

/* Allocate new spline_link. */
plCurve_spline *linklib_spline_link_new(int components, 
                                               const int *ns, 
                                               const int *open, 
                                               const int *cc);

/* Free memory for spline_link. */
void linklib_spline_link_free(plCurve_spline *L);

/* Convert conventional link to spline_link. */
plCurve_spline *convert_to_spline_link(plCurve *L);

/* Convert spline_link to conventional link (with resampling). */
plCurve *convert_spline_to_link(plCurve_spline *spL,int *nv);

/* Evaluate a spline_link at a particular s value. */
plcl_vector evaluate_spline_link(plCurve_spline *spL,int cmp,double s);

/* Define the error codes */
#define PLCL_E_BAD_RANDOM     1
#define PLCL_E_TOO_FEW_VERTS  2
#define PLCL_E_CANT_ALLOC     3
#define PLCL_E_TOO_FEW_COMPS  4
#define PLCL_E_NULL_PTR       5
#define PLCL_E_NEG_CST        6
#define PLCL_E_BAD_CST_KIND   7
#define PLCL_E_TOO_FEW_DBLS   8
#define PLCL_E_TOO_FEW_INTS   9
#define PLCL_E_NO_VECT       10
#define PLCL_E_BAD_CVC_LINE  11
#define PLCL_E_BAD_CVRT_LINE 12
#define PLCL_E_BAD_CLR_LINE  13
#define PLCL_E_NEW_FAILED    14
#define PLCL_E_BAD_COLOR     15
#define PLCL_E_BAD_VERT_LINE 16
#define PLCL_E_INF_KAPPA     17
#define PLCL_E_BAD_COMPONENT 18
#define PLCL_E_BAD_VERTEX    19
#define PLCL_E_ZERO_VECTOR   20
#define PLCL_E_TOO_FEW_SAMPS 21
#define PLCL_E_SMP_TOO_CLOSE 22
#define PLCL_E_CANT_FIND_POS 23
#define PLCL_E_TOO_MANY_VRTS 24

#if (__cplusplus || c_plusplus)
};
#endif
#endif
