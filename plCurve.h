/*
 *
 * Data structures and prototypes for the plCurve library
 *
 *  $Id: plCurve.h,v 1.9 2006-02-05 04:18:54 ashted Exp $
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
int  plCurve_error_num;
char plCurve_error_str[80];

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
  int    cp;      /* Component */
  int    vt;      /* Vertex */
  int    kind;    /* What kind of constraint */
  double coef[6]; /* Coefficients for defining plane or line */
} plCurve_constraint;
  
typedef struct plCurve_type {	
  int nc;			/* Number of components */
  plCurve_pline *cp;            /* Components */
  int ncst;                     /* Number of constraints */
  plCurve_constraint *cst;      /* The constraints themselves */
} plCurve;

/*
 * Prototypes for vector routines. 
 *
 */

inline plcl_vector linklib_vplus(plcl_vector A,plcl_vector B);
inline plcl_vector linklib_vminus(plcl_vector A,plcl_vector B);
inline plcl_vector linklib_cross(plcl_vector A,plcl_vector B);
inline plcl_vector linklib_scalarmult(double x,plcl_vector A);
inline plcl_vector linklib_vdivide(plcl_vector A,plcl_vector B);
inline plcl_vector plcl_vector_random();

inline double  linklib_vdist(plcl_vector A,plcl_vector B);

inline void plcl_vector_normalize(plcl_vector *V);

/*
 * Macros for vector work (requires careful programming)
 *
 */

#define linklib_dot(A,B)    ((A).c[0]*(B).c[0] + (A).c[1]*(B).c[1] + (A).c[2]*(B).c[2])
#define linklib_norm(A)     sqrt(linklib_dot((A),(A)))
#define linklib_vadd(A,B)  \
  (A).c[0] += (B).c[0]; (A).c[1] += (B).c[1]; (A).c[2] += (B).c[2];
#define linklib_vsub(A,B)  \
  (A).c[0] -= (B).c[0]; (A).c[1] -= (B).c[1]; (A).c[2] -= (B).c[2];
#define linklib_vsmult(s,V) (V).c[0] *= s; (V).c[1] *= s; (V).c[2] *= s;
  /* Add a multiple of B to A */
#define linklib_vmadd(A,s,B)  (A).c[0] += (s)*(B).c[0]; \
                              (A).c[1] += (s)*(B).c[1]; \
                              (A).c[2] += (s)*(B).c[2];
  /* A = B + s(C-B)                               *
   * equivalent to                                *
   *   A = C; vsub(A,B); vsmult(s,A); vsadd(A,B); */
#define linklib_vweighted(A,s,B,C)  \
    (A).c[0] = (B).c[0] + s*((C).c[0] - (B).c[0]); \
    (A).c[1] = (B).c[1] + s*((C).c[1] - (B).c[1]); \
    (A).c[2] = (B).c[2] + s*((C).c[2] - (B).c[2]); 

#define linklib_vlincombine(A,s,B,t,C) \
    (C).c[0] = s*(A).c[0] + t*(B).c[0]; \
    (C).c[1] = s*(A).c[1] + t*(B).c[1]; \
    (C).c[2] = s*(A).c[2] + t*(B).c[2]; 

#define linklib_vdiv(A,B) \
    (A).c[0] /= (B).c[0]; \
    (A).c[1] /= (B).c[1]; \
    (A).c[2] /= (B).c[2];

#define linklib_vmul(A,B) \
    (A).c[0] *= (B).c[0]; \
    (A).c[1] *= (B).c[1]; \
    (A).c[2] *= (B).c[2];

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

/* Read link data from a file */
plCurve *plCurve_read(FILE *infile);

/* Write link data to a file */
int plCurve_write(FILE *outfile, const plCurve *L);

/* Fix the "hidden vertices" for easy handling of closed components */
void plCurve_fix_wrap(const plCurve *L);

/* Count the edges in a link (correctly handling open/closed) */
int plCurve_edges(const plCurve *L);

/* Compute the (minrad-based) curvature of L at vertex vt of component cp */
double plCurve_curvature(const plCurve *L, 
                              const int cp, 
                              const int vt);

/* Copy a link */
plCurve *plCurve_copy(const plCurve *L);

/* Compute tangent vector */
plcl_vector plCurve_tangent_vector(plCurve *link,int cp, int vt);

/* Find the arclength of a link. */
double plCurve_length(plCurve *L,double *component_lengths);

/* Find the arclength position of a point on a link. */
double plCurve_parameter(plCurve *L,int cmp,int vertnum);

/* Force a plCurve to close as gently as possible */
void plCurve_force_closed(plCurve *link);

/* Return a value for how far a constraint is from being satisfied */
double plCurve_cst_check(const plCurve L, const plCurve_constraint cst);

/* 
 * Force a constraint to be satisfied and return a value for how far from being
 * satisfied it was 
 */
double plCurve_cst_fix(const plCurve L, const plCurve_constraint cst);

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

#if (__cplusplus || c_plusplus)
};
#endif
#endif
