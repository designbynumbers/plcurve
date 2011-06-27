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

/* We now introduce a data type encoding a symmetry of a link. Such a symmetry has to include
   both a geometric transformation of space AND a corresponding map from each vertex of the 
   plcurve to a target vertex. */

typedef double (plc_matrix)[3][3];

struct plc_vertex_loc {

  int cp;
  int vt;

};

typedef struct plc_type plCurve; /* We need to forward declare the plCurve type. */

typedef struct plc_symmetry_type {

  plc_matrix              *transform;
  plCurve                 *curve;
  struct plc_vertex_loc  **target; /* Array of curve->nc arrays of curve->cp[cp].nv arrays of plc_vertex_loc */ 

} plc_symmetry;

typedef struct plc_symmetry_group_type { /* This requires a little bit of the group structure to be specified. */

  int n;
  plc_symmetry **sym; /* This is just a list of the group elements. */
  int *inverse; /* the (index of) the inverse of symmetry i is given by G->inverse[i] */

} plc_symmetry_group;


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

struct plc_type {
  int         nc;                              /* Number of components */
  plc_strand *cp;                              /* Components */
  /*@only@*/ /*@null@*/ plc_constraint *cst;   /* Constraints */
  /*@only@*/ /*@null@*/ plc_vert_quant *quant; /* per-vertex quantifiers */
  plc_symmetry_group *G;                       /* Symmetry group (may be null) */
};

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
  plc_spline_strand  *cp;     /* Components */
  plc_constraint     *cst;    /* Constraints */
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

/* Returns a vector randomly distributed on the surface of the unit sphere */
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

#define plc_M_vect_diff(A,B,C)						\
  (A).c[0] = (B).c[0] - (C).c[0];  (A).c[1] = (B).c[1] - (C).c[1];  (A).c[2] = (B).c[2] - (C).c[2];

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
/* The arrays nv, open, and cc are expected to be of length components. */

/* nv[i] is the number of vertices of the ith component */
/* open[i] is true if the ith component is open */
/* cc[i] is the number of colors for that component (0,1, or nv[i]) */

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

/* Set the color of a plCurve to a single color */
void plc_set_color(plCurve *L, const plc_color inColor);

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

/* Test whether a vertex is constrained. If constraint is non-null, set it to the active constraint. */
bool plc_is_constrained(plCurve * const L,int cp, int vt,plc_constraint **constraint);

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

/* Compute an index between 0 and plc_num_verts(L) - 1 for a (cp,vt)
   pair in a plCurve, using full wraparound addressing for closed
   components and repeating the last or first vertex for open ones. We
   guarantee that these numbers occur consecutively in dictionary
   order on the pairs (cp,vt). */

int plc_vertex_num(const plCurve * const L, const int cp, const int vt);

/* Convert back from a vertex number given by plc_vertex_num to a (cp,vt) pair. */
int plc_cp_num(const plCurve * const L, int wrapVt);
int plc_vt_num(const plCurve * const L, int wrapVt);

/* Compute the turning angle at a vertex. Uses wraparound addressing if needed. */
double plc_turning_angle(plCurve * const L, const int cmp, const int vert, 
			 bool *ok);

/* Compute the MinRad-based curvature of L at vertex vt of component cp */
double plc_MR_curvature(plCurve * const L, const int cmp, const int vert);

/* Find total curvature for plCurve, including value for each component */
/* if component_tc is non-null. */
double plc_totalcurvature(const plCurve * const L,
/*@null@*/ /*@out@*/ double *component_tc);

/* Copy a plCurve */
plCurve *plc_copy(const plCurve * const L);

/*
 * Compute a (unit) tangent vector to L at vertex vert of component cmp by
 * taking the incoming tangent and outgoing tangent and averaging them *with
 * their lengths taken into account*.  That is, the average of (0.0,0.0,6.0)
 * and (8.0,0.0,0.0) is (4.0,0.0,3.0) which normalizes to (0.8,0.0,0.6)
 *
 */
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

/* Find the arclength position of a vertex on the a plCurve. */
/* On a multicomponent curve, s values add from 0 (0th vert, component 0) */
/* to the total arclength of the curve (last vert, last component) */

double plc_s(const plCurve * const L, const int cmp, const int vert);


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

/* Calculate the center of mass of the plCurve, thinking of the vertices as
     equal mass points. */
  plc_vector plc_center_of_mass(const plCurve * const L);
  
/* Find the nearest vertex on a plCurve to a query point. Precomputes some
   data which can be saved to speed up future queries on the same plCurve if pc_data is non-NULL. 
   
   If we are passed a non-null pointer in pc_data, we check *pc_data. If it is NULL, we assume this 
   is the first call on this curve, precompute data, and return a pointer to that data in *pc_data.
   If *pc_data is non-NULL, we assume that this is precomputed data from a previous call for the 
   same plCurve and use it. (If this is NOT true, hilarity is likely to ensue.)
   
   If precomputed data is saved, it is the callers responsibility to free the pointer later using 
   plc_nearest_vertex_pc_data_free. The data structure is not exposed to the user because it may 
   change without warning in future versions of plCurve. The precomputed data stores the pointer
   *plCurve, so it will not work on a copy of the curve. 

   Returns both the space location of the vertex and its cp, vt index. */

  struct plc_nearest_neighbor_pc_data { /* THIS MAY CHANGE WITHOUT WARNING IN FUTURE VERSIONS! */
    
    plc_vector *check_buffer;
    int search_dimension;
    int *sorted_buffer;

  };

  struct plc_nearest_vertex_pc_data { /* THIS MAY CHANGE WITHOUT WARNING IN FUTURE VERSIONS! */

    plCurve *check_curve; /* Note that we will LOSE memory if this pointer goes bad before free. */
    struct plc_nearest_neighbor_pc_data **component_data;

  };
  
  plc_vector plc_nearest_vertex(const plc_vector pt,plCurve *L,int *cp, int *vt, 
				struct plc_nearest_vertex_pc_data **pc_data, int *plc_error);

  void plc_nearest_vertex_pc_data_free(struct plc_nearest_vertex_pc_data **pc_data);
  
  /* Finds the nearest point in a buffer of plc_vectors to a query point. Precomputes
     some data which can be saved to speed up future queries on the same buffer if pc_data is non-NULL. 
     
     If we are passed a non-null pointer in pc_data, we check *pc_data. If it is NULL, we assume this 
     is the first call on this curve, precompute data, and return a pointer to that data in *pc_data.
     If *pc_data is non-NULL, we assume that this is precomputed data from a previous call for the 
     SAME buffer and use it. (If this is NOT true, hilarity is likely to ensue.)
     
     Returns the index of the closest point in the query buffer. */

  int plc_nearest_neighbor(const plc_vector pt,const int n, plc_vector *buffer,
			   struct plc_nearest_neighbor_pc_data **pc_data, int *plc_error);
  void plc_nearest_neighbor_pc_data_free(struct plc_nearest_neighbor_pc_data **pc_data);


/* Scale a plCurve (and its' constraints!) by a factor. */
void plc_scale( plCurve *L, const double alpha);  

/* Perform a "whitten group" operation on L, mirroring, reversing and 
   permuting components of L. The syntax is 
   
   mirror = +1 or -1, with -1 to mirror entire link over xy plane.
   eps    = array of L->nc integers, each +1 or -1, with -1 to reverse 
   perm   = array of 2*L->nc integers so that if i is the first instance
            of index j in the array, then perm[i+1] = p(j).

   In the output link, the components are in the order

   eps[0] K_p(0), ... , eps[L->nc-1] K_p(L->nc-1)

   The function operates in place on L and respects constraints and colors,
   but not quantifiers. */
              
void plc_whitten(plCurve *L, int mirror, int *eps, int *perm);      

/* Perform a ``fold'' move on a plCurve */
void plc_pfm( plCurve *L, int cp, int vt0, int vt1, double angle);

/* Rotate a plCurve around an axis. */
void plc_rotate( plCurve *L, plc_vector axis, double angle);

/* Rotate a plCurve so the given axis points in the direction (0,0,1). */
void plc_random_rotate(plCurve *link, plc_vector axis);

/* Perform a random perturbation on a plCurve. Does not perturb
   constrained vertices. */
void plc_perturb( plCurve *L, double radius); 

/****************************** plCurve Symmetry Functions ********************/

void plc_identity_matrix(plc_matrix *A);
void plc_rotation_matrix(plc_vector axis, double angle,plc_matrix *A);
void plc_reflection_matrix(plc_vector axis,plc_matrix *A);

/* We now define a high level interface for dealing with symmetries. */

plc_symmetry *plc_symmetry_new(plCurve *model);
void plc_symmetry_free(plc_symmetry **A);
/* Make a new-memory copy of A */
plc_symmetry *plc_symmetry_copy(plc_symmetry *A);

/* This creates a plc_symmetry from a transform by searching to try to figure
   out the "intended" target of each vertex under the transform A. */
plc_symmetry *plc_build_symmetry(plc_matrix *A,plCurve *L);

/* This is a combination of matrix multiplication and applying the permutation
   of vertices in the symmetries to build a new symmetry (matrix product BA). 
   Returns NULL on fail. */
plc_symmetry *plc_compose_symmetries(plc_symmetry *A,plc_symmetry *B);

/* We now need constructors and destructors for the symmetry group */
plc_symmetry_group *plc_symmetry_group_new(int n);
void plc_symmetry_group_free(plc_symmetry_group **G);
plc_symmetry_group *plc_symmetry_group_copy(plc_symmetry_group *G);

/* We define a couple of standard groups as well. Return NULL if the build fails. */
/* Remember that the curves have to basically have the desired symmetry to start. */
plc_symmetry_group *plc_rotation_group(plCurve *L,plc_vector axis, int n);
plc_symmetry_group *plc_reflection_group(plCurve *L,plc_vector axis);

/* This symmetrizes a plCurve over the group L->G. */
void plc_symmetrize(plCurve *L);

  /* This symmetrizes a variation (a buffer of vectors of length plc_num_verts), assumed to 
     represent vectors located at the vertices of L over the symmetry group L->G. */

void plc_symmetrize_variation(plCurve *L,plc_vector *buffer);

/* Checks the distance between the position of each vertex and it's target after the 
   symmetry transform and returns the maximum. This serves as a check on the quality 
   of a symmetry possessed by a curve. */

double plc_symmetry_check(plCurve *L,plc_symmetry *A);

  /* To check an entire group, use plc_symmetry_group_check, which checks the entire 
     group L->G and returns the maximum error. */

double plc_symmetry_group_check(plCurve *L);


  /************************ plCurve Topology Library ********************/

/* This contains some functionality designed to work with plCurves as knots,
   including converting them to an abstract ``crossing'' representation, 
   computing their HOMFLY polynomials (using lmpoly) and identifying their
   knot types (by HOMFLY). */

/* The Millett/Ewing representation of a knot diagram numbers the
   crossings from 1 to ncrossings and then stores for each crossing
   the crossing connected to each arc coming from the crossing in the
   order

        a
	|
	|
    b---|-->d
        |
	V
        c

   So a crossing code representation of a plCurve is a char buffer 
   containing lines of the form

   17+2b10c11c31a

   meaning that crossing 17 is a positive crossing 

   connected in the a position to the b position of crossing 2,
   connected in the b position to the c position of crossing 10,
   connected in the c position to the c position of crossing 11 and
   connected in the d position to the a position of crossing 31.

   In order to simplify communication with the lmpoly code of Ewing
   and Millett, we store the crossing code as a standard (0
   terminated) string, including newlines. We will read from 
   that string using a replacement version of the "read" primitive.

*/

#define MAXPRIMEFACTORS 10
#define MAXHOMFLY       1024

typedef struct knottypestruct {

    int  nf;                            /* Number of prime factors */
    int  cr[MAXPRIMEFACTORS];           /* Crossing number of each prime factor */
    int  ind[MAXPRIMEFACTORS];           /* Index (in Rolfsen or Cerf) of each prime factor */
    char sym[MAXPRIMEFACTORS][128];     /* Symmetry tag (Whitten group element) for each prime factor */
    char homfly[MAXHOMFLY];             /* Homfly polynomial (as plc_lmpoly output) */

} plc_knottype;

/* Convert a plCurve to Millett/Ewing crossing code. */
char *plc_ccode( plCurve *L);

/* Compute the HOMFLY polynomial of a plCurve (returned as string) */
char *plc_homfly( plCurve *L);

/* Find the knot type of a plCurve */
/* Sets nposs to the number of possible knottypes found for the curve. If we cannot
   classify the knot, return 0 for nposs and NULL for the buffer of knot types. */
plc_knottype *plc_classify( plCurve *L, int *nposs);
 
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
#define PLC_E_STALE_PCDATA  13

#if (__cplusplus || c_plusplus)
};
#endif
#endif /* PLCURVE_H */
/*@=exportlocal@*/
