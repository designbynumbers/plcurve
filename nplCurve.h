/* nplCurve.h.  Generated by configure.  */
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
  
  /*  This header file includes prototypes for n-dimensional plCurves. 
      Not all of the functionality of the standard plCurve library has
      been ported, and we retain the original code for performance reasons. 
      
  */
  
  
  /*@-exportlocal@*/
#ifndef NPLCURVE_H
#define NPLCURVE_H
  
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
    
    typedef struct nplc_vector_type {
  double *c;
  int n;
} nplc_vector;
    
    /* There's nothing n-dimensional here, but we want to make
       nplCurve independent from plCurve, so we duplicate the
       functionality. Unfortunately, to have the option of using both
       together, this means either installing a separate plColor
       library and header, or changing the name and duplicating the
       functionality. We go with the latter. */

    typedef struct nplc_color_type {  
      
      double r;
      double g;
      double b;
      double alpha;
    
    } nplc_color;
    
    typedef struct nplc_strand_type {
      int            nv;     /* Number of vertices */
      bool           open;   /* This is an "open" strand (with distinct ends) */
      int            cc;     /* Color count (number of colors) */
      nplc_vector   *vt;     /* Actual vertices */
      nplc_color    *clr;    /* Colors */
    } nplc_strand;
    
    /* Curve constraint kind */
    typedef enum nplc_cst_kind_type {
      nunconstrained = 0,
      nfixed,
      nline,
      nplane,
    } nplc_cst_kind;

    typedef struct nplc_constraint_type {
      nplc_cst_kind kind;    /* What kind of constraint */
      nplc_vector   vect[2]; /* Vectors to define plane, line or fixed point */
      int           cmp;     /* Component */
      int           vert;    /* Starting vertex */
      int           num_verts; /* Length of run */
      /*@only@*/ /*@null@*/ struct nplc_constraint_type *next;
    } nplc_constraint;
    
    /* The vectors in vect change meaning depending on the type of constraint:

  fixed : vect[0] is the fixed point, vect[1] is ignored 
  line:   vect[0] is the tangent vector, vect[1] is a point on the line
  plane:  vect[0] is the normal vector, vect[1].x is the distance from origin.

  */

typedef struct nplc_vert_quant_type { /* Vertex quantifiers */
  int              cmp;    /* Component */
  int              vert;   /* Vertex */
  char             tag[4]; /* 3-character tag */
  double           quant;  /* Quantifier */
  /*@only@*/ /*@null@*/ struct nplc_vert_quant_type *next;
} nplc_vert_quant;

typedef struct nplc_type {
  int             nc;                           /* Number of components */
  nplc_strand *cp;                           /* Components */
  /*@only@*/ /*@null@*/ nplc_constraint *cst;   /* Constraints */
  /*@only@*/ /*@null@*/ nplc_vert_quant *quant; /* per-vertex quantifiers */
} nplCurve;

/* PlCurve_spline types */
typedef struct nplc_spline_strand_type {
  bool         open;     /* This is an "open" strand (with distinct ends) */
  int          ns;       /* Number of samples used to build spline. */
  double      *svals;    /* s values at samples */
  nplc_vector *vt;       /* positions at these s values */
  nplc_vector *vt2;      /* _second_ derivatives at these s values */
  int          cc;       /* Number of colors */
  nplc_color  *clr;      /* color values at samples */
} nplc_spline_strand;

typedef struct nplc_spline_type {
  int                 nc;     /* Number of components */
  nplc_spline_strand *cp;     /* Components */
} nplc_spline;
  
/*
 * Prototypes for vector routines.
 *
 */
  
nplc_vector nplc_vect_new(int n);
void nplc_vect_free(nplc_vector *nv);
void nplc_vect_copy(nplc_vector to,nplc_vector from);

nplc_vector *nplc_vect_buf_new(int dim,int num_vects);
void nplc_vect_buf_free(int bufsize,nplc_vector *buf);

nplc_vector nplc_vect_sum(nplc_vector A,nplc_vector B);   /* A + B */
nplc_vector nplc_vect_diff(nplc_vector A,nplc_vector B);  /* A - B */
  /* There is no general cross product operation for n-vectors */
nplc_vector nplc_scale_vect(double s,nplc_vector A);      /* sA */
nplc_vector nplc_normalize_vect(const nplc_vector V,
                                /*@null@*/ bool *ok);     /* V/|V| */

nplc_vector nplc_random_vect(int n);

/* Translate three doubles into a vector */
nplc_vector nplc_build_vect(int n,...);


/* Multiply or divide two ordered triples componetwise */
nplc_vector nplc_component_mult(nplc_vector A, nplc_vector B);
nplc_vector nplc_component_div(nplc_vector A, nplc_vector B,
                               /*@null@*/ bool *ok);

/* Return a linear combination: a*A + b*B */
nplc_vector nplc_vlincomb(double a,nplc_vector A, double b,nplc_vector B);
nplc_vector nplc_vmadd(nplc_vector A, double s, nplc_vector B); /* A + sB */
nplc_vector nplc_vweighted(double s, nplc_vector A, nplc_vector B); /* (1-s)A + sB */
  
/* Different vector measurements */
double nplc_dot_prod(nplc_vector A,nplc_vector B);
double nplc_norm(nplc_vector A);
double nplc_distance(nplc_vector A, nplc_vector B);
double nplc_sq_dist(nplc_vector A, nplc_vector B);
double nplc_angle(nplc_vector A, nplc_vector B, bool *ok);

/* Do two vectors match ? */
bool nplc_vecteq(nplc_vector A, nplc_vector B);

  /* We can't have macros for the various operations, since 
     we don't know the dimension of the vectors. */

char *nplc_vect_print(nplc_vector A);

    /* Returns a pointer to a static string printing the vector. */

/*
 * Prototypes for routines to deal with nplCurves. INCOMPLETE!
 *
 */

int nplc_dim(nplCurve *L);
    /* Returns the dimension of the curve. */

/* Build a new nplCurve (with associated strands) */
/*@only@*/ nplCurve *nplc_new(const int dim,
			      const int components,
			      const int * const nv,
			      const bool * const open,
			      const int * const cc);

/* Free the nplCurve (and strands) */
void nplc_free(/*@only@*/ /*@null@*/ nplCurve *L);


#ifdef CONVERTED

/* Add a component to the curve which will become component number add_as. */
void nplc_add_component(nplCurve *L, const int add_as, const int nv, 
                       const bool open, const int cc,
                       const nplc_vector * const vt,
            /*@null@*/ const nplc_color  * const clr);

/* And remove one */
void nplc_drop_component(nplCurve *L, const int cmp);

/* Set a constraint on a vertex or run of vertices */
void nplc_set_fixed(nplCurve * const L,
                   const int cmp,
                   const int vert,
                   const nplc_vector point);

void nplc_constrain_to_line(nplCurve * const L,
                           const int cmp,
                           const int vert,
                           const int num_verts,
                           const nplc_vector tangent,
                           const nplc_vector point_on_line);

void nplc_constrain_to_plane(nplCurve * const L,
                            const int cmp,
                            const int vert,
                            const int num_verts,
                            const nplc_vector normal,
                            const double dist_from_origin);

void nplc_unconstrain(nplCurve * const L, const int cmp,
                     const int vert, const int num_verts);

/* Remove a constraint from the list of constraints returning the number of
 * vertices thus set nunconstrained.  */
int nplc_remove_constraint(nplCurve * const L,
                          const nplc_cst_kind kind,
                          const nplc_vector vect[]);

/* Remove all constraints */
void nplc_remove_all_constraints(nplCurve * const L);

#endif

/* Read nplCurve data from a file */
/*@only@*/ /*@null@*/ nplCurve *nplc_read(FILE *file,
                              /*@out@*/ int *error_num,
                              /*@out@*/ char error_str[],
                                        size_t error_str_len);

/* Write nplCurve data to a file */
void nplc_write(FILE *outfile, nplCurve * const L);

/* Fix the "hidden vertices" for easy handling of closed components */
void nplc_fix_wrap(nplCurve * const L);

/* Count the edges in a nplCurve (correctly handling open/closed) */
/* Deprecated in versions > 1.3 in favor of nplc_edges call below. */
int nplc_num_edges(const nplCurve * const L);

/* Count edges in nplCurve, returning total and storing #edges for */
/* each component in component_edges if this is non-NULL. */
int nplc_edges(const nplCurve * const L,
/*@null@*/ /*@out@*/ int *component_edges);

/* Count the vertices in a nplCurve */
int nplc_num_verts(const nplCurve * const L);

#ifdef CONVERTED

/* Compute the MinRad-based curvature of L at vertex vt of component cp */
double nplc_MR_curvature(nplCurve * const L, const int cmp, const int vert);

/* Copy a nplCurve */
nplCurve *nplc_copy(const nplCurve * const L);

/* Compute average of inward and outward tangents (and normalize) */
nplc_vector nplc_mean_tangent(const nplCurve * const L, const int cmp,
                            const int vert, bool *ok);

/* Find the arclength of a nplCurve. Total arclength is returned, arclength */
/* of individual strands stored in component_lengths if this pointer is non-NULL. */
double nplc_arclength(const nplCurve * const L,
/*@null@*/ /*@out@*/ double *component_lengths);

/* Find the arclength distance from one vertex to another.  On closed
 * strands, give the shortest of the two options.  */
double nplc_subarc_length(const nplCurve * const L, const int cmp,
                         const int vert1, const int vert2);

/* Return how far a constraint is from being satisfied (sup norm). */
double nplc_check_cst(const nplCurve * const L);

/* Fix all the vertices which are out of compliance with their constraints. */
void nplc_fix_cst(nplCurve * const L);

/* Either return (if given a char *) or print out the library version number */
void nplc_version(/*@null@*/ char *version, size_t strlen);

/* Put 4 doubles together into a color */
nplc_color plc_build_color(const double r, const double g,
			   const double b, const double alpha);

/* Allocate new spline. */
nplc_spline *nplc_spline_new(const int          components,
                           const int  * const ns,
                           const bool * const open,
                           const int  * const cc);

/* Free memory for spline. */
void nplc_spline_free(/*@only@*/ /*@null@*/ nplc_spline *L);

/* Convert nplCurve to spline representation. */
nplc_spline *nplc_convert_to_spline(nplCurve * const L, bool *ok);

/* Convert splined curve to nplCurve (with resampling). */
nplCurve *nplc_convert_from_spline(const nplc_spline * const spL,
                                 const int * const nv);

/* Samples a spline at a particular s value. */
nplc_vector nplc_sample_spline(const nplc_spline * const spL,
                             const int cmp,
                             double s);

/* Calculate the diameter of the nplCurve, thinking of the vertices as 
   a set of points in R^3 */
double nplc_pointset_diameter(const nplCurve * const L);

/* Scale a nplCurve (and its' constraints!) by a factor. */
void nplc_scale( nplCurve *link, const double alpha);  

#endif

/* Define the error codes */
#define NPLC_E_NO_VECT       1
#define NPLC_E_BAD_CVC_LINE  2
#define NPLC_E_BAD_CMP_NUM   3
#define NPLC_E_BAD_CVRT_LINE 4
#define NPLC_E_BAD_CLR_LINE  5
#define NPLC_E_BAD_COLOR     6
#define NPLC_E_BAD_VERT_LINE 7
#define NPLC_E_BAD_CST_NUM   8
#define NPLC_E_BAD_CST_LINE  9
#define NPLC_E_BAD_CST_KIND  10
#define NPLC_E_BAD_CST_NUMS  11
#if (__cplusplus || c_plusplus)
};
#endif
#endif /* NPLCURVE_H */
/*@=exportlocal@*/
