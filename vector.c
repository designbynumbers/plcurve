/*
 * @COPYRIGHT@
 *
 * Routines for working with vectors.
 *
 * $Id: vector.c,v 1.29 2006-03-02 05:32:57 ashted Exp $
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

#ifdef HAVE_MATH_H
  #include <math.h>
#endif
#ifdef HAVE_STRING_H
  #include <string.h>
#endif
#ifdef HAVE_STDLIB_H
  #include <stdlib.h>
#endif
#ifdef HAVE_FLOAT_H
  #include <float.h>
#endif
#ifdef HAVE_ASSERT_H
  #include <assert.h>
#endif


/**************************************************************/
/*   Basic Linear Algebra Operations                          */
/**************************************************************/

/* Returns A + B. */
inline plcl_vector plcl_vect_sum(plcl_vector A,plcl_vector B) {
  plcl_M_add_vect(A,B);
  return A;
}

/* Returns A - B. */
inline plcl_vector plcl_vect_diff(plcl_vector A,plcl_vector B) {
  plcl_M_sub_vect(A,B);
  return A;
}

/* Returns A x B. */
inline plcl_vector plcl_cross_prod(plcl_vector A,plcl_vector B) {
  plcl_vector C;
  plcl_M_cross(C,A,B);
  return C;
}

/* Returns sA. */
inline plcl_vector plcl_scale_vect(double s,plcl_vector A) {
  plcl_M_scale_vect(s,A);
  return A;
}

inline plcl_vector plcl_component_mult(plcl_vector A,plcl_vector B) {
  plcl_M_component_mult(A,B);
  return A;
}

/* Should we add an "ok" parameter here as in _normalize_vect? */
inline plcl_vector plcl_component_div(plcl_vector A,plcl_vector B,
                                      /*@null@*/ bool *ok) {
  if (ok != NULL) {
    *ok = ((fabs(B.c[0]) > DBL_EPSILON) &&
           (fabs(B.c[1]) > DBL_EPSILON) &&
           (fabs(B.c[2]) > DBL_EPSILON));
    if (*ok) {
      plcl_M_component_div(A,B);
    }
  } else {
    if ((fabs(B.c[0]) > DBL_EPSILON) && (fabs(B.c[1]) > DBL_EPSILON) &&
      (fabs(B.c[2]) > DBL_EPSILON)) {
      plcl_M_component_div(A,B);
    } else {
      fprintf(stderr,"plcl_component_div: Divisor has zero component.\n");
      exit(EXIT_FAILURE);
    }
  }
  return A;
}

/* Returns the dot product of A and B */
inline double plcl_dot_prod(plcl_vector A, plcl_vector B)
/*@modifies nothing@*/ {
  return plcl_M_dot(A,B);
}

inline double plcl_norm(plcl_vector A) {
  return plcl_M_norm(A);
}

inline double plcl_distance(plcl_vector A, plcl_vector B) {
  return plcl_M_norm(plcl_vect_diff(A,B));
}

/* The square of the distance between A and B (faster than _distance) */
inline double plcl_sq_dist(plcl_vector A, plcl_vector B) {
  return plcl_M_sq_dist(A,B);
}

inline bool plcl_vecteq(plcl_vector A, plcl_vector B) /*@modifies nothing@*/ {
  return plcl_M_vecteq(A,B);
}

/* Procedure returns a vector which points in the same direction as V but has
 * length 1.  It sets *ok to false if the norm is too small. */
inline plcl_vector plcl_normalize_vect(const plcl_vector V,
                                       /*@null@*/ bool *ok) {
  double vnrm;

  vnrm = plcl_M_norm(V);
  if (vnrm < DBL_EPSILON && -vnrm < DBL_EPSILON) {
    if (ok != NULL) {
      *ok = false;
    } else {
      fprintf(stderr,
        "plcl_normalize_vect: Attempted to normalize zero vector.\n");
      exit(EXIT_FAILURE);
    }
  }
  return plcl_scale_vect(1.0/vnrm,V);
}

/*
 * George Masaglia's "new method" for finding a random point on a 3-sphere,
 * from his short article "Choosing a Point from the Surface of a Sphere"
 * in The Annals of Mathematical Statistics, v. 43, no. 2 (Apr, 1972) 645-646.
 *
 */
plcl_vector plcl_random_vect()
{
  int i;
  plcl_vector R;
  double V1 = 0.0, V2 = 0.0;
  double S = 0.0;
  double sqt;

  /*@+loopexec@*/
  for (i = 0; i < 1000 &&
              (S - 1.0 > DBL_EPSILON ||
               S - 0.01 < DBL_EPSILON); i++) {
    V1 = 2*(double)(rand())/(double)(RAND_MAX) - 1;
    V2 = 2*(double)(rand())/(double)(RAND_MAX) - 1;
    S = V1*V1 + V2*V2;
  }
  /*@=loopexec@*/
  assert(S - 0.01 >= DBL_EPSILON && S - 1.0 <= DBL_EPSILON);
  sqt = sqrt(1-S);
  R.c[0] = 2*V1*sqt;
  R.c[1] = 2*V2*sqt;
  R.c[2] = 1-2*S;

  return R;
}

/* Return a linear combination: a*A + b*B */
inline plcl_vector plcl_vlincomb(double a,plcl_vector A,
                                 double b,plcl_vector B) {
  plcl_vector R;

  plcl_M_vlincomb(R,a,A,b,B);
  return R;
}

inline plcl_vector plcl_vmadd(plcl_vector A, double s, plcl_vector B) {
  plcl_M_vmadd(A,s,B);
  return A;
}

inline plcl_vector plcl_vweighted(double s, plcl_vector A, plcl_vector B) {
  plcl_vector R;
  plcl_M_vweighted(R,s,A,B);
  return R;
}

/* Put together a vector from 3 doubles */
inline plcl_vector plcl_build_vect(const double x,
                                   const double y,
                                   const double z) /*@modifies nothing@*/ {
  plcl_vector V = { { x, y, z } };
  return V;
}
