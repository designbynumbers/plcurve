/*
 * @COPYRIGHT@
 *
 * Routines for working with vectors.
 *
 * $Id: vector.c,v 1.37 2007-07-12 15:27:52 cantarel Exp $
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
inline plc_vector plc_vect_sum(plc_vector A,plc_vector B) {
  plc_M_add_vect(A,B);
  return A;
}

/* Returns A - B. */
inline plc_vector plc_vect_diff(plc_vector A,plc_vector B) {
  plc_M_sub_vect(A,B);
  return A;
}

/* Returns A x B. */
inline plc_vector plc_cross_prod(plc_vector A,plc_vector B) {
  plc_vector C;
  plc_M_cross(C,A,B);
  return C;
}

/* Returns sA. */
inline plc_vector plc_scale_vect(double s,plc_vector A) {
  plc_M_scale_vect(s,A);
  return A;
}

inline plc_vector plc_component_mult(plc_vector A,plc_vector B) {
  plc_M_component_mult(A,B);
  return A;
}

/* Should we add an "ok" parameter here as in _normalize_vect? */
inline plc_vector plc_component_div(plc_vector A,plc_vector B,
                                      /*@null@*/ bool *ok) {

  if (ok != NULL) {
    *ok = ((fabs(B.c[0]) > DBL_EPSILON) &&
           (fabs(B.c[1]) > DBL_EPSILON) &&
           (fabs(B.c[2]) > DBL_EPSILON));
    if (*ok) {
      plc_M_component_div(A,B);
    }
  } else {
    if ((fabs(B.c[0]) > DBL_EPSILON) && (fabs(B.c[1]) > DBL_EPSILON) &&
      (fabs(B.c[2]) > DBL_EPSILON)) {
      plc_M_component_div(A,B);
    } else {
      fprintf(stderr,"plc_component_div: Divisor has zero component.\n");
      exit(EXIT_FAILURE);
    }
  }
  return A;
}

/* Returns the dot product of A and B */
inline double plc_dot_prod(plc_vector A, plc_vector B)
/*@modifies nothing@*/ {
  return plc_M_dot(A,B);
}

inline double plc_norm(plc_vector A) {
  return plc_M_norm(A);
}

inline double plc_distance(plc_vector A, plc_vector B) {
  return plc_M_norm(plc_vect_diff(A,B));
}

/* The square of the distance between A and B (faster than _distance) */
inline double plc_sq_dist(plc_vector A, plc_vector B) {
  return plc_M_sq_dist(A,B);
}

#define plc_M_sqr(A) \
   ((A)*(A))

/* Computes the angle between two vectors. Can fail if one or the other 
   has norm zero. */
double plc_angle(plc_vector A, plc_vector B, bool *ok) 
 
  /* We use the algorithm suggested by Schatte. 
     
  @article{312261,
           author = {Peter Schatte},
	   title = {Computing the angle between vectors},
           journal = {Computing},
           volume = {63},
           number = {1},
           year = {1999},
           issn = {0010-485X},
           pages = {93--96},
           doi = {http://dx.doi.org/10.1007/s006070050052},
           publisher = {Springer-Verlag New York, Inc.},
	   address = {New York, NY, USA}
	   }

      which has better numerical stability than the standard algorithm,
      especially when dealing with angles that are very small or large. */

  {
    double L = 0.0; /* L = |A|^2|B|^2 - (A.B)^2 = \sum_{i<j} (A_i B_j - A_j B_i)^2 */
    double D = 0.0; /* D = A.B */
    double angle;

    L += plc_M_sqr(A.c[0]*B.c[1] - A.c[1]*B.c[0]);
    L += plc_M_sqr(A.c[0]*B.c[2] - A.c[2]*B.c[0]);
    L += plc_M_sqr(A.c[1]*B.c[2] - A.c[2]*B.c[1]);
    L = sqrt(L);

    D = plc_M_dot(A,B);
    angle = atan2(L,D);

    /* This can go wrong only if the numerator and denominator are _both_
       very small (assuming that the system atan2 is fairly robust). */

    *ok = (fabs(L) > 1e-12 || fabs(D) > 1e-12);

    return angle;
  }

inline bool plc_vecteq(plc_vector A, plc_vector B) /*@modifies nothing@*/ {
  return plc_M_vecteq(A,B);
}

/* Procedure returns a vector which points in the same direction as V but has
 * length 1.  It sets *ok to false if the norm is too small. */
inline plc_vector plc_normalize_vect(const plc_vector V,
                                       /*@null@*//*@out@*/ bool *ok) {
  double vnrm;

  vnrm = plc_M_norm(V);
  if (vnrm < DBL_EPSILON && -vnrm < DBL_EPSILON) {
    if (ok != NULL) {
      *ok = false;
    } else {
      fprintf(stderr,
        "plc_normalize_vect: Attempted to normalize zero vector.\n");
      exit(EXIT_FAILURE);
    }
  } else {
    if (ok != NULL) { *ok = true; }
  }
    
  return plc_scale_vect(1.0/vnrm,V);
}

/*
 * George Masaglia's "new method" for finding a random point on a 3-sphere,
 * from his short article "Choosing a Point from the Surface of a Sphere"
 * in The Annals of Mathematical Statistics, v. 43, no. 2 (Apr, 1972) 645-646.
 *
 */
plc_vector plc_random_vect()
{
  int i;
  plc_vector R;
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
inline plc_vector plc_vlincomb(double a,plc_vector A,
                                 double b,plc_vector B) {
  plc_vector R;

  plc_M_vlincomb(R,a,A,b,B);
  return R;
}

inline plc_vector plc_vmadd(plc_vector A, double s, plc_vector B) {
  plc_M_vmadd(A,s,B);
  return A;
}

inline plc_vector plc_vweighted(double s, plc_vector A, plc_vector B) {
  plc_vector R;
  plc_M_vweighted(R,s,A,B);
  return R;
}

/* Put together a vector from 3 doubles */
inline plc_vector plc_build_vect(const double x,
                                   const double y,
                                   const double z) /*@modifies nothing@*/ {
  plc_vector V = { { x, y, z } };
  return V;
}
