/*
 * @COPYRIGHT@
 * 
 * Routines for working with vectors.
 *
 * $Id: vector.c,v 1.22 2006-02-22 17:08:24 cantarel Exp $
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

#include "plCurve.h"

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
  plcl_vector C;

  C.c[0] = A.c[0] + B.c[0];
  C.c[1] = A.c[1] + B.c[1];
  C.c[2] = A.c[2] + B.c[2];

  return C;
}
  
/* Returns A - B. */ 
inline plcl_vector plcl_vect_diff(plcl_vector A,plcl_vector B) { 
  plcl_vector C;

  C.c[0] = A.c[0] - B.c[0];
  C.c[1] = A.c[1] - B.c[1];
  C.c[2] = A.c[2] - B.c[2];

  return C;
}
  
/* Returns A x B. */
inline plcl_vector plcl_cross_prod(plcl_vector A,plcl_vector B) { 
  plcl_vector C;

  C.c[0] = A.c[1] * B.c[2] - A.c[2] * B.c[1];
  C.c[1] = A.c[2] * B.c[0] - A.c[0] * B.c[2];
  C.c[2] = A.c[0] * B.c[1] - A.c[1] * B.c[0];

  return C;
}

/* Returns sA. */
inline plcl_vector plcl_scale_vect(double s,plcl_vector A) { 
  plcl_vector C;

  C.c[0] = s * A.c[0];
  C.c[1] = s * A.c[1];
  C.c[2] = s * A.c[2];

  return C;
}

inline plcl_vector plcl_component_mult(plcl_vector A,plcl_vector B) {
  plcl_vector C;

  C.c[0] = A.c[0]*B.c[0];
  C.c[1] = A.c[1]*B.c[1];
  C.c[2] = A.c[2]*B.c[2];

  return C;
}

inline plcl_vector plcl_component_div(plcl_vector A,plcl_vector B) {
  plcl_vector C;

  C.c[0] = A.c[0]/B.c[0];
  C.c[1] = A.c[1]/B.c[1];
  C.c[2] = A.c[2]/B.c[2];

  return C;
}

/* Returns the dot product of A and B */
inline double plcl_dot_prod(plcl_vector A, plcl_vector B) {
  return A.c[0]*B.c[0] + A.c[1]*B.c[1] + A.c[2]*B.c[2];
}

inline double plcl_norm(plcl_vector A) {
  return sqrt(plcl_dot_prod(A,A));
}

/* Procedure replaces V with a unit vector (if possible). */
inline plcl_vector plcl_normalize_vect(const plcl_vector V) {
  double vnrm;

  plcl_error_num = 0;

  vnrm = plcl_M_norm(V);
  if (vnrm < DBL_EPSILON && -vnrm < DBL_EPSILON) {
    plcl_error_num = PLCL_E_ZERO_VECTOR;
    snprintf(plcl_error_str,sizeof(plcl_error_str),
      "plcl_normalize_vect: Can't normalize zero vector.\n");
    return V;
  }
  return plcl_scale_vect(1.0/vnrm,V);
}

/*
 * George Masaglia's "new method" for finding a random point on a 3-sphere,
 * from his short article "Choosing a Point from the Surface of a Sphere"
 * in The Annals of Mathematical Statistics, v. 43, no. 2 (Apr, 1972) 645-646.
 *
 */
plcl_vector plcl_random_2()
{
  int i;
  plcl_vector R;
  double V1 = 0.0, V2 = 0.0;
  double S = 0.0;
  double sqt;

  for (i = 0; i < 1000 && 
              (S - 1.0 > DBL_EPSILON ||
               S - 0.01 < DBL_EPSILON); i++) {
    V1 = 2*(double)(rand())/(double)(RAND_MAX) - 1;
    V2 = 2*(double)(rand())/(double)(RAND_MAX) - 1;
    S = V1*V1 + V2*V2;
  }
  assert(S - 0.01 >= DBL_EPSILON && S - 1.0 <= DBL_EPSILON);
  sqt = sqrt(1-S);
  R.c[0] = 2*V1*sqt;
  R.c[1] = 2*V2*sqt;
  R.c[2] = 1-2*S;

  return R;
}

/*
 * Procedure assumes that the "rand" random number generator exists on this
 * system and has been seeded. It generates a random unit vector. 
 */
plcl_vector plcl_random_vect()
{
  int i;
  plcl_vector R;

  for(i=0;i<1000;i++) {
    R.c[0] = 2*(double)(rand())/(double)(RAND_MAX) - 1;
    R.c[1] = 2*(double)(rand())/(double)(RAND_MAX) - 1;
    R.c[2] = 2*(double)(rand())/(double)(RAND_MAX) - 1;

    if (plcl_M_norm(R) > 0.1 && plcl_M_norm(R) < 0.9) {
      return plcl_normalize_vect(R);
    }
  }
  /* If we make it to here, there is evidently a problem with the random 
     number generator, we haven't had any acceptable vectors in 1000 tries. */
  assert(false);
  return plcl_build_vect(1,0,0);
}

/* Put together a vector from 3 doubles */
inline plcl_vector plcl_build_vect(const double x, 
                                   const double y,
                                   const double z) {
  plcl_vector V = { { x, y, z } };
  return V;
}
