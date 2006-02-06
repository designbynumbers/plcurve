/*
 * @COPYRIGHT@
 * 
 * Routines for working with vectors.
 *
 * $Id: vector.c,v 1.13 2006-02-06 00:10:12 ashted Exp $
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

#ifdef HAVE_MATH_H
#include <math.h>
#endif

#ifdef HAVE_STRING_H
#include <string.h>
#endif 

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif 

#include "plCurve.h"

/**************************************************************/
/*   Basic Linear Algebra Operations                          */
/**************************************************************/

/* Returns A + B. */
inline plcl_vector linklib_vplus(plcl_vector A,plcl_vector B) { 
  plcl_vector C;

  plcl_error_num = plcl_error_str[0] = 0;

  C.c[0] = A.c[0] + B.c[0];
  C.c[1] = A.c[1] + B.c[1];
  C.c[2] = A.c[2] + B.c[2];

  return C;
}
  
/* Returns A - B. */ 
inline plcl_vector linklib_vminus(plcl_vector A,plcl_vector B) { 
  plcl_vector C;

  plcl_error_num = plcl_error_str[0] = 0;

  C.c[0] = A.c[0] - B.c[0];
  C.c[1] = A.c[1] - B.c[1];
  C.c[2] = A.c[2] - B.c[2];

  return C;
}
  
/* Returns A x B. */
inline plcl_vector linklib_cross(plcl_vector A,plcl_vector B) { 
  plcl_vector C;

  plcl_error_num = plcl_error_str[0] = 0;

  C.c[0] = A.c[1] * B.c[2] - A.c[2] * B.c[1];
  C.c[1] = A.c[2] * B.c[0] - A.c[0] * B.c[2];
  C.c[2] = A.c[0] * B.c[1] - A.c[1] * B.c[0];

  return C;
}

/* Returns xA. */
inline plcl_vector linklib_scalarmult(double x,plcl_vector A) { 
  plcl_vector C;

  plcl_error_num = plcl_error_str[0] = 0;

  C.c[0] = x * A.c[0];
  C.c[1] = x * A.c[1];
  C.c[2] = x * A.c[2];

  return C;
}

inline plcl_vector linklib_vdivide(plcl_vector A,plcl_vector B) {
  plcl_vector C;

  plcl_error_num = plcl_error_str[0] = 0;

  C.c[0] = A.c[0]/B.c[0];
  C.c[1] = A.c[1]/B.c[1];
  C.c[2] = A.c[2]/B.c[2];

  return C;
}

/* Returns the dot product of A and B */
inline double plcl_dot(plcl_vector A, plcl_vector B) {
  plcl_error_num = plcl_error_str[0] = 0;
  return A.c[0]*B.c[0] + A.c[1]*B.c[1] + A.c[2]*B.c[2];
}

inline double plcl_norm(plcl_vector A) {
  plcl_error_num = plcl_error_str[0] = 0;
  return sqrt(plcl_dot(A,A));
}

inline double linklib_vdist(plcl_vector A,plcl_vector B) 
{
  plcl_vector diff;

  plcl_error_num = plcl_error_str[0] = 0;

  diff = A;
  linklib_vsub(diff,B);
  
  return plcl_M_norm(diff);
}

/* Procedure replaces V with a unit vector (if possible). */
void plcl_vector_normalize(plcl_vector *V)
{
  double vnrm;

  plcl_error_num = plcl_error_str[0] = 0;

  vnrm = plcl_M_norm((*V));
  if (vnrm == 0) {
    plcl_error_num = PLCL_E_ZERO_VECTOR;
#ifdef HAVE_STRLCPY
    strlcpy(plcl_error_str,
      "plcl_vector_normalize: Can't normalize zero vector.\n",
      sizeof(plcl_error_str));
#else
    strncpy(plcl_error_str,
      "plcl_vector_normalize: Can't normalize zero vector.\n",
      sizeof(plcl_error_str)-1);
    plcl_error_str[sizeof(plcl_error_str)-1] = '\0';
#endif
  }
  linklib_vsmult(1/vnrm,(*V));
}

/*
 * Procedure assumes that the "rand" random number generator exists on this
 * system and has been seeded. It generates a random unit vector. 
 */
plcl_vector plcl_vector_random()
{
  int i;
  plcl_vector R;

  plcl_error_num = plcl_error_str[0] = 0;

  R.c[0] = 1;            /* The nonrandom vector */
  R.c[1] = R.c[2] = 0;

  for(i=0;i<1000;i++) {
    R.c[0] = 2*(double)(rand())/(double)(RAND_MAX) - 1;
    R.c[1] = 2*(double)(rand())/(double)(RAND_MAX) - 1;
    R.c[2] = 2*(double)(rand())/(double)(RAND_MAX) - 1;

    if (plcl_M_norm(R) > 0.1 && plcl_M_norm(R) < 0.9) {
      plcl_vector_normalize(&R);
      return R;
    }
  }

  plcl_error_num = PLCL_E_BAD_RANDOM;
  sprintf(plcl_error_str,"plcl_vector_random: Apparent error in rand().\n");
  return R;
}
