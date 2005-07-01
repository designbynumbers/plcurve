/*
 * 
 * Prototypes for routines in vector.c.
 *
 * $Id: vector.h,v 1.3 2005-07-01 01:56:33 cantarel Exp $
 *
 */

/* Copyright 2004 The University of Georgia */

/* This file is part of vecttools.
   
vecttools is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

vecttools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with vecttools; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#ifndef __LINKLIB_VECTOR_H
#define __LINKLIB_VECTOR_H


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_ASSERT_H
#include <assert.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

typedef struct linklib_vector_type {    /* A point in 3-space. */
  double c[3];
} linklib_vector;

/*
 * Prototypes for vector routines. 
 *
 */

inline linklib_vector linklib_vplus(linklib_vector A,linklib_vector B);
inline linklib_vector linklib_vminus(linklib_vector A,linklib_vector B);
inline linklib_vector linklib_cross(linklib_vector A,linklib_vector B);
inline linklib_vector linklib_scalarmult(double x,linklib_vector A);
inline linklib_vector linklib_vdivide(linklib_vector A,linklib_vector B);

inline double  linklib_vdist(linklib_vector A,linklib_vector B);

inline void linklib_vector_normalize(linklib_vector *V);

/*
 * Macro replacements (requires some reprogramming)
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


linklib_vector linklib_vector_random();

#endif
