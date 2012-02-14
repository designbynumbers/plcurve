/*
 *  Routines to generate random polygons as part of plCurve.
 *
 */

/* Copyright 2009 The University of Georgia. */

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

#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#ifdef HAVE_MATH_H
  #include <math.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef HAVE_STDARG_H
#include <stdarg.h>
#endif
#ifdef HAVE_CTYPE_H
#include <ctype.h>
#endif
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#ifdef HAVE_ASSERT_H
#include <assert.h>
#endif
#ifdef HAVE_FLOAT_H
#include <float.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

double *gaussian_array(int n)
/* Returns 2n independent standard Gaussians, generated using the Box-Muller method and drand48 */
{
  double *out,*step;
  double pi = 3.14159265358979323846264338327;
  int i;
  double r,theta;
  out = malloc(2*n*sizeof(double));
  if (out == NULL) { return NULL; }
  
  for(step=out;step<out+2*n;) {

    r = sqrt(-2*log(drand48())); theta = 2*pi*drand48();
    *step = r*cos(theta); step++;
    *step = r*sin(theta); step++;

  }

  return out;
}

plc_vector hopfImap(double q0,double q1,double q2,double q3)
{
  plc_vector ret;
  ret.c[0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
  ret.c[1] = 2*q1*q2 - 2*q0*q3;
  ret.c[2] = 2*q0*q2 + 2*q1*q3;
}

complex double HermitianDot(complex double *A,complex double *B,int n) 
{
  complex double ret = 0;  
  for(i=0;i<n;i++,A++,B++) { ret += (*A) * conj(*B); }
  return ret;
} 

void ComplexScalarMultiply(complex double s,complex double *A,int n) 
{
  int i;
  for(i=0;i<n;i++,A++) { (*A) *= s; }
}

plCurve *plc_random_closed_polygon(int nEdges) 

{

  /* 1. Generate vectors of 2n independent Gaussians. */

  double *Araw,*Braw;
  Araw = gaussian_array(2*nEdges);
  Braw = gaussian_array(2*nEdges);

  /* 2. Convert to complex. */

  complex double *A,*B;

  A = malloc(nEdges,sizeof(complex double));
  B = malloc(nEdges,sizeof(complex double));

  int i;
  for(i=0;i<nEdges;i++) {
    A[i] = Araw[2*i] + I*Araw[2*i + 1];
    B[i] = Braw[2*i] + I*Araw[2*i + 1];
  }

  /* 3. Normalize A. */

  complex double norm = 0;
  norm = creal(sqrt(HermitianDot(A,A,nEdges)));
  ComplexScalarMultiply(1/norm,A,nEdges);

  /* 4. Set B to B - conj(<A,B>) A */

  complex double s;
  s = conj(HermitianDot(A,B,nEdges));
  for(i=0;i<nEdges;i++) { B[i] -= s*A[i]; }

  /* 5. Normalize B. */

  norm = creal(sqrt(HermitianDot(B,B,nEdges)));
  ComplexScalarMultiply(1/norm,B,nEdges);

  /* 6. Assemble Polygon. */

  bool open={false};
  int cc=0,nv = nEdges;
  plCurve *L;

  L = plc_new(1,&nv,&open,&cc);
  L->cp[0].vt[0] = plc_build_vect(0,0,0);
  
  for(i=1;i<nEdges;i++) {
    L->cp[0].vt[i] = plc_vect_sum(L->cp[0].vt[i-1],hopfImap(creal(A[i]),cimag(A[i]),creal(B[i]),cimag(B[i])));
  }

  /* 7. Cleanup memory. */

  free(Araw); free(Braw);
  free(A); free(B);

  return L;

}
    
				 

  
  


  
  



  


    

  

     

