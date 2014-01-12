/*
 * @COPYRIGHT@
 *
 * Routines for working with n-dimensional vectors.
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
#include <nplCurve.h>
  
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

#ifdef HAVE_STDARG_H
#include <stdarg.h>
#endif
  
/**************************************************************/
/*   Basic Linear Algebra Operations                          */
/**************************************************************/
  
  nplc_vector nplc_vect_new(int n)
{
  nplc_vector nv;
  nv.n = n;
  nv.c = calloc(n,sizeof(double));
  assert(nv.c != NULL);
  return nv;
}

void nplc_vect_free(nplc_vector *nv)
{
  nv->n = 0;
  if (nv->c != NULL) { free(nv->c); }
  nv->c = NULL;
}

void nplc_vect_copy(nplc_vector to,nplc_vector from)
{
  assert(from.n == to.n);
  assert(from.c != NULL);

  int i;

  for(i=0;i<from.n;i++) { to.c[i] = from.c[i]; }
}

nplc_vector *nplc_vect_buf_new(int dim,int num_vects)
     /* Allocate a buffer of nplc vectors */
{
  nplc_vector *retbuf;
  int i;

  retbuf = calloc(num_vects,sizeof(nplc_vector));
  for (i=0;i<num_vects;i++) {
    
    retbuf[i] = nplc_vect_new(dim);

  }
  return retbuf;
}

void nplc_vect_buf_free(int bufsize,nplc_vector *buf)
{
  int i;

  for(i=0;i<bufsize;i++) {

    nplc_vect_free(&(buf[i]));

  }

}

int nplc_vect_dim(nplc_vector A,nplc_vector B) 
  /* Checks dimensions and returns common dimension. */
{
  assert(A.n == B.n);
  return A.n;
}

/* Returns A + B. */
inline nplc_vector nplc_vect_sum(nplc_vector A,nplc_vector B) {
  int i,n;
  nplc_vector ret;
  n = nplc_vect_dim(A,B);
  ret = nplc_vect_new(n);
  for(i=0;i<n;i++) { ret.c[i] = A.c[i] + B.c[i]; }
  return ret;
}

/* Returns A - B. */
inline nplc_vector nplc_vect_diff(nplc_vector A,nplc_vector B) {
  int i,n;
  nplc_vector ret;
  n = nplc_vect_dim(A,B);
  ret = nplc_vect_new(n);
  for(i=0;i<n;i++) { ret.c[i] = A.c[i] - B.c[i]; }
  return ret;
}


/* Returns sA. */
inline nplc_vector nplc_scale_vect(double s,nplc_vector A) {
  int i,n;
  nplc_vector ret;
  n = nplc_vect_dim(A,A);
  ret = nplc_vect_new(n);
  for(i=0;i<n;i++) { ret.c[i] = s*A.c[i]; }
  return ret;
}

inline nplc_vector nplc_component_mult(nplc_vector A,nplc_vector B) {
  int i,n;
  nplc_vector ret;
  n = nplc_vect_dim(A,B);
  ret = nplc_vect_new(n);
  for(i=0;i<n;i++) { ret.c[i] = A.c[i] * B.c[i]; }
  return ret;
}

/* Should we add an "ok" parameter here as in _normalize_vect? */
inline nplc_vector nplc_component_div(nplc_vector A,nplc_vector B,
                                      /*@null@*/ bool *ok) {

  int i,n;
  nplc_vector ret;

  n = nplc_vect_dim(A,B);
  ret = nplc_vect_new(n);

  for(i=0;i<n;i++) { 
    
    if (fabs(B.c[i]) < DBL_EPSILON) { 

      if (ok != NULL) { *ok = false; return ret; }
      else { fprintf(stderr,"nplc_component_div: Divisor has zero component.\n"); exit(EXIT_FAILURE); }
      
    }

  }

  for (i=0;i<n;i++) {

    ret.c[i] = A.c[i] / B.c[i]; 

  }

  if (ok != NULL) { *ok = true; }

  return ret;

}

/* Returns the dot product of A and B */
inline double nplc_dot_prod(nplc_vector A, nplc_vector B)
{
  int i,n;
  double ret=0;
  n = nplc_vect_dim(A,B);
  for(i=0;i<n;i++) { ret += A.c[i] * B.c[i]; }
  return ret;
}


inline double nplc_norm(nplc_vector A) {
  return sqrt(nplc_dot_prod(A,A));
}

inline double nplc_distance(nplc_vector A, nplc_vector B) {
  nplc_vector diff;
  double dist;
  diff = nplc_vect_diff(A,B);
  dist = nplc_norm(diff);
  nplc_vect_free(&diff);
  return dist;
}

/* The square of the distance between A and B (faster than _distance) */
inline double nplc_sq_dist(nplc_vector A, nplc_vector B) {

  int i,n;
  double ret=0;
  n = nplc_vect_dim(A,B);
  for(i=0;i<n;i++) { ret += A.c[i] * B.c[i]; }
  return ret;
}

#define nplc_M_sqr(A) \
   ((A)*(A))

/* Computes the angle between two vectors. Can fail if one or the other 
   has norm zero. */
double nplc_angle(nplc_vector A, nplc_vector B, bool *ok) 
 
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

    int i,j,n;

    n = nplc_vect_dim(A,B);

    for(j=1;j<n;j++) {

      for(i=0;i<j;i++) {

	L += nplc_M_sqr(A.c[i]*B.c[j] - A.c[j]*B.c[i]);

      }

    }

    L = sqrt(L);
    D = nplc_dot_prod(A,B);
   
    angle = atan2(L,D);

    /* This can go wrong only if the numerator and denominator are _both_
       very small (assuming that the system atan2 is fairly robust). */

    *ok = (fabs(L) > 1e-12 || fabs(D) > 1e-12);

    return angle;
  }

inline bool nplc_vecteq(nplc_vector A, nplc_vector B) /*@modifies nothing@*/ {
  int i,n;

  n = nplc_vect_dim(A,B);

  for(i=0;i<n;i++) { 
    
    if ((A.c[i] - B.c[i] <= -2*DBL_EPSILON) ||
      (A.c[i] - B.c[i] >=  2*DBL_EPSILON)) { return false; }

  }

  return true;
}
	
/* Procedure returns a vector which points in the same direction as V but has
 * length 1.  It sets *ok to false if the norm is too small. */
inline nplc_vector nplc_normalize_vect(const nplc_vector V,
                                       /*@null@*//*@out@*/ bool *ok) {
  double vnrm;

  vnrm = nplc_norm(V);
  if (vnrm < DBL_EPSILON && -vnrm < DBL_EPSILON) {
    if (ok != NULL) {
      *ok = false;
    } else {
      fprintf(stderr,
        "nplc_normalize_vect: Attempted to normalize zero vector.\n");
      exit(EXIT_FAILURE);
    }
  } else {
    if (ok != NULL) { *ok = true; }
  }
    
  return nplc_scale_vect(1.0/vnrm,V);
}

/* Random unit n-vector */

nplc_vector nplc_random_vect(int n)
{
  int i,j;
  nplc_vector ret,nrm;

  ret = nplc_vect_new(n);
  
  for(j=0;j<1000;j++) {

    for(i=0;i<n;i++) {
      
      ret.c[i] = 2*(double)(rand())/(double)(RAND_MAX) - 1;
    
    }

    if (nplc_norm(ret) < 1.0 && nplc_norm(ret) > 0.1) {

      nrm = nplc_normalize_vect(ret,NULL);
      nplc_vect_free(&ret);

      return nrm;

    }

  }

  fprintf(stderr,
	  "nplc_random_vect: Library error in generating random vect.\n"
	  "                  rand() may be broken on this machine.\n");
  exit(EXIT_FAILURE);
}

/* Return a linear combination: a*A + b*B */
inline nplc_vector nplc_vlincomb(double a,nplc_vector A,
                                 double b,nplc_vector B) {

  int i,n;
  nplc_vector ret;
  n = nplc_vect_dim(A,B);
  ret = nplc_vect_new(n);
  for(i=0;i<n;i++) { ret.c[i] += a*A.c[i] + b*B.c[i]; }
  return ret;

}

inline nplc_vector nplc_vmadd(nplc_vector A, double s, nplc_vector B) {
 
  int i,n;
  n = nplc_vect_dim(A,B);
  for(i=0;i<n;i++) { A.c[i] += s*B.c[i]; }
  return A;

}

inline nplc_vector nplc_vweighted(double s, nplc_vector A, nplc_vector B) {

  nplc_vector ret;
  ret = nplc_vlincomb((1-s),A,s,B);
  return ret;

}

/* Put together a vector from 3 doubles */
inline nplc_vector nplc_build_vect(int n,...) {

  int i;
  nplc_vector ret;
  va_list ap;

  ret = nplc_vect_new(n);

  va_start(ap,n);
  for(i=0;i<n;i++) {
    ret.c[i] = va_arg(ap,double);
  }
  va_end(ap);

  return ret;
}

char *nplc_vect_print(nplc_vector A)
{
  static char pbuf[16000];
  int i,used=0;
  
  used += snprintf(&(pbuf[used]),sizeof(pbuf),"(");
  
  for(i=0;i<A.n;i++) {
    
    used += snprintf(&(pbuf[used]),sizeof(pbuf)-used,"%g,",A.c[i]);
    
  }
  
  used += snprintf(&(pbuf[used-1]),sizeof(pbuf)-used,")");
  return pbuf;
}

char *nplc_vect_clist(nplc_vector A)
{
  static char pbuf[16000];
  int i,used=0;
  
  for(i=0;i<A.n;i++) {
    
    used += snprintf(&(pbuf[used]),sizeof(pbuf)-used,"%.16g ",A.c[i]);
    
  }
  
  return pbuf;
}

  
