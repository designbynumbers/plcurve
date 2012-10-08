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
  double r,theta;
  out = calloc(2*n,sizeof(double));
  if (out == NULL) { return NULL; }
  
  for(step=out;step<out+2*n;) {

    r = sqrt(-2*log(drand48())); 
    theta = 2*pi*drand48();
    
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
  return ret;
}

complex double HermitianDot(complex double *A,complex double *B,int n) 
{
  complex double ret = 0;  int i;
  for(i=0;i<n;i++,A++,B++) { ret += (*A) * conj(*B); }
  return ret;
} 

void ComplexScalarMultiply(complex double s,complex double *A,int n) 
{
  int i;
  for(i=0;i<n;i++,A++) { (*A) *= s; }
}

plCurve *plc_random_closed_polygon_internal(int nEdges, bool selfcheck)
{

  /* 1. Generate vectors of 2n independent Gaussians. */

  double *Araw,*Braw;
  Araw = gaussian_array(nEdges);
  Braw = gaussian_array(nEdges);

  /* 2. Convert to complex. */

  complex double *A,*B;

  A = malloc(nEdges*sizeof(complex double));
  B = malloc(nEdges*sizeof(complex double));

  int i;
  for(i=0;i<nEdges;i++) {
    A[i] = Araw[2*i] + I*Araw[2*i + 1];
    B[i] = Braw[2*i] + I*Braw[2*i + 1];
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

  /* 5a. Selfcheck, if needed. */

  if (selfcheck) { 

    complex double aa, ab, bb;

    aa = HermitianDot(A,A,nEdges);
    bb = HermitianDot(B,B,nEdges);
    ab = HermitianDot(A,B,nEdges);

    if (fabs(creal(aa) - 1.0) > 1e-10 || fabs(cimag(aa)) > 1e-10) {

      fprintf(stderr,"plc_closed_polygon_selfcheck: <A,A> = %g + %g i != 1.0\n",creal(aa),cimag(aa));
      exit(1);

    }

    if (fabs(creal(bb) - 1.0) > 1e-10 || fabs(cimag(bb)) > 1e-10) {

      fprintf(stderr,"plc_closed_polygon_selfcheck: <A,A> = %g + %g i != 1.0\n",creal(bb),cimag(bb));
      exit(1);

    }

    if (fabs(creal(ab)) > 1e-10 || fabs(cimag(ab)) > 1e-10) {

      fprintf(stderr,"plc_closed_polygon_selfcheck: <A,B> = %g + %g i != 0.0\n",creal(ab),cimag(ab));
      exit(1);

    }

    plc_vector edgesum = {{0,0,0}};

    for(i=0;i<nEdges;i++) {

      plc_M_add_vect(edgesum,hopfImap(creal(A[i]),cimag(A[i]),creal(B[i]),cimag(B[i])));

    }

    if (plc_M_norm(edgesum) > 1e-10) { 

      fprintf(stderr,"plc_closed_polygon_selfcheck: Sum of edges is (%g,%g,%g) with norm %g != 0.0\n",
	      plc_M_clist(edgesum),plc_M_norm(edgesum));
      exit(1);

    }

  } 

  /* 6. Assemble Polygon. */

  bool open={false};
  int cc=0,nv = nEdges;
  plCurve *L;

  L = plc_new(1,&nv,&open,&cc);
  L->cp[0].vt[0] = hopfImap(creal(A[0]),cimag(A[0]),creal(B[0]),cimag(B[0]));
  
  for(i=1;i<nEdges;i++) {
    L->cp[0].vt[i] = plc_vect_sum(L->cp[0].vt[i-1],hopfImap(creal(A[i]),cimag(A[i]),creal(B[i]),cimag(B[i])));
  }

  plc_fix_wrap(L);

  /* 7. Cleanup memory. */

  free(Araw); free(Braw);
  free(A); free(B);

  return L;

}
    
plCurve *plc_random_closed_polygon(int nEdges) 
{
  return plc_random_closed_polygon_internal(nEdges,false);
}

plCurve *plc_random_closed_polygon_selfcheck(int nEdges) 
{
  return plc_random_closed_polygon_internal(nEdges,true);
}

/******************** Random Polygon in Plane *********************/

plCurve *plc_random_closed_plane_polygon_internal(int nEdges, bool selfcheck)
{

  /* 1. Generate vector of 2n independent Gaussians. */

  double *Raw;
  Raw = gaussian_array(nEdges);

  /* 2. Convert to real frame vectors (keep complex type for consistency with space version ) */

  complex double *A,*B;

  A = malloc(nEdges*sizeof(complex double));
  B = malloc(nEdges*sizeof(complex double));

  int i;
  for(i=0;i<nEdges;i++) {
    A[i] = Raw[2*i];
    B[i] = Raw[2*i+1];
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

  /* 5a. Selfcheck, if needed. */

  if (selfcheck) { 

    complex double aa, ab, bb;

    aa = HermitianDot(A,A,nEdges);
    bb = HermitianDot(B,B,nEdges);
    ab = HermitianDot(A,B,nEdges);

    if (fabs(creal(aa) - 1.0) > 1e-10 || fabs(cimag(aa)) > 1e-10) {

      fprintf(stderr,"plc_closed_plane_polygon_selfcheck: <A,A> = %g + %g i != 1.0\n",creal(aa),cimag(aa));
      exit(1);

    }

    if (fabs(creal(bb) - 1.0) > 1e-10 || fabs(cimag(bb)) > 1e-10) {

      fprintf(stderr,"plc_closed_plane_polygon_selfcheck: <A,A> = %g + %g i != 1.0\n",creal(bb),cimag(bb));
      exit(1);

    }

    if (fabs(creal(ab)) > 1e-10 || fabs(cimag(ab)) > 1e-10) {

      fprintf(stderr,"plc_closed_plane_polygon_selfcheck: <A,B> = %g + %g i != 0.0\n",creal(ab),cimag(ab));
      exit(1);

    }

    plc_vector edgesum = {{0,0,0}};

    for(i=0;i<nEdges;i++) {

      plc_M_add_vect(edgesum,hopfImap(creal(A[i]),cimag(A[i]),creal(B[i]),cimag(B[i])));

    }

    if (plc_M_norm(edgesum) > 1e-10) { 

      fprintf(stderr,"plc_closed_plane_polygon_selfcheck: Sum of edges is (%g,%g,%g) with norm %g != 0.0\n",
	      plc_M_clist(edgesum),plc_M_norm(edgesum));
      exit(1);

    }

  } 

  /* 6. Assemble Polygon. */

  bool open={false};
  int cc=0,nv = nEdges;
  plCurve *L;

  L = plc_new(1,&nv,&open,&cc);
  L->cp[0].vt[0] = hopfImap(creal(A[0]),cimag(A[0]),creal(B[0]),cimag(B[0]));
  
  for(i=1;i<nEdges;i++) {
    L->cp[0].vt[i] = plc_vect_sum(L->cp[0].vt[i-1],hopfImap(creal(A[i]),cimag(A[i]),creal(B[i]),cimag(B[i])));
  }

  plc_fix_wrap(L);

  /* 7. Cleanup memory. */

  free(Raw);
  free(A); free(B);

  return L;

}
    
plCurve *plc_random_closed_plane_polygon(int nEdges) 
{
  return plc_random_closed_plane_polygon_internal(nEdges,false);
}
 
plCurve *plc_random_closed_plane_polygon_selfcheck(int nEdges) 
{
  return plc_random_closed_plane_polygon_internal(nEdges,true);
}


 
/******************** Random Arm in Space *************************/

plCurve *plc_random_open_polygon_internal(int nEdges, bool selfcheck)
{

  /* 1. Generate vectors of 2n independent Gaussians. */

  double *Araw,*Braw;
  Araw = gaussian_array(nEdges);
  Braw = gaussian_array(nEdges);

  /* 2. Convert to complex. */

  complex double *A,*B;

  A = malloc(nEdges*sizeof(complex double));
  B = malloc(nEdges*sizeof(complex double));

  int i;
  for(i=0;i<nEdges;i++) {
    A[i] = Araw[2*i] + I*Araw[2*i + 1];
    B[i] = Braw[2*i] + I*Braw[2*i + 1];
  }

  /* 3. Normalize A and B so that |A|^2 + |B|^2 = 2. */

  complex double norm = 0;
  norm = sqrt(creal(HermitianDot(A,A,nEdges)) + creal(HermitianDot(B,B,nEdges)));
 
  ComplexScalarMultiply(sqrt(2.0)/norm,A,nEdges);
  ComplexScalarMultiply(sqrt(2.0)/norm,B,nEdges);

  /* 3a. Selfcheck, if needed. */

  if (selfcheck) { 

    complex double aa, bb;

    aa = HermitianDot(A,A,nEdges);
    bb = HermitianDot(B,B,nEdges);

    if (fabs(creal(aa) + creal(bb) - 2.0) > 1e-10 || fabs(cimag(aa)) > 1e-10 || fabs(cimag(bb)) > 1e-10) {

      fprintf(stderr,"plc_open_polygon_selfcheck: <A,A> + <B,B> = %g + %g i != 2.0\n",creal(aa) + creal(bb),cimag(aa) + cimag(bb));
      exit(1);

    }

  } 

  /* 6. Assemble Polygon. */

  bool open={true};
  int cc=0,nv = nEdges+1;
  plCurve *L;

  L = plc_new(1,&nv,&open,&cc);
  L->cp[0].vt[0] = plc_build_vect(0,0,0);
  
  for(i=1;i<nEdges+1;i++) {
    L->cp[0].vt[i] = plc_vect_sum(L->cp[0].vt[i-1],hopfImap(creal(A[i-1]),cimag(A[i-1]),creal(B[i-1]),cimag(B[i-1])));
  }

  plc_fix_wrap(L);

  /* 7. Cleanup memory. */

  free(Araw); free(Braw);
  free(A); free(B);

  return L;

}
    
plCurve *plc_random_open_polygon(int nEdges) 
{
  return plc_random_open_polygon_internal(nEdges,false);
}

plCurve *plc_random_open_polygon_selfcheck(int nEdges)
{
  return plc_random_open_polygon_internal(nEdges,true);
}

/***************** Random Arm in Plane *****************************/

plCurve *plc_random_open_plane_polygon_internal(int nEdges, bool selfcheck)
{

  /* 1. Generate vector of 2n independent Gaussians. */

  double *Raw;
  Raw = gaussian_array(nEdges);

  /* 2. Convert to reals (keep complex type for consistency with space polygon version of code). */

  complex double *A,*B;

  A = malloc(nEdges*sizeof(complex double));
  B = malloc(nEdges*sizeof(complex double));

  int i;
  for(i=0;i<nEdges;i++) {
    A[i] = Raw[2*i];
    B[i] = Raw[2*i+1];
  }

  /* 3. Normalize A and B so that |A|^2 + |B|^2 = 2. */

  complex double norm = 0;
  norm = sqrt(creal(HermitianDot(A,A,nEdges)) + creal(HermitianDot(B,B,nEdges)));
 
  ComplexScalarMultiply(sqrt(2.0)/norm,A,nEdges);
  ComplexScalarMultiply(sqrt(2.0)/norm,B,nEdges);

  /* 3a. Selfcheck, if needed. */

  if (selfcheck) { 

    complex double aa, bb;

    aa = HermitianDot(A,A,nEdges);
    bb = HermitianDot(B,B,nEdges);

    if (fabs(creal(aa) + creal(bb) - 2.0) > 1e-10 || fabs(cimag(aa)) > 1e-10 || fabs(cimag(bb)) > 1e-10) {

      fprintf(stderr,"plc_open_plane_polygon_selfcheck: <A,A> + <B,B> = %g + %g i != 2.0\n",creal(aa) + creal(bb),cimag(aa) + cimag(bb));
      exit(1);

    }

  } 

  /* 6. Assemble Polygon. */

  bool open={true};
  int cc=0,nv = nEdges+1;
  plCurve *L;

  L = plc_new(1,&nv,&open,&cc);
  L->cp[0].vt[0] = plc_build_vect(0,0,0);
  
  for(i=1;i<nEdges+1;i++) {
    L->cp[0].vt[i] = plc_vect_sum(L->cp[0].vt[i-1],hopfImap(creal(A[i-1]),cimag(A[i-1]),creal(B[i-1]),cimag(B[i-1])));
  }

  plc_fix_wrap(L);

  /* 7. Cleanup memory. */

  free(Raw);
  free(A); free(B);

  return L;

}
    
plCurve *plc_random_open_plane_polygon(int nEdges) 
{
  return plc_random_open_plane_polygon_internal(nEdges,false);
}
    
plCurve *plc_random_open_plane_polygon_selfcheck(int nEdges)
{
  return plc_random_open_plane_polygon_internal(nEdges,true);
}

/*************************** PseudoEquilateral Polygons ********************/


plCurve *plc_random_closed_polygon_PE_internal(int nEdges,double LOWER, double UPPER, int *n_attempts,bool selfcheck)

{
  int attempt = 0;
  plCurve *L;
  double longedge,shortedge,meanedge,moment2edge;
  int MAX_ATTEMPTS = 10*nEdges;

  for(attempt=0;attempt < MAX_ATTEMPTS;attempt++ ) {

    L = plc_random_closed_polygon_internal(nEdges,selfcheck);
    plc_edgelength_stats(L,&longedge,&shortedge,&meanedge,&moment2edge);

    if (longedge < UPPER && shortedge > LOWER) {

      if (n_attempts != NULL) { *n_attempts = attempt+1; }
      return L;

    }

    plc_free(L);

  }

  /* We should only get here if the number of attempts is equal to MAX_ATTEMPTS. */
  *n_attempts = attempt;
  fprintf(stderr,"plc_random_closed_polygon_PE: After %d attempts, failed to generate a %d-gon with edges in [%g,%g].\n",attempt,nEdges,LOWER,UPPER);
  return NULL;

}
  
plCurve *plc_random_closed_plane_polygon_PE_internal(int nEdges,double LOWER, double UPPER, int *n_attempts,bool selfcheck)

{
  int attempt = 0;
  plCurve *L;
  double longedge,shortedge,meanedge,moment2edge;
  int MAX_ATTEMPTS = 10*nEdges;

  for(attempt=0;attempt < MAX_ATTEMPTS;attempt++ ) {

    L = plc_random_closed_plane_polygon_internal(nEdges,selfcheck);
    plc_edgelength_stats(L,&longedge,&shortedge,&meanedge,&moment2edge);

    if (longedge < UPPER && shortedge > LOWER) {

      if (n_attempts != NULL) { *n_attempts = attempt+1; }
      return L;

    }

    plc_free(L);

  }

  /* We should only get here if the number of attempts is equal to MAX_ATTEMPTS. */

  *n_attempts = attempt;  
  fprintf(stderr,"plc_random_closed_plane_polygon_PE: After %d attempts, failed to generate a %d-gon with edges in [%g,%g].\n",attempt,nEdges,LOWER,UPPER);
  return NULL;

}


plCurve *plc_random_open_polygon_PE_internal(int nEdges,double LOWER, double UPPER, int *n_attempts,bool selfcheck)

{
  int attempt = 0;
  plCurve *L;
  double longedge,shortedge,meanedge,moment2edge;
  int MAX_ATTEMPTS = 10*nEdges;

  for(attempt=0;attempt < MAX_ATTEMPTS;attempt++ ) {

    L = plc_random_open_polygon_internal(nEdges,selfcheck);
    plc_edgelength_stats(L,&longedge,&shortedge,&meanedge,&moment2edge);

    if (longedge < UPPER && shortedge > LOWER) {

      if (n_attempts != NULL) { *n_attempts = attempt+1; }
      return L;

    }

    plc_free(L);

  }

  /* We should only get here if the number of attempts is equal to MAXATTEMPTS. */

  *n_attempts = attempt;
  fprintf(stderr,"plc_random_open_polygon_PE: After %d attempts, failed to generate a %d-gon with edges in [%g,%g].\n",attempt,nEdges,LOWER,UPPER);
  return NULL;

}


plCurve *plc_random_open_plane_polygon_PE_internal(int nEdges,double LOWER, double UPPER, int *n_attempts,bool selfcheck)

{
  int attempt = 0;
  plCurve *L;
  double longedge,shortedge,meanedge,moment2edge;
  int MAX_ATTEMPTS = 10*nEdges;

  for(attempt=0;attempt < MAX_ATTEMPTS;attempt++ ) {

    L = plc_random_open_plane_polygon_internal(nEdges,selfcheck);
    plc_edgelength_stats(L,&longedge,&shortedge,&meanedge,&moment2edge);

    if (longedge < UPPER && shortedge > LOWER) {

      if (n_attempts != NULL) { *n_attempts = attempt+1; }
      return L;

    }

    plc_free(L);

  }

  /* We should only get here if the number of attempts is equal to MAXATTEMPTS. */

  *n_attempts = attempt;
  fprintf(stderr,"plc_random_open_plane_polygon_PE: After %d attempts, failed to generate a %d-gon with edges in [%g,%g].\n",attempt,nEdges,LOWER,UPPER);
  return NULL;

}  


plCurve *plc_random_closed_polygon_PE(int nEdges,double LOWER,double UPPER) {

  return plc_random_closed_polygon_PE_internal(nEdges,LOWER,UPPER,NULL,false);

}

plCurve *plc_random_open_polygon_PE(int nEdges,double LOWER, double UPPER) {

 return plc_random_open_polygon_PE_internal(nEdges,LOWER,UPPER,NULL,false);

}

plCurve *plc_random_closed_plane_polygon_PE(int nEdges,double LOWER, double UPPER) {
 
  return plc_random_closed_plane_polygon_PE_internal(nEdges,LOWER,UPPER,NULL,false);

}

plCurve *plc_random_open_plane_polygon_PE(int nEdges, double LOWER, double UPPER) {

  return plc_random_open_plane_polygon_PE_internal(nEdges,LOWER,UPPER,NULL,false);

}

plCurve *plc_random_closed_polygon_PE_selfcheck(int nEdges,double LOWER,double UPPER,int *n_attempts) {

  return plc_random_closed_polygon_PE_internal(nEdges,LOWER,UPPER,n_attempts,true);

}

plCurve *plc_random_open_polygon_PE_selfcheck(int nEdges,double LOWER, double UPPER, int *n_attempts)  {

  return plc_random_open_polygon_PE_internal(nEdges,LOWER,UPPER,n_attempts,true);

}

plCurve *plc_random_closed_plane_polygon_PE_selfcheck(int nEdges,double LOWER, double UPPER, int *n_attempts)  {

  return plc_random_closed_plane_polygon_PE_internal(nEdges,LOWER,UPPER,n_attempts,true);

}

plCurve *plc_random_open_plane_polygon_PE_selfcheck(int nEdges, double LOWER, double UPPER, int *n_attempts)  {

  return plc_random_open_plane_polygon_PE_internal(nEdges,LOWER,UPPER,n_attempts,true);

}

/************ Random Space Arm with specified failure to close *************/

#ifdef FTC_IN

plCurve *plc_random_closed_polygon_internal(int nEdges, bool selfcheck)
{

  /* 1. Generate vectors of 2n independent Gaussians. */

  double *Araw,*Braw;
  Araw = gaussian_array(nEdges);
  Braw = gaussian_array(nEdges);

  /* 2. Convert to complex. */

  complex double *A,*B;

  A = malloc(nEdges*sizeof(complex double));
  B = malloc(nEdges*sizeof(complex double));

  int i;
  for(i=0;i<nEdges;i++) {
    A[i] = Araw[2*i] + I*Araw[2*i + 1];
    B[i] = Braw[2*i] + I*Braw[2*i + 1];
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

  /* 5a. Selfcheck, if needed. */

  if (selfcheck) { 

    complex double aa, ab, bb;

    aa = HermitianDot(A,A,nEdges);
    bb = HermitianDot(B,B,nEdges);
    ab = HermitianDot(A,B,nEdges);

    if (fabs(creal(aa) - 1.0) > 1e-10 || fabs(cimag(aa)) > 1e-10) {

      fprintf(stderr,"plc_closed_polygon_selfcheck: <A,A> = %g + %g i != 1.0\n",creal(aa),cimag(aa));
      exit(1);

    }

    if (fabs(creal(bb) - 1.0) > 1e-10 || fabs(cimag(bb)) > 1e-10) {

      fprintf(stderr,"plc_closed_polygon_selfcheck: <A,A> = %g + %g i != 1.0\n",creal(bb),cimag(bb));
      exit(1);

    }

    if (fabs(creal(ab)) > 1e-10 || fabs(cimag(ab)) > 1e-10) {

      fprintf(stderr,"plc_closed_polygon_selfcheck: <A,B> = %g + %g i != 0.0\n",creal(ab),cimag(ab));
      exit(1);

    }

    plc_vector edgesum = {{0,0,0}};

    for(i=0;i<nEdges;i++) {

      plc_M_add_vect(edgesum,hopfImap(creal(A[i]),cimag(A[i]),creal(B[i]),cimag(B[i])));

    }

    if (plc_M_norm(edgesum) > 1e-10) { 

      fprintf(stderr,"plc_closed_polygon_selfcheck: Sum of edges is (%g,%g,%g) with norm %g != 0.0\n",
	      plc_M_clist(edgesum),plc_M_norm(edgesum));
      exit(1);

    }

  } 

  /* 6. Assemble Polygon. */

  bool open={false};
  int cc=0,nv = nEdges;
  plCurve *L;

  L = plc_new(1,&nv,&open,&cc);
  L->cp[0].vt[0] = hopfImap(creal(A[0]),cimag(A[0]),creal(B[0]),cimag(B[0]));
  
  for(i=1;i<nEdges;i++) {
    L->cp[0].vt[i] = plc_vect_sum(L->cp[0].vt[i-1],hopfImap(creal(A[i]),cimag(A[i]),creal(B[i]),cimag(B[i])));
  }

  plc_fix_wrap(L);

  /* 7. Cleanup memory. */

  free(Araw); free(Braw);
  free(A); free(B);

  return L;

}    

#endif  

     

