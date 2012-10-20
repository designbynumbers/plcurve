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

#ifdef HAVE_GSL_GSL_RNG_H
#include <gsl/gsl_rng.h>
#endif

#ifdef HAVE_GSL_GSL_RANDIST_H
#include <gsl/gsl_randist.h>
#endif

double *gaussian_array(gsl_rng *r, int n)
/* Returns 2n independent standard Gaussians, generated using gsl_ran_ugaussian_ratio_method() */
{
  double *out,*step;
 
  out = calloc(2*n,sizeof(double));
  if (out == NULL) { return NULL; }
  
  for(step=out;step<out+2*n;step++) {

    *step = gsl_ran_ugaussian_ratio_method(r);

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

plCurve *plc_random_closed_polygon_internal(gsl_rng *r, int nEdges, bool selfcheck)
{

  /* 1. Generate vectors of 2n independent Gaussians. */

  double *Araw,*Braw;
  Araw = gaussian_array(r,nEdges);
  Braw = gaussian_array(r,nEdges);

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
    
plCurve *plc_random_closed_polygon(gsl_rng *r, int nEdges) 
{
  return plc_random_closed_polygon_internal(r,nEdges,false);
}

plCurve *plc_random_closed_polygon_selfcheck(gsl_rng *r, int nEdges) 
{
  return plc_random_closed_polygon_internal(r,nEdges,true);
}

/******************** Random Polygon in Plane *********************/

plCurve *plc_random_closed_plane_polygon_internal(gsl_rng *r, int nEdges, bool selfcheck)
{

  /* 1. Generate vector of 2n independent Gaussians. */

  double *Raw;
  Raw = gaussian_array(r,nEdges);

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
    
plCurve *plc_random_closed_plane_polygon(gsl_rng *r,int nEdges) 
{
  return plc_random_closed_plane_polygon_internal(r,nEdges,false);
}
 
plCurve *plc_random_closed_plane_polygon_selfcheck(gsl_rng *r,int nEdges) 
{
  return plc_random_closed_plane_polygon_internal(r,nEdges,true);
}


 
/******************** Random Arm in Space *************************/

plCurve *plc_random_open_polygon_internal(gsl_rng *r, int nEdges, bool selfcheck)
{

  /* 1. Generate vectors of 2n independent Gaussians. */

  double *Araw,*Braw;
  Araw = gaussian_array(r,nEdges);
  Braw = gaussian_array(r,nEdges);

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
    
plCurve *plc_random_open_polygon(gsl_rng *r,int nEdges) 
{
  return plc_random_open_polygon_internal(r,nEdges,false);
}

plCurve *plc_random_open_polygon_selfcheck(gsl_rng *r,int nEdges)
{
  return plc_random_open_polygon_internal(r,nEdges,true);
}

/***************** Random Arm in Plane *****************************/

plCurve *plc_random_open_plane_polygon_internal(gsl_rng *r,int nEdges, bool selfcheck)
{

  /* 1. Generate vector of 2n independent Gaussians. */

  double *Raw;
  Raw = gaussian_array(r,nEdges);

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
    
plCurve *plc_random_open_plane_polygon(gsl_rng *r,int nEdges) 
{
  return plc_random_open_plane_polygon_internal(r,nEdges,false);
}
    
plCurve *plc_random_open_plane_polygon_selfcheck(gsl_rng *r,int nEdges)
{
  return plc_random_open_plane_polygon_internal(r,nEdges,true);
}


/************ Random Equilateral Polygons (and Arms) ***********************/

plCurve *plc_random_equilateral_closed_polygon(gsl_rng *r,int nEdges)

{

  if (nEdges < 2 || nEdges > 100000000) {

    fprintf(stderr,"plc_random_equilateral_closed_polygon: Can't generate polygon with %d edges.\n",nEdges);
    exit(1);
    
  }

  /* Step 0. Generate nEdges unit vectors which sum to zero */

  double pi = 3.14159265358979;
  plc_vector *edgeset;
  int i;

  edgeset = calloc(nEdges,sizeof(plc_vector));
  if (edgeset == NULL) { 

    fprintf(stderr,"plc_random_equilateral_closed_polygon: Can't allocate %d edge vectors.\n",nEdges);
    exit(1);
    
  } 
  
  if (nEdges % 2 == 1) { /* The number of edges is odd; add a triple */

    edgeset[0] = plc_build_vect(1,0,0);
    edgeset[1] = plc_build_vect(cos(2.0*pi/3.0),sin(2.0*pi/3.0),0);
    edgeset[2] = plc_build_vect(cos(4.0*pi/3.0),sin(4.0*pi/3.0),0);
    
    i = 3;
  
  } else {

    i = 0;

  }

  for(;i<nEdges;i+=2) {

    gsl_ran_dir_3d (r, &(edgeset[i].c[0]), &(edgeset[i].c[1]), &(edgeset[i].c[2]));
    edgeset[i+1] = plc_scale_vect(-1.0,edgeset[i]);
    
  }

  /* Step 1. Randomly permute them. */

  gsl_ran_shuffle(r,edgeset,nEdges,sizeof(plc_vector));

  /* Step 2. Apply ``hedgehog moves'' */

  int *indexlist;
  indexlist = calloc(nEdges,sizeof(int));
  if (indexlist == NULL) { 

    fprintf(stderr,"plc_random_equilateral_closed_polygon: Can't allocate %d integers.\n",nEdges);
    exit(1);
    
  } 

  for(i=0;i<nEdges;i++) { indexlist[i] = i; }

  for(i=0;i<nEdges*nEdges;i++) {
    
    int edges[2];
    gsl_ran_choose(r,edges,2,indexlist,nEdges,sizeof(int));  /* Choose a pair of edges to twist */

    plc_vector axis;
    bool ok;

    axis = plc_normalize_vect(plc_vect_sum(edgeset[edges[0]],edgeset[edges[1]]),&ok);

    if (ok) {

      double theta;
      theta = gsl_ran_flat(r,0,2*pi);
      edgeset[edges[0]] = plc_rotate_vect(edgeset[edges[0]],axis,theta);
      edgeset[edges[1]] = plc_rotate_vect(edgeset[edges[1]],axis,theta);
	    
    }
  
  }
      
  /* Step 3. Assemble the edgeset into vertices for a new polygon, scaling as we go. */

  plCurve *L;
  bool open = { false };
  int  cc = { 0 };
  int  nv = nEdges;

  L = plc_new(1,&nv,&open,&cc);
  L->cp[0].vt[0] = plc_build_vect(0,0,0);

  for(i=1;i<nEdges;i++) {
    
    L->cp[0].vt[i] = plc_vect_sum(L->cp[0].vt[i-1],plc_scale_vect(2.0/nEdges,edgeset[i-1]));

  }

  /* Step 4. Housekeeping and return. */

  free(edgeset);
  free(indexlist);

  return L;

}   
    

plCurve *plc_random_equilateral_open_polygon(gsl_rng *r,int nEdges) 
{

  plCurve *L;
  bool open = { true };
  int  cc = { 0 };
  int  nv = nEdges+1;
  int i;

  L = plc_new(1,&nv,&open,&cc);
  L->cp[0].vt[0] = plc_build_vect(0,0,0);

  for(i=1;i<nv;i++) {

    plc_vector dir; 
    gsl_ran_dir_3d(r, &(dir.c[0]), &(dir.c[1]), &(dir.c[2]));
    
    L->cp[0].vt[i] = plc_vect_sum(L->cp[0].vt[i-1],plc_scale_vect(2.0/nEdges,dir));

  }

  return L;

}

/************ Random Space Arm with specified failure to close *************/

#ifdef FTC_IN

plCurve *plc_random_closed_polygon_internal(gsl_rng *r, int nEdges, bool selfcheck)
{

  /* 1. Generate vectors of 2n independent Gaussians. */

  double *Araw,*Braw;
  Araw = gaussian_array(r,nEdges);
  Braw = gaussian_array(r,nEdges);

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

     

