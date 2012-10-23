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

/* We are going to implement the "triple rotation" algorithm of Varela, Hinson,
   Arsuaga, and Diao. */

plc_vector plc_subsphere_sample(gsl_rng *r,plc_vector axis,double height) 

/* Sample uniformly from the portion of a sphere with "axis" coordinate < height. */

{
  plc_vector z_sample;
  plc_vector n_axis;
  plc_vector z_axis = { {0,0,1} };
  double samp_z;
  double samp_theta;
  double pi = 3.14159265358979;
  bool ok;


  n_axis = plc_normalize_vect(axis,&ok);
  if (!ok) { 
    fprintf(stderr,"plc_subsphere_sample: Axis too short.\n");
    exit(1);
  }

  samp_z = gsl_ran_flat(r,-1,height);
  samp_theta = gsl_ran_flat(r,0,2*pi);

  z_sample = plc_build_vect(cos(samp_theta)*sqrt(1 - samp_z*samp_z),
			    sin(samp_theta)*sqrt(1 - samp_z*samp_z),
			    samp_z);

  /* Now we rotate to make an "axis sample" instead of a z sample. */

  plc_vector raxis;
  raxis = plc_normalize_vect(plc_vect_sum(n_axis,z_axis),&ok);

  if (!ok) { /* The n_axis is - the z-axis. No problem; just invert the sample in z and we're done. */

    z_sample.c[2] *= -1.0;
    return z_sample;

  } else {

    return plc_rotate_vect(z_sample,raxis,pi);

  }

}

plc_vector plc_sphere_intersection_sample(gsl_rng *r,plc_vector A, plc_vector B) 

/* Sample uniformly from the sphere of intersection of the unit radius
   spheres centered at A and B. */

{
  double pi = 3.14159265358979;

  plc_vector diff_axis,z_axis = {{0,0,1.0}},minusz = {{0,0,-1.0}};
  plc_vector frameA, frameB, frameC;

  double     theta_sample;
  double     ABdist;

  plc_vector sample;
  bool       ok;

  /* Make sure that there _is_ an intersection */ 
  
  diff_axis = plc_vect_diff(A,B); 
  ABdist = plc_norm(diff_axis);
  
  if (ABdist > 2.0) {
    return plc_vweighted(0.5,A,B);  /* Return the midpoint: spheres just touch */
  }

  if (ABdist < 1e-8) { /* Return a random point on the sphere: spheres coincide */
    gsl_ran_dir_3d(r,&sample.c[0],&sample.c[1],&sample.c[2]);
    return plc_vect_sum(A,sample);
  }

  diff_axis = plc_normalize_vect(diff_axis,&ok);

  /* Now sample on an appropriate circle */

  theta_sample = gsl_ran_flat(r,0,4*pi);

  if (plc_distance(diff_axis,z_axis) < 1e-5) {

    frameA = plc_build_vect(1,0,0);
    frameB = plc_build_vect(0,1,0);
    frameC = plc_build_vect(0,0,1);

  } else if (plc_distance(diff_axis,minusz) < 1e-5) {

    frameA = plc_build_vect(1,0,0);
    frameB = plc_build_vect(0,1,0);
    frameC = plc_build_vect(0,0,-1);

  } else {
    
    bool ok;
    frameA = plc_normalize_vect(plc_cross_prod(diff_axis,z_axis),&ok);
    frameB = plc_normalize_vect(plc_cross_prod(diff_axis,frameA),&ok);
    frameC = diff_axis;

  }

  double radius = sqrt(1 - (ABdist*ABdist)/4.0);

  sample = plc_vect_sum(
			plc_vlincomb(cos(theta_sample)*radius,frameA,
				     sin(theta_sample)*radius,frameB),
			plc_scale_vect(ABdist/2.0,frameC)
			);

  sample = plc_vect_sum(sample,B);

  return sample;
  
}

void plc_double_rotation(gsl_rng *r, 
			 plc_vector r1,plc_vector r2,plc_vector r3,
			 plc_vector *r1p, plc_vector *r2p, plc_vector *r3p)

/* Implementation of the randomized "double rotation" operation of
   Hinson, Varela, Diao, Arsuaga */

{
  double nR,height;
  plc_vector R,Xp,Yp,O = {{0,0,0}};
  
  R = plc_vect_sum(r1,plc_vect_sum(r2,r3)); /* r1 + r2 + r3 */
  nR = plc_norm(R);
  height = (nR < 1.0) ? 1.0 : (3 - nR*nR)/(2*nR);
  
  Yp = plc_vect_sum(plc_subsphere_sample(r,R,height),R);
  Xp = plc_sphere_intersection_sample(r,O,Yp);
  
  *r1p = Xp;
  *r2p = plc_vect_diff(Yp,Xp);
  *r3p = plc_vect_diff(R,Yp);
 
}

void plc_single_rotation(gsl_rng *r,
			 plc_vector r1,plc_vector r2,
			 plc_vector *r1p, plc_vector *r2p) {

  /* An implementation of the "hedgehog move" of spinning r1, r2 around r1 + r2 */

  double theta_samp;
  plc_vector r_axis;
  bool ok;
  double pi = 3.14159265358979;

  r_axis = plc_normalize_vect(plc_vect_sum(r1,r2),&ok);

  if (!ok) { /* r2 = -r1, no need to rotate */

    *r1p = r1; *r2p = r2;

  } else {

    theta_samp = gsl_ran_flat(r,0,2*pi);
    
    *r1p = plc_rotate_vect(r1,r_axis,theta_samp);
    *r2p = plc_rotate_vect(r2,r_axis,theta_samp);
    
  }
  
}

plCurve *plc_random_equilateral_closed_polygon(gsl_rng *r,int nEdges)

/* An implementation of the "fast ergodic double rotation" algorithm of 
   Varela, Hinson, Diao, Arsuaga, J. Phys. A. 42 (2009) */

{
  if (nEdges < 2 || nEdges > 100000000) {

    fprintf(stderr,
	    "plc_random_equilateral_closed_polygon: Can't "
	    "generate polygon with %d edges.\n",nEdges);
    exit(1);
    
  }
  
  plc_vector *edgeset;
  int i;
  int *indexset;

  edgeset = calloc(nEdges,sizeof(plc_vector));
  indexset = calloc(nEdges,sizeof(int));

  if (edgeset == NULL || indexset == NULL) { 

    fprintf(stderr,"plc_random_equilateral_closed_polygon: "
	    "Can't allocate %d edge vectors or integers.\n",nEdges);
    exit(1);
    
  } 
  
  for(i=0;i<nEdges;i++) { indexset[i] = i; }
  
  /* Now we start work. */

  if (nEdges % 2 == 1) { /* Odd case; start by generating triangle. */

    plc_vector Y;

    gsl_ran_dir_3d(r,&(edgeset[0].c[0]),&(edgeset[0].c[1]),&(edgeset[0].c[2]));
    Y = plc_sphere_intersection_sample(r,plc_build_vect(0,0,0),edgeset[0]);
    edgeset[1] = plc_vect_diff(Y,edgeset[0]);
    edgeset[2] = plc_scale_vect(-1.0,Y);
    i = 3;

  } else { /* Even case; start by generating 4-gon */
    
    plc_vector r1,r2;

    gsl_ran_dir_3d(r,&(r1.c[0]),&(r1.c[1]),&(r1.c[2]));
    gsl_ran_dir_3d(r,&(r2.c[0]),&(r2.c[1]),&(r2.c[2]));
    plc_single_rotation(r,r1,r2,&edgeset[0],&edgeset[1]);

    edgeset[2] = plc_scale_vect(-1.0,r1);
    edgeset[3] = plc_scale_vect(-1.0,r2);
    i = 4;

  }

  for(;i<nEdges;i+=2) {  /* Now generate pairs of edges */

    int rc[2];
    plc_vector rV;

    gsl_ran_dir_3d(r,&(rV.c[0]),&(rV.c[1]),&(rV.c[2]));
    gsl_ran_choose(r,rc,2,indexset,i,sizeof(int)); /* Choose 2 previous */

    /* Now double-rotate those two and the new rV, replacing the old 2 and 
       adding the new guy to the edgelist. This preserves sum of these 3. */

    plc_double_rotation(r,rV,edgeset[rc[0]],edgeset[rc[1]],
			&edgeset[i],&edgeset[rc[0]],&edgeset[rc[1]]);

    /* To keep the sum of the whole thing equal to zero, add -rV to the list */

    edgeset[i+1] = plc_scale_vect(-1.0,rV);

  }

  gsl_ran_shuffle(r,edgeset,nEdges,sizeof(plc_vector));

  /* Now turn these edges into a polygon. */

  plCurve *L;
  bool open = { false };
  int  cc = { 0 };
  int  nv = nEdges;

  L = plc_new(1,&nv,&open,&cc);
  L->cp[0].vt[0] = plc_build_vect(0,0,0);

  for(i=1;i<nEdges;i++) {
    
    L->cp[0].vt[i] = plc_vect_sum(L->cp[0].vt[i-1],
				  plc_scale_vect(2.0/nEdges,edgeset[i-1]));

  }

  plc_fix_wrap(L);  /* As we've modified vertex buffer! */

  /* Step 4. Housekeeping and return. */

  free(edgeset);
  free(indexset);

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

  plc_fix_wrap(L); /* Because we modified the vertex buffer! */

  return L;

}

/************ Loop closure algorithm *************/


plCurve *plc_random_ftc_internal(gsl_rng *r, int nEdges, double length, double ftc, bool selfcheck)

/* Generate open arcs with a specified failure to close (along the x-axis) and length. */

{
  if (length < ftc) {

    return NULL;  /* Can't generate any arc whose ftc > length! */

  } 

  if (length < ftc + 1e-8) { length = ftc; } /* The almost straight case; make straight. */

  if (nEdges == 1 && fabs(length-ftc) < 1e-10) {  
    
    /* There is a straight arm with this ftc, but the method below won't generate it. */
    /* In this case, we do it manually. */

    plCurve *closedL;
    int nv = 2, cc = 0;
    bool open = true;

    closedL = plc_new(1,&nv,&open,&cc);
    closedL->cp[0].vt[0] = plc_build_vect(0,0,0);
    closedL->cp[0].vt[1] = plc_build_vect(length,0,0);

    plc_fix_wrap(closedL);
    return closedL;

  }

  if (nEdges == 1 && fabs(length - ftc) > 1e-10) {

    return NULL;

  }
    
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

      fprintf(stderr,"plc_random_ftc_internal_selfcheck: <A,A> = %g + %g i != 1.0\n",creal(aa),cimag(aa));
      exit(1);

    }

    if (fabs(creal(bb) - 1.0) > 1e-10 || fabs(cimag(bb)) > 1e-10) {

      fprintf(stderr,"plc_random_ftc_internal_selfcheck: <A,A> = %g + %g i != 1.0\n",creal(bb),cimag(bb));
      exit(1);

    }

    if (fabs(creal(ab)) > 1e-10 || fabs(cimag(ab)) > 1e-10) {

      fprintf(stderr,"plc_random_ftc_internal_selfcheck: <A,B> = %g + %g i != 0.0\n",creal(ab),cimag(ab));
      exit(1);

    }

  }

  /* Now apply the failure-to-close scaling. */
  
  double ell;
  double aScale,bScale;
  
  ell = 2.0*(ftc/length);  /* This is the failure to close for a polygon of length 2.0 */
  aScale = sqrt(1 + ell/2.0);
  bScale = sqrt(1 - ell/2.0);
  
  for(i=0;i<nEdges;i++) {  A[i] *= aScale; B[i] *= bScale; }
  
  /* Now we can check ftc, edge sum */

  if (selfcheck) {

    plc_vector edgesum = {{0,0,0}};

    for(i=0;i<nEdges;i++) {

      plc_M_add_vect(edgesum,hopfImap(creal(A[i]),cimag(A[i]),creal(B[i]),cimag(B[i])));

    }

    if (fabs(plc_M_norm(edgesum) - ell) > 1e-8) { 

      fprintf(stderr,"plc_random_ftc_internal_selfcheck: Sum of edges is (%g,%g,%g) with norm %g != %g\n",
	      plc_M_clist(edgesum),plc_M_norm(edgesum),ell);
      exit(1);

    }

    if (fabs(edgesum.c[0] - ell) > 1e-8 || fabs(edgesum.c[1]) > 1e-8 || fabs(edgesum.c[2]) > 1e-8) {

       fprintf(stderr,"plc_random_ftc_internal_selfcheck: Sum of edges is (%g,%g,%g) not in positive x direction.\n",
	      plc_M_clist(edgesum));
      exit(1);

    } 

  }

  /* 6. Assemble Polygon. */

  bool open={true};
  int cc=0,nv = nEdges+1;
  plCurve *L;
  double scale = length/2.0;

  L = plc_new(1,&nv,&open,&cc);
  L->cp[0].vt[0] = plc_build_vect(0,0,0);
  
  for(i=1;i<nv;i++) {
    L->cp[0].vt[i] = plc_vect_sum(L->cp[0].vt[i-1],
				  plc_scale_vect(scale,hopfImap(creal(A[i-1]),cimag(A[i-1]),creal(B[i-1]),cimag(B[i-1]))));
  }

  plc_fix_wrap(L);

  /* 7. Cleanup memory. */

  free(Araw); free(Braw);
  free(A); free(B);

  return L;

} 

plCurve *plc_loop_closure(gsl_rng *r,int cp,plCurve *openL,int nEdges)

/* Generate a random closure for the component cp of the open curve openL with nEdges (total). */

{
  plCurve *closureArc;
  plCurve *closedL;
  double   closureLength,ftc;
  double   pi = 3.141592653589793;

  if (!openL->cp[cp].open) { 

    fprintf(stderr,"plc_loop_closure: Can't close a closed curve.\n");
    return NULL;

  }
  
  closureLength = 2.0 - plc_arclength(openL,NULL);
  ftc = plc_distance(openL->cp[cp].vt[0],openL->cp[cp].vt[openL->cp[cp].nv-1]);
  
  if (closureLength < 0) {

    fprintf(stderr,"plc_loop_closure: Length of open curve is >= 2.0. Can't generate length 2 closure.\n");
    return NULL;

  }

  if (closureLength < ftc) {

    fprintf(stderr,"plc_loop_closure: Closure length %g too small to close gap of %g.\n",
	    closureLength,ftc);
    return NULL;

  }

  if (nEdges - (openL->cp[cp].nv - 1) < 1) { 

    fprintf(stderr,"plc_loop_closure: Total number of edges in closed loop %d smaller than edges in open section (%d) + 1.\n",
	    nEdges,openL->cp[cp].nv);
    return NULL;

  }

  closureArc = plc_random_ftc_internal(r,nEdges - (openL->cp[cp].nv-1),closureLength,ftc,false);

  /* There are two cases. If the original curve is already (numerically) closed, there is no operation here. */
  /* Otherwise, we'll need to rotate and translate the closureArc into place. */

  bool ok;
  plc_vector desiredFTCdir = plc_normalize_vect(plc_vect_diff(openL->cp[cp].vt[0],openL->cp[cp].vt[openL->cp[cp].nv-1]),&ok);

  if (ok) {

    plc_vector r_axis,x_axis = {{1,0,0}};
    r_axis = plc_normalize_vect(plc_vect_sum(desiredFTCdir,x_axis),&ok);
    if (!ok) { r_axis = plc_build_vect(0,1,0); }
    plc_rotate(closureArc,r_axis,pi);

  }

  plc_translate(closureArc,openL->cp[cp].vt[openL->cp[cp].nv-1]);

  /* Now we need to splice the two curves together. */

  closedL = plc_copy(openL);
  plc_drop_component(closedL,cp);

  plc_vector *newVerts;
  int vt,ca_vt,clr;

  newVerts = calloc(nEdges,sizeof(plc_vector));
  assert(newVerts != NULL);

  for(vt=0;vt<openL->cp[cp].nv;vt++) { newVerts[vt] = openL->cp[cp].vt[vt]; }
  for(ca_vt=1;ca_vt<closureArc->cp[0].nv-1 && vt < nEdges;ca_vt++,vt++) { newVerts[vt] = closureArc->cp[cp].vt[ca_vt]; }

  assert(ca_vt == closureArc->cp[0].nv-1 && vt == nEdges); /* These should finish at the same time */

  plc_color *cbuf = NULL;

  if (openL->cp[cp].cc != 0) {
  
    cbuf = calloc(nEdges,sizeof(plc_color));
    for(clr=0;clr<openL->cp[cp].cc;clr++) { 
      cbuf[clr] = openL->cp[cp].clr[clr];
    }

  }

  int cc = 0;
  if (openL->cp[cp].cc > 0) {
    cc = (openL->cp[cp].cc == 1) ? 1 : nEdges;
  }

  plc_add_component(closedL,cp,nEdges,false,cc,newVerts,cbuf);
  plc_fix_wrap(closedL);

  /* Now we do some housekeeping */

  free(newVerts); if (cbuf != NULL) { free(cbuf); } 
  plc_free(closureArc);

  return closedL;
   
}

    
    
  

   


     

