/*
 *  Routines to generate random self-avoiding polygons as part of plCurve.
 *
 */

/* Copyright 2023 The University of Georgia. */

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
#include"plCurve.h"

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

bool plc_is_sap_internal(plCurve *L, bool verbose) {

  /* This is an internal debugging function which double-checks that all
     vertices of L are at least (almost) distance 1.0 away from each other.
     It's written for the general case where L has many components, even
     though we don't have code to generate such SAPs yet. */

  int cp1,cp2;
  int vt1,vt2;

  for(cp1=0;cp1<L->nc;cp1++) {

    for(cp2=0;cp2<=cp1;cp2++) {

      if (cp1 == cp2) {
	
	for(vt1=1;vt1<L->cp[cp1].nv;vt1++) {
	  
	  for(vt2=0;vt2<vt1;vt2++) {
	    
	    if (plc_sq_dist(L->cp[cp1].vt[vt1],L->cp[cp2].vt[vt2]) < 1.0-1e-10) {

	      if (verbose) {

		printf("plc_is_sap: Squared distance between L->cp[%d].vt[%d] and L->cp[%d].vt[%d] is %g < %g\n",
		       cp1,vt1,cp2,vt2,plc_sq_dist(L->cp[cp1].vt[vt1],L->cp[cp2].vt[vt2]),1.0-1e-10);
		       
		return false;

	      }
	      
	    }
	    
	  }
	  
	}
	
      } else {
	
	for(vt1=0;vt1<L->cp[cp1].nv;vt1++) {
	  
	  for(vt2=0;vt2<L->cp[cp2].nv;vt2++) {
	    
	     if (verbose) {

		printf("plc_is_sap: Squared distance between L->cp[%d].vt[%d] and L->cp[%d].vt[%d] is %g < %g\n",
		       cp1,vt1,cp2,vt2,plc_sq_dist(L->cp[cp1].vt[vt1],L->cp[cp2].vt[vt2]),1.0-1e-10);
		       
		return false;

	     }

	  }

	}

      }

    }

  }

  return true;

}

bool plc_is_sap(plCurve *L) {

  plc_is_sap_internal(L,false);

}


plCurve *plc_random_equilateral_closed_self_avoiding_polygon(gsl_rng *rng,int n)
/* 
   Generates random closed polygons where vertices are surrounded by a 
   disjoint spheres of radius 1/2 (the "string of pearls" model) by 
   rejection sampling. This is exponentially slow, so it's limited both 
   in number of attempts and by the max and min n which will be attempted. 

   We also do a fair amount of one-off optimization here in order to increase 
   performance. 

*/
{
  int nv = n,cc = {0};
  bool open = {false};
  plCurve *L = plc_new(1,&nv,&open,&cc);

  /* We'll copy over our sap to the buffer in L once we generate it. */

  double *X[3];  /* This will be our working buffer of vertex positions. */
 
  assert(n >= 5);
  if (n < 5) {

    printf("Error. Can't call plc_random_sap\n"
	   "with n < 5\n");
    exit(1);

  }

  assert(n <= 25);
  if (n > 25) {

    printf("Error. Can't call plc_random_sap\n"
	   "with n > 25\n");
    exit(1);
  }

  /* The design of this code is basically straight from the Mathematica notebook spaam-sampling.nb */
  /* These are buffers of vertex coordinates. In order to address the individual vectors, we use */
  /* the **seemingly backwards** notation X[k][i] to get the k-th entry of the i-th vector. */

  X[0] = malloc(n * sizeof(double));
  X[1] = malloc(n * sizeof(double));
  X[2] = malloc(n * sizeof(double));

  assert(X[0] != NULL); assert(X[1] != NULL); assert(X[2] != NULL);

  X[0][0] = 0.; X[1][0] = 0.; X[2][0] = 0.;
  X[0][1] = 1.; X[1][1] = 0.; X[2][1] = 0.;

  double frame[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};

  /* To match the memory usage in the X buffers, we'll address these "backwards" too: frame[k][i] is */
  /* the k-th coordinate of the i-th vector. */
  
  double last_diag, this_diag = 0.;
  double normal[3];
  int i,j;
  int trials=0,max_trials = 50000000;

  double TWOPI = 6.2831853071795864769;
  double cos_phi, sin_phi;
  double theta, cos_theta, sin_theta;
  double nm;

  restart_trial:

  trials++;
  assert(trials < max_trials);
  
  normal[0] = 0.; normal[1] = 0.; normal[2] = 1.;
  last_diag = 1.0;

  for(i=2;i<(n-1); i++, last_diag = this_diag) {

    this_diag = last_diag + 2.0*(gsl_rng_uniform(rng)-0.5);
    if ((this_diag + last_diag < 1.0) || (this_diag <= 0.0)) { goto restart_trial; };
    
    /* frame[-][0] = X[-][i-1]/||X[-][i-1]|| */
    
    nm = sqrt(X[0][i-1]*X[0][i-1] + X[1][i-1]*X[1][i-1] + X[2][i-1]*X[2][i-1]);
    frame[0][0] = X[0][i-1]/nm; frame[1][0] = X[1][i-1]/nm; frame[2][0] = X[2][i-1]/nm;
    
    /* frame[-][1] = normal x frame[-][0] */
      
    frame[0][1] = normal[1]*frame[2][0] - normal[2]*frame[1][0];  
    frame[1][1] = normal[2]*frame[0][0] - normal[0]*frame[2][0];
    frame[2][1] = normal[0]*frame[1][0] - normal[1]*frame[0][0];
    
    /* frame[-][1] *= 1/||f[-][1]|| */
    
    nm = sqrt(frame[0][1]*frame[0][1] + frame[1][1]*frame[1][1] + frame[2][1]*frame[2][1]);
    frame[0][1] /= nm;
    frame[1][1] /= nm;
    frame[2][1] /= nm;
    
    cos_phi = (this_diag * this_diag + last_diag * last_diag - 1.0)/(2.0 * this_diag * last_diag);
    sin_phi = sqrt(1 - cos_phi * cos_phi);
      
    /* X[-][i] = this_diag (cos_phi * frame[-][0] + sin_phi * frame[-][1]); */
    
    X[0][i] = this_diag * (cos_phi * frame[0][0] + sin_phi * frame[0][1]);
    X[1][i] = this_diag * (cos_phi * frame[1][0] + sin_phi * frame[1][1]);
    X[2][i] = this_diag * (cos_phi * frame[2][0] + sin_phi * frame[2][1]);
    
    /* frame[-][2] = -sin_phi frame[-][0] + cos_phi frame[-][1] */
    
    frame[0][2] = -sin_phi * frame[0][0] + cos_phi * frame[0][1];
    frame[1][2] = -sin_phi * frame[1][0] + cos_phi * frame[1][1];
    frame[2][2] = -sin_phi * frame[2][0] + cos_phi * frame[2][1];
    
    /* Now pick a dihedral angle theta to try. */
    
    theta = TWOPI*gsl_rng_uniform(rng);
    cos_theta = cos(theta);
    sin_theta = sin(theta);
      
    /* normal = cos(theta) normal + sin(theta) frame[-][2]; */
      
    normal[0] = cos_theta * normal[0] + sin_theta * frame[0][2];
    normal[1] = cos_theta * normal[1] + sin_theta * frame[1][2];
    normal[2] = cos_theta * normal[2] + sin_theta * frame[2][2];
    
    /* Now check distances. We go backwards because the new sphere is most likely */
    /* to collide with nearby spheres. */
    
    for(j=i-1;j>=0;j--) {

      if (((X[0][i]-X[0][j])*(X[0][i]-X[0][j]) +
	   (X[1][i]-X[1][j])*(X[1][i]-X[1][j]) +
	   (X[2][i]-X[2][j])*(X[2][i]-X[2][j])) < 1.0-(1e-9)) { goto restart_trial; }
	
    }

  }

  /* If we got to this point, we're at i = n - 1; the last vertex! */

  assert(i==n-1); 
  last_diag = this_diag;
  
  if (last_diag <= 0.0 || last_diag >= 2.0) { goto restart_trial; }
    
  /* frame[-][0] = X[-][i-1]/||X[-][i-1]|| */
  
  nm = sqrt(X[0][i-1]*X[0][i-1] + X[1][i-1]*X[1][i-1] + X[2][i-1]*X[2][i-1]);
  frame[0][0] = X[0][i-1]/nm; frame[1][0] = X[1][i-1]/nm; frame[2][0] = X[2][i-1]/nm;
  
  /* frame[-][1] = normal x frame[-][0] */
  
  frame[0][1] = normal[1]*frame[2][0] - normal[2]*frame[1][0];  
  frame[1][1] = normal[2]*frame[0][0] - normal[0]*frame[2][0];
  frame[2][1] = normal[0]*frame[1][0] - normal[1]*frame[0][0];
  
  /* frame[-][1] *= 1/||f[-][1]|| */
  
  nm = sqrt(frame[0][1]*frame[0][1] + frame[1][1]*frame[1][1] + frame[2][1]*frame[2][1]);
  frame[0][1] /= nm;
  frame[1][1] /= nm;
  frame[2][1] /= nm;
  
  /* Since we know that this_diag is 1.0, this formula simplifies... */
	
  cos_phi = last_diag/2.0;
  sin_phi = sqrt(1.0 - cos_phi * cos_phi);
  
  /* X[-][i] = this_diag (cos_phi * frame[-][0] + sin_phi * frame[-][1]); */
	
  X[0][i] = cos_phi * frame[0][0] + sin_phi * frame[0][1];
  X[1][i] = cos_phi * frame[1][0] + sin_phi * frame[1][1];
  X[2][i] = cos_phi * frame[2][0] + sin_phi * frame[2][1];
  
  /* This is the last vertex, so there's no need to update the frame, but we still */
  /* check distances. We go backwards because the new sphere is most likely */
  /* to collide with nearby spheres. */
	
  for(j=i-1;j>=0;j--) {
	  
    if (((X[0][i]-X[0][j])*(X[0][i]-X[0][j]) +
	 (X[1][i]-X[1][j])*(X[1][i]-X[1][j]) +
	 (X[2][i]-X[2][j])*(X[2][i]-X[2][j])) < 1.0-(1e-9)) { goto restart_trial; };
	
  }
	
  /* If we make it to here, it's time to copy the contents of the buffers into the plCurve data structure */

  for(i=0;i<n;i++) {
    
    L->cp[0].vt[i].c[0] = X[0][i];
    L->cp[0].vt[i].c[1] = X[1][i];
    L->cp[0].vt[i].c[2] = X[2][i];
    
  }

  plc_fix_wrap(L);

  /* Now we need to clean up the temporary buffers */

  free(X[0]); free(X[1]); free(X[2]);

  return L;

} 
