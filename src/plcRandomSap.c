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

plCurve *plc_random_equilateral_closed_self_avoiding_polygon(gsl_rng *rng,int n);
/* 
   Generates random closed polygons where vertices are surrounded by a 
   disjoint spheres of radius 1/2 (the "string of pearls" model) by 
   rejection sampling. This is exponentially slow, so it's limited both 
   in number of attempts and by the max and min n which will be attempted. 
*/
{
  int nv = n,cc = {0};
  bool open = {false};
  plCurve *L = plc_new(1,&nv,&open,&cc);

  L->cp[0].vt[0] = plc_build_vect(0,0,0);
  L->cp[0].vt[1] = plc_build_vect(1,0,0);

  double last_diag, this_diag;
  int i;
  int safety_check;
  double TWOPI = 6.2831853071795864769;
  double theta;

  plc_vector normal;

  assert(n >= 5);
  if (n < 5) {

    printf("Error. Can't call plc_random_sap\n"
	   "with n < 5\n");
    exit(1);

  }

  assert(n <= 30);
  if (n > 30) {

    printf("Error. Can't call plc_random_sap\n"
	   "with n > 30\n");
    exit(1);
  }
  
  do {

    for(i=2,distances_ok = true,last_diag=1.0,normal=plc_build_vect(0,0,1);
	i<n && distances_ok;
	i++,last_diag=this_diag) {
      
      if (i < n-1) {
	
	/* Generate the next diagonal by adding a uniform random */
	/* variate in [-1,1]. Remember that i is the 0-indexed vertex */
	/* number, so "this_diag" will be used to build vertex i. */
	
	this_diag = last_diag + (2.0*gsl_rng_uniform(rng)-1.0);
	
	/* In this case, |this_diag - last_diag| < 1.0 by construction, 
	   but we have to check that this_diag + last_diag > 1.0, and
	   that this_diag > 0. */
	
	distances_ok = (this_diag + last_diag > 1.0 && this_diag > 1e-12);
	
      } else {
	
	/* The last diagonal is always 1.0. */
	
	this_diag = 1.0;
	
	/* In this case, we already know that this_diag + last_diag > 1.0 and 
	   that this_diag > 0, but we have to explicitly check
	   |this_diag - last_diag| < 1.0 */
	
	distances_ok = (fabs(this_diag - last_diag) < 1.0);
	
      }
      
      if (distances_ok) {
	
	/* The diagonal distance is ok, so we can build the next vertex.
	   We now do that, and add it to the polygon-in-progress L. 
	   
	   At this point, we're building vt[i] using the current value 
	   of normal, the position of vt[i], and the diagonal d[i-1]. 
	   
	   Define the angle of the triangles at the vertex across from
	   the edge vt[i-1]->vt[i] to be alpha. */
	
	double cos_alpha = (last_diag*last_diag + this_diag*this_diag - 1.0)/(2.0*last_diag*this_diag);
	double sin_alpha = sqrt(1 - cos_alpha*cos_alpha);
	
	/*
	  
	    this_diag 
              |   <=== new dihedral angle applied HERE
           f2 |   vt[i]  
      f3\   ^ v /\
         \  |  /  \  1 
          \ | /    \ 
           \|/a  f1 \
            *----->--* vt[i-1]
       vt[0]    ^ 	    
                |
              last_diag

	*/

	plc_vector f1, f2, f3;
	bool ok;
	
	f1 = plc_normalize_vect(L->cp[0].vt[i-1],&ok);
	assert(ok);
	f2 = plc_normalize_vect(plc_cross_prod(normal,f1),&ok);
	assert(ok);
	
	L->cp[0].vt[i] = plc_vlincomb(this_diag*cos_alpha,f1,
				      this_diag*sin_alpha,f2);
	
	/* Now we have to rotate the normal to encode the new
	   dihedral angle. Technically, we could skip this part for
	   the last diagonal, but it would take longer to check for
	   the special case than to just do it. */
	
	f3 = plc_cross_prod(normal,plc_normalize_vect(L->cp[0].vt[i],&ok));
	assert(ok);
	f3 = plc_normalize_vect(f3,&ok);  /* This shouldn't do anything */
	
	theta = TWOPI*gsl_rng_uniform(rng);
	
	normal = plc_vlincomb(cos(theta),normal,
			      sin(theta),f3);
	normal = plc_normalize_vect(normal,&ok);
	
      }

      /* We now have to check the self-avoiding condition, which means 
         checking the distances between the new vertex and each of the 
         existing ones. */

      int j;
      for(j=0;j<i && distances_ok;j++) {

	distances_ok =
	  distances_ok && (plc_sq_distance(L->cp[0].vt[j],L->cp[0].vt[i]) >= 1.0);

      }	
      
    }
    
    safety_check++;
    
  } while (!distances_ok && safety_check < 50000);
  
  assert(safety_check < 50000); 
  assert(fabs(plc_distance(L->cp[0].vt[0],L->cp[0].vt[n-1]) - 1.0) < 1e-8);
  
  plc_fix_wrap(L);
  
  return L;
  
}


  

plCurve *fantriangulation_action_angle(int n,double *theta, double *d)
/* 
   Now we need to assemble the action-angle coordinates into a
   polygon.  We've written this code before, but in a very general
   way. Now we're going to use the fact that we know this is the fan
   triangulation in order to run this as fast as possible. 
*/
{
  
  plc_vector normal = {{0,0,1}};
  int i;

  /* The way this is going to work is that we'll retain a normal vector 
     for every triangle as we build it, updating as we apply the d[i]
     and theta[i] in sequence. There are n-1 of the d's, starting with 
     1.0 and ending with 1.0, but only n-3 of the thetas. */

  for(i=1;i<n-1;i++) {

    /* At this point, we're building vt[i+1] using the current value 
       of normal, the position of vt[i], and the diagonal d[i-1]. 
    
       Define the angles of the triangles at the vertex across from
       the edge vt[i]->vt[i+1] to be the alpha[i]; we won't store
       the history of them, but we will need the current one. */

    double cos_alpha = (d[i-1]*d[i-1] + d[i]*d[i] - 1.0)/(2.0*d[i-1]*d[i]);
    double sin_alpha = sqrt(1 - cos_alpha*cos_alpha);

    /*

             d[i] 
              |   <=== dihedral angle theta[i-1] applied HERE
           f2 |   vt[i+1]  
      f3\   ^ v /\
         \  |  /  \  1 
          \ | /    \ 
           \|/a  f1 \
            *----->--* vt[i]
       vt[0]    ^ 	    
                |
              d[i-1]

    */

    plc_vector f1, f2, f3;
    bool ok;
    
    f1 = plc_normalize_vect(L->cp[0].vt[i],&ok);
    assert(ok);
    f2 = plc_normalize_vect(plc_cross_prod(normal,f1),&ok);
    assert(ok);

    L->cp[0].vt[i+1] = plc_vlincomb(d[i]*cos_alpha,f1,
				    d[i]*sin_alpha,f2);

    /* If we're going to have another triangle, we have to 
       rotate the normal. (By definition, we're going to let 
       the dihedral angles be the angles between these normals.) */

    if (i<n-2) { 

      f3 = plc_cross_prod(normal,plc_normalize_vect(L->cp[0].vt[i+1],&ok));
      assert(ok);
      f3 = plc_normalize_vect(f3,&ok);  /* This shouldn't do anything */
    
      normal = plc_vlincomb(cos(theta[i-1]),normal,
			    sin(theta[i-1]),f3);
      normal = plc_normalize_vect(normal,&ok);

    }

  }

  assert(fabs(plc_distance(L->cp[0].vt[0],L->cp[0].vt[n-1]) - 1.0) < 1e-8);
  plc_fix_wrap(L);
  return L;
}
