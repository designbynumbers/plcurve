/* 
 * test_findpoca: This program is part of the liboctrope test suite. In
 *                particular, it attempts to exercise all the cases in the
 *                find_poca routines. 
 */

/* Copyright 2004 The University of Georgia. */

/* This file is part of liboctrope.
   
liboctrope is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

liboctrope is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with liboctrope; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#ifdef HAVE_MATH_H
#include <math.h>
#endif
#ifdef HAVE_TIME_H
#include <time.h>
#endif
#ifdef HAVE_ASSERT_H
#include <assert.h>
#endif


#include "octrope.h"

int NUM_PASSED = 0;
int gPar = 0;

/* These functions are in the octrope library, but not visible in octrope.h */

void find_poca(plc_vector edgeA[2], plc_vector edgeB[2],
	       double *s, double *t, double *dist);  

bool octrope_solveMatrixLU(double A[2][2],double b[2],double x[2]);


#define max(A,B) (A) > (B) ? (A) : (B)
#define min(A,B) (A) < (B) ? (A) : (B)

#define __DISTANCE(s,t,n2v1v2,b0,b1,a00,a01,a11)  (n2v1v2 - 2*(s)*(b0) - 2*(t)*(b1) + (s)*(s)*a00 + 2*(s)*(t)*a01 + (t)*(t)*a11)

void original_find_poca(plc_vector edgeA[2], plc_vector edgeB[2],
               double *s, double *t, double *dist)

{
  plc_vector d[2],v1v2;
  double b[2],A[2][2];
  double det,n2v1v2;
  double sprime;
  double Aends[2],Bends[2],iends[2],avg;
  plc_vector temp_vect1,temp_vect2;
  
  /*
   * We first prepare and solve the matrix equation to find the critical point
   * of the distance squared functional. 
   *
   */ 
  
  d[0].c[0] = edgeA[1].c[0] - edgeA[0].c[0];
  d[1].c[0] = edgeB[1].c[0] - edgeB[0].c[0];  
  v1v2.c[0] = edgeA[0].c[0] - edgeB[0].c[0];

  d[0].c[1] = edgeA[1].c[1] - edgeA[0].c[1];
  d[1].c[1] = edgeB[1].c[1] - edgeB[0].c[1];  
  v1v2.c[1] = edgeA[0].c[1] - edgeB[0].c[1];

  d[0].c[2] = edgeA[1].c[2] - edgeA[0].c[2];
  d[1].c[2] = edgeB[1].c[2] - edgeB[0].c[2];  
  v1v2.c[2] = edgeA[0].c[2] - edgeB[0].c[2];

  n2v1v2 = plc_M_dot(v1v2,v1v2);
  
  b[0] = -plc_M_dot(v1v2,d[0]);
  b[1] =  plc_M_dot(v1v2,d[1]);

  A[0][0] = plc_M_dot(d[0],d[0]);
  A[1][0] = A[0][1] = -plc_M_dot(d[0],d[1]);
  A[1][1] = plc_M_dot(d[1],d[1]);
  
  det = A[0][0] * A[1][1] - A[0][1]*A[1][0];
  
  if (det > 1e-10 || det < -1e-10) { /* The simple solution will work fine */
    
    det = 1/det;
    *t = det * (-A[1][0]*b[0] + A[0][0]*b[1]);
    *s = det * (A[1][1]*b[0] - A[0][1]*b[1]);
    
  } else {  /* Use LU decomposition to solve the matrix (more stable) */

    double x[2];

    if (octrope_solveMatrixLU(A,b,x)) {
      
      /* Now set s and t from the matrix solution */
      *s = x[0]; *t = x[1];
    
    } else {
      
      /* The matrix is singular, meaning that the lines are parallel. */
      /* We write the line segments as intervals in the parallel direction. */
      
      Aends[0] = 0;
      Aends[1] = sqrt(A[0][0]);     /* Length of edgeA */
      
      Bends[0] = b[0]/Aends[1];
      Bends[1] = Bends[0] - A[0][1]/Aends[1];
      /* A[0][1] = - <d[0],d[1]>, Aends[1] = |d[0]| */
      
      /* Now we compute the intersection and take the midpoint */
      
      iends[0] = max(Aends[0],min(Bends[0],Bends[1]));
      iends[1] = min(Aends[1],max(Bends[0],Bends[1]));
      
      avg = (iends[0] + iends[1])/2;
      
      *t = (avg - Bends[0])/(Bends[1] - Bends[0]);
      *s = avg/Aends[1];  /* Expressed in s coordinates. */
      /* Deal with no overlap conditions */
      if (*s > 1) { *s = 1; } else if (*s < 0) { *s = 0; }
      if (*t > 1) { *t = 1; } else if (*t < 0) { *t = 0; }
      
      temp_vect1 = d[0];
      plc_M_scale_vect(*s,temp_vect1);
      temp_vect2 = d[1];
      plc_M_scale_vect(*t,temp_vect2);
      plc_M_sub_vect(temp_vect1,temp_vect2);
      plc_M_add_vect(temp_vect1,v1v2);
      *dist = plc_M_norm(temp_vect1);
            
    }

  }   

  /* Now we check to see whether s and t are in range. */
  
  if (*s < 0 || *s > 1) { /* s is out of range */ 

    if (*t < 0 || *t > 1) { /* Both are out of range */

      if (*s < 0) { *s = 0; } else { *s = 1; }
      if (*t < 0) { *t = 0; } else { *t = 1; }

      sprime = (b[0] - (*t) * A[0][1])/A[0][0];  

      if (sprime < 0) { sprime = 0; } else if (sprime > 1) { sprime = 1; }
      if (sprime != *s) { 
	/*
	 * The minimum along the S-side is not at the corner where it 
	 * meets the T-side in question, so we have found the actual
	 * minimum. 
	 *
	 */ 
	*s = sprime;
      } else {
	/* 
	 * The minimum along the S-side is at the corner.  Thus the minimum
	 * will be found at the minimum along the T-side.
	 *
	 */
	*t = (b[1] - (*s) * A[0][1])/A[1][1];
	if (*t < 0) { *t = 0; } else if (*t > 1) { *t = 1; }
      }
    } else { /* t is in but s is out */
      if (*s < 0) { 
	*s = 0; 
	*t = b[1]/A[1][1];
	if (*t < 0) { *t = 0; } else if (*t > 1) { *t = 1; }
      } else {  /* we know that s > 1 range if we reach this spot */
	*s = 1; 
	*t = (b[1] - A[0][1])/A[1][1]; 
	if (*t < 0) { *t = 0; } else if (*t > 1) { *t = 1; }
      }
    }
  } else { /* s is in range */
    if (*t < 0) { 
      *t = 0; 
      *s = b[0]/A[0][0];
      if (*s < 0) { *s = 0; } else if (*s > 1) { *s = 1; }
    } else if (*t > 1) { 
      *t = 1; 
      *s = (b[0] - A[0][1])/A[0][0];  
      if (*s < 0) { *s = 0; } else if (*s > 1) { *s = 1; }
    } /* otherwise, both were in range, leave them alone */
  }
  
  *dist = sqrt(__DISTANCE(*s,*t,n2v1v2,b[0],b[1],A[0][0],A[0][1],A[1][1]));
  
} /* find_poca */




void test_case(plc_vector edgeA[2], plc_vector edgeB[2], 
	       double s_expected, double t_expected, double dist_expected)

     /* Tests this case, exits with code 1 and prints error message if fails,
	returns silently if passed. */

{
  double s,t,dist;

  find_poca(edgeA,edgeB,&s,&t,&dist);

  if (octrope_error_num != 0) { /* should never happen but will stop the
                                   compiler warnings */
    fprintf(stderr,"test_findpoca: %s",octrope_error_str);
    exit(-1);
  }

  if (s_expected > 0 && t_expected > 0) {

    if (fabs(s - s_expected) > 1e-10 || 
	fabs(t - t_expected) > 1e-10 || 
	fabs(dist - dist_expected) > 1e-10) {
      
      fprintf(stderr,
	      "FAILED test_findpoca case.\n"
	      "EdgeA                               EdgeB            \n"
	      "(% 02.5f,% 02.5f,%02.5f)      (% 02.5f,% 02.5f,% 02.5f)    \n"
	      "(% 02.5f,% 02.5f,%02.5f)      (% 02.5f,% 02.5f,% 02.5f)    \n"
	      "\n"
	      "s returned: %5g\t t returned: %5g\t dist returned: %5g\n"
	      "s expected: %5g\t t expected: %5g\t dist expected: %5g\n",
	      edgeA[0].c[0],edgeA[0].c[1],edgeA[0].c[2],
	      edgeB[0].c[0],edgeB[0].c[1],edgeB[0].c[2],
	      edgeA[1].c[0],edgeA[1].c[1],edgeA[1].c[2],
	      edgeB[1].c[0],edgeB[1].c[1],edgeB[1].c[2],
	      s,t,dist,s_expected,t_expected,dist_expected);
      exit(1);
      
    }
    
  } else {

    /* Just test distance. */
    
    if (fabs(dist - dist_expected) > 1e-9) {
      
      fprintf(stderr,
	      "FAILED parallel test_findpoca case.\n"
	      "EdgeA                                 EdgeB            \n"
	      "(% 02.5f,% 02.5f,% 02.5f)      (% 02.5f,% 02.5f,% 02.5f)    \n"
	      "(% 02.5f,% 02.5f,% 02.5f)      (% 02.5f,% 02.5f,% 02.5f)    \n"
	      "\n"
	      "dist returned: %.16g\n"
	      "dist expected: %.16g\n",
	      edgeA[0].c[0],edgeA[0].c[1],edgeA[0].c[2],
	      edgeB[0].c[0],edgeB[0].c[1],edgeB[0].c[2],
	      edgeA[1].c[0],edgeA[1].c[1],edgeA[1].c[2],
	      edgeB[1].c[0],edgeB[1].c[1],edgeB[1].c[2],
	      dist,dist_expected);
      exit(1);
      
    }
   
  }

  NUM_PASSED++;

}
  
void poca_ordinary_correct(int times_to_run) {

  plc_vector dirA = {{1,0,0}}, dirB = {{1,1,0}};
  plc_vector startA = {{0,0,0}}, startB = {{0,0,1}};
  plc_vector edgeA[2], edgeB[2];
  
  int i;
  

  /* There's nothing particularly intelligent about this code. 
     We list a bunch of cases for which we know the solutions, 
     and check the results of find_poca. */

  for(i=0;i<times_to_run;i++) {
    
    /* s,t in range. */
    
    edgeA[0] = plc_vect_sum(startA,plc_scale_vect(-1,dirA));
    edgeA[1] = plc_vect_sum(startA,plc_scale_vect(1,dirA));
    
    edgeB[0] = plc_vect_sum(startB,plc_scale_vect(-1,dirB));
    edgeB[1] = plc_vect_sum(startB,plc_scale_vect(1,dirB));
    
    /* We expect the closest point to be at s=0.5, t=0.5, dist = 1.*/
    
    test_case(edgeA,edgeB,0.5,0.5,1);
    
    /* s in range, t = 0. */
    
    edgeB[0] = plc_vect_sum(startB,plc_scale_vect(0.5,dirB));
    edgeB[1] = plc_vect_sum(startB,plc_scale_vect(3,dirB));
    
    test_case(edgeA,edgeB,0.75,0,sqrt(1.25)); 
    
    /* s in range, t = 1. */
    
    edgeB[0] = plc_vect_sum(startB,plc_scale_vect(-3,dirB));
    edgeB[1] = plc_vect_sum(startB,plc_scale_vect(-0.5,dirB));
    
    test_case(edgeA,edgeB,0.25,1,sqrt(1.25)); 
    
    /* t in range, s = 0. */
    
    edgeB[0] = plc_vect_sum(startB,plc_scale_vect(-1,dirB));
    edgeB[1] = plc_vect_sum(startB,plc_scale_vect(1,dirB));
    
    edgeA[0] = plc_vect_sum(startA,plc_scale_vect(0.5,dirA));
    edgeA[1] = plc_vect_sum(startA,plc_scale_vect(3,dirA));
    
    test_case(edgeA,edgeB,0,0.625,sqrt(9.0/8.0)); 
    
    /* t in range, s = 1. */
    
    edgeA[0] = plc_vect_sum(startA,plc_scale_vect(-3,dirA));
    edgeA[1] = plc_vect_sum(startA,plc_scale_vect(-0.5,dirA));
    
    test_case(edgeA,edgeB,1,0.375,sqrt(9.0/8.0)); 
    
    /* We now test the corners. */
    
    /* s = 0, t = 0 */
    
    edgeB[0] = plc_vect_sum(startB,plc_scale_vect(0.5,dirB));
    edgeB[1] = plc_vect_sum(startB,plc_scale_vect(3,dirB));
    
    edgeA[0] = plc_vect_sum(startA,plc_scale_vect(0.5,dirA));
    edgeA[1] = plc_vect_sum(startA,plc_scale_vect(3,dirA));
    
    test_case(edgeA,edgeB,0,0,sqrt(1.25)); 
    
    /* s = 1, t = 1 */
    
    edgeB[0] = plc_vect_sum(startB,plc_scale_vect(-3,dirB));
    edgeB[1] = plc_vect_sum(startB,plc_scale_vect(-0.5,dirB));
    
    edgeA[0] = plc_vect_sum(startA,plc_scale_vect(-3,dirA));
    edgeA[1] = plc_vect_sum(startA,plc_scale_vect(-0.5,dirA));
    
    test_case(edgeA,edgeB,1,1,sqrt(1.25)); 
    
    /* s = 0, t = 1 */
    
    /* These cases require a little different geometry. */
    
    edgeB[0] = plc_vect_sum(startB,plc_scale_vect(-3,dirB));
    edgeB[1] = plc_vect_sum(startB,plc_scale_vect(-0.5,dirB));
    
    edgeA[0] = plc_vect_sum(startA,plc_scale_vect(0,dirA));
    edgeA[1] = plc_vect_sum(startA,plc_scale_vect(1,dirA));
    
    test_case(edgeA,edgeB,0,1,sqrt(1.5)); 
    
    /* s = 1, t = 0 */
    
    edgeB[0] = plc_vect_sum(startB,plc_scale_vect(0.5,dirB));
    edgeB[1] = plc_vect_sum(startB,plc_scale_vect(3,dirB));
    
    edgeA[0] = plc_vect_sum(startA,plc_scale_vect(-1,dirA));
    edgeA[1] = plc_vect_sum(startA,plc_scale_vect(0,dirA));
    
    test_case(edgeA,edgeB,1,0,sqrt(1.5)); 
    
    /* (-8,0,2)--(-1,0,2), (0,-8,0)--(0,-1,0), expect s=1, t=1 */
    edgeA[0].c[0] = -8; edgeA[0].c[1] = 0; edgeA[0].c[2] = 2;
    edgeA[1].c[0] = -1; edgeA[1].c[1] = 0; edgeA[1].c[2] = 2;
    edgeB[0].c[0] = 0; edgeB[0].c[1] = -8; edgeB[0].c[2] = 0;
    edgeB[1].c[0] = 0; edgeB[1].c[1] = -1; edgeB[1].c[2] = 0;
    
    test_case(edgeA,edgeB,1,1,sqrt(6));

    /* Planar test case */

    edgeA[0] = plc_build_vect(1,0,0);
    edgeA[1] = plc_build_vect(2,0,0);
    edgeB[0] = plc_build_vect(0,1,0);
    edgeB[1] = plc_build_vect(0,2,0);

    test_case(edgeA,edgeB,0,0,sqrt(2));

    /* Check from randomness */

    edgeA[0] = plc_build_vect( 0.00000, 0.00000, 0.00000); edgeB[0] = plc_build_vect( 10.00000, 0.00000, 0.00000);   
    edgeA[1] = plc_build_vect( 0.53385, 0.21041, 0.81898); edgeB[1] = plc_build_vect( 10.53385, 0.21041, 0.81898);

    double s,t,dist;
    original_find_poca(edgeA,edgeB,&s,&t,&dist);
    test_case(edgeA,edgeB,-1,-1,dist);

  }
}

void poca_parallel_correct() {

  plc_vector dirA = {{1,0,0}};
  plc_vector startA = {{0,0,0}}, startB = {{0,0,1}};
  plc_vector edgeA[2], edgeB[2];

   /* The two segments overlap. */

  edgeB[0] = plc_vect_sum(startB,plc_scale_vect(0,dirA));
  edgeB[1] = plc_vect_sum(startB,plc_scale_vect(2,dirA));

  edgeA[0] = plc_vect_sum(startA,plc_scale_vect(-1,dirA));
  edgeA[1] = plc_vect_sum(startA,plc_scale_vect(1,dirA));
    
  test_case(edgeA,edgeB,0.75,0.25,1.0); 

  /* Edge B is to the right of edge A. */

  edgeB[0] = plc_vect_sum(startB,plc_scale_vect(2,dirA));
  edgeB[1] = plc_vect_sum(startB,plc_scale_vect(3,dirA));

  edgeA[0] = plc_vect_sum(startA,plc_scale_vect(0,dirA));
  edgeA[1] = plc_vect_sum(startA,plc_scale_vect(1,dirA));
    
  test_case(edgeA,edgeB,1,0,sqrt(2)); 

  /* Edge B is to the left of edge A. */

  edgeB[0] = plc_vect_sum(startB,plc_scale_vect(0,dirA));
  edgeB[1] = plc_vect_sum(startB,plc_scale_vect(1,dirA));

  edgeA[0] = plc_vect_sum(startA,plc_scale_vect(2,dirA));
  edgeA[1] = plc_vect_sum(startA,plc_scale_vect(3,dirA));
    
  test_case(edgeA,edgeB,0,1,sqrt(2)); 

  printf("\t all precoded segment pair poca tests PASSED.\n");

}

void   poca_random_correct(int ntests)

/* Tests the findpoca code against the reference original_find_poca implementation over ntests pairs of random vectors. */

{
  plc_vector A[2],B[2];
  int i;
  
  A[0] = plc_build_vect(0,0,0);
  B[0] = plc_build_vect(10,0,0);

  double s,t,dist;
  double ts,tt,tdist;

  srand(177); // seed the random number generator

  for(i=0;i<ntests;i++) {

    A[1] = plc_vect_sum(A[0],plc_random_vect());
    B[1] = plc_vect_sum(B[0],plc_random_vect());

    original_find_poca(A,B,&s,&t,&dist);
    test_case(A,B,s,t,dist);
  
  }
  
  printf("%d random segment pair poca cases PASSED.\n",ntests);
  
  for(i=0;i<ntests;i++) {

    A[1] = plc_vect_sum(A[0],plc_random_vect());
    B[1] = plc_vect_sum(B[0],A[1]);

    original_find_poca(A,B,&s,&t,&dist);
    test_case(A,B,-1,-1,dist); // we don't check the s and t values, because there can be a range of options
  
  }
  
  printf("%d random parallel segment pair poca cases PASSED.\n",ntests);

  int findpoca_wins=0,original_find_poca_wins=0;
  double max_fpwin = 0, max_ofpwin = 0;

  for(i=0;i<ntests;i++) {

    A[1] = plc_vect_sum(A[0],plc_build_vect(0,1,0));
    B[1] = plc_vect_sum(B[0],plc_vect_sum(A[1],plc_scale_vect(1e-7,plc_random_vect())));

    original_find_poca(A,B,&s,&t,&dist);
    find_poca(A,B,&ts,&tt,&tdist);
    
    if (fabs(tdist-dist) > 1e-11) {
      if (tdist < dist) { 
	if (fabs(tdist - dist) > max_fpwin) { max_fpwin = fabs(tdist - dist);}
	findpoca_wins++; 
      } else {
	if (fabs(tdist - dist) > max_ofpwin) { max_ofpwin = fabs(tdist - dist);}
	original_find_poca_wins++;
      }

    }
  }
  
  printf("find_poca wins %d times (by <= %g) and original_find_poca wins %d times (by <= %g).\n",
	 findpoca_wins,max_fpwin,original_find_poca_wins,max_ofpwin);
  printf("%d random parallel + noise segment pair poca cases PASSED.\n",ntests);

  int contcheck_fails = 0;

  A[1] = plc_vect_sum(A[0],plc_build_vect(0,1,0));
  B[1] = plc_vect_sum(B[0],plc_vect_sum(A[1],plc_scale_vect(1e-7,plc_random_vect())));
  
  find_poca(A,B,&ts,&tt,&tdist);

  plc_vector oldB = B[1];

  for(i=0;i<ntests;i++) {

    B[1] = plc_vect_sum(B[1],plc_scale_vect(1e-10,plc_random_vect()));
    find_poca(A,B,&s,&t,&dist);
    
    if (fabs(tdist-dist) > 4e-10) {

      contcheck_fails++;

    }
    
    tdist = dist;
    oldB = B[1];
    ts = s;
    tt = t;

  }

  if (contcheck_fails > 0) {

    printf("%d distance function continuity violations on small variations of nearly parallel edges.\n",
	   contcheck_fails);
    
    exit(1);
    
  } else {
    
    printf("%d continuity checks PASSED.\n",ntests);
    
  } 

}

void performance_test() 
{

  int ntests = 10000; 
  int nreps = 1000;
  int i,j;
  clock_t timeused,starttime;
  double  performance;
  double s,t,dist;
  plc_vector A[2] = {{{0,0,0}},{{0,1,0}}},B[2] = {{{0,0,10}},{{0,1,0}}};

  fprintf(stderr,"Now making a performance test of find_poca.\n");

  srand(177); // seed the randomnumber generator

  starttime = clock();

  for(i=0;i<ntests;i++) {

    A[1] = plc_vect_sum(A[0],plc_random_vect());
    B[1] = plc_vect_sum(B[0],plc_random_vect());

    for(j=0;j<nreps;j++) {
      
      original_find_poca(A,B,&s,&t,&dist);
      
    }

  }
 
  timeused = clock();

  performance = (double)(ntests*nreps)*(double)(CLOCKS_PER_SEC)/
    (double)(timeused - starttime);

  fprintf(stderr,"original_find_poca yields approximately %g pocas per second on this machine.\n"
	  ,performance);

  starttime = clock();

  for(i=0;i<ntests;i++) {

    A[1] = plc_vect_sum(A[0],plc_random_vect());
    B[1] = plc_vect_sum(B[0],plc_random_vect());

    for(j=0;j<nreps;j++) {
      
      find_poca(A,B,&s,&t,&dist);
      
    }

  }
 
  timeused = clock();

  performance = (double)(ntests*nreps)*(double)(CLOCKS_PER_SEC)/
    (double)(timeused - starttime);

  fprintf(stderr,"find_poca yields approximately %g pocas per second on this machine.\n"
	  ,performance);
  
}


int main() {

  poca_ordinary_correct(1);
  //poca_parallel_correct();
  poca_random_correct(1000000);
  performance_test();

  exit(0);

}


