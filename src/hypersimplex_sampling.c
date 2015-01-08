/*

  hypersimplex_sampling.c : This file provides an efficient way to
                            sample the (k,n)th hypersimplex using
                            Stanley's triangulation and the methods in
                            kascentpermutation.c, plus the GSL random
                            number generators.

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

/* We will link this in from kascentpermutation.c */

void random_k_ascent_permutation(int n,int k,
				 int *perm,
				 gsl_rng *rng);

double *psi_inverse(int n,double *y)
/* Given a vector of points  y (strictly) in (0,1)^{n-1},
   perform the Stanley psi^{-1} map to generate a (new)
   vector x of (n-1) doubles (allocated internally). */
{
  double *x = malloc((n-1)*sizeof(double));
  assert(x != NULL);
  int i;

  x[0] = y[0];
  
  for(i=1;i<n-1;i++) {

    if (y[i-1] < y[i]) {

      x[i] = y[i] - y[i-1];

    } else {

      x[i] = y[i] - y[i-1] + 1.0;

    }

  }
  return x;
}

double *hypersimplex_sample(int k, int n,gsl_rng *rng)
/* Use Stanley triangulation to sample randomly from the 
   (k,n)th hypersimplex: the subpolytope of [0,1]^{n-1}
   defined by the inequalities

   k-1 <= Sum x_i <= k.

   Allocates a buffer of n-1 doubles to hold the result,
   this is the caller's responsibility to free. */
{
  int *perm = malloc(n*sizeof(int));
  assert(perm != NULL);

  /* Our first job is to generate a sample from the standard
     (permutation-ordered) simplex in [0,1]^{n-1}. Then we'll
     transform it by permutation and by psi_inverse. 

     To do this, we start by sampling the standard simplex
     in R^n using the Dirichlet distribution sampler in GSL;
     then we'll take a linear transform of this point to the
     standard simplex, then map that by psi_inverse. 

     Consulting the docs, we see that setting all the alpha
     to 1.0 should make the distribution flat on the simplex.*/

  double *theta = malloc(n*sizeof(double));
  double *alpha = malloc(n*sizeof(double));
  assert(alpha != NULL); assert(theta != NULL);

  int i;
  for(i=0;i<n;i++) { alpha[i] = 1.0; }
  gsl_ran_dirichlet(rng,n,alpha,theta);

  /* Now these theta are the coefficients in our point as a linear
     combination of the standard basis vectors e[i] in R^n. But we
     want them to be rewritten as a linear combination of the vectors
     (0,...,0), (0,0,...,1), (0,...,1,1), (1,...,1) which define the
     sorted simplex in R^{n-1}.

     It's worth doing some thinking about this before coding it.

     If we let theta[0] be the coefficient of (1,...,1) and theta[n]
     be the coefficient of (0,...,0) for convenience, note that when
     we sum the sorted vectors, the coordinates of the resulting sum
     are just partial sums of the original theta_i, from the first
     coordinate:

     (theta[0], theta[0]+theta[1], ..., theta[0] +...+ theta[n-1])

     Note that theta[n] doesn't come into it anywhere, but you don't
     lose any information since theta[n] = 1 - Sum theta[i] anyway.
     Note also that these are obviously sorted, which is a good sign.*/

  double *sorted_x = malloc((n-1)*sizeof(double));
  assert(sorted_x != NULL);

  sorted_x[0] = theta[0];
  
  for(i=1;i<n-1;i++) {
    
    sorted_x[i] = sorted_x[i-1] + theta[i];

  }

  /* Now we need to resort these guys according to a (k-1)-ascent
     permutation of (n-1) letters. */

  random_k_ascent_permutation(n-1,k-1,perm,rng);

  double *x = malloc((n-1)*sizeof(double));
  assert(x != NULL);

  for(i=0;i<n-1;i++) {
  
    x[i] = sorted_x[perm[i]];
  
  }

  /* Now we apply the psi_inverse transformation. */

  double *y = psi_inverse(n,x);

  /* Everything else is housecleaning. */

  free(theta);
  free(alpha);
  free(sorted_x);
  free(x);

  return y;
}

double *hypercube_slice_sample(int n,gsl_rng *rng)
  
/* Use hypersimplex sampling to sample the slice of [-1,1]^n where 
   
   (z_1, ... ,z_n) = Sum z_i = 0. 
   
   Returns a newly allocated buffer; it's the caller's responsibility
   to free it.

   Proof of algorithm.

   Let's call coordinates on [0,1]^{n-1} the y_i, with sum Y, and
   coordinates on [-1,1]^n the z_i with sum z_i = "sum" = 0.

   We are going to map [0,1]^{n-1} -> [-1,1]^n in the following way:

   (y_1,...,y_{n-1}) -> (2y_1 - 1, ... , 2y_{n-1} - 1) 
                     =  (z_1, ..., z_{n-1}) 

   These coordinates have sum 2Y - (n-1) = 2Y - n + 1. 

   Now if the y_i were sampled from the (k,n)th hypersimplex, we know
   the sum Y obeyed k-1 <= Y <= k.

   So the sum of the new coordinates obeys

   2(k-1) - n + 1 <= Sum z_i <= 2k - n + 1

       2k - n - 1 <= Sum z_i <= 2k - n + 1

   This means that by setting z_n = 2k - n - Sum_{i=1}^{n-1} z_i,
   we get some z_n in [-1,1] so that 

   0 = Sum_{i=1}^n z_i = 2k - n.

   That is, we want to solve this k using "sum" and "n". At this point,
   the procedure bifurcates. 

   If n is even, we can solve 

             k = n/2

   and sample from the (n/2,n) hypersimplex.    
*/     
{
  if (n%2 == 0) {

    int k = n/2;
    double *y = hypersimplex_sample(k,n,rng);
    double *z = malloc(n*sizeof(double));
    assert(z != NULL);

    int i;
    z[n-1] = 0;
    
    for(i=0;i<n;i++) {
      
      z[i] = 2*y[i] - 1.0;
      z[n-1] -= z[i];

    }

    free(y);
    return z;

  } else {

    /* If n is odd, we take one of the two hypersimplices

       (floor(n/2),n) or (ceil(n/2),n)

       and perform the same linear transformation. We're
       rejection sampling here, so we may have to repeat
       this process. 
    */

    double *z = malloc(n*sizeof(double));
    assert(z != NULL);
    
    do {

      int k = (gsl_rng_uniform_int(rng,2) == 0) ? floor(n/2) : ceil(n/2);
      double *y = hypersimplex_sample(k,n,rng);
     
      int i;
      z[n-1] = 0;
    
      for(i=0;i<n;i++) {
      
	z[i] = 2*y[i] - 1.0;
	z[n-1] -= z[i];

      }

      free(y);

    } while (z[n-1] < -1.0 || z[n-1] > 1.0);
    
    return z;

  }

}
