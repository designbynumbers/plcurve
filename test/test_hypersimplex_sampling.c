/*
  
  test_hypersimplex_sampling: Unit tests for the code in hypersimplex_sampling.c.

*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
#endif

#ifdef HAVE_STDIO_H
   #include<stdio.h>
#endif

#ifdef HAVE_MATH_H
   #include<math.h>
#endif

#ifdef HAVE_ASSERT_H
   #include<assert.h>
#endif

#ifdef HAVE_STRING_H
   #include<string.h>
#endif

#ifdef HAVE_STDINT_H
   #include<stdint.h>
#endif

#ifdef HAVE_STDLIB_H
   #include<stdlib.h>
#endif

#ifdef HAVE_STDBOOL_H
   #include<stdbool.h>
#endif

#ifdef HAVE_ARGTABLE2_H
  #include<argtable2.h>
#endif

#ifdef HAVE_GSL_GSL_RNG_H
#include <gsl/gsl_rng.h>
#endif

#ifdef HAVE_GSL_GSL_RANDIST_H
#include <gsl/gsl_randist.h>
#endif


/* Most of the functions in hypersimplex_sampling.c aren't exposed, so we 
   include prototypes here. */

double *psi_inverse(int n,double *y);
double *hypersimplex_sample(int k, int n,gsl_rng *rng);
double *hypercube_slice_sample(int n,gsl_rng *rng);

int PD_VERBOSE=50;

double *psi(int n, double *x)
{
  double *y = malloc((n-1)*sizeof(double));
  double accum = 0;
  int i;
  
  for(i=0;i<n-1;i++) {

    accum += x[i];
    y[i] = accum - floor(accum);

  }

  return y;
}

bool psi_inverse_test()
{ 

  printf("----------------------------------------------\n"
	 "testing psi_inverse code                      \n"
	 "----------------------------------------------\n");

  printf("testing the case (1/3,5/7,5/13,2/10,3/5)...");

  double test1[5] = {1.0/3.0,5.0/7.0,5.0/13.0,2.0/10.0,3.0/5.0};
  double *x, *y;

  y = psi(6,test1);
  x = psi_inverse(6,y);

  int i;
  for(i=0;i<5;i++) {

    if (fabs(x[i] - test1[i]) > 1e-7) {

      printf("fail at x[%d] == %g != %g\n",i,x[i],test1[i]);
      printf("original values   : %1.8f %1.8f %1.8f %1.8f %1.8f\n",
	     test1[0],test1[1],test1[2],test1[3],test1[4]);
      printf("psi         values: %1.8f %1.8f %1.8f %1.8f %1.8f\n",
	     y[0],y[1],y[2],y[3],y[4]);
      printf("psi_inv(psi) vals : %1.8f %1.8f %1.8f %1.8f %1.8f\n",
	     x[0],x[1],x[2],x[3],x[4]);

      return false;
    }
  }

  free(x);
  free(y);

  x = psi_inverse(6,test1);
  y = psi(6,x);

  for(i=0;i<5;i++) {

    if (fabs(y[i] - test1[i]) > 1e-7) {

      printf("fail\n");
      printf("original values   : %1.8f %1.8f %1.8f %1.8f %1.8f\n",
	     test1[0],test1[1],test1[2],test1[3],test1[4]);
      printf("psi_inv     values: %1.8f %1.8f %1.8f %1.8f %1.8f\n",
	     x[0],x[1],x[2],x[3],x[4]);	 
      printf("psi(psi_inv) vals : %1.8f %1.8f %1.8f %1.8f %1.8f\n",
	     y[0],y[1],y[2],y[3],y[4]);

      return false;
    }
  }

  free(x);
  free(y);
  
  printf("pass (psi_inv(psi) == psi(psi_inv) == orig)\n");

  const gsl_rng_type * rng_T;
  gsl_rng *rng;
  
  gsl_rng_env_setup();  
  rng_T = gsl_rng_default;
  rng = gsl_rng_alloc(rng_T);
  
  int seedi;
  seedi = time(0); 
  gsl_rng_set(rng,seedi);

  printf("testing 100 random inputs of length 367\n");   
  printf("with %s random number gen, seeded with %d.\n",
	 gsl_rng_name(rng),seedi);
  
  double *xin = malloc(367*sizeof(double));
  for(i=0;i<367;i++) {
    xin[i] = gsl_rng_uniform_pos(rng);
  }
  
  y = psi(368,xin);
  x = psi_inverse(368,y);

  for(i=0;i<367;i++) {

    if (fabs(x[i] - xin[i]) > 1e-7) {

      printf("fail\n");
      int idx[5];
      if (i >= 2 && i < 365) {
	idx[0] = i-2; idx[1] = i-1; idx[2] = i; idx[3] = i+1; idx[4] = i+2; 
      } else if (i < 2) {
	idx[0] = i; idx[1] = i+1; idx[2] = i+2; idx[3] = i+3; idx[4] = i+4; 
      } else {
	idx[0] = i-4; idx[1] = i-3; idx[2] = i-2; idx[3] = i-1; idx[4] = i; 
      }

      printf("values differ at index %d\n",i);

      printf("indices           : %8d %8d %8d %8d %8d\n",
	     idx[0],idx[1],idx[2],idx[3],idx[4]);      
      printf("original values   : %1.8f %1.8f %1.8f %1.8f %1.8f\n",
	     xin[idx[0]],xin[idx[1]],xin[idx[2]],
	     xin[idx[3]],xin[idx[4]]);
      printf("psi         values: %1.8f %1.8f %1.8f %1.8f %1.8f\n",
	     y[idx[0]],y[idx[1]],y[idx[2]],y[idx[3]],y[idx[4]]);
      printf("psi_inv(psi) vals : %1.8f %1.8f %1.8f %1.8f %1.8f\n",
	     x[idx[0]],x[idx[1]],x[idx[2]],x[idx[3]],x[idx[4]]);

      return false;
      
    }
  }

  free(x); free(y);

  x = psi_inverse(368,xin);
  y = psi(368,x);

  for(i=0;i<367;i++) {

    if (fabs(y[i] - xin[i]) > 1e-7) {

      printf("fail\n");
      int idx[5];
      if (i >= 2 && i < 365) {
	idx[0] = i-2; idx[1] = i-1; idx[2] = i; idx[3] = i+1; idx[4] = i+2; 
      } else if (i < 2) {
	idx[0] = i; idx[1] = i+1; idx[2] = i+2; idx[3] = i+3; idx[4] = i+4; 
      } else {
	idx[0] = i-4; idx[1] = i-3; idx[2] = i-2; idx[3] = i-1; idx[4] = i; 
      }

      printf("values differ at index %d\n",i);

      printf("indices           : %8d %8d %8d %8d %8d\n",
	     idx[0],idx[1],idx[2],idx[3],idx[4]);      
      printf("original values   : %1.8f %1.8f %1.8f %1.8f %1.8f\n",
	     xin[idx[0]],xin[idx[1]],xin[idx[2]],
	     xin[idx[3]],xin[idx[4]]);
      printf("psi         values: %1.8f %1.8f %1.8f %1.8f %1.8f\n",
	     x[idx[0]],x[idx[1]],x[idx[2]],x[idx[3]],x[idx[4]]);
      printf("psi_inv(psi) vals : %1.8f %1.8f %1.8f %1.8f %1.8f\n",
	     y[idx[0]],y[idx[1]],y[idx[2]],y[idx[3]],y[idx[4]]);

      return false;
      
    }
  }

  printf("pass (100 random vectors checked)\n");
  
 
  printf("----------------------------------------------\n"
	 "testing psi_inverse_code: pass  \n"
	 "----------------------------------------------\n");

  gsl_rng_free(rng);

  return true;
 
}


bool hypersimplex_test() {

  printf("---------------------------------------------\n"
	 "testing sample_hypersimplex code             \n"
	 "---------------------------------------------\n");

  clock_t start,end;
  double cpu_time_used;

  const gsl_rng_type * rng_T;
  gsl_rng *rng;
  
  gsl_rng_env_setup();  
  rng_T = gsl_rng_default;
  rng = gsl_rng_alloc(rng_T);
  
  int seedi;
  seedi = time(0); 
  gsl_rng_set(rng,seedi);

  printf("testing sampling from hypersimplices     \n");   
  printf("with %s random number gen, seeded with %d.\n",
	 gsl_rng_name(rng),seedi);

  printf("generating 10 (5,10) hypersimplex samples for inspection...\n");

  int i,j;
  for(i=0;i<10;i++) {

    double *x;
    x = hypersimplex_sample(5,10,rng);

    printf("\t(");

    for(j=0;j<9;j++) {

      printf("%1.5g%s",x[j],
	     j < 8 ? "," : ")\n");

    }

    double total = 0;
    for(j=0;j<9;j++) {
      total += x[j];
    }

    for(j=0;j<9;j++) {
      if (x[j] < 0 || x[j] > 1.0) {
	printf("FAIL (element x[%d] = %g, not in [0,1])\n",
	       j,x[j]);
	return false;
      }
    }

    if (total < 4.0 || total > 5.0) {

      printf("FAIL (sum %g is not in (4,5))\n",total);
      return false;

    }

    free(x);

  }

  printf("generating samples in all (k,39) hypersimplices...");
  
  int k;
  for(k=2;k<37;k++) {

    int i;
    for(i=0;i<1000;i++) {

      double *x;
      x = hypersimplex_sample(k,39,rng);

      double total = 0;
      for(j=0;j<38;j++) {
	total += x[j];
      }

      for(j=0;j<38;j++) {
	if (x[j] < 0 || x[j] > 1.0) {
	  
	  printf("FAIL on sample %d from (%d,%d)th hypersimplex\n"
		 "element x[%d] = %g, not in [0,1]\n",
		 i,k,39,
		 j,x[j]);
	  return false;
	}
      }
      
      if (total < k-1 || total > k) {

	printf("FAIL on sample %d from (%d,%d)th hypersimplex\n"
	       "FAIL (sum %g is not in (%d,%d))\n",
	       i,k,39,total,k-1,k);
	       
	return false;
	
      }

      free(x);

    }

  }

  printf("pass (entries in-range, sums correct)\n");
  
  printf("performance testing for the (500,1000) hypersimplex...");

  start = clock();

  for(i=0;i<1000;i++) {

    double *x;
    x = hypersimplex_sample(500,1000,rng);
    free(x);

  }

  end = clock();
  cpu_time_used = (((double)(end - start))/CLOCKS_PER_SEC)/1000.0;

  printf("%g sec (for this system)\n",cpu_time_used);    
  gsl_rng_free(rng);

  printf("---------------------------------------------\n"
	 "testing hypersimplex_sample: pass            \n"
	 "---------------------------------------------\n");

  return true;

}

/* double *hypercube_slice_sample(int n,gsl_rng *rng) */


bool hypercube_test() {

  printf("---------------------------------------\n"
	 "test hypercube_slice_sample            \n"
	 "---------------------------------------\n");
  
  clock_t start,end;
  double cpu_time_used;

  const gsl_rng_type * rng_T;
  gsl_rng *rng;
  
  gsl_rng_env_setup();  
  rng_T = gsl_rng_default;
  rng = gsl_rng_alloc(rng_T);
  
  int seedi;
  seedi = time(0); 
  gsl_rng_set(rng,seedi);

  printf("testing generation hypercube slice samples\n");   
  printf("with %s random number gen, seeded with %d.\n",
	 gsl_rng_name(rng),seedi);

  printf("generating 10 slice samples from the 10-cube for inspection...\n");

  int i,j;
  for(i=0;i<10;i++) {

    double *x;
    x = hypercube_slice_sample(10,rng);

    printf("\t(");

    for(j=0;j<10;j++) {

      printf("%1.5g%s",x[j],
	     j < 9 ? "," : ")\n");

    }

    double total = 0;
    for(j=0;j<10;j++) {
      total += x[j];
    }

    for(j=0;j<10;j++) {
      if (x[j] < -1.0 || x[j] > 1.0) {
	printf("FAIL (element x[%d] = %g, not in [-1.0,1.0])\n",
	       j,x[j]);
	return false;
      }
    }

    if (fabs(total) > 1e-8) {

      printf("FAIL (sum %g != 0.0)\n",total);
      return false;

    }

    free(x);

  }

  printf("generating samples in slices of cubes from n=3 to n=39...");
  
  int n;
  for(n=3;n<39;n++) {

    int i;
    for(i=0;i<1000;i++) {

      double *x;
      x = hypercube_slice_sample(n,rng);

      double total = 0;
      for(j=0;j<n;j++) {
	total += x[j];
      }

      for(j=0;j<n;j++) {
	if (x[j] < -1.0 || x[j] > 1.0) {
	  
	  printf("FAIL on sample %d from hypercube %d\n"
		 "element x[%d] = %g, not in [-1.0,1.0]\n",
		 i,n,
		 j,x[j]);
	  return false;
	}
      }
      
      if (fabs(total) > 1e-8) {

	printf("FAIL on sample %d from hypercube %d\n"
	       "sum %g != 0.0\n",
	       i,n,total);
	       
	return false;
	
      }

      free(x);

    }

  }

  printf("pass (entries in-range, sums correct)\n");
  
  printf("performance testing for the 1000-cube...");

  start = clock();

  for(i=0;i<100;i++) {

    double *x;
    x = hypercube_slice_sample(1000,rng);
    free(x);

  }

  end = clock();
  cpu_time_used = (((double)(end - start))/CLOCKS_PER_SEC)/100.0;

  printf("%g sec (for this system)\n",cpu_time_used);    
  printf("performance testing for the 1001-cube...");

  start = clock();

  for(i=0;i<100;i++) {

    double *x;
    x = hypercube_slice_sample(1001,rng);
    free(x);

  }

  end = clock();
  cpu_time_used = (((double)(end - start))/CLOCKS_PER_SEC)/100.0;

  printf("%g sec (for this system)\n",cpu_time_used);    
  gsl_rng_free(rng);

  printf("---------------------------------------------\n"
	 "testing hypercube_slice_sample: pass            \n"
	 "---------------------------------------------\n");

  return true;

}
  

int main() {

  printf("test_hypersimplex_sampling (%s)\n",PACKAGE_STRING);
  printf("--------------------------------------------------------\n"
	 "Unit tests for sampling hypersimplex and hypercube slice. \n"
	 "========================================================\n");

  if (!psi_inverse_test() || !hypersimplex_test() || !hypercube_test()) {

    printf("=======================================================\n");
    printf("test_hypersimplex_sampling:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================================\n");
    printf("test_hypersimplex_sampling:  PASS.\n");
    exit(0);

  }

  return 0;

}
