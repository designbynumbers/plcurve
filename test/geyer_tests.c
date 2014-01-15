/* 

   These are unit tests for the ips estimator code, just to make 
   sure that everything is working as expected. 

*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <assert.h>
#include <stdbool.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "geyer_estimators.h"

#define PASS 1
#define FAIL 0

double simple_sample_mean(double *data,int dLen)
{
  double result = 0.0;
  int i;

  for(i=0;i<dLen;i++) { 
    
    result += data[i];

  }

  return result/((double) dLen);
}
  

int sample_mean_test() {

  double *data;
  int    dLen,i;
  double result;
  
  clock_t start,end;
  double cpu_time_used;
  double simple_result;
  double simple_cpu_time_used;

  data = calloc(10000000,sizeof(double));
  assert(data != NULL);
  dLen = 10000000;

  printf("\n"
	 "sample_mean Test Suite\n"
	 "----------------------\n");

  printf("\tInitializing array...");
  for(i=0;i<dLen;i++) { data[i] = 1.0; }
  printf("done.\n");


  printf("\t10 sample identical value test...");
  result = sample_mean(data,10);
  if (fabs(result-1.0) > 1e-12) { 
    printf("FAIL\n");
    return FAIL;
  }
  printf("pass\n");

  printf("\t%d sample identical value test...",dLen);
 
  start = clock();
  result = sample_mean(data,dLen);
  if (fabs(result-1.0) > 1e-8) { 
    printf("FAIL\n");
    return FAIL;
  }
  end = clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  start = clock();
  simple_result = sample_mean(data,dLen);
  end = clock();
  simple_cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("pass. \n\t\t%4.10f (%3.3g sec) vs %4.10f (%3.3g sec) for simple version.\n",
	 result,cpu_time_used,simple_result,simple_cpu_time_used);

  
  printf("\t%d sample variable value test...",dLen);

  double correct_value;

  for(i=0;i<dLen;i++) { 
    data[i] = pow(10.0,i%5);
  }

  correct_value = (1.0 + 10.0 + 100.0 + 1000.0 + 10000.0)/5.0;
 
  start = clock();
  result = sample_mean(data,dLen);
  if (fabs(result-correct_value) > 1e-8) { 
    printf("FAIL\n");
    return FAIL;
  }
  end = clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  start = clock();
  simple_result = sample_mean(data,dLen);
  end = clock();
  simple_cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("pass.\n\t\t%4.10f (%3.3g sec) vs %4.10f (%3.3g sec) for simple version.\n",
	 result,cpu_time_used,simple_result,simple_cpu_time_used);

  free(data);
  return PASS;
}

typedef struct point_struct {

  double c[3];

} point;

void swap(double buf[2]) { 

  double hold;
  hold = buf[0]; buf[0] = buf[1]; buf[1] = hold;

}

double min(double a,double b) { 

  if (a < b) { return a; } else { return b; }

} 

double max(double a,double b) {

  if (a > b) { return a; } else { return b; }

}

point markov_step(const gsl_rng *rng,point x) {

  /* Does a hit and run step on the unit square using the random number generator T. */

  double theta;
  double PI = 3.1415926;
  point v;
  
  theta = gsl_ran_flat(rng,0,2*PI);
  v.c[0] = cos(theta); v.c[1] = sin(theta);

  double tupper,tlower;

  /* We now work out the limits on t which keep the point x + t v 
     inside the unit square. We need to solve the equations: 

     x.c[0] + t v.c[0] = 1
     x.c[0] + t v.c[0] = 0
     
     x.c[1] + t v.c[1] = 1
     x.c[1] + t v.c[1] = 0

  */

  double intervalA[2], intervalB[2];
  
  intervalA[0] = (1 - x.c[0])/v.c[0];
  intervalA[1] = (-x.c[0])/v.c[0];

  intervalB[0] = (1 - x.c[1])/v.c[1];
  intervalB[1] = (-x.c[1])/v.c[1];

  /* Now we need to take the intersection of the intervals. It will help to sort endpoints. */

  if (intervalA[1] < intervalA[0]) { swap(intervalA); }
  if (intervalB[1] < intervalB[0]) { swap(intervalB); }

  /* We don't need to worry about the disjoint case, because 0 is always in both intervals. */

  tupper = min(intervalA[1],intervalB[1]);
  tlower = max(intervalA[0],intervalB[0]);

  /* Now we just choose a point along the segment. */

  double t;
  t = gsl_ran_flat(rng,tlower,tupper);

  point ret;

  ret.c[0] = x.c[0] + t*v.c[0];
  ret.c[1] = x.c[1] + t*v.c[1];

  return ret;

}
  

int ips_sigma_estimator_test() 
{

  printf("ips_sigma_estimator Test Suite\n"
	 "------------------------------\n");

  /* The first thing we do is create a zero-meaned data set 
     where we know the lagged covariances. */

  int m = 5;
  int i;
  double *data;

  data = malloc(m*sizeof(double));
  assert(data != NULL);

  for(i=0;i<m;i++) {

    if (i%5 == 0)      { data[i] = 1.0;  }
    else if (i%5 == 1) { data[i] = 1.0;  }
    else if (i%5 == 2) { data[i] = 0.0;  }
    else if (i%5 == 3) { data[i] = -1.0; }
    else if (i%5 == 4) { data[i] = -1.0; } 

  }

  /* So the data consists of 1, 1, 0, -1, -1. 
     This is clearly zero-mean. Let's work out the covariances. 

     lag 0 (variance): 1^2 + 1^2 + 0^2 + (-1)^2 + (-1)^2 = 4.0 (and then it repeats), so we should get 4/5

     lag 1 : 1*1 + 1*0 + 0*-1 + -1*-1 = 2.0, so we should get 2/5.

     lag 2 : 1*0 + 1*-1 + 0*-1 = -1.0, so we should get -1/5

     lag 3 : 1*-1 + 1*-1 = -2.0, so we should get -2/5.

     This means that the first big gamma should be -2/5 - 4/5 = -6/5,
     and should NOT be used to compute the ips estimator.  The ips
     estimator should then be:

     sqrt(gamma0 + 2 gamma1) = sqrt(4/5 + 2 * 2/5) = sqrt(8/5).

     Let's test this.
 
  */

  double result;
  double *Gammas;
  int    N;
  double gamma0,gamma1;

  clock_t start,end;
  double cpu_time_used;

  printf("\tdata set 1,1,0,-1,-1 tests...");

  start = clock();
  result = ips_sigma_estimator_with_gammas(data,m,&gamma0,&gamma1,&Gammas,&N);
  end = clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("(%3.3g sec)\n",cpu_time_used);

  printf("\t\tgamma0 == 4/5 test...");

  if (fabs(gamma0 - 4.0/5.0) < 1e-12) { 

    printf("pass\n");

  } else {
    
    printf("FAIL (gamma0 = %4.10f) \n",gamma0);
    return FAIL;

  }

  printf("\t\tgamma1 == 2/5 test...");
  
  if (fabs(gamma1 - 2.0/5.0) < 1e-12) { 

    printf("pass\n");

  } else {
    
    printf("FAIL (gamma1 = %4.10f) \n",gamma1);
    return FAIL;

  }
    
  printf("\t\tN == 0 test...");

  if (N == 0) { 

    printf("pass\n");

  } else {

    printf("FAIL (N = %d) \n",N);
    return FAIL;

  }

  if (Gammas != NULL) { free(Gammas); Gammas = NULL; }

  printf("\t\tsigma estimate == sqrt(8/5) == %g test...",sqrt(8.0/5.0));

  if (fabs(result - sqrt(8.0/5.0)) < 1e-12) { 

     printf("pass\n");

  } else {

    printf("FAIL (sigma esimate = %4.10f) \n",result);
    return FAIL;

  }

  printf("\tdata set 1,1,0,-1,-1 tests...PASS\n\n");

  /* For this data set, we should get:

     lag 0: 1^2 + 1^2 + 0^2 + 1^2 + 1^2 + (-1)^2 + (-1)^2 + 0^2 + (-1)^2 + (-1)^2 = 8/10.
     lag 1: 1*1 + 1*0 + 0*1 + 1*1 + 1*-1 + -1*-1 + -1*0 + 0*-1 + -1*-1 = 3/10.
     
     so gamma0 = 8/10 and gamma1 = 3/10. For the Gammas,

     lag 2: 1*0 + 1*1 + 0*1 + 1*-1 + 1*-1 + -1*0 + -1*-1 + 0*-1 = 0.
     lag 3: 1*1 + 1*1 + 0*-1 + 1*-1 + 1*0 + -1*-1 + -1*-1 = 3/10.
     lag 4: 1*1 + 1*-1 + 0*-1 + 1*0 + 1*-1 + -1*-1 = 0.
     lag 5: 1*-1 + 1*-1 + 0*0 + 1*-1 * 1*-1 = -4/10.

     This means that N should equal 1, and we should get Gamma[0] =
     3/10 (Gamma[1] = -4/10 won't be included in the ips estimator).

     The overall estimate should then be sqrt of 

     gamma0 + 2 gamma1 + 2 Gamma[0] = 8/10 + 2 * 3/10 + 2 * 3/10 = 20/10 = 2.0.

     or sqrt(2).

  */

  m = 10;
  free(data);
  data = malloc(m*sizeof(double));
  assert(data != NULL);
  
  data[0] = 1; data[1] = 1; data[2] = 0; data[3] = 1; data[4] = 1;
  data[5] = -1; data[6] = -1; data[7] = 0; data[8] = -1; data[9] = -1;

  printf("\tdata set 1,1,0,1,1,-1,-1,0,-1,-1 tests...");

  start = clock();
  result = ips_sigma_estimator_with_gammas(data,m,&gamma0,&gamma1,&Gammas,&N);
  end = clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("(%3.3g sec)\n",cpu_time_used);

  printf("\t\tgamma0 == 8/10 test...");

  if (fabs(gamma0 - 8.0/10.0) < 1e-12) { 

    printf("pass\n");

  } else {
    
    printf("FAIL (gamma0 = %4.10f) \n",gamma0);
    return FAIL;

  }

  printf("\t\tgamma1 == 3/10 test...");
  
  if (fabs(gamma1 - 3.0/10.0) < 1e-12) { 

    printf("pass\n");

  } else {
    
    printf("FAIL (gamma1 = %4.10f) \n",gamma1);
    return FAIL;

  }
    
  printf("\t\tN == 1 test...");

  if (N == 1) { 

    printf("pass\n");

  } else {

    printf("FAIL (N = %d) \n",N);
    return FAIL;

  }

  printf("\t\tGamma[0] == 3/10 test...");

  if (fabs(Gammas[0] - 3.0/10.0) < 1e-12) { 

    printf("pass\n");

  } else {

    printf("FAIL (Gammas[0] == %g)\n",Gammas[0]);
    return FAIL;

  }

  printf("\t\tsigma estimate == sqrt(2) == %g test...",sqrt(2.0));

  if (fabs(result - sqrt(2.0)) < 1e-12) { 

     printf("pass\n");

  } else {

    printf("FAIL (sigma esimate = %4.10f) \n",result);
    return FAIL;

  }

  if (Gammas != NULL) { free(Gammas); Gammas = NULL; }

  printf("\tdata set 1,1,0,1,1,-1,-1,0,-1,-1 tests...PASS\n\n");

  const gsl_rng_type * T;
  gsl_rng * rng;
  
  gsl_rng_env_setup();
  
  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  
  int seed = time(0);
  gsl_rng_set(rng,seed);
    
  printf("\tSimple Markov chain tests\n"
	 "\t\t%s random number gen initialized, seeded with %d...done\n"
	 ,gsl_rng_name(rng),seed);

  /* We are now going to try to compute sample variance for some actual data. */

  m = 1000000;
  free(data);
  data = malloc(m*sizeof(double));
  assert(data != NULL);
  
  point markov_pt;
  printf("\t\tgenerating %d point data set...",m);

  start = clock();
 
  int k;
  for(k=0,markov_pt.c[0] = 0.5,markov_pt.c[1] = 0.5;
      k<m;k++) {

    markov_pt = markov_step(rng,markov_pt);
    
    double x,y;
    x = markov_pt.c[0]; y = markov_pt.c[1];

    data[k] = sin(x*x*y - 7.0*x + 5.0*y + 2.0*x*x); /* A weird function, just for testing. */

  }

  double mean = sample_mean(data,m);
  for(k=0;k<m;k++) { data[k] -= mean; }
  
  end = clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("done (%3.3g sec)\n",cpu_time_used);

  printf("\t\tcomputing ips_sigma_estimator...");

  start = clock();
  result = ips_sigma_estimator_with_gammas(data,m,&gamma0,&gamma1,&Gammas,&N);
  end = clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("complete (%3.3g sec)\n",cpu_time_used);
  printf("\t\tgamma0 = %4.4g\n",gamma0);
  printf("\t\tgamma1 = %4.4g\n",gamma1);
  printf("\t\tN = %d\n",N);
  
  printf("\t\tGammas are ");

  for(k=0;k<N;k++) { 

    printf("%3.3g ",Gammas[k]);

  }
  printf("\n");

  double expected;
  expected = gamma0 + 2*gamma1;
  for(k=0;k<N;k++) { expected += 2*Gammas[k]; }
  expected = sqrt(expected);

  printf("\t\tsigma_estimator test (should be %g)...",expected);

  if(fabs(result - expected) < 1e-12) { 

    printf("pass\n");

  } else {

    printf("FAIL (result == %g)\n",result);
    return FAIL;

  }

  printf("\tSimple Markov chain tests...pass\n");

  if (Gammas != NULL) { free(Gammas); }
  free(data);
  gsl_rng_free(rng);
       
  return PASS;

}

int sin4xy_test(int m,const gsl_rng *rng) {
  
  double *data;

  data = malloc(m*sizeof(double));
  assert(data != NULL);
  
  point markov_pt;
  printf("\ttesting at %d points of Markov chain...\n",m);
  printf("\t\tgenerating %d point sin4xy data set...",m);

  clock_t start,end;
  double  cpu_time_used;

  start = clock();
 
  int k;
  for(k=0,markov_pt.c[0] = 0.5,markov_pt.c[1] = 0.5;
      k<m;k++) {

    markov_pt = markov_step(rng,markov_pt);
    
    double x,y;
    x = markov_pt.c[0]; y = markov_pt.c[1];

    data[k] = sin(4*x*y); 

  }
  
  end = clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("done (%3.3g sec)\n",cpu_time_used);

  printf("\t\tcomputing ips_value and 95%% confidence error bars...",m);
  double ipsval,ips_error; 

  start = clock();
  ipsval = ips_value(data,m,&ips_error);
  end = clock();

  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;
  printf("done (%3.3g sec)\n",cpu_time_used);
  
  /* Do housecleaning now because we're going to terminate in the next conditional. */
  free(data);

  double correct_answer = 0.5261229310;

  printf("\t\tactual error = %3.5g, error bound = %3.5g\n",fabs(ipsval-correct_answer),ips_error);
  printf("\t\tcomparing %g +- %g to correct value of %g...",ipsval,ips_error,correct_answer);

  if (fabs(ipsval - correct_answer) < ips_error) { 

    printf("pass\n");
    printf("\ttesting at %d points of Markov chain...pass\n",m);
    return PASS;

  } else {

    printf("FAIL\n");
    printf("\ttesting at %d points of Markov chain...FAIL\n",m);
    return FAIL;

  }

}


double spherical_distance(point a,point b) {

  /* Computes the (shorter) spherical distance between a and b */
  /* Doesn't require the vectors a, b to be unit length, but should
     not be zero length */

  double d  = a.c[0]*b.c[0] + a.c[1]*b.c[1] + a.c[2]*b.c[2];
  double aa = a.c[0]*a.c[0] + a.c[1]*a.c[1] + a.c[2]*a.c[2];
  double bb = b.c[0]*b.c[0] + b.c[1]*b.c[1] + b.c[2]*b.c[2];
  double theta = acos(d/sqrt(aa*bb)); 
  
  return theta;

}

int spherical_distance_test() {

  /* Just a couple of quick checks to make sure that the 
     arithmetic above isn't doing something stupid. */

  point a,b,c,d;

  a.c[0] =  1; a.c[1] =  1; a.c[2] =  1;
  b.c[0] = -1; b.c[1] = -1; b.c[2] =  1;
  c.c[0] = -1; c.c[1] =  1; c.c[2] = -1;
  d.c[0] =  1; d.c[1] = -1; d.c[2] = -1;

  /* These are vertices of a regular tetrahedron, so the spherical
     distances should all be 2 arctan(sqrt(2)) ~= 1.910633236 */

  printf("spherical_distance Test Suite\n"
	 "-----------------------------\n");
  printf("\ttesting distances in regular tetrahedron...");
  
  double correct_answer = 1.910633236249019;

  if (fabs(spherical_distance(a,b) - correct_answer) > 1e-12) { 

    printf("FAIL. (%g vs correct %g) for points (%g,%g,%g) and (%g,%g,%g)\n",
	   spherical_distance(a,b), correct_answer,
	   a.c[0],a.c[1],a.c[2],
	   b.c[0],b.c[1],b.c[2]);
    return FAIL;

  }

  if (fabs(spherical_distance(a,c) - correct_answer) > 1e-12) { 

    printf("FAIL. (%g vs correct %g) for points (%g,%g,%g) and (%g,%g,%g)\n",
	   spherical_distance(a,c), correct_answer,
	   a.c[0],a.c[1],a.c[2],
	   c.c[0],c.c[1],c.c[2]);
    return FAIL;

  }

  if (fabs(spherical_distance(a,d) - correct_answer) > 1e-12) { 

    printf("FAIL. (%g vs correct %g) for points (%g,%g,%g) and (%g,%g,%g)\n",
	   spherical_distance(a,d), correct_answer,
	   a.c[0],a.c[1],a.c[2],
	   d.c[0],d.c[1],d.c[2]);
    return FAIL;

  }

  if (fabs(spherical_distance(b,c) - correct_answer) > 1e-12) { 

    printf("FAIL. (%g vs correct %g) for points (%g,%g,%g) and (%g,%g,%g)\n",
	   spherical_distance(b,c), correct_answer,
	   b.c[0],b.c[1],b.c[2],
	   c.c[0],c.c[1],c.c[2]);
    return FAIL;

  }

  if (fabs(spherical_distance(c,d) - correct_answer) > 1e-12) { 

    printf("FAIL. (%g vs correct %g) for points (%g,%g,%g) and (%g,%g,%g)\n",
	   spherical_distance(c,d), correct_answer,
	   c.c[0],c.c[1],c.c[2],
	   d.c[0],d.c[1],d.c[2]);
    return FAIL;

  }

  if (fabs(spherical_distance(d,b) - correct_answer) > 1e-12) { 

    printf("FAIL. (%g vs correct %g) for points (%g,%g,%g) and (%g,%g,%g)\n",
	   spherical_distance(d,b), correct_answer,
	   d.c[0],d.c[1],d.c[2],
	   b.c[0],b.c[1],b.c[2]);
    return FAIL;

  }

  printf("pass\n");

  return PASS;
}

double spherical_area(point a,point b,point c) {

  /* Computes the spherical area of the triangle given by a, b, c 
     using L'Huiller's formula: 

     tan(E/4)^2 = tan(s/2) tan((s-A)/2) tan((s-B)/2) tan((s-C)/2).

     where S = (A + B + C)/2 and A, B, C are the lengths of the sides
     opposite the points a, b, and c.
  */

  double s,A,B,C;

  A = spherical_distance(b,c);
  B = spherical_distance(a,c);
  C = spherical_distance(a,b);
  s = (A + B + C)/2.0;
  
  double E;
  E = 4*atan(sqrt(tan(s/2.0)*tan((s-A)/2.0)*tan((s-B)/2.0)*tan((s-C)/2.0)));

  return E;
}

int spherical_area_test() {

  /* Just a couple of quick checks to make sure that the 
     arithmetic above isn't doing something stupid. */

  point a,b,c,d;

  a.c[0] =  1; a.c[1] =  1; a.c[2] =  1;
  b.c[0] = -1; b.c[1] = -1; b.c[2] =  1;
  c.c[0] = -1; c.c[1] =  1; c.c[2] = -1;
  d.c[0] =  1; d.c[1] = -1; d.c[2] = -1;

  /* These are vertices of a regular tetrahedron, so the spherical
     areas should all be pi. */

  printf("spherical_area Test Suite\n"
	 "-----------------------------\n");
  printf("\ttesting face areas in regular tetrahedron...");
  
  double correct_answer = 3.141592653589793;

  if (fabs(spherical_area(a,b,c) - correct_answer) > 1e-12) { 

    printf("FAIL. (%g vs correct %g) for points (%g,%g,%g), (%g,%g,%g), (%g,%g,%g)\n",
	   spherical_area(a,b,c), correct_answer,
	   a.c[0],a.c[1],a.c[2],
	   b.c[0],b.c[1],b.c[2],
	   c.c[0],c.c[1],c.c[2]);
    return FAIL;

  }

  if (fabs(spherical_area(a,b,d) - correct_answer) > 1e-12) { 

    printf("FAIL. (%g vs correct %g) for points (%g,%g,%g), (%g,%g,%g), (%g,%g,%g)\n",
	   spherical_area(a,b,c), correct_answer,
	   a.c[0],a.c[1],a.c[2],
	   b.c[0],b.c[1],b.c[2],
	   d.c[0],d.c[1],d.c[2]);
    return FAIL;

  }
  
  if (fabs(spherical_area(a,c,d) - correct_answer) > 1e-12) { 

    printf("FAIL. (%g vs correct %g) for points (%g,%g,%g), (%g,%g,%g), (%g,%g,%g)\n",
	   spherical_area(a,c,d), correct_answer,
	   a.c[0],a.c[1],a.c[2],
	   c.c[0],c.c[1],c.c[2],
	   d.c[0],d.c[1],d.c[2]);
    return FAIL;

  }

  if (fabs(spherical_area(b,c,d) - correct_answer) > 1e-12) { 

    printf("FAIL. (%g vs correct %g) for points (%g,%g,%g), (%g,%g,%g), (%g,%g,%g)\n",
	   spherical_area(b,c,d), correct_answer,
	   b.c[0],b.c[1],b.c[2],
	   c.c[0],c.c[1],c.c[2],
	   d.c[0],d.c[1],d.c[2]);
    return FAIL;

  }

  if (fabs(spherical_area(c,b,d) - correct_answer) > 1e-12) { 

    printf("FAIL. (%g vs correct %g) for points (%g,%g,%g), (%g,%g,%g), (%g,%g,%g)\n",
	   spherical_area(c,b,d), correct_answer,
	   c.c[0],c.c[1],c.c[2],
	   b.c[0],b.c[1],b.c[2],
	   d.c[0],d.c[1],d.c[2]);
    return FAIL;

  }
  
  printf("pass\n");

  return PASS;
}

typedef struct point_s23_type {

  point p[3];

} point_s23;

point_s23 tsmcmc_s23_step(const gsl_rng *rng,point_s23 p) {

  /* Does a tsmcmc step on the product of three 2-spheres using the
     random number generator rng. 

     First, we need to assemble the moment polytope point x and do a 
     hit and run step. */

  point x;
  x.c[0] = p.p[0].c[2]; x.c[1] = p.p[1].c[2]; x.c[2] = p.p[2].c[2];

  point v;
  
  v.c[0] = gsl_ran_ugaussian(rng);
  v.c[1] = gsl_ran_ugaussian(rng);
  v.c[2] = gsl_ran_ugaussian(rng);

  double tupper,tlower;

  /* We now work out the limits on t which keep the point x + t v 
     inside the unit cube. We need to solve the equations: 

     x.c[0] + t v.c[0] = 1
     x.c[0] + t v.c[0] = 0
     
     x.c[1] + t v.c[1] = 1
     x.c[1] + t v.c[1] = 0

     x.c[2] + t v.c[2] = 1
     x.c[2] + t v.c[2] = 0

  */

  double intervalA[2], intervalB[2], intervalC[2];
  
  if (fabs(v.c[0]) > 1e-12) {
    intervalA[0] = (1 - x.c[0])/v.c[0];
    intervalA[1] = (-1 -x.c[0])/v.c[0];
  } else {
    intervalA[0] = -1e16; intervalA[1] = 1e16;
  }

  if (fabs(v.c[1]) > 1e-12) { 
    intervalB[0] = (1 - x.c[1])/v.c[1];
    intervalB[1] = (-1 -x.c[1])/v.c[1];
  } else {
    intervalB[0] = -1e16; intervalB[1] = 1e16;
  }

  if (fabs(v.c[2]) > 1e-12) { 
    intervalC[0] = (1 - x.c[2])/v.c[2];
    intervalC[1] = (-1 -x.c[2])/v.c[2];
  } else {
    intervalC[0] = -1e16; intervalC[1] = 1e16;
  }

  /* Now we need to take the intersection of the intervals. It will
     help to sort endpoints. */

  if (intervalA[1] < intervalA[0]) { swap(intervalA); }
  if (intervalB[1] < intervalB[0]) { swap(intervalB); }
  if (intervalC[1] < intervalC[0]) { swap(intervalC); }

  /* We don't need to worry about the disjoint case, because 0 is
     always in both intervals. */

  tupper = min(min(intervalA[1],intervalB[1]),intervalC[1]);
  tlower = max(max(intervalA[0],intervalB[0]),intervalC[0]);

  /* Now we just choose a point along the segment. */

  double t;
  t = gsl_ran_flat(rng,tlower,tupper);

  point xret;

  xret.c[0] = x.c[0] + t*v.c[0];
  xret.c[1] = x.c[1] + t*v.c[1];
  xret.c[2] = x.c[2] + t*v.c[2];

  /* Now that we have the new moment polytope point, we need to build
     the corresponding triple of sphere points. */

  point_s23 ret;
  double theta[3];
  int k;
  double PI = 3.141592653589793;

  for(k=0;k<3;k++){ 

    double r;
    theta[k] = gsl_ran_flat(rng,0,2*PI);
    ret.p[k].c[2] = xret.c[k]; /* The z coordinate comes from the moment polytope */

    r = sqrt(1 - ret.p[k].c[2]*ret.p[k].c[2]);
    ret.p[k].c[0] = cos(theta[k])*r;
    ret.p[k].c[1] = sin(theta[k])*r;

  }

  return ret;

}

int tsmcmc_spherical_area_test(int m,const gsl_rng *rng) {
  
  double *data;

  data = malloc(m*sizeof(double));
  assert(data != NULL);
  
  point_s23 markov_pt;
  printf("\ttesting at %d points of tsmcmc Markov chain on (S^2)^3...\n",m);
  printf("\t\tgenerating %d point spherical_area data set...",m);

  clock_t start,end;
  double  cpu_time_used;

  start = clock();
 
  int k;

  /* We start with a triangle from the regular tetrahedron */

  point a,b,c;
    
  a.c[0] =  1/sqrt(3.0); a.c[1] =  1/sqrt(3.0); a.c[2] =  1/sqrt(3.0);
  b.c[0] = -1/sqrt(3.0); b.c[1] = -1/sqrt(3.0); b.c[2] =  1/sqrt(3.0);
  c.c[0] = -1/sqrt(3.0); c.c[1] =  1/sqrt(3.0); c.c[2] = -1/sqrt(3.0);

  markov_pt.p[0] = a; markov_pt.p[1] = b; markov_pt.p[2] = c;


  for(k=0;k<m;k++) {

    markov_pt = tsmcmc_s23_step(rng,markov_pt);
    data[k] = spherical_area(markov_pt.p[0],markov_pt.p[1],markov_pt.p[2]); 

  }
  
  end = clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("done (%3.3g sec)\n",cpu_time_used);

  printf("\t\tcomputing ips_value and 95%% confidence error bars...",m);
  double ipsval,ips_error; 

  start = clock();
  ipsval = ips_value(data,m,&ips_error);
  end = clock();

  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;
  printf("done (%3.3g sec)\n",cpu_time_used);
  
  /* Do housecleaning now because we're going to terminate in the next conditional. */
  free(data);

  double correct_answer = 1.570796326794897;

  printf("\t\tactual error = %3.5g, error bound = %3.5g\n",fabs(ipsval-correct_answer),ips_error);
  printf("\t\tcomparing %g +- %g to correct value of %g...",ipsval,ips_error,correct_answer);

  if (fabs(ipsval - correct_answer) < ips_error) { 

    printf("pass\n");
    printf("\ttesting spherical area at %d points of Markov chain...pass\n",m);
    return PASS;

  } else {

    printf("FAIL\n");
    printf("\ttesting spherical area at %d points of Markov chain...FAIL\n",m);
    return FAIL;

  }

}


int ips_value_tests() {

  const gsl_rng_type * T;
  gsl_rng * rng;
  
  gsl_rng_env_setup();
  
  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  
  int seed = 1794912;
  gsl_rng_set(rng,seed);
    
  printf("ips_value tests\n"
	 "---------------\n"
	 "integrating sin(4xy) over unit square with hit-and-run...\n"
	 "%s random number gen initialized, seeded with %d...done\n"
	 ,gsl_rng_name(rng),seed);

  double correct_answer = 0.5261229310;
  printf("correct value according to Mathematica...1/4 (EulerGamma - CosIntegral[4] + Log[4]) (%0.10g)\n",correct_answer);
  
  if (sin4xy_test(1000,rng) == FAIL) {return FAIL;}
  if (sin4xy_test(10000,rng) == FAIL) {return FAIL;}
  if (sin4xy_test(100000,rng) == FAIL) {return FAIL;}
  if (sin4xy_test(1000000,rng) == FAIL) {return FAIL;}
  if (sin4xy_test(10000000,rng) == FAIL) {return FAIL;}

  printf("\nintegrating average area of spherical triangle...\n"
	 "%s random number gen initialized, seeded with %d...done\n"
	 ,gsl_rng_name(rng),seed);

  correct_answer = 1.570796326794897;

  printf("correct value...pi/2) (%0.10g)\n",correct_answer);

  if (tsmcmc_spherical_area_test(1000,rng) == FAIL) {return FAIL;}
  if (tsmcmc_spherical_area_test(10000,rng) == FAIL) {return FAIL;}
  if (tsmcmc_spherical_area_test(100000,rng) == FAIL) {return FAIL;}
  if (tsmcmc_spherical_area_test(1000000,rng) == FAIL) {return FAIL;}
  if (tsmcmc_spherical_area_test(10000000,rng) == FAIL) {return FAIL;}

  gsl_rng_free(rng);

  return PASS;
}

bool gips_mean_tests() { 

  printf("\nGeyer Initial Positive Sequence (gips) Sample Mean Tests\n"
	 "--------------------------------------------------------\n");

  int datasizes[6] = {3,999,12456,25000,500000,5000037};
  int k;

  for(k=0;k<6;k++) {
    
    printf("\n\ttesting mean value of constant sequence of %d 1.0...\n",datasizes[k]);
    
    int bufsizes[8] = {5,10,25,50,500,1000,10000,1000000};
    int size;
    
    for(size=0;size<8;size++) { 
      
      clock_t start,end;

      start = clock();
      geyer_ips_t *gips;
      gips = geyer_ips_new(bufsizes[size]); 
      
      int i;
      for(i=0;i<datasizes[k];i++) { 
	
	geyer_ips_add(gips,1.0);
	
      }
      
      bool ok;
      double val = geyer_ips_value(gips,NULL,&ok);
      
      if(fabs(val - 1.0) > 1e-12) { 
	
	printf("\t\tbufsize %d test... FAIL (value %g != 1.0)\n",gips->buffer_size,val);
	printf("FAIL\n");
	return false;
	
      } 
      
      end = clock();
      double cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

      printf("\t\tbufsize %d test... pass (%g sec, %g dot prods, %g ms/add) \n",
	     gips->buffer_size,cpu_time_used,gips->dot_product_time,1000*cpu_time_used/(double)(datasizes[k]));

      geyer_ips_free(&gips);
      geyer_ips_free(&gips);
      
    }

  }

  int harmonic_sizes[2] = {1000,1738};
  double harmonic_answers[2] = {0.007485470860550345,0.004624852491447343};

  for(k=0;k<2;k++) {

    printf("\n\ttesting average value of (1/n) on n in 1..%d...\n",harmonic_sizes[k]);

    geyer_ips_t *gips;
    gips = geyer_ips_new(134);
    int i;

    for(i=1;i<=harmonic_sizes[k];i++) {

      geyer_ips_add(gips,1.0/(double)(i));
      
    }
    
    bool ok;
    double val = geyer_ips_value(gips,NULL,&ok);

    if (fabs(val - harmonic_answers[k]) > 1e-12) { 
    
      printf("\t\tvalue %0.13f != correct answer (from Mathematica) %0.13f\n",
	     val,harmonic_answers[k]);
      printf("\tFAIL\n");
      return false;

    }

    printf("\t\tcomputed %0.17f == correct %0.17f\n\t\tpass\n",val,harmonic_answers[k]);
    geyer_ips_free(&gips);

  } 
  
  printf("\tedge case testing - compute mean of 0 data points...");
  
  geyer_ips_t *gips = geyer_ips_new(12);
  bool ok;
  double val = geyer_ips_value(gips,NULL,&ok);
  geyer_ips_free(&gips);

  if (fabs(val) < 1e-12) { 
    printf("pass\n");
  } else {
    printf("fail (mean = %g != 0) \n",val);
    return false;
  }
  
  return true;
}

double split_lagged_dot(int lag,int sizeA,double *A,int sizeB,double *B);

int split_lagged_dot_tests() {

  printf("split_lagged_dot test suite\n"
	 "-----------------------------------------\n");

  printf("vecA = (1,2,3,4,5) vecB = (1,2,3,4,5)...\n");
  double vecA[5] = {1,2,3,4,5}, vecB[5] = {1,2,3,4,5};

  printf("\tlag 1...");
  double val = split_lagged_dot(1,5,vecA,5,vecB);
  if (fabs(val-5.0) <1e-10) { 
    printf("pass (result = %g == 1*5)\n",val);
  } else { 
    printf("FAIL (result = %g != 1*5)\n",val);
    return 0;
  }

  printf("\tlag 2...");

  val = split_lagged_dot(2,5,vecA,5,vecB);
  if (fabs(val-14.0) <1e-10) { 
    printf("pass (result = %g == 4*1 + 5*2)\n",val);
  } else { 
    printf("FAIL (result = %g != 4*1 + 5*2)\n",val);
    return 0;
  }
    
  printf("\tlag 3...");

  val = split_lagged_dot(3,5,vecA,5,vecB);
  if (fabs(val-26.0) <1e-10) { 
    printf("pass (result = %g == 3*1 + 4*2 + 5*3)\n",val);
  } else { 
    printf("FAIL (result = %g != 3*1 + 4*2 + 5*3)\n",val);
    return 0;
  }

  printf("\tlag 4...");

  val = split_lagged_dot(4,5,vecA,5,vecB);
  if (fabs(val-40.0) <1e-10) { 
    printf("pass (result = %g == 2*1 + 3*2 + 4*3 + 5*4)\n",val);
  } else { 
    printf("FAIL (result = %g != 2*1 + 3*2 + 4*3 + 5*4)\n",val);
    return 0;
  }

  printf("\tlag 5...");

  val = split_lagged_dot(5,5,vecA,5,vecB);
  if (fabs(val-55.0) <1e-10) { 
    printf("pass (result = %g == 1*1 + 2*2 + 3*3 + 4*4 + 5*5)\n",val);
  } else { 
    printf("FAIL (result = %g != 1*1 + 2*2 + 3*3 + 4*4 + 5*5)\n",val);
    return 0;
  }

  printf("\tlag 5, bsize = 3...");
 
  val = split_lagged_dot(5,5,vecA,3,vecB);
  if (fabs(val-14.0) <1e-10) { 
    printf("pass (result = %g == 1*1 + 2*2 + 3*3)\n",val);
  } else { 
    printf("FAIL (result = %g != 1*1 + 2*2 + 3*3)\n",val);
    return 0;
  }

  printf("\tlag 5, asize = 4, bsize = 3...");
 
  val = split_lagged_dot(5,4,vecA,3,vecB);
  if (fabs(val-8.0) <1e-10) { 
    printf("pass (result = %g == 1*2 + 2*3)\n",val);
  } else { 
    printf("FAIL (result = %g != 1*2 + 2*3)\n",val);
    return 0;
  }

  printf("-----------------------------------------\n");

  return 1;

}

double lagged_dot(int lag,int n,double *vec);

int geyer_ips_integer_datatests(int bufsize, int DATAPOINTS) { 

  printf("\tinteger sequence test %d data points, buffer size %d...\n",DATAPOINTS,bufsize);

  geyer_ips_t *gips = geyer_ips_new(bufsize);
  double *data = calloc(DATAPOINTS,sizeof(double));

  int i;

  clock_t add_start,add_end;
  double  add_cpu_seconds_used;

  add_start = clock();

  for(i=0;i<DATAPOINTS;i++) { 

    data[i] = i;
    geyer_ips_add(gips,data[i]);

  }

  add_end = clock();
  add_cpu_seconds_used = ((double)(add_end-add_start))/CLOCKS_PER_SEC;

  double correct_val, correct_error;
  correct_val = ips_value(data,DATAPOINTS,&correct_error);

  double *zero_meaned_data;
  double correct_gamma0, correct_gamma1;
  double *correct_Gammas;
  int    correct_N;
  double correct_sigma;

  zero_meaned_data = calloc(DATAPOINTS,sizeof(double));
  for(i=0;i<DATAPOINTS;i++) { zero_meaned_data[i] = data[i] - correct_val; }

  double *lagged_dots;
  int nldots = 1000 < DATAPOINTS ? 1000 : DATAPOINTS;
  lagged_dots = calloc(nldots,sizeof(double)); /* Compute the first thousand lagged dots */
  for(i=0;i<nldots;i++) { lagged_dots[i] = lagged_dot(i,DATAPOINTS,zero_meaned_data); }
  
  correct_sigma = ips_sigma_estimator_with_gammas(zero_meaned_data,DATAPOINTS,
						  &correct_gamma0, &correct_gamma1,
						  &correct_Gammas, &correct_N);

  double ips_val,ips_error;
  bool ok;
  clock_t start,end;
  double  cpu_seconds_used;

  start = clock();
  ips_val = geyer_ips_value(gips,&ips_error,&ok);
  end = clock();
  cpu_seconds_used = ((double)(end - start))/CLOCKS_PER_SEC;

  if (fabs(ips_val - correct_val) > 1e-10) { 
    
    printf("\t\tips sample mean = %.13g != sample_mean = %.13g\n\tFAIL\n",ips_val,correct_val);
    return FAIL;

  } else {

    printf("\t\tchecking sample mean...pass (%.13g ~ %.13g)\n",ips_val,correct_val);

  }

  if (ok) {

    /* Now check some internals. */
    
    printf("\t\tused N = %d Gammas (of %d ldots) to compute sigma = %.13g \n",
	   gips->N,gips->nldots,gips->sigma);
    
    /* There's a subtle issue here: ips_sigma_estimator_with_gammas doesn't compute Gamma[0],
       or count it in correct_Gammas. So actually correct_N is one LESS than the correct value
       of gips->N. */
    
    if (gips->N != correct_N+1) { 
      
      printf("\t\tnumber of Gammas used = %d != correct value = %d\nFAIL\n",gips->N,correct_N+1);
      return FAIL;
      
    } else {
      
      printf("\t\tchecking number of Gammas used...pass (%d == %d)\n",gips->N,correct_N+1);
      
    }
    
    if (gips->N > 0) {
      
      printf("\t\tchecking Gammas...");
      
      if (fabs(correct_gamma0 + correct_gamma1 - gips->Gamma[0]) > 1e-8) {
	
	printf("Gammas[0] = %0.13g != correct Gamma[0] = %0.13g\nFAIL\n",
	       gips->Gamma[0],correct_gamma0 + correct_gamma1);
	return FAIL;
	
      } 
      
      for(i=1;i<gips->N;i++) {
	
	if (fabs(gips->Gamma[i]-correct_Gammas[i-1]) > 1e-8) { 
	  
	  printf("Gammas[%d] = %0.13g != correct Gamma[%d] = %0.13g\nFAIL\n",
		 i,gips->Gamma[i],i,correct_Gammas[i-1]);
	  return FAIL;
	  
	} 
	
      }
      
      printf("pass (all %d Gammas match)\n",gips->N);
      
    }
    
    if (fabs(gips->sigma - correct_sigma) > 1e-8 || isnan(gips->sigma)) { 
      
      printf("\t\tsigma = %0.13g != correct sigma = %0.13g\n",gips->sigma,correct_sigma);
      return FAIL;
      
    } else {
      
      printf("\t\tchecking sigma...pass (%0.13g == %0.13g)\n",gips->sigma,correct_sigma);
      
    }
    
    if (fabs(ips_error - correct_error) > 1e-8 || isnan(ips_error)) { 
      
      printf("\t\tips error = %.13g != correct error = %.13g\n\tFAIL\n",ips_error,correct_error);
      return FAIL;
      
    } else {
      
      printf("\t\tchecking error estimate...pass (%0.13g == %0.13g)\n",ips_error,correct_error);
      
    }

  } else {

    printf("\t\tran out of ldots while computing initial positive sequence...pass\n");

  }

  free(data); free(zero_meaned_data); free(correct_Gammas); free(lagged_dots);
  geyer_ips_free(&gips);

  printf("\t...pass (%g sec to add data, %g ms/item, %g sec to compute estimator)\n",
	 add_cpu_seconds_used,1000*add_cpu_seconds_used/DATAPOINTS,cpu_seconds_used);
  return PASS;

}
  

int geyer_ips_random_datatests(const gsl_rng *rng,int bufsize,int DATAPOINTS) {

  printf("\trandom data test %d data points, buffer size %d...\n",DATAPOINTS,bufsize);

  geyer_ips_t *gips = geyer_ips_new(bufsize);
  
  double *data = calloc(DATAPOINTS,sizeof(double));
  int i;
  
  for(i=0;i<DATAPOINTS;i++) { 

    data[i] = gsl_ran_flat(rng,0,1);
    geyer_ips_add(gips,data[i]);

  }

  double correct_val, correct_error;
  correct_val = ips_value(data,DATAPOINTS,&correct_error);

  double *zero_meaned_data;
  double correct_gamma0, correct_gamma1;
  double *correct_Gammas;
  int    correct_N;
  double correct_sigma;

  zero_meaned_data = calloc(DATAPOINTS,sizeof(double));
  for(i=0;i<DATAPOINTS;i++) { zero_meaned_data[i] = data[i] - correct_val; }

  double *lagged_dots;
  int nldots = 1000 < DATAPOINTS ? 1000 : DATAPOINTS;
  lagged_dots = calloc(nldots,sizeof(double)); /* Compute the first thousand lagged dots */
  for(i=0;i<nldots;i++) { lagged_dots[i] = lagged_dot(i,DATAPOINTS,zero_meaned_data); }
  
  correct_sigma = ips_sigma_estimator_with_gammas(zero_meaned_data,DATAPOINTS,
						  &correct_gamma0, &correct_gamma1,
						  &correct_Gammas, &correct_N);

  double ips_val,ips_error;
  bool ok;
  ips_val = geyer_ips_value(gips,&ips_error,&ok);

  if (fabs(ips_val - correct_val) > 1e-10) { 
    
    printf("\t\tips sample mean = %.13g != sample_mean = %.13g\n\tFAIL\n",ips_val,correct_val);
    return FAIL;

  } else {

    printf("\t\tchecking sample mean...pass (%.13g ~ %.13g)\n",ips_val,correct_val);

  }

  if (ok) { 

    /* Now check some internals. */

    printf("\t\tused N = %d Gammas (of %d ldots) to compute sigma = %.13g \n",
	   gips->N,gips->nldots,gips->sigma);
    
    
    
    /* There's a subtle issue here: ips_sigma_estimator_with_gammas doesn't compute Gamma[0],
       or count it in correct_Gammas. So actually correct_N is one LESS than the correct value
       of gips->N. */
    
    if (gips->N != correct_N+1) { 
      
      printf("\t\tnumber of Gammas used = %d != correct value = %d\nFAIL\n",gips->N,correct_N+1);
      return FAIL;
      
    } else {
      
      printf("\t\tchecking number of Gammas used...pass (%d == %d)\n",gips->N,correct_N+1);
      
    }
    
    if (gips->N > 0) {
      
      printf("\t\tchecking Gammas...");
      
      if (fabs(correct_gamma0 + correct_gamma1 - gips->Gamma[0]) > 1e-12) {
	
	printf("Gammas[0] = %0.13g != correct Gamma[0] = %0.13g\nFAIL\n",
	       gips->Gamma[0],correct_gamma0 + correct_gamma1);
	return FAIL;
	
      } 
      
      for(i=1;i<gips->N;i++) {
	
	if (fabs(gips->Gamma[i]-correct_Gammas[i-1]) > 1e-8) { 
	  
	  printf("Gammas[%d] = %0.13g != correct Gamma[%d] = %0.13g\nFAIL\n",
		 i,gips->Gamma[i],i,correct_Gammas[i-1]);
	  return FAIL;
	  
	} 
	
      }
      
      printf("pass (all %d Gammas match)\n",gips->N);
      
    }
    
    if (fabs(gips->sigma - correct_sigma) > 1e-8 || isnan(gips->sigma)) { 
      
      printf("\t\tsigma = %0.13g != correct sigma = %0.13g\n",gips->sigma,correct_sigma);
      return FAIL;
      
    } else {
      
      printf("\t\tchecking sigma...pass (%0.13g == %0.13g)\n",gips->sigma,correct_sigma);
      
    }
    
    if (fabs(ips_error - correct_error) > 1e-8 || isnan(ips_error)) { 
      
      printf("\t\tips error = %.13g != correct error = %.13g\n\tFAIL\n",ips_error,correct_error);
      return FAIL;
      
    } else {
      
      printf("\t\tchecking error estimate...pass (%0.13g == %0.13g)\n",ips_error,correct_error);
      
    }

  } else {

    printf("\t\tran out of ldots before the end of the initial positive sequence...pass\n");
    
  }
    
  free(data); free(zero_meaned_data); free(correct_Gammas); free(lagged_dots);
  geyer_ips_free(&gips);

  printf("\t...pass\n");
  return PASS;

}

int geyer_ips_forced_correlation_datatests(const gsl_rng *rng,int bufsize,int DATAPOINTS) {

  printf("\tforced correlation data test %d data points (bufsize %d)...\n",DATAPOINTS,bufsize);

  geyer_ips_t *gips = geyer_ips_new(bufsize);
  
  double *data = calloc(DATAPOINTS,sizeof(double));
  int i,j;
  clock_t add_start,add_end;
  double add_cpu_seconds_used;

  add_start = clock();
  
  for(i=0;i<DATAPOINTS;i++) { 

    data[i] = gsl_ran_flat(rng,0,1);
    if (i > 20) { for(j=1;j<10;j++) { data[i] += (1/10.0)*data[i-j]; } }
    geyer_ips_add(gips,data[i]);

  }

  add_end = clock();
  add_cpu_seconds_used = ((double)(add_end - add_start))/CLOCKS_PER_SEC;

  double correct_val, correct_error;
  correct_val = ips_value(data,DATAPOINTS,&correct_error);

  double ips_val,ips_error;
  bool ok;
  clock_t val_start,val_end;
  double val_cpu_seconds_used;

  val_start = clock();
  ips_val = geyer_ips_value(gips,&ips_error,&ok);
  val_end = clock();
  val_cpu_seconds_used = ((double)(val_end-val_start))/CLOCKS_PER_SEC;

  if (fabs(ips_val - correct_val) > 1e-10) { 
    
    printf("\t\tips sample mean = %.13g != sample_mean = %.13g\n\tFAIL\n",ips_val,correct_val);
    return FAIL;

  } else {

    printf("\t\tchecking sample mean...pass (%.13g =~ %.13g)\n",ips_val,correct_val);

  }

  if (ok) { 

    /* Now check some internals. */
    
    printf("\t\tused N = %d Gammas (%d of %d ldots) to compute sigma = %.13g \n",
	   gips->N,2*gips->N,gips->nldots,gips->sigma);
    
    double *zero_meaned_data;
    double correct_gamma0, correct_gamma1;
    double *correct_Gammas;
    int    correct_N;
    double correct_sigma;
    
    zero_meaned_data = calloc(DATAPOINTS,sizeof(double));
    for(i=0;i<DATAPOINTS;i++) { zero_meaned_data[i] = data[i] - correct_val; }
    
    correct_sigma = ips_sigma_estimator_with_gammas(zero_meaned_data,DATAPOINTS,
						    &correct_gamma0, &correct_gamma1,
						    &correct_Gammas, &correct_N);
    
    if (gips->N != correct_N+1) { 
      
      printf("\t\tnumber of Gammas used = %d != correct value = %d\nFAIL\n",gips->N,correct_N+1);
      return FAIL;
      
    } else {
      
      printf("\t\tchecking number of Gammas used...pass (%d == %d)\n",gips->N,correct_N+1);
      
    }
    
    if (gips->N > 0) {

      printf("\t\tchecking Gammas...");
      
      for(i=1;i<gips->N;i++) {
	
	if (fabs(gips->Gamma[i]-correct_Gammas[i-1]) > 1e-8) { 
	  
	  printf("Gammas[%d] = %0.13g != correct Gamma[%d] = %0.13g\nFAIL\n",
		 i,gips->Gamma[i],i,correct_Gammas[i-1]);
	  return FAIL;
	  
	} 
	
      }
      
      printf("pass (all %d Gammas match)\n",gips->N);
      
    }
    
    if (fabs(gips->sigma - correct_sigma) > 1e-8) { 
      
      printf("\t\tsigma = %0.13g != correct sigma = %0.13g\n",gips->sigma,correct_sigma);
      return FAIL;
      
    } else {
      
      printf("\t\tchecking sigma...pass (%0.13g == %0.13g)\n",gips->sigma,correct_sigma);
      
    }
    
    if (fabs(ips_error - correct_error) > 1e-8) { 
      
      printf("\t\t ips error = %.13g != correct error = %.13g\n\tFAIL\n",ips_error,correct_error);
      return FAIL;
      
    } else {

      printf("\t\tchecking error estimate...pass (%0.13g == %0.13g)\n",ips_error,correct_error);
    
    }

    free(zero_meaned_data); free(correct_Gammas);

  } else {

    printf("\t\tran out of ldots before computing estimate...inconclusive\n");

  }

  geyer_ips_free(&gips);

  printf("\t...pass (%g sec for add, %g ms/item, %g sec for val)\n",
	 add_cpu_seconds_used,1000.0*add_cpu_seconds_used/DATAPOINTS,
	 val_cpu_seconds_used);
  return PASS;

}

int geyer_ips_actual_markov_chain_datatests(const gsl_rng *rng,int bufsize,int DATAPOINTS) {
  
  double *data;
  clock_t start,end;
  double cpu_seconds_used;

  data = malloc(DATAPOINTS*sizeof(double));
  assert(data != NULL);
  
  point markov_pt;
  printf("\tactual Markov chain test with gips (bufsize %d) and %d points of Markov chain...\n",bufsize,DATAPOINTS);

  geyer_ips_t *gips = geyer_ips_new(bufsize);

  start = clock();

  int k,i;
  for(k=0,markov_pt.c[0] = 0.5,markov_pt.c[1] = 0.5;
      k<DATAPOINTS;k++) {

    markov_pt = markov_step(rng,markov_pt);
    
    double x,y;
    x = markov_pt.c[0]; y = markov_pt.c[1];

    data[k] = sin(4*x*y); 
    geyer_ips_add(gips,data[k]);

  }

  end=clock();
  cpu_seconds_used = ((double)(end-start))/CLOCKS_PER_SEC;

  printf("\t\ttime to run geyer_ips_add... %g ms/datapoint\n",1000*cpu_seconds_used/(double)(DATAPOINTS));
  
  double correct_val, correct_error;
  correct_val = ips_value(data,DATAPOINTS,&correct_error);

  double ips_val,ips_error;
  bool ok;
  ips_val = geyer_ips_value(gips,&ips_error,&ok);

  if (fabs(ips_val - correct_val) > 1e-10) { 
    
    printf("\t\tips sample mean = %.13g != sample_mean = %.13g\n\tFAIL\n",ips_val,correct_val);
    return FAIL;

  } else {

    printf("\t\tchecking sample mean...pass (%.13g =~ %.13g)\n",ips_val,correct_val);

  }

  if (ok) { 

    /* Now check some internals. */
    
    printf("\t\tused N = %d Gammas (%d of %d ldots) to compute sigma = %.13g \n",
	   gips->N,2*gips->N,gips->nldots,gips->sigma);
    
    double *zero_meaned_data;
    double correct_gamma0, correct_gamma1;
    double *correct_Gammas;
    int    correct_N;
    double correct_sigma;
    
    zero_meaned_data = calloc(DATAPOINTS,sizeof(double));
    for(i=0;i<DATAPOINTS;i++) { zero_meaned_data[i] = data[i] - correct_val; }
    
    correct_sigma = ips_sigma_estimator_with_gammas(zero_meaned_data,DATAPOINTS,
						    &correct_gamma0, &correct_gamma1,
						    &correct_Gammas, &correct_N);
     
    if (gips->N != correct_N+1) { 
      
      printf("\t\tnumber of Gammas used = %d != correct value = %d\nFAIL\n",gips->N,correct_N+1);
      return FAIL;
      
    } else {
      
      printf("\t\tchecking number of Gammas used...pass (%d == %d)\n",gips->N,correct_N+1);
      
    }
    
    if (gips->N > 0) {
      
      printf("\t\tchecking Gammas...");
      
      for(i=1;i<gips->N;i++) {
	
	if (fabs(gips->Gamma[i]-correct_Gammas[i-1]) > 1e-8) { 
	  
	  printf("Gammas[%d] = %0.13g != correct Gamma[%d] = %0.13g\nFAIL\n",
		 i,gips->Gamma[i],i,correct_Gammas[i-1]);
	  return FAIL;
	  
	} 
	
      }
      
      printf("pass (all %d Gammas match)\n",gips->N);
      
    }
    
    if (fabs(gips->sigma - correct_sigma) > 1e-8) { 
      
      printf("\t\tsigma = %0.13g != correct sigma = %0.13g\n",gips->sigma,correct_sigma);
      return FAIL;
      
    } else {
      
      printf("\t\tchecking sigma...pass (%0.13g == %0.13g)\n",gips->sigma,correct_sigma);
      
    }
    
    if (fabs(ips_error - correct_error) > 1e-8) { 
    
      printf("\t\t ips error = %.13g != correct error = %.13g\n\tFAIL\n",ips_error,correct_error);
      return FAIL;
    
    } else {

      printf("\t\tchecking error estimate...pass (%0.13g == %0.13g)\n",ips_error,correct_error);
      
    }

    free(data); free(zero_meaned_data); free(correct_Gammas);

  } else {

    printf("\t\tran out of ldots before ips ended...inconclusive\n");

  }

  geyer_ips_free(&gips);

  printf("\t...pass\n");
  return PASS;

}

int geyer_ips_tests(const gsl_rng *rng) {
  
  printf("Geyer ips error calculation test suite\n"
	 "--------------------------------------\n");

  int bufsizes[6] = {10,300,5,500,5000,500000};

  int sizes[8] = {10,450,10,25,50,75,513,
		   12768,15346};
  int i,j;

  for(j=0;j<6;j++) {

    for(i=0;i<8;i++) {

      if (geyer_ips_integer_datatests(bufsizes[j],sizes[i]) == FAIL) { 
	
	printf("\tFAIL\n");
	return FAIL;
	
      } else {
	
	printf("\tpass\n");
	
      }
      
    }

  }

  for(j=0;j<6;j++) {
  
    for(i=0;i<8;i++) {
  
      if (geyer_ips_random_datatests(rng,bufsizes[j],sizes[i]) == FAIL) {

	printf("\tFAIL\n");
	return FAIL;
	
      } else {
	
	printf("\tpass\n");
	
      }
      
    }   

  }

  return PASS;

  for(j=0;j<6;j++) {
   
    for(i=0;i<8;i++) {
  
      if (geyer_ips_forced_correlation_datatests(rng,bufsizes[j],sizes[i]) == FAIL) {

	printf("\tFAIL\n");
	return FAIL;
	
      } else {
	
	printf("\tpass\n");
	
      }
      
    }   

  }

  
  if (geyer_ips_actual_markov_chain_datatests(rng,500000,2300000) == FAIL) {
    
    printf("\tFAIL\n");
    return FAIL;
    
  } else {
    
    printf("\tpass\n");
    
  }
  
  printf("-------------------------------------\n");
  return PASS;
    
}

int main()
{
  const gsl_rng_type * T;
  gsl_rng * rng;
  
  gsl_rng_env_setup();
  
  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  
  int seed = 1794912;
  gsl_rng_set(rng,seed);

  printf("Geyer IPS test suite\n");
  printf("-----------------------------------------------------\n");
  printf("with %s random number gen, seeded with %d.\n",gsl_rng_name(rng),seed);

  if (split_lagged_dot_tests() == PASS) {
  
    printf("split_lagged_dot tests... PASS\n\n");

  } else {

    printf("split_lagged_dot tests... FAIL\n\n");
    exit(1);
    
  }

  if (geyer_ips_tests(rng) == PASS) {
  
    printf("geyer_ips error calc tests... PASS\n\n");

  } else {

    printf("geyer_ips error calc tests... FAIL\n\n");
    exit(1);
    
  }

  if (gips_mean_tests() == PASS) {
  
    printf("gips mean tests... PASS\n\n");

  } else {

    printf("gips mean tests... FAIL\n\n");
    exit(1);
    
  }

  
  
  if (sample_mean_test() == PASS) {

    printf("sample_mean tests... PASS\n\n");

  } else {

    printf("sample_mean tests... FAIL\n\n");
    exit(1);

  }

  if (spherical_distance_test() == PASS) {

    printf("spherical_distance tests...PASS\n\n");

  } else {

    printf("spherical_distance tests...FAIL\n\n");
    exit(1);

  }

  if (spherical_area_test() == PASS) {

    printf("spherical_area tests...PASS\n\n");

  } else {

    printf("spherical_area tests...FAIL\n\n");
    exit(1);

  }


  if (ips_sigma_estimator_test() == PASS) {

    printf("ips_sigma_estimator tests...PASS\n\n");

  } else {

    printf("ips_sigma_estimator tests... FAIL\n\n");
    exit(1);

  }

  if (ips_value_tests() == PASS) {

    printf("ips_value tests... PASS\n\n");

  } else {

    printf("ips_value tests... FAIL\n\n");
    exit(1);

  }

  gsl_rng_free(rng);
  exit(0);

}
  
