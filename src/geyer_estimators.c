/* 

   This is an implementation of the initial positive sequence
   estimator for the variance of errors in a Markov chain integration
   (following Geyer). We use BLAS calls to do this, since all the
   covariance computations are basically about dot products.

   Of course, we could do a much better job if we used OpenCL or
   something, but we're not going that far in the implementation
   direction at the moment.

   This is meant as library code, so it doesn't include a main.
   Run estimator_tests.c to try some simple tests on the method.

*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>

#define TRUE  (1==1)
#define FALSE (1==0)

// Turn asserts ON.
#define DEBUG 1 

#ifdef DARWIN
  #include<Accelerate/Accelerate.h>
#else
  #include<gsl/gsl_cblas.h>
#endif

#include "geyer_estimators.h"

/* Now we enter the basic functions. */

double ips_sigma_estimator(double *data,int m) {

  return ips_sigma_estimator_with_gammas(data,m,NULL,NULL,NULL,NULL);

}

double ips_sigma_estimator_with_gammas(double *zero_meaned_data, int m, 
				       double *gamma0,
				       double *gamma1,
				       double **Gammas, 
				       int *N) {

  double *internal_Gammas;
  int     internal_Gammas_length = 256;
  int     internal_Gammas_used = 0;
  double  internal_Gamma_sum = 0;
  
  int     k;

  internal_Gammas = calloc(internal_Gammas_length,sizeof(double));
  assert(internal_Gammas != NULL);

  assert(zero_meaned_data != NULL);
  assert(m > 0);
  assert((Gammas == NULL && N == NULL) || (Gammas != NULL && N != NULL));

  bool nan_flag = false;
  for(k=0;k<m;k++) { if (isnan(zero_meaned_data[k])) { nan_flag = true; } }
  assert(nan_flag == false);

  /* We start by computing little gamma_0 (the variance) and little gamma_1
     (the lag-1 autocovariance) */

  double internal_gamma0, internal_gamma1;

  internal_gamma0 = (1.0/((double) m)) * cblas_ddot(m,zero_meaned_data,1,zero_meaned_data,1);
  internal_gamma1 = (1.0/((double) m)) * cblas_ddot(m-1,zero_meaned_data,1,zero_meaned_data+1,1);

  /* We are now going to compute big gammas we run out or find a negative one */
  /* We have to be a little careful about control flow in the for loop here,  */
  /* so we use a flag to make sure we're testing the right gamma. */

  bool done_flag = FALSE;

  for(k=1;k<m/2 && !done_flag;k++) { 

    /* 
       Notes. 
       
       1. Geyer defines Gamma_k to be the sum of the lag-2k and lag-2k+1
       autocovariance.

       2. The buffer of internal_Gammas is offset by one
       since there really isn't a Gamma_0.  

       3. Following Geyer, we are using the "biased" estimator for
       autocovariance given by dividing by the full number of samples m,
       regardless of lag. */
    
    internal_Gammas[k-1] = (1/((double) m)) * (cblas_ddot(m-2*k,zero_meaned_data,1,zero_meaned_data+2*k,1) + 
					     cblas_ddot(m-(2*k+1),zero_meaned_data,1,zero_meaned_data+2*k+1,1));

    /* Worry about running out of buffer space. */
    
    if (k > internal_Gammas_length - 10) {   

      internal_Gammas_length *= 2;
      internal_Gammas = realloc(internal_Gammas,internal_Gammas_length * sizeof(double));
      assert(internal_Gammas != NULL);
      
    }	

    if (internal_Gammas[k-1] <= 0) { /* Remember that the CURRENT Gamma is stored in k-1 */

      done_flag = TRUE;

    } else {

      internal_Gamma_sum += internal_Gammas[k-1];

    }
    
  }
  
  if (gamma0 != NULL) { *gamma0 = internal_gamma0; }
  if (gamma1 != NULL) { *gamma1 = internal_gamma1; }
  
  if (Gammas != NULL) { 
    
    *Gammas = internal_Gammas;
    *N = k-2; 
    /* Remember that k got incremented on the way out of the loop when we detected the first negative,
       AND that we wanted to throw that one out. */
    
  } else { 
    
    free(internal_Gammas);
    
  }

  /* It can happen that the sum here is negative if you add 2*internal_gamma1. In that case, 
     we are really far from convergence, but we still need to make sure that a result is 
     returned (if only for testing). We know that internal_Gamma_sum is positive, so the 
     only case where this could happen is if 2*internal_gamma1 is negative and 

     internal_gamma0 + 2 * internal_gamma1 < 0.

  */

  assert(internal_Gamma_sum >= 0);
  assert(internal_gamma0 >= 0);
  
  double ips_est = internal_gamma0 + 2*internal_gamma1 + 2*internal_Gamma_sum;

  if (ips_est >= 0) { 

    return sqrt(ips_est);

  } else {

    return sqrt(internal_gamma0);
   
  }

}

double sample_mean(double *data,int dLen)
/* This is a quick halving adder for the numbers in data. */
{
  int newLen; /* The largest power of 2 which is smaller than dLen */
  newLen = (int)(pow(2.0,floor(log2((double)(dLen))))); 
  double *newdata;
  
  newdata = malloc(newLen*sizeof(double));
  assert(newdata != NULL);

  /* First, we fold down the data into newdata */

  memcpy(newdata,data,newLen*sizeof(double));

  int i;

  for(i=0;i<dLen - newLen;i++) {

    newdata[i] += data[i+newLen];

  }
  
  /* Now we go ahead and start halving */

  for(;newLen > 1;newLen /= 2) {

    for(i=0;i<newLen/2;i++) {

      newdata[i] += newdata[i+newLen/2];

    }

  }

  /* The last value is stored in newdata[0] */

  double mean;
  mean = newdata[0]/((double) dLen);
  free(newdata);

  return mean;

}
  
double ips_value(double *raw_data,int dLen,double *error)
/* Sample mean, with 95% confidence interval computed using ips estimator. */
{
  /* The first thing to do is copy and zero-mean the data. */

  double *zero_meaned_data;
  zero_meaned_data = malloc(dLen*sizeof(double));
  double smean = sample_mean(raw_data,dLen);
  int k;
  for(k=0;k<dLen;k++) { zero_meaned_data[k] = raw_data[k] - smean; }
 
  double sigma = ips_sigma_estimator(zero_meaned_data,dLen);

  assert(error != NULL);
  *error = 1.96 * sigma / sqrt(dLen);
  free(zero_meaned_data);

  return smean;
}

/*****************************************************************************/

/* gips version 2.0 

   We learned from the last incremental strategy that it's incredibly,
   painfully slow to compute things in this way. So we come up with a
   new strategy which we call "buffered ips" estimators.

   Here's the basic idea.

   1.) If we knew in advance how many lagged covariances we were going
   to need, we could simply keep a buffer of that many data points
   around and update the covariances as we went. This would
   (potentially, eventually) pose "adding small to large" problems for
   very long runs. I don't think those errors are likely to happen in
   practice, but the real problem is that it's the least efficient way
   to do the actual covariance calculation.

   2.) Doing the actual dot products is something that can be productively
   hardware- or GPU-accelerated if only we make the data buffer long
   enough to make it worthwhile.

   So the new strategy is to take a long data buffer, and update as
   infrequently as possible. That is,

   Allocate a long data buffer (say, 1 million data points). 

   When it fills (or when we are called on to evaluate the ips
   estimator), compute as many lagged covariances as we need using the
   fastest linear algebra library available to us, keeping track of the 
   number of entries in the initial positive sequence (ipslength). 

   Double the number of lagged covariances in the initial positive
   sequence, and store this many data points when we reset the buffer,
   along with the corresponding collection of (current) lagged covariances.
   
   When the buffer refills, update the current collection of lagged
   covariances. (Whether we need them or not.) The buffer size should
   be invisible to the API because (in principle) any number of lagged
   covariances can be dealt with gracefully by the code.

   The actual details of the blas calls can be encapsulated so that if 
   (for instance) cuda can be linked with, but doesn't work, we fall
   back to regular blas (or in the worst case to handwritten dot products).

***************************************************************************/

#define GIPS_MIN(a,b) (((a)<(b))?(a):(b))
#define GIPS_MAX(a,b) (((a)<(b))?(b):(a))

/* Our linear algebra interface comes down to two functions */

double lagged_dot(int lag,int n,double *vec) {

  /* Computes the (lagged) dot product of vec and its' shift by "lag" positions */

  assert(vec != NULL);
  if (lag > n) { return 0; }
  else {
    return cblas_ddot(n-lag,vec,1,vec+lag,1);
  }
}

double split_lagged_dot(int lag,int sizeA,double *A,int sizeB,double *B) {

  /* 
     Computes lagged dots between pairs at the END of vector A and the
     START of vector B (with given sizes). There's no particular
     relationship between lag, sizeA, and sizeB, so we have to be
     ready for anything.
  */
  
  double val = 0;

  /* Indexing is a tricky issue here. We have

     ...A[sizeA-3], A[sizeA-2], A[sizeA-1], B[0], B[1], B[2] ...

     at the junction. Consider the lag-3 case as an example. We should
     start with aidx = sizeA - 3, and bidx = 0, and continue until
     aidx = sizeA - 1 or bidx = sizeB - 1. */

  int aidx,bidx;  /* These vectors will index into A and B. */

  for(aidx=sizeA-lag,bidx=0;
      aidx<sizeA && bidx<sizeB;aidx++,bidx++) { 

    if (aidx >= 0) { /* It could be negative if we haven't gotten to the start yet */

      val += A[aidx]*B[bidx];

    }

  }

  return val;
}

  
double dot(int n,double *vecA,double *vecB) { 

  /* Compute an (ordinary) dot product. */
  return cblas_ddot(n,vecA,1,vecB,1);

}

double geyer_ips_data_val(geyer_ips_t *gips,int i) { 

  /* Look up a value from gips using databuf, initbuf, and finalbuf,
     if we have it. */

  assert(i >= 0 && i < gips->samples); 

  if (gips->samples < gips->buffer_size) { /* We have everything! */

    return gips->databuf[i];

  } 

  if (gips->ifsize != 0) { /* We've initialized initbuf and finalbuf */

    if (i < gips->ifsize) { 

      return gips->initbuf[i];

    } else if (i >= gips->samples - gips->ifsize) {

      return gips->finalbuf[i % gips->ifsize];

    } else {

      printf("geyer_ips_data_val: Requested index %d is not within ifsize == %d of start or end (%d) of data\n",
	     i,gips->ifsize,gips->samples);
      return 0;

    }

  }

}
   

double zero_meaned_covariance(geyer_ips_t *gips,int lag,double lagged_dot,double sample_mean)
  
/* 
   Lagged_dot holds the lagged dot of the original data. But keep in
   mind that this data ISN'T zero-meaned. Assuming that 

   gips->ifsize > lag = k

   we can compensate, because we have access to the initial and final
   elements of the data.

   Here's how:

   If we write the data as d[i], and write sample_mean = mu, the
   lagged dot we want (of the zero-meaned data) is
   
   sum_{i=1}^{m-k} (d[i] - mu)(d[i+k] - mu) = 
   
   (a)    sum_{i=1}^{m-k} d[i] d[i+k] 
   (b)   -sum_{i=1}^{m-k} d[i] mu
   (c)   -sum_{i=1}^{m-k} d[i+k] mu
   (d)   +sum_{i=1}^{m-k} mu^2.
   
   We already have (a), which is lagged_dot.
   Now (b) and (c) are roughly 
   
   mu sum_{i=1}^m d[i] = mu * mu * m = mu^2 m,
   
   but we're missing some terms. In fact, 
   
   -(b) = sum_{i=1}^m d[i] mu - sum_{i=m-k+1}^m d[i] mu 
   = mu^2 m - mu sum_{i=m-k+1}^m d[i].
   
   and
   
   -(c) = sum_{i=1}^m d[i] mu - sum_{i=1}^k d[i] mu
   = mu^2 m - mu sum_{i=1}^k d[i].
   
   Similarly, (d) is (m-k) mu^2, proving that we want
   
   (a) - 2 mu^2 m + (m-k) mu^2 = (a) + (m - k - 2 m) mu^2 = (a) - (k+m) mu^2
   
   + correction terms from (b) and (c).  */

{
  int i;
  double result = lagged_dot;
  assert(lag < gips->ifsize);
  double initial_correction = 0, final_correction = 0;

  for(i=0;i<lag;i++) { 

    initial_correction += geyer_ips_data_val(gips,i);
    final_correction   += geyer_ips_data_val(gips,gips->samples-1-i);

  }

  result += sample_mean*(initial_correction + final_correction);
  result -= sample_mean*sample_mean*(double)((gips->samples + lag));
  result /= (double)(gips->samples); 
      
  /* Note: Even though there are fewer terms in the actual
     computation of lagged covariance, and hence this is a biased
     estimator, Geyer recommends division by the full number of
     samples (see p. 475 of "Practical MCMC"). */

  return result;
  
}

/*********************************************************/

geyer_ips_t *geyer_ips_new(int gips_buffer_size) {

  geyer_ips_t *gips;
  int i;

  gips = calloc(1,sizeof(geyer_ips_t));
  assert(gips != NULL);

  /* Allocate a new storage buffer. NaN it to make sure we don't use
     uninitialized values later. */

  gips->buffer_size = gips_buffer_size;
  gips->databuf = malloc(gips->buffer_size*sizeof(double));
  assert(gips->databuf != NULL);
  for(i=0;i<gips->buffer_size;i++) { gips->databuf[i] = 1.0/0.0; } 
  gips->databuf_used = 0;

  gips->samples = 0;

  /* We can't allocate initbuf and finalbuf until we know how big they
     are supposed to be, which won't be computed until the first
     buffer fill. */

  gips->initbuf  = NULL;
  gips->finalbuf = NULL;
  gips->ifsize = 0;

  /* We won't store (or compute!) any lagged dot product until forced
     to. */

  gips->ldot  = NULL;
  gips->nldots = 0;
  gips->overlapbuf = NULL;

  /* We use these to compute the actual ips error estimate */

  gips->sample_sum = 0;
  gips->Gamma = NULL;
  gips->N = 0;
  gips->sigma = 0;

  /* And include some self-profiling code. */

  gips->dot_product_time = 0;

  return gips;

}

void geyer_ips_add(geyer_ips_t *gips,double data) 

{
  if (gips->databuf_used < gips->buffer_size) { /* There's room */

    gips->databuf[(gips->samples % gips->buffer_size)] = data;
    gips->databuf_used++;
    gips->sample_sum += data;
    if (gips->finalbuf != NULL & gips->ifsize > 0) { gips->finalbuf[gips->samples % gips->ifsize] = data; }
    gips->samples++;

  } else if (gips->samples == gips->buffer_size) { /* The buffer is full for the first time. */

    /* The first buffer fill requires us to compute the number of
       (lagged) covariances to keep.  At this point, ALL the data is
       in the buffer. So the smartest thing to do is simply to perform
       the computation in full using the old code and record the
       number of Gammas. */

    double *zero_meaned_data = calloc(gips->buffer_size,sizeof(double));
    assert(zero_meaned_data != NULL);
    double sample_mean = gips->sample_sum /(double)(gips->buffer_size);
    int i,N;
    for(i=0;i<gips->buffer_size;i++) { zero_meaned_data[i] = gips->databuf[i] - sample_mean; }
    double gamma0, gamma1, *Gammas, sigma;
    sigma = ips_sigma_estimator_with_gammas(zero_meaned_data,gips->buffer_size,&gamma0,&gamma1,&Gammas,&N);

    /* We want to keep at least three times as many lagged dots as we expect to use, 
       (and at least 50 in any event), but we have to cap this at buffer size to 
       deal with the artifically small buffers which occur in testing. */

    gips->nldots = GIPS_MIN(GIPS_MAX(6*(N+1),50),gips->buffer_size); 
    free(Gammas); free(zero_meaned_data);

    /* Now that we know we can't simply save everything (we ran out of buffer space!) 
       we get ready to transition to a different strategy. We'll cache the first and 
       last nldots samples first. */

    gips->ifsize = gips->nldots;
    gips->initbuf = malloc(gips->ifsize*sizeof(double));
    gips->finalbuf = malloc(gips->ifsize*sizeof(double));
    memcpy(gips->initbuf,gips->databuf,gips->ifsize*sizeof(double));

    /* Filling the finalbuf isn't so easy, because we expect the data
       to be in a particular PLACE. Further, we have to deal with the
       fact that ifsize and buffer_size have no particular relation to
       one another. 

       So we fake the last "ifsize" values of gips->samples and copy
       the corresponding stuff from databuf to finalbuf, as if we'd
       been recording the data all along. */

    for(i=gips->samples-gips->ifsize;i<gips->samples;i++) { 
      
      gips->finalbuf[i % gips->ifsize] = gips->databuf[i % gips->buffer_size];

    }
      
    /* Now we go ahead and store lagged dots for a rainy day. */

    gips->ldot = calloc(gips->nldots,sizeof(double));
    assert(gips->ldot != NULL);

    clock_t start,end;

    start = clock();

    for(i=0;i<gips->nldots;i++) { 

      gips->ldot[i] = lagged_dot(i,gips->buffer_size,gips->databuf);
     
    }
  
    end = clock();
    gips->dot_product_time += ((double)(end - start))/CLOCKS_PER_SEC;
   
    /* 
       We now know the maximum number of lagged dot products to
       compute. We need to store this many data points to in overlapbuf to 
       compute future lagged dots the next time the buffer fills. 
    */

    gips->overlapbuf = calloc(gips->nldots,sizeof(double));
    assert(gips->overlapbuf != NULL);

    memcpy(gips->overlapbuf,
	   gips->databuf + gips->buffer_size - gips->nldots,
	   gips->nldots * sizeof(double));

    /* We have just copied max_ipslength data values to the overlap
       buffer for safekeeping. As a bug defense, we now NaN out the
       data buffer and store the new value. */
    
    for(i=0;i<gips->buffer_size;i++) { gips->databuf[i] = 1.0/0.0; }

    /* Finally, we store the new value. */

    gips->databuf[0] = data;
    gips->databuf_used = 1;
    gips->sample_sum += data;
    gips->finalbuf[gips->samples % gips->ifsize] = data;
    gips->samples++;

  } else { /* The buffer is full, but not for the first time. */

    int i;
    clock_t start,end;

    /* The first step is to update the lagged covariances. */

    start = clock();

    for(i=0;i<gips->nldots;i++) { 
     
      gips->ldot[i] += lagged_dot(i,gips->buffer_size,gips->databuf);

    }

    /* We now need to update the missing terms in the lagged
       dots. Here's the idea: the skip-k pairs which are fully
       contained within overlapbuf were already computed (when
       overlapbuf was part of the last databuf), and the skip-k pairs
       which are fully contained within (the new) databuf were just
       computed. But the others are "split" lagged dots. */

    for(i=0;i<gips->nldots;i++) {

      gips->ldot[i] += split_lagged_dot(i,gips->nldots,gips->overlapbuf,gips->databuf_used,gips->databuf);

    }

    end = clock();
    gips->dot_product_time += ((double)(end - start))/CLOCKS_PER_SEC;

    /* Now that we've updated the lagged covariances, we shift the
       data to overlapbuf (nuking overlapbuf in the process) and NaN
       the data buffer again. */

    memcpy(gips->overlapbuf,
	   gips->databuf + gips->buffer_size - gips->nldots,
	   gips->nldots * sizeof(double));

    /* We have just copied max_ipslength data values to the overlap
       buffer for safekeeping. As a bug defense, we now NaN out the
       data buffer and store the new value. */
    
    for(i=0;i<gips->buffer_size;i++) { gips->databuf[i] = 1.0/0.0; }

    /* Finally, we store the new value. */

    gips->databuf[0] = data;
    gips->databuf_used = 1;
    gips->sample_sum += data;

    gips->finalbuf[gips->samples % gips->ifsize] = data; 
    /* We know we HAVE finalbuf because we set it up at the first buffer fill. */

    gips->samples++;
 
  }

}


double geyer_ips_value(geyer_ips_t *gips,double *error,bool *ok)

{
  /* Our cards have been called, and we need to report a sample 
     mean and error. This needs to be done a bit carefully because
     we may call this repeatedly (and keep adding data in between).
     So we can't mess up the internal elements (which should be 
     opaque to the caller, anyway). 

     Sets ok to false if we ran out of lagged dots (in this case, 
     the buffer size is generally too small).
  */

  if (gips->samples == 0) { 

    if (error != NULL) { *error = 0; } 
    return 0;

  }

  /* Step 0. Compute the current sample mean. */

  double mu; /* The sample mean. All data was recorded in sample_sum as it came in. */
  int i;

  mu = gips->sample_sum/(double)(gips->samples);

  /* Step 0.1. Check to see if we need to compute error. */

  if (error == NULL) { *ok = true; return mu; } /* If error is not desired, this is easy! */

  /* Step 0.2. Check to see if we've got everything in the buffer and can use existing code. */

  if (gips->samples <= gips->buffer_size) { 

    double *zero_meaned_data = malloc(gips->samples*sizeof(double));
    assert(zero_meaned_data != NULL);

    int i;
    for(i=0;i<gips->samples;i++) { zero_meaned_data[i] = gips->databuf[i] - mu; }

    double gamma0,gamma1,*Gammas;
    int N;
    gips->sigma = ips_sigma_estimator_with_gammas(zero_meaned_data,gips->samples,&gamma0,&gamma1,&Gammas,&N);
    
    /* Now store the results. */
    if (gips->Gamma != NULL) { free(gips->Gamma); gips->Gamma = NULL; }

    gips->N = N+1;
    gips->Gamma = calloc(gips->N,sizeof(double));
    assert(gips->Gamma != NULL);
   
    gips->Gamma[0] = gamma0 + gamma1;
    for(i=1;i<gips->N;i++) { gips->Gamma[i] = Gammas[i-1]; }
   
    free(zero_meaned_data);
    free(Gammas); 

    *error = 1.96 * gips->sigma / sqrt(gips->samples);
    *ok = true;
    return mu;

  }
    
  /* Step 0.3. We need to compute sigma and not everything is in the buffer. 
     Get ready to rumble-- this is the hard one. We want to assemble a complete 
     list of lagged dots first. This is probably not the most efficient strategy
     since, we'll get "extra" ldots, but it's conceptually the cleanest. */

  double *ldots; 
  ldots = malloc(gips->nldots*sizeof(double));
  assert(ldots != NULL); 

  /* First, copy the old cache of ldots. */

  memcpy(ldots,gips->ldot,gips->nldots*sizeof(double)); 

  /* Second, add contributions from the current contents of databuf */
  
  int max_lag = GIPS_MIN(gips->nldots,gips->databuf_used); 
  for(i=0;i<max_lag;i++) { ldots[i] += lagged_dot(i,gips->databuf_used,gips->databuf); }

  /* Third, add contributions from the "straddle" pairs where one
     element has been cached to overlapbuf and the other lies in
     databuf. */

  for(i=0;i<gips->nldots;i++) { 

    ldots[i] += split_lagged_dot(i,gips->nldots,gips->overlapbuf,gips->databuf_used,gips->databuf);

  }

  /* Now we need to turn the ldots into lagged covariances of zero-meaned data (littlegammas). */

  double *littlegamma;
  littlegamma = calloc(gips->nldots,sizeof(double));
  assert(littlegamma != NULL); 

  for(i=0;i<gips->nldots;i++) {

    littlegamma[i] = zero_meaned_covariance(gips,i,ldots[i],mu);

  }
    
  /* These are all of the littlegammas we're going to get to work with. 
     We go ahead and assemble the initial positive sequence of Gammas from them. */

  /* Make sure we have a fresh buffer of Gammas to work with. */
  
  if (gips->Gamma != NULL) { free(gips->Gamma); }
  gips->Gamma = calloc(gips->nldots,sizeof(double)); 
  assert(gips->Gamma != NULL);
 
  /* Fill the Gamma array, working our way through until Gammas are
     no longer positive or until we run out of lagged dots to use. */
    
  bool doneflag = false;

  for(gips->Gamma[0] = littlegamma[0] + littlegamma[1],gips->N=1;
      gips->Gamma[gips->N-1] > 0 && (2*gips->N)+1 < gips->nldots;
      gips->N++) { 

    gips->Gamma[gips->N] = littlegamma[2*gips->N] + littlegamma[(2*gips->N)+1];
    
  }

  if (gips->Gamma[gips->N-1] < 0) { gips->N--; *ok = true; } /* The last Gamma was negative. */
  else { *ok = false; } /* We computed all of the Gammas that we could, but they are still positive! */
  /* Now we can actually compute the estimator. */

  gips->sigma = -littlegamma[0]; 
  for(i=0;i<gips->N;i++) { gips->sigma += 2*gips->Gamma[i]; }
  gips->sigma = sqrt(gips->sigma); 
  
  *error = 1.96 * gips->sigma / sqrt(gips->samples);

  free(littlegamma); free(ldots);

  return mu;

}

void geyer_ips_free(geyer_ips_t **gips) { 

  if (*gips != NULL) { 

    (*gips)->buffer_size = 0;

    if ((*gips)->databuf != NULL) { free((*gips)->databuf); (*gips)->databuf=NULL; }
    (*gips)->samples = 0;
    
    if ((*gips)->overlapbuf != NULL) { free((*gips)->overlapbuf); (*gips)->overlapbuf=NULL; }
    (*gips)->sample_sum = 0;

    if ((*gips)->initbuf != NULL) { free((*gips)->initbuf); (*gips)->initbuf=NULL; }
    if ((*gips)->finalbuf != NULL) { free((*gips)->finalbuf); (*gips)->finalbuf=NULL; }
    (*gips)->ifsize = 0;

    if ((*gips)->ldot != NULL) { free((*gips)->ldot); (*gips)->ldot=NULL; }
    (*gips)->nldots = 0;

    (*gips)->dot_product_time = 0;
    
    if ((*gips)->Gamma != NULL) { free((*gips)->Gamma); (*gips)->Gamma=NULL; }
    (*gips)->N = 0;

    free(*gips);
    *gips = NULL;

  }

}
