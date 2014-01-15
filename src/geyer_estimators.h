/* 

   Geyer gives various estimators for the variance in a Markov chain
   central limit theorem application. This is an implementation of the 
   various estimators.

   initial positive sequence estimator: ips_estimator (and friends)
   initial monotone sequence estimator: ims_estimator (and friends)
   initial convex sequence estimator:   ics_estimator (and friends)

   These estimators are meant to be applied to a (large) buffer of
   doubles coming from a Markov chain applied to some sample space.
   We can use them to compute approximate error bars for Markov chain
   integrations.

   We use BLAS calls to do this, since all the covariance computations
   are basically about dot products. Of course, we could do a much
   better job if we used OpenCL or something, but we're not going that
   far in the implementation direction at the moment.

*/

#ifndef __GEYER_ESTIMATORS_H__
#define __GEYER_ESTIMATORS_H__ 1

double ips_sigma_estimator(double *zero_meaned_data,int dLen);
/* Estimate of variance from zero-meaned sequence of values at Markov chain points. */

double ips_sigma_estimator_with_gammas(double *zero_meaned_data, int m, 
				       double *gamma0,
				       double *gamma1,
				       double **Gammas, 
				       int *N);
/* Estimate of variance together with lagged autocovariances used in computation, assumes zero-meaned data. */

double sample_mean(double *data,int dLen);
/* Sample mean for data. */

double ips_value(double *raw_data,int dLen,double *error);
/* Sample mean, with 95% confidence interval computed using ips estimator. Expects raw (i.e. NOT zero-meaned) data*/

/* This is meant to be an opaque type. */

typedef struct geyer_incremental_mean_struct {
  
  int nbuffers;         /* Number of levels */
  double **buffers;     /* Pointers to the actual buffers. */
  int *buf_used;        /* Number of entries in buffer[i] (often 0 if we haven't reached this level yet) */

  int buf_size;      /* Number of elements stored in each buffer. */
  
} geyer_incremental_mean_t;

geyer_incremental_mean_t *geyer_incremental_mean_new(int buffer_size);
void geyer_incremental_mean_add(geyer_incremental_mean_t *gim,double data); 
double geyer_incremental_mean_compute(geyer_incremental_mean_t *gim);
void geyer_incremental_mean_free(geyer_incremental_mean_t **gim); 

/* Again, we have an opaque type. */

typedef struct geyer_incremental_ips_struct {

  unsigned long int stepcount; 
  geyer_incremental_mean_t *sample_mean;

  int nac;
  geyer_incremental_mean_t **raw_lagged_cov;

  /* We need to store only the first nac and the last nac values */
  
  double *initial_buf;
  double *databuf;

  /* We also store data from the last time we
     computed a sigma estimate */

  double mu;          /* Sample mean */
  double sigma;       /* Last sigma estimate */

  double gamma0;      /* Sample variance */
  double gamma1;      /* Lag-1 sample variance */
  double *gammas;     /* gamma[k] = lag-k autocovariance */

  int    N;           /* Number of Gammas used in last ips sigma estimate */
  int    nGammas;     /* Total number of Gammas available */
  double *Gammas;     /* Gamma[k] = gamma[2k] + gamma[2k+1] */  

} geyer_incremental_ips_t;

/* The idea of the incremental ips type is that we'll buffer "nac" 
   data points in a circular buffer and then use it to keep track of  
   "nac" raw autocovariances. (We can't zero-mean the data since we
   don't know what the mean is going to be.) We can then use the
   raw autocovariances to compute an ips estimator when called for. 
   We need to know "nac" upfront, since we won't buffer extra data. 

   We return our Markov chain estimate for the value, together with the ips_sigma 
   error estimate. */

geyer_incremental_ips_t *geyer_incremental_ips_new(int nac);
void   geyer_incremental_ips_add(geyer_incremental_ips_t *giips,double data);
double geyer_incremental_ips_value(geyer_incremental_ips_t *giips,double *error,bool *ok);
void   geyer_incremental_ips_free(geyer_incremental_ips_t **giips);

typedef struct geyer_ips_struct {
  
  int     buffer_size;         /* the size of databuf */
  double *databuf;
  int     databuf_used;        /* number of elements in databuf currently filled */

  int     samples;             /* total number of samples recorded */
  double *overlapbuf; 

  double *initbuf;             /* inital buffer stores the first data points recorded */
  double *finalbuf;            /* finalbuf stores the last data points */
  int     ifsize;              /* Both of these are size "ifsize" */

  long double sample_sum;      /* sum of all samples */
  
  double *ldot;                /* ldot = ith lagged dot product: sum_j data[j]*data[i+j] */
  int     nldots;              /* Number of lagged dots stored in buffer. */

  double  dot_product_time;    /* Self-profiling code */ 
  
  double  sigma;               /* The actual error estimate (as of last call to geyer_ips_value */
  double  *Gamma;              /* The initial positive sequence of Gammas. */
  int     N;                   /* The number of positive Gammas. */

} geyer_ips_t;

geyer_ips_t *geyer_ips_new(int gips_buffer_size);
void         geyer_ips_add(geyer_ips_t *gips,double data); 
double       geyer_ips_value(geyer_ips_t *gips,double *error_estimate,bool *ok);
void         geyer_ips_free(geyer_ips_t **gips);

#endif
