/*

  random_sap_test.c : Test code for the random self-avoiding polygon generation functions in plCurve.

*/

#include"plCurve.h"
#include"plcRandomPolygon.h"

#include<config.h>

#ifdef HAVE_MATH_H
 #include<math.h>
#endif

#ifdef HAVE_STDLIB_H
 #include<stdlib.h>
#endif

#ifdef HAVE_TIME_H
#include<time.h>
#endif

#ifdef HAVE_COMPLEX_H
#include<complex.h>
#endif

#ifdef HAVE_STDIO_H
#include<stdio.h>
#endif

#ifdef HAVE_GSL_GSL_RNG_H
#include<gsl/gsl_rng.h>
#endif

gsl_rng *rng; /* The global random number generator */

bool plc_is_sap_internal(plCurve *L,bool verbose); /* An internal debugging function in plcRandomSap.c */

bool PAPERMODE;
FILE *outfile;

void timestamp(FILE *outfile) {

  fprintf(outfile,"# Run generated:");

  time_t curtime;
  struct tm *loctime;

  /* Get the current time. */
  curtime = time (NULL);

  /* Convert it to local time representation. */
  loctime = localtime (&curtime);

  /* Print out the date and time in the standard format. */
  fputs (asctime (loctime), outfile);

}

bool test_equilateral_prediction(gsl_rng *rng,int n,double integrand(plCurve *L, void *args),double expected,char *prediction_name)
{
  printf("---------------------------------------"
	 "-----------------------------------\n"
	 "Testing %s for equilateral %d-gons.\n"
	 "Running for 10 seconds...",prediction_name,n);

  tsmcmc_run_parameters run_params = tsmcmc_default_unconfined_parameters();
  tsmcmc_run_stats run_stats;
  double error,result;

  tsmcmc_triangulation_t T = tsmcmc_spiral_triangulation(n);
  result = tsmcmc_equilateral_expectation(rng,integrand,NULL,50000,2,T,run_params,&run_stats,&error);

  printf("done.\n"
         "Run statistics: \n"
	 "\tstep counts (dihedral/momentpoly/permute): %d/%d/%d\n"
	 "\ttotal time                               : %2.2g sec\n"
	 "\tlagged covariances used in computing err : %d (of %d)\n",
	 run_stats.dihedral_steps,
	 run_stats.mp_steps,
	 run_stats.permute_steps,
	 run_stats.total_seconds,
	 run_stats.max_lagged_covariance_used,
	 run_stats.lagged_covariances_available);

  printf("checking expected value (%g) within error bounds...",expected);

  if (fabs(expected - result) < error) {

    printf("pass\n"
	   "Details:\n"
	   "\tComputed expectation: %4.4g\n"
	   "\tExact expectation   : %4.4g\n"
	   "\tActual Error        : %4.2g\n"
	   "\tError Bound         : %4.2g\n",
	   result,
	   expected,
	   fabs(expected-result),error);

  } else if (fabs(expected - result) < 4*error) {

    printf("conditional pass\n"
	   "Details:\n"
	   "\tComputed expectation: %4.4g\n"
	   "\tExact expectation   : %4.4g\n"
	   "\tActual Error        : %4.2g\n"
	   "\t95%% Error Bound    : %4.2g\n"
	   "\tNever Error Bound   : %4.2g\n",
	   result,
	   expected,
	   fabs(expected-result),error,4*error);

    printf("Test will continue. The 95%% confidence interval fails 5%% of the time\n"
	   "so you are likely to observe this from time to time. However, we \n"
	   "should NEVER exceed the error bound above.\n");
    return true;

  } else {

    printf("FAIL (actual error %g >= 4*error bound %g)\n",fabs(expected-result),error);
    printf("This failure means the Markov chain has not converged (or there\n"
	   "is a serious bug somewhere in the library). In either case, it's\n"
	   "a reportable failure, and I'd like to know your platform so that\n"
	   "I can try to reproduce the failure.\n");
    return false;

  }

  printf("result...pass\n");
  printf("------------------------------------------"
	 "-----------------------"
	 "--------\n");
  return true;

}

int glob_skip;
double skip_squared_chordlength(plCurve *L, void *args) {

  /* Uses the global glob_skip to set the skip length,
     then computes the average squared length of chords skipping
     glob_skip edges in the plCurve L. */

  double total = 0;
  int vt;

  for(vt=0;vt<L->cp[0].nv;vt++) {

    total += plc_sq_dist(L->cp[0].vt[vt],L->cp[0].vt[(vt+glob_skip)%L->cp[0].nv]);

  }

  total /= (double)(L->cp[0].nv);
  return total;

}

double eq_pol_prediction(int n,int k) {
  return (double)((n-k)*k)/(double)((n-1));
}
char eq_pol_predstring[256] = "((n-k)/(n-1)) k";

bool equilateral_unconfined_chordlength_tests(gsl_rng *rng)
{
  int nvals[10] = {50,100,200,250};
  int kvals[10] = {10,20,30};
  int num_n = 2, num_k = 1;
  int n,k;
  char predname[2048];

  for(n=0;n<num_n;n++) {

    for(k=0;k<num_k;k++) {

      glob_skip = kvals[k];
      sprintf(predname,"avg squared length of skip %d chords",kvals[k]);

      if (!test_equilateral_prediction(rng,nvals[n],skip_squared_chordlength,eq_pol_prediction(nvals[n],kvals[k]),predname)) {

	return false;

      }

    }

  }

  return true;

}


bool test_ftc_prediction(gsl_rng *rng,double ftc,int n,double integrand(plCurve *L, void *args),double expected,char *prediction_name)
{
  printf("--------------------------------------"
	 "-----------------------------------\n"
	 "Testing %s for %d-gons with edge 0 length %g\n"
	 "all other edges equilateral (fixed failure-to-close)\n"
	 "Running for 10 seconds...",prediction_name,n,ftc);

  tsmcmc_run_parameters run_params = tsmcmc_default_unconfined_parameters();
  tsmcmc_run_stats run_stats;
  double error,result;

  tsmcmc_triangulation_t T = tsmcmc_fan_triangulation(n);
  result = tsmcmc_fixed_ftc_expectation(rng,integrand,NULL,ftc,50000,2,T,run_params,&run_stats,&error);

  printf("done.\n"
         "Run statistics: \n"
	 "\tstep counts (dihedral/momentpoly/permute): %d/%d/%d\n"
	 "\ttotal time                               : %2.2g sec\n"
	 "\tlagged covariances used in computing err : %d (of %d)\n",
	 run_stats.dihedral_steps,
	 run_stats.mp_steps,
	 run_stats.permute_steps,
	 run_stats.total_seconds,
	 run_stats.max_lagged_covariance_used,
	 run_stats.lagged_covariances_available);

  printf("checking expected value (%g) within error bounds...",expected);

  if (fabs(expected - result) < error) {

    printf("pass\n"
	   "Details:\n"
	   "\tComputed expectation: %4.4g\n"
	   "\tExact expectation   : %4.4g\n"
	   "\tActual Error        : %4.2g\n"
	   "\tError Bound         : %4.2g\n",
	   result,
	   expected,
	   fabs(expected-result),error);

  } else if (fabs(expected - result) < 4*error) {

    printf("conditional pass\n"
	   "Details:\n"
	   "\tComputed expectation: %4.4g\n"
	   "\tExact expectation   : %4.4g\n"
	   "\tActual Error        : %4.2g\n"
	   "\t95%% Error Bound     : %4.2g\n"
	   "\tNever Error Bound   : %4.2g\n",
	   result,
	   expected,
	   fabs(expected-result),error,4*error);

    printf("Test will continue. The 95%% confidence interval fails 5%% of the time\n"
	   "so you are likely to observe this much error during an average test a \n"
	   "couple of times.\n");
    return true;

  } else {

    printf("FAIL\n"
	   "Details:\n"
	   "\tComputed expectation: %4.4g\n"
	   "\tExact expectation   : %4.4g\n"
	   "\tActual Error        : %4.2g\n"
	   "\t95%% Error Bound     : %4.2g\n"
	   "\tNever Error Bound   : %4.2g\n",
	   result,
	   expected,
	   fabs(expected-result),error,4*error);

    printf("This failure means the Markov chain has not converged (or there\n"
	   "is a serious bug somewhere in the library). In either case, it's\n"
	   "a reportable failure, and I'd like to know your platform so that\n"
	   "I can try to reproduce the failure.\n");
    return false;

  }

  printf("result...pass\n");
  printf("------------------------------------------"
	 "-----------------------"
	 "--------\n");
  return true;

}

int glob_v0;
int glob_v1;

double chordlength(plCurve *L, void *args) {

  return plc_distance(L->cp[0].vt[glob_v0],L->cp[0].vt[glob_v1]);

}



int sample_quality_worker(gsl_rng *rng,int n,int s)
{

  plCurve *L;
  int samp;
  clock_t start, end;
  double cpu_seconds_used;

  printf("Testing sample-quality for %d-gons over %d samples...",n,s);
  start = clock();

  for(samp=0;samp<s;samp++) {

    L = plc_random_equilateral_closed_self_avoiding_polygon(rng,n);

    if (!plc_is_sap_internal(L,true)) {

      printf("fail\n\tTest failed at sample %d, which is not self-avoiding.\n",
	     samp);

      plc_write(stdout,L);
      plc_free(L);
      return false;

    }

    double longest,shortest,mean,var;
    plc_edgelength_stats(L,&longest,&shortest,&mean,&var);

    if (longest > 1.0 + 1e-9 || shortest < 1.0 - 1e-9) {

      printf("fail.\n\t Test failed at sample %d, which is not equilateral.\n\tLongest edge: %g, Shortest edge: %g\n",samp,longest,shortest);

      plc_write(stdout,L);
      
      plc_free(L);
      return false;

    }

    plc_free(L);
    
  }

  end = clock();
  cpu_seconds_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("pass (%g polygons/sec)\n",(double)(s)/cpu_seconds_used);
  return true;

}


bool sample_quality_tests(gsl_rng *rng)
{
  printf("--------------------------------------"
	 "-----------------------------------\n"
	 "Testing sample quality for saps.\n\n");

  if (!sample_quality_worker(rng,5,100000) ||
      !sample_quality_worker(rng,6,100000) ||
      !sample_quality_worker(rng,7,100000) ||
      !sample_quality_worker(rng,8,50000)  ||
      !sample_quality_worker(rng,9,50000)  ||
      !sample_quality_worker(rng,10,10000)  ||
      !sample_quality_worker(rng,11,1000)  ||
      !sample_quality_worker(rng,12,1000)  ||
      !sample_quality_worker(rng,13,1000)  ||
      !sample_quality_worker(rng,14,500)  ||
      !sample_quality_worker(rng,15,500)  ||
      !sample_quality_worker(rng,16,500) ||
      !sample_quality_worker(rng,17,100)  ||
      !sample_quality_worker(rng,18,100) ||
      !sample_quality_worker(rng,19,100)  ||
      !sample_quality_worker(rng,20,100)) {

    printf("Sample Quality Test: FAIL\n");
    printf("--------------------------------------"\
	   "-----------------------------------\n");
	   
    return false;
    
  }

  printf("Sample Quality Test: PASS\n");
  printf("--------------------------------------"
	 "-----------------------------------\n");

  return true;

}
  
double eq_pol_gyradius_prediction(int n) {
  return (double)(n+1)/(double)(3.0*n*n);
}
char eq_pol_gyradius_predstring[256] = "(1/3) (n+1)/n^2";

int main(int argc, char *argv[]) {

  bool PASS = {true};

  const gsl_rng_type * rng_T;

  gsl_rng_env_setup();
  rng_T = gsl_rng_default;
  rng = gsl_rng_alloc(rng_T);

  int seedi = time(0);

  //if (seed->count > 0) { seedi = seed->ival[0]; }
  //else { seedi = time(0); }

  gsl_rng_set(rng,seedi);

  if (argc > 1) { PAPERMODE = true; }

  printf("plcRandomSap test suite \n"
	 "------------------------------- \n"
	 "plCurve can generate random self-avoiding polygons either using the \n"
	 "folding-and-rejection Markov chain or by direct rejection sampling.\n"
	 "Burnin and skip steps for the folding chain are currently unknown,\n"\
	 "but the rejection sampler should produce perfect samples (eventually).\n"\
	 "\n"
	 "At the moment, all we can test is sample quality (are the polygons \n"\
	 "really closed, equilateral, and self-avoiding) and timing to convergence\n"
	 "for several integrands, as we don't know any theoretical expectations\n"\
	 "for off-lattice SAPS at this point. All polygons are self-avoiding\n"\
	 "in the hard-sphere model.\n"\
	 "\n"\
	 "This test suite code is using the random number generator %s \n"
	 "with seed %d.\n"
	 "\n"
	 "===========================================================\n"
	 ,gsl_rng_name(rng),seedi);

  if (!sample_quality_tests(rng)) { PASS = false; }
  gsl_rng_free(rng);

  printf("=========================================================================\n");

  if (PASS) {

    printf("Random SAP Test Suite: PASS\n");
    exit(0);

  } else {

    printf("Random SAP Test Suite: FAIL\n");
    exit(1);

  }

}
