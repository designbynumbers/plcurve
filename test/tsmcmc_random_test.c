/* 

   tsmcmc_random_test.c : Test code for the random polygon generation functions in plCurve. 

*/

#include<plCurve.h>
#include<plcRandomPolygon.h>

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

bool test_equilateral_prediction(gsl_rng *rng,int n,double integrand(plCurve *L),double expected,char *prediction_name)
{
  printf("---------------------------------------"
	 "-----------------------------------\n"
	 "Testing %s for equilateral %d-gons.\n"
	 "Running for 10 seconds...",prediction_name,n);

  tsmcmc_run_parameters run_params = tsmcmc_default_unconfined_parameters();
  tsmcmc_run_stats run_stats;
  double error,result;

  tsmcmc_triangulation_t T = tsmcmc_spiral_triangulation(n);
  result = tsmcmc_equilateral_expectation(rng,integrand,500000,10,T,run_params,&run_stats,&error);

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
double skip_squared_chordlength(plCurve *L) {

  /* Uses the global glob_skip to set the skip length,
     then computes the average sqaured length of chords skipping
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
  int num_n = 4, num_k = 3;
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


bool test_ftc_prediction(gsl_rng *rng,double ftc,int n,double integrand(plCurve *L),double expected,char *prediction_name)
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
  result = tsmcmc_fixed_ftc_expectation(rng,integrand,ftc,500000,10,T,run_params,&run_stats,&error);

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

double chordlength(plCurve *L) {

  return plc_distance(L->cp[0].vt[glob_v0],L->cp[0].vt[glob_v1]);

}

bool fixed_failure_to_close_chordlength_tests(gsl_rng *rng) 

/* We got the data for this on an ad-hoc basis from polymake, so 
   only a few tests are possible. */

{

  double ftc5gon25[2] = {2.0,1.5};
  int n = 10;

  double ftc2[7] = {1.95118, 1.90235, 1.8472, 1.76726, 1.65451, 1.50066, 1.27735};
  double ftc3[7] = {3967.0/1428.0,1825.0/714.0,1111.0/476.0,144701.0/68544.0,8051.0/4284.0,74423.0/45696.0,635.0/476.0};
  double ftc5[7] = {2143.0/476.0,953.0/238.0,1669.0/476.0,358.0/119.0,1195.0/476.0,30655.0/15232.0,719.0/476.0};

  char pred_name[2048];
  int i;

  /* First, run some very simple checks on a 5-gon with failure-to-close = 2.5*/

  for(i=0;i<2;i++) {

    glob_v0 = 0; glob_v1 = i+2;
    sprintf(pred_name,"length of %d-%d chord",glob_v0,glob_v1);
    if (!test_ftc_prediction(rng,2.5,5,chordlength,ftc5gon25[i],pred_name)) {

	return false;
	
      }
  }	

  /* Now test equilateral 10-gons. */

  for(i=2;i<8;i++) {

    glob_skip = i;
    sprintf(pred_name,"squared length of skip %d chord",i);
    if (!test_ftc_prediction(rng,1.0,10,skip_squared_chordlength,
			     eq_pol_prediction(10,i),pred_name)) { 

      return false;

    }

  }


  for(i=0;i<7;i++) { 

    glob_v0 = 0; glob_v1 = i+2;
    sprintf(pred_name,"length of %d-%d chord",glob_v0,glob_v1);
    if(!test_ftc_prediction(rng,3.0,10,chordlength,ftc3[i],pred_name)) { 

    }

  }

  for(i=0;i<7;i++) { 

    glob_v0 = 0; glob_v1 = i+2;
    sprintf(pred_name,"length of %d-%d chord",glob_v0,glob_v1);
    if(!test_ftc_prediction(rng,2.0,10,chordlength,ftc2[i],pred_name)) { 

      return false;

    }

  }

 

  for(i=0;i<7;i++) { 

    glob_v0 = 0; glob_v1 = i+2;
    sprintf(pred_name,"length of %d-%d chord",glob_v0,glob_v1);
    if(!test_ftc_prediction(rng,5.0,10,chordlength,ftc5[i],pred_name)) { 

      return false;

    }

  }

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

  printf("plcRandomPolygon test suite \n"
	 "------------------------------- \n"
	 "plCurve generates Markov chains of equilateral polygons using the\n"
	 "toric symplectic algorithm of Cantarella and Shonkwiler [arxiv,2013].\n"
	 "\n"
	 "This test suite code is using the random number generator %s \n"
	 "with seed %d.\n"
	 "\n"
	 "This program tests the polygons generated against theoretical calculations\n"
	 "==========================================================================\n"
	 ,gsl_rng_name(rng),seedi);

  if (!fixed_failure_to_close_chordlength_tests(rng)) { PASS = false; }
  if (!equilateral_unconfined_chordlength_tests(rng)) { PASS = false; }
       
  gsl_rng_free(rng);

  printf("=========================================================================\n");
 
  if (PASS) { 

    printf("Random Polygon Test Suite: PASS\n");
    exit(0); 

  } else { 

    printf("Random Polygon Test Suite: FAIL\n");
    exit(1); 

  }

}
