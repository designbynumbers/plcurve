#include "config.h"
#include "plCurve.h"
#include "plcTopology.h"

#ifdef HAVE_MATH_H
#include<math.h>
#endif
#ifdef HAVE_STDLIB_H
#include<stdlib.h>
#endif
#ifdef HAVE_TIME_H
#include<time.h>
#endif
#ifdef HAVE_STDBOOL_H
#include<stdbool.h>
#endif

#define VERBOSE 1

bool randomwalk_test(gsl_rng *rng,int nedges,bool verbose,int reps);
int randomwalk_test_instance(gsl_rng *rng, int nedges, bool verbose);
void set_pd_code_from_plCurve_verbose(bool val);

int main() {
    const gsl_rng_type *T;
    gsl_rng *rng;

    gsl_rng_env_setup();

    T = gsl_rng_default;
    rng = gsl_rng_alloc(T);

    int seedi;

    seedi = time(0);
    gsl_rng_set(rng, seedi);

    randomwalk_test(rng, 250, false, 1000000);
}

int randomwalk_test_instance(gsl_rng *rng, int nedges, bool verbose) {
    plCurve *L = plc_random_equilateral_closed_polygon(rng,nedges);
    pd_code_t *projected_pd = pd_code_from_plCurve(rng,L);
    double ncross;

    if (projected_pd == NULL) {

        printf("FAIL\n");
        printf("unable to generate pd code from plcurve\n");
        return -1;

    }

    ncross = (double)(projected_pd->ncross);

    if (!pd_ok(projected_pd)) {

        printf("FAIL\n");
        printf("pd code generated from plCurve was:\n");
        printf("-----------------------------------\n");
        pd_write(stdout,projected_pd);
        pd_code_free(&projected_pd);
        
        return -1;
    }
    
    pd_code_free(&projected_pd);
    plc_free(L);
    return ncross;
}

bool randomwalk_test(gsl_rng *rng,int nedges,bool verbose,int reps) {

  clock_t start,end;
  double cpu_time_used;

  printf("------------------------------------------------\n"
	 "random equilateral %d-gon test\n"
	 "------------------------------------------------\n",nedges);

  start = clock();
  printf("computing %d pd_codes (random rotation enabled)...",reps);
  set_pd_code_from_plCurve_verbose(verbose);

  if (verbose) { printf("\n\n"); }
  
  int i;
  int n_success = 0;
  double totalcross = 0;
  double ncross;
  for(i=0;i<reps;i++) {
      ncross = randomwalk_test_instance(rng, nedges, verbose);
      
      if (ncross < 0) {
          // test failure
      } else {
          // test success
          n_success += 1;
          totalcross += ncross;
      }
  }
  
  end = clock();
  cpu_time_used = ((double)(end - start)/CLOCKS_PER_SEC);
  printf("done (%2.4g sec, %i of %i (=%f) passed pd_ok)\n",
         cpu_time_used, n_success, reps, (1.0*n_success)/(1.0*reps));

  printf("average #crossings of produced pd codes...%3.4g\n",totalcross/(double)(reps));
  
  return true;
 
}
