/* randomwalk_test_test.c
author: Harrison Chapman

This is a file which tests `randomwalk_test()`, which is a part of the
ccode_test test suite.

We know that there will be some cases in which `randomwalk_test()`
will fail to project to a pdcode object; principle among reasons are
those in which edges of the random polygon generated are
self-intersecting (in which case we would be unable to determine the
sign of that given crossing in constructing the pdcode object).

Based on a prior run of this program, the expected success rate of a
given loop run in randomwalk_test() run on random polygons with 250
edges is about 999974 in 1000000 (99.9974%). Hence, the expected
failure probability of a given random 250-gon is 2.6e-5.

In ccode_test, randomwalk_test() is called to process 1001 random
250-gons into pdcode objects. If a single loop 'fails', the entire
test 'fails'. Based on the prior run of this program, we estimate (as
of March 2016) that this step should fail about 97.431% of the time,
or about 3-in-100 times.
*/

#include "config.h"
#include "plCurve.h"
#include "plcTopology.h"

#include<math.h>
#include<stdlib.h>
#include<time.h>

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
