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

uint64_t *xos; /* The global Xoshiro random number state */

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

bool equilateral_unconfined_chordlength_tests(uint64_t *xos)
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

      if (!0) {

	return false;

      }

    }

  }

  return true;

}

int glob_v0;
int glob_v1;

double chordlength(plCurve *L, void *args) {

  return plc_distance(L->cp[0].vt[glob_v0],L->cp[0].vt[glob_v1]);

}

int sample_quality_worker(uint64_t *xos,int n,int s)
{

  plCurve *L;
  int samp;
  clock_t start, end;
  double cpu_seconds_used;

  printf("Testing sample-quality for %d-gons over %d samples...",n,s);
  start = clock();

  for(samp=0;samp<s;samp++) {

    L = plc_random_equilateral_closed_self_avoiding_polygon(xos,n);

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


bool sample_quality_tests(uint64_t *xos)
{
  printf("--------------------------------------"
	 "-----------------------------------\n"
	 "Testing sample quality for saps.\n\n");

  if (!sample_quality_worker(xos,5,100000) ||
      !sample_quality_worker(xos,6,100000) ||
      !sample_quality_worker(xos,7,100000) ||
      !sample_quality_worker(xos,8,50000)  ||
      !sample_quality_worker(xos,9,50000)  ||
      !sample_quality_worker(xos,10,10000)  ||
      !sample_quality_worker(xos,11,1000)  ||
      !sample_quality_worker(xos,12,1000)  ||
      !sample_quality_worker(xos,13,1000)  ||
      !sample_quality_worker(xos,14,500)  ||
      !sample_quality_worker(xos,15,500)  ||
      !sample_quality_worker(xos,16,500) ||
      !sample_quality_worker(xos,17,100)  ||
      !sample_quality_worker(xos,18,100) ||
      !sample_quality_worker(xos,19,100)  ||
      !sample_quality_worker(xos,20,100)) {

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
  int seedi = time(0);

  xos = plc_xoshiro_init((uint64_t)(time(0)));
    
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
	 ,"xoshiro256+",seedi);

  if (!sample_quality_tests(xos)) { PASS = false; }

  plc_xoshiro_free(xos);
			 
			 
  printf("=========================================================================\n");

  if (PASS) {

    printf("Random SAP Test Suite: PASS\n");
    exit(0);

  } else {

    printf("Random SAP Test Suite: FAIL\n");
    exit(1);

  }

}
