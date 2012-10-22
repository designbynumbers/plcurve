/* 

   loopclosure_test.c : Test code for the loop closure functions in plCurve. 

*/

#include<plCurve.h>
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

gsl_rng *r;  /* The global random number generator */

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

plCurve *plc_random_ftc_internal(gsl_rng *r, int nEdges, double length, double ftc, bool selfcheck);

bool timing_and_selftest(int nPolygons,int nSizes,int *Sizes,double *lengths,double *ftcs) {

  clock_t start,end;
  double cpu_time_used;
  int i,j;
  plCurve *test;
  bool PASS = true;
  double length = 0, ftc = 0;
  bool localPASS;

  printf("\n"
	 "Timing and length, with selftest, for fixed failure-to-close polygons\n");

  printf("Verts   Samples    Length (desired)   FTC (desired)   Timing          Result\n");
  printf("----------------------------------------------------------------------------\n");

  if (PAPERMODE) {
    fprintf(outfile,"TimingDataFixedFailureToCloseSpace ");
  }

  for(i=0;i<nSizes;i++) {

    localPASS = true;

    start = clock();

    for(j=0;j<nPolygons;j++) {

      test = plc_random_ftc_internal(r,Sizes[i],lengths[i],ftcs[i],true);

      if (test == NULL) {

	length = 0; ftc = 0;
	localPASS = false;

      } else {
   
	length = plc_arclength(test,NULL);
	ftc    = plc_distance(test->cp[0].vt[0],test->cp[0].vt[test->cp[0].nv-1]);
	
	if (fabs(length - lengths[i])  > 1e-10) { localPASS = false; }
	if (fabs(ftc - ftcs[i]) > 1e-10) {localPASS = false; }
	
	plc_free(test);

      }

    }

    end = clock();
    cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;
    cpu_time_used /= (double)(nPolygons);
    
    printf("%-5d   %-5d      %-5g (%-5g)      %-3.2f (%-3.2f)     %-#8.8g  ",Sizes[i],nPolygons,length,lengths[i],ftc,ftcs[i],cpu_time_used);

    if (PAPERMODE) {
      fprintf(outfile,", %d , %g ",Sizes[i],cpu_time_used);
    }
    
    if (localPASS) {
      printf(" pass \n");
    } else {
      printf(" FAIL \n"); PASS = false;
    }
    
  }

  printf("----------------------------------------------------------------------------\n\n");

  if (PAPERMODE) {
    fprintf(outfile,"\n");
  }
  
  return PASS;

}




int main(int argc, char *argv[]) {

  bool PASS = {true};

  const gsl_rng_type * T;
     
  gsl_rng_env_setup();
  
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  int seed = time(0);
  gsl_rng_set(r,seed);

  if (argc > 1) { PAPERMODE = true; }

  printf("Loop closure tests \n"
	 "------------------------------- \n"
	 "plCurve generates loop closures for open polygons using the algorithm\n"
	 "of Cantarella, Rawdon, Shonkwiler (et. al.) [2013]. This program tests\n"
	 "the loop closure functionality.\n"
	 "\n"
	 "This test suite code is using the random number generator %s \n"
	 "with seed %d.\n"
	 ,gsl_rng_name(r),seed);

  double lengths[10] = {1.0,0.25,0.33,0.05,1.9};
  double ftcs[10]    = {1.0,0.1 ,0.0 ,0.02,0.1};
  int    sizes[10]   = {1  ,2   ,100 ,500 ,1000};
  int    nSizes = 5;

  int    nPolygons = 5000;

   if (PAPERMODE) { 
    outfile = fopen("loop_closure_timing.csv","w");
    timestamp(outfile);
  }

   if (!timing_and_selftest(nPolygons,nSizes,sizes,lengths,ftcs)) {
    PASS = false;
  }

  if (PAPERMODE) { fclose(outfile); }

  gsl_rng_free(r);

  printf("======================================\n");

  if (PASS) {

    printf("LOOPCLOSURE_TEST: PASS\n");
    exit(0);

  } else {

    printf("LOOPCLOSURE_TEST: FAIL\n");
    exit(1);

  }

}
