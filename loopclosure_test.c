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

bool test_plc_delete_arc() {

  plCurve *testL;
  int     nc = 3, nv[3] = {92,189,341}, cc[3] = {0,0,0};
  bool    open[3] = { false, true, false };
  int     vt,cp;
  plCurve *deletedL;

  testL = plc_new(nc,nv,open,cc); 

  for(cp=0;cp<3;cp++) {

    for(vt=0;vt<testL->cp[cp].nv;vt++) {

      testL->cp[cp].vt[vt] = plc_build_vect(cp,vt,vt);
      
    }

  }

  plc_fix_wrap(testL);

  /* Now start the tests. */

  printf("plc_delete_arc test suite\n"
	 "Starting with a 3 arc curve, closed, open, closed,\n" 
	 "with %d, %d, %d vertices.\n"
	 "--------------------------------------------------\n"
	 ,testL->cp[0].nv,testL->cp[1].nv,testL->cp[2].nv);

  
  printf("Closed curve operations\n"
	 "Delete verts 0..%d (entire component)...",testL->cp[0].nv-1);

  deletedL = plc_delete_arc(testL,0,0,testL->cp[0].nv-1);
  
  if (deletedL->nc != 2) { printf("FAIL (didn't delete cp 0)\n"); return false; }
  if (deletedL->cp[0].nv != testL->cp[1].nv || deletedL->cp[1].nv != testL->cp[2].nv ||
      deletedL->cp[0].open != testL->cp[1].open || deletedL->cp[1].open != testL->cp[2].open) {
    printf("FAIL (didn't shift other components after deleting 0).\n"); return false; }

  plc_free(deletedL);

  printf("PASS \n");
  
  printf("Delete verts 5..4 (entire component)...");

  deletedL = plc_delete_arc(testL,0,5,4);
  
  if (deletedL->nc != 2) { printf("FAIL (didn't delete cp 0)\n"); return false; }
  if (deletedL->cp[0].nv != testL->cp[1].nv || deletedL->cp[1].nv != testL->cp[2].nv ||
      deletedL->cp[0].open != testL->cp[1].open || deletedL->cp[1].open != testL->cp[2].open) {
    printf("FAIL (didn't shift other components after deleting 0).\n"); return false; }

  plc_free(deletedL);

  printf("PASS \n");

  printf("Delete verts -10..512 (entire component)...");

  deletedL = plc_delete_arc(testL,0,-10,512);
  
  if (deletedL->nc != 2) { printf("FAIL (didn't delete cp 0)\n"); return false; }
  if (deletedL->cp[0].nv != testL->cp[1].nv || deletedL->cp[1].nv != testL->cp[2].nv ||
      deletedL->cp[0].open != testL->cp[1].open || deletedL->cp[1].open != testL->cp[2].open) {
    printf("FAIL (didn't shift other components after deleting 0).\n"); return false; }

  plc_free(deletedL);

  printf("PASS \n");

  printf("Delete verts 17..17 (single vertex)...");

  deletedL = plc_delete_arc(testL,0,17,17);
  
  if (deletedL->nc != 3) { printf("FAIL (deleted cp 0)\n"); return false; }
  if (deletedL->cp[0].nv   != testL->cp[0].nv-1 || deletedL->cp[1].nv != testL->cp[1].nv || deletedL->cp[2].nv != testL->cp[2].nv ||
      deletedL->cp[0].open != true              || deletedL->cp[1].open != testL->cp[1].open || deletedL->cp[2].open != testL->cp[2].open) {
    printf("FAIL (component 0 is not open, has incorrect # of verts, or other comps changed).\n"); return false; }

  plc_free(deletedL);

  printf("PASS \n");

  printf("Delete verts 32..40 (conventional arc, 9 vertices lost, %d should survive)...",testL->cp[0].nv-9);

  deletedL = plc_delete_arc(testL,0,32,40);
  
  if (deletedL->nc != 3) { printf("FAIL (deleted cp 0)\n"); return false; }
  if (deletedL->cp[0].nv   != testL->cp[0].nv-9 || deletedL->cp[1].nv != testL->cp[1].nv || deletedL->cp[2].nv != testL->cp[2].nv ||
      deletedL->cp[0].open != true              || deletedL->cp[1].open != testL->cp[1].open || deletedL->cp[2].open != testL->cp[2].open) {
    printf("FAIL (component 0 is not open, has incorrect # of verts, or other comps changed).\n"); return false; }

  plc_free(deletedL);

  printf("PASS \n");

  printf("Delete verts 40..32 (wraparound arc, 7 vertices should survive)...");

  deletedL = plc_delete_arc(testL,0,40,32);
  
  if (deletedL->nc != 3) { printf("FAIL (deleted cp 0)\n"); return false; }
  if (deletedL->cp[0].nv   != 7                 || deletedL->cp[1].nv != testL->cp[1].nv || deletedL->cp[2].nv != testL->cp[2].nv ||
      deletedL->cp[0].open != true              || deletedL->cp[1].open != testL->cp[1].open || deletedL->cp[2].open != testL->cp[2].open) {
    printf("FAIL (component 0 is not open, has incorrect # of verts, or other comps changed).\n"); return false; }

  plc_free(deletedL);

  printf("PASS \n");

  printf("Delete verts 40..38 (wraparound arc, 1 vertex should survive)...");

  deletedL = plc_delete_arc(testL,0,40,38);
  
  if (deletedL->nc != 3) { printf("FAIL (deleted cp 0)\n"); return false; }
  if (deletedL->cp[0].nv   != 1                 || deletedL->cp[1].nv != testL->cp[1].nv || deletedL->cp[2].nv != testL->cp[2].nv ||
      deletedL->cp[0].open != true              || deletedL->cp[1].open != testL->cp[1].open || deletedL->cp[2].open != testL->cp[2].open) {
    printf("FAIL (component 0 is not open, has incorrect # of verts, or other comps changed).\n"); return false; }

  plc_free(deletedL);

  printf("PASS \n");

  printf("Delete verts 40..%d (conventional arc to end of component, 40 verts survive)...",testL->cp[0].nv-1);

  deletedL = plc_delete_arc(testL,0,40,testL->cp[0].nv-1);
  
  if (deletedL->nc != 3) { printf("FAIL (deleted cp 0)\n"); return false; }
  if (deletedL->cp[0].nv   != 40                 || deletedL->cp[1].nv != testL->cp[1].nv || deletedL->cp[2].nv != testL->cp[2].nv ||
      deletedL->cp[0].open != true              || deletedL->cp[1].open != testL->cp[1].open || deletedL->cp[2].open != testL->cp[2].open) {
    printf("FAIL (component 0 is not open, has incorrect # of verts, or other comps changed).\n"); return false; }

  plc_free(deletedL);

  printf("PASS \n");

  printf("Delete verts 0..53 (conventional arc, deleting head of component, %d verts survive)...",testL->cp[0].nv-54);

  deletedL = plc_delete_arc(testL,0,0,53);
  
  if (deletedL->nc != 3) { printf("FAIL (deleted cp 0)\n"); return false; }
  if (deletedL->cp[0].nv   != testL->cp[0].nv-54|| deletedL->cp[1].nv != testL->cp[1].nv || deletedL->cp[2].nv != testL->cp[2].nv ||
      deletedL->cp[0].open != true  || deletedL->cp[1].open != testL->cp[1].open || deletedL->cp[2].open != testL->cp[2].open) {
    printf("FAIL (component 0 is not open, has incorrect # of verts, or other comps changed).\n"); return false; }

  plc_free(deletedL);

  printf("PASS \n");

  printf("\n"
	 "Open component operations\n"
	 "Delete vertices -10..512 (entire component) ... ");

  deletedL = plc_delete_arc(testL,1,-10,512);
  
  if (deletedL->nc != 2) { printf("FAIL (didn't delete cp 1)\n"); return false; }
  if (deletedL->cp[0].nv   != testL->cp[0].nv  || deletedL->cp[1].nv   != testL->cp[2].nv  ||
      deletedL->cp[0].open != testL->cp[0].open|| deletedL->cp[1].open != testL->cp[2].open ) {

    printf("FAIL (other comps changed).\n"); return false; }

  plc_free(deletedL);

  printf("PASS \n");

  printf("Delete vertices 0..512 (entire component) ... ");

  deletedL = plc_delete_arc(testL,1,0,512);
  
  if (deletedL->nc != 2) { printf("FAIL (didn't delete cp 1)\n"); return false; }
  if (deletedL->cp[0].nv   != testL->cp[0].nv  || deletedL->cp[1].nv   != testL->cp[2].nv  ||
      deletedL->cp[0].open != testL->cp[0].open|| deletedL->cp[1].open != testL->cp[2].open ) {

    printf("FAIL (other comps changed).\n"); return false; }

  plc_free(deletedL);

  printf("PASS \n");

  printf("Delete vertices 0..%d (entire component) ... ",testL->cp[1].nv-1);

  deletedL = plc_delete_arc(testL,1,0,testL->cp[1].nv-1);
  
  if (deletedL->nc != 2) { printf("FAIL (didn't delete cp 1)\n"); return false; }
  if (deletedL->cp[0].nv   != testL->cp[0].nv  || deletedL->cp[1].nv   != testL->cp[2].nv  ||
      deletedL->cp[0].open != testL->cp[0].open|| deletedL->cp[1].open != testL->cp[2].open ) {

    printf("FAIL (other comps changed).\n"); return false; }

  plc_free(deletedL);

  printf("PASS \n");

  printf("Delete vertices 0..17 (shorten component 1, %d verts should survive) ... ", testL->cp[1].nv - 18);

  deletedL = plc_delete_arc(testL,1,0,17);
  
  if (deletedL->nc != 3) { printf("FAIL (deleted cp 1)\n"); return false; }
  if (deletedL->cp[0].nv   != testL->cp[0].nv      || deletedL->cp[0].open != testL->cp[0].open ||
      deletedL->cp[1].nv   != testL->cp[1].nv - 18 || deletedL->cp[1].open != true ||
      deletedL->cp[2].nv   != testL->cp[2].nv      || deletedL->cp[2].open != testL->cp[2].open ) {

    printf("FAIL (other comps changed).\n"); return false; }

  plc_free(deletedL);

  printf("PASS \n");

  printf("Delete vertices -15..17 (shorten component 1, %d verts should survive) ... ", testL->cp[1].nv - 18);

  deletedL = plc_delete_arc(testL,1,-15,17);
  
  if (deletedL->nc != 3) { printf("FAIL (deleted cp 1)\n"); return false; }
  if (deletedL->cp[0].nv   != testL->cp[0].nv      || deletedL->cp[0].open != testL->cp[0].open ||
      deletedL->cp[1].nv   != testL->cp[1].nv - 18 || deletedL->cp[1].open != true ||
      deletedL->cp[2].nv   != testL->cp[2].nv      || deletedL->cp[2].open != testL->cp[2].open ) {

    printf("FAIL (other comps changed).\n"); return false; }

  plc_free(deletedL);

  printf("PASS \n");

  printf("Delete vertices 15..512 (shorten component 1, %d vertices should remain) ... ",15);

  deletedL = plc_delete_arc(testL,1,15,512);
  
  if (deletedL->nc != 3) { printf("FAIL (deleted cp 1)\n"); return false; }
  if (deletedL->cp[0].nv   != testL->cp[0].nv      || deletedL->cp[0].open != testL->cp[0].open ||
      deletedL->cp[1].nv   != 15 || deletedL->cp[1].open != true ||
      deletedL->cp[2].nv   != testL->cp[2].nv      || deletedL->cp[2].open != testL->cp[2].open ) {

    printf("FAIL (other comps changed).\n"); return false; }

  plc_free(deletedL);

  printf("PASS \n");

  printf("Delete vertices 42..%d (shorten component 1, %d vertices should remain) ... ",testL->cp[1].nv,42);

  deletedL = plc_delete_arc(testL,1,42,testL->cp[1].nv);
  
  if (deletedL->nc != 3) { printf("FAIL (deleted cp 1)\n"); return false; }
  if (deletedL->cp[0].nv   != testL->cp[0].nv      || deletedL->cp[0].open != testL->cp[0].open ||
      deletedL->cp[1].nv   != 42                   || deletedL->cp[1].open != true ||
      deletedL->cp[2].nv   != testL->cp[2].nv      || deletedL->cp[2].open != testL->cp[2].open ) {

    printf("FAIL (other comps changed).\n"); return false; }

  plc_free(deletedL);

  printf("PASS \n");

  printf("Delete vertices 10..42 (split component 1, should have 10 and %d vertices in cps 1 and 2) ... ",testL->cp[1].nv-43);

  deletedL = plc_delete_arc(testL,1,10,42);
  
  if (deletedL->nc != 4) { printf("FAIL (didn't split cp 1)\n"); return false; }
  
  if (deletedL->cp[0].nv   != testL->cp[0].nv      || deletedL->cp[0].open != testL->cp[0].open ||
      deletedL->cp[1].nv   != 10                   || deletedL->cp[1].open != true ||
      deletedL->cp[2].nv   != testL->cp[1].nv-43   || deletedL->cp[2].open != true || 
      deletedL->cp[3].nv   != testL->cp[2].nv      || deletedL->cp[3].open != testL->cp[2].open)
    {

    printf("FAIL (incorrect # of verts in split, or other components changed).\n"); return false; 
    
    }

  plc_free(deletedL);
  plc_free(testL);

  printf("PASS \n");
  
  printf("-------------------------------------------\n");
  printf("plc_delete_arc: PASS \n\n");

  fflush(stdout);

  return true;

}

bool test_plc_loop_closure(gsl_rng *r, plCurve *L, int cp, int nVerts, int nPolygons,bool diskout) {

  clock_t start,end;
  double cpu_time_used;
  int j;
  plCurve *test;
  bool PASS = true;
  bool localPASS;

  double *lengths;
  lengths = calloc(L->nc,sizeof(double));

  printf("\n"
	 "Timing, length, and vertex count check for loop closures\n\n");

  printf("Verts   Samples    Length (desired)    Timing          Result\n");
  printf("-------------------------------------------------------------\n");

  if (PAPERMODE) {
    fprintf(outfile,"LoopClosureTimingData ");
  }

  localPASS = true;

  start = clock();
  
  for(j=0;j<nPolygons;j++) {
    
    test = plc_loop_closure(r,cp,L,nVerts);
    
    if (test == NULL) {
      
      localPASS = false;
      
    } else {
      
      plc_arclength(test,lengths);
      
      if (fabs(lengths[cp] - 2.0)  > 1e-10) { localPASS = false; }
      if (test->cp[cp].nv != nVerts) {localPASS = false; }

      if (diskout && j < 10) {

	char fname[256];
	FILE *outfile;
	
	sprintf(fname,"random_closure_%d.vect",j);
	outfile = fopen(fname,"w");
	plc_write(outfile,test);
	fclose(outfile);
	
      }
      
      plc_free(test);
      
    }
    
  }
  
  end = clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;
  cpu_time_used /= (double)(nPolygons);
  
  printf("%-5d   %-5d      %-5g  (%-5g)      %-#8.8g  ",nVerts,nPolygons,lengths[cp],2.0,cpu_time_used);
  
  if (PAPERMODE) {
    fprintf(outfile,", %d , %g ",nVerts,cpu_time_used);
  }
  
  if (localPASS) {
    printf("   pass \n");
  } else {
    printf("   FAIL \n"); PASS = false;
  }

  if (diskout) { printf("\n"
			"Wrote files to random_closure_0--%d.vect\n\n",(j < 9) ? j : 9); }
          
  printf("-------------------------------------------------------------\n\n");
  
  if (PAPERMODE) {
    fprintf(outfile,"\n");
  }

  free(lengths);
  
  return PASS;
  
}

bool circle_closure_test(gsl_rng *r,int nEdges,int nTrials,bool diskout) {
  
  double pi = 3.141592653589793;
  plCurve *circle;
  int nc = 1,nv = {nEdges},cc = {0};
  bool open = false;
  double theta = 2.0*pi/(double)(nEdges);
  int vt;

  printf("Circle closure test with %d edges and %d trials...",nEdges,nTrials);

  circle = plc_new(nc,&nv,&open,&cc);
  for(vt=0;vt<nEdges;vt++) {
    circle->cp[0].vt[vt] = plc_build_vect(cos(vt*theta),sin(vt*theta),0);
  }
  plc_fix_wrap(circle);

  double len;

  len = plc_arclength(circle,NULL);
  plc_scale(circle,1.0/len);

  if (fabs(plc_arclength(circle,NULL) - 1.0) > 1e-10) {

    printf("FAIL: Couldn't construct length 2.0 circle with %d edges.\n",nEdges);
    return false;

  }

  plCurve *subarc;
  subarc = plc_delete_arc(circle,0,(int)(nEdges/2.0),(int)((3.0/4.0)*nEdges));

  plCurve *circleclosure;
  int trial;

  for(trial=0;trial<nTrials;trial++) {

    circleclosure = plc_loop_closure(r,0,subarc,nEdges);

    if (fabs(plc_arclength(circleclosure,NULL) - 2.0) > 1e-10 ||
	circleclosure->cp[0].open ||
	circleclosure->cp[0].nv != nEdges) {

      printf("FAIL: Closure arclength, open, or verts error at trial %d.\n",trial);
    
    }

    if (diskout) {

      FILE *outfile;
      char  circname[256];
      
      sprintf(circname,"circle_closure_%d.vect",trial);
      outfile = fopen(circname,"w");
      plc_write(outfile,circleclosure);
      fclose(outfile);
      
    }

    plc_free(circleclosure);

  }

  if (diskout) {

    printf("(wrote files to circle_closure_0--%d.vect)...",nTrials-1);
  
  }
  
  plc_free(subarc);
  plc_free(circle);

  printf("pass\n");
  return true;

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

  int    nPolygons = 50;

  if (PAPERMODE) { 
    outfile = fopen("loop_closure_timing.csv","w");
    timestamp(outfile);
  }

  if (!timing_and_selftest(nPolygons,nSizes,sizes,lengths,ftcs)) {
    PASS = false;
  }
  
  if (PAPERMODE) { fclose(outfile); }

  if (!test_plc_delete_arc()) { PASS = false; }

  plCurve *random,*subarc;
  random = plc_random_closed_polygon(r,513);
  subarc = plc_delete_arc(random,0,0,147);

  if (!test_plc_loop_closure(r,subarc,0,513,5000,false)) { PASS = false; }

  plc_free(random);
  plc_free(subarc);

  random = plc_random_closed_polygon(r,10);
  subarc = plc_delete_arc(random,0,5,10);

  if (!test_plc_loop_closure(r,subarc,0,10,5000,false)) { PASS = false; }

  plc_free(random);
  plc_free(subarc);

  random = plc_random_closed_polygon(r,129);
  subarc = plc_delete_arc(random,0,0,100);

  if (!test_plc_loop_closure(r,subarc,0,129,5000,false)) { PASS = false; }

  plc_free(random);
  plc_free(subarc);


  if (!circle_closure_test(r,8,1,false) ||
      !circle_closure_test(r,192,1000,false) ||
      !circle_closure_test(r,593,5000,false)) { PASS = false; }
  
  gsl_rng_free(r);

  printf("======================================\n");

  if (PASS) {

    printf("LOOPCLOSURE_TEST: PASS\n");
    fflush(stdout);
    exit(0);

  } else {

    printf("LOOPCLOSURE_TEST: FAIL\n");
    fflush(stdout);
    exit(1);

  }

}
