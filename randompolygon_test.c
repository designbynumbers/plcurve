/* 

   randompolygon_test.c : Test code for the symmetry functions in plCurve. 

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

plCurve *plc_random_polygon_selfcheck(int nEdges,bool selfcheck); /* A private version which turns on debugging code. */

double *gaussian_array(int n);

void test_gaussianarray()
{
  int N = 500000;
  double *testvar;
  clock_t start,end;
  double cpu_time_used;

  int SampNum;
  
  printf("Testing Box-Muller gaussian distribution code.\n\n");

  
  for(SampNum = 0; SampNum < 10; SampNum++) {

    
    printf("Generating array of %d samples...",2*N);

    start = clock();
    testvar = gaussian_array(N);
    end = clock();
    cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;
    
    printf("done. (%3.2g sec). ",cpu_time_used);
    
    double mean=0,moments[4];
    int i;
    
    for(i=0;i<2*N;i++) { mean += testvar[i]; }
    mean /= 2*N;
    
    for(i=0;i<2*N;i++) {
      
      moments[2] += pow(testvar[i] - mean,2.0);
      moments[3] += pow(testvar[i] - mean,3.0);
      moments[4] += pow(testvar[i] - mean,4.0);
      
    }
    
    moments[2] /= 2*N; moments[3] /= 2*N; moments[4] /= 2*N;
    
    double skewness,kurtosis;
    
    skewness = moments[3]/pow(moments[2],1.5);
    kurtosis = moments[4]/pow(moments[2],2.0);
    
    double JB; /* Jarque-Bera statistic, from Wikipedia */
    
    JB = (1.0*N/6.0)*(pow(skewness,2.0) + (1.0/4.0)*pow(kurtosis - 3.0,2.0));
    
    printf(" S = %+7.5f. K = %+7.5f. JB statistic = %7.5g.\n",skewness,kurtosis,JB);
    
    free(testvar);

  }
    
  printf("\n");

}

complex double HermitianDot(complex double *A,complex double *B,int n);

void test_hdot()
{
  complex double A[3] = {1.0 + 3.0*I,-2.0 - 1.9*I,1.0 - 5.0*I};
  complex double B[3] = {2.0 - 2.0*I,5.0 - 2.3*I, 3 + 9*I};
  complex double hDot;

  printf("Testing hermitian dot product\n\n"
	 "<(1.0 + 3.0*I,-2.0 - 1.9*I,1.0 - 5.0*I),(2.0 - 2.0*I,5.0 - 2.3*I, 3 + 9*I)> = -51.63 - 30.1 I...");

  hDot = HermitianDot(A,B,3);   /* Should be -51.63 - 30.1 I according to Mathematica */

  if (fabs(creal(hDot) - (-51.63)) > 1e-10 || fabs(cimag(hDot) - (-30.1)) > 1e-10) { 

    printf("FAIL. Result is %g + %g I. \n",creal(hDot),cimag(hDot));
    exit(1);

  } else {

    printf("pass. Result is %g + %g I. \n",creal(hDot),cimag(hDot));

  }

  printf("\n");

} 

 /* Generate chordlength data for N samples of NEDGE polygons at skips J, 2J, 3J, .... */

double *mean_squared_chordlengths(int N,int NEDGE,int J)

 {
   plCurve *test;
   double *chordlengthdata;
   chordlengthdata = calloc(NEDGE,sizeof(double));

   int skip;
   int i;
   
   for(i=0;i<N;i++) {

    test = plc_random_closed_polygon(NEDGE);

    for(skip=J;skip<NEDGE;skip+=J) {

	chordlengthdata[skip] += plc_M_sq_dist(test->cp[0].vt[0],test->cp[0].vt[skip]);

    }

    plc_free(test);

   }
   
   for(i=0;i<NEDGE;i++) { chordlengthdata[i] /= N; }
   
   return chordlengthdata;

 }

 /* Generate chordlength data for N samples of NEDGE polygons at skips J, 2J, 3J, .... */

double *open_mean_squared_chordlengths(int N,int NEDGE,int J)

 {
   plCurve *test;
   double *chordlengthdata;
   chordlengthdata = calloc(NEDGE,sizeof(double));

   int skip;
   int i;
   
   for(i=0;i<N;i++) {

    test = plc_random_open_polygon(NEDGE);

    for(skip=J;skip<NEDGE;skip+=J) {

	chordlengthdata[skip] += plc_M_sq_dist(test->cp[0].vt[0],test->cp[0].vt[skip]);

    }

    plc_free(test);

   }
   
   for(i=0;i<NEDGE;i++) { chordlengthdata[i] /= N; }
   
   return chordlengthdata;

 }

bool timing_and_length2_test(int nPolygons,int nSizes,int *Sizes,plCurve *polygen(int nEdges),char *typestring) {

  clock_t start,end;
  double cpu_time_used;
  int i,j;
  plCurve *test;
  bool PASS = true;
  double length = 0;
  bool localPASS;

  printf("Timing and length test for %s polygons\n",typestring);

  printf("Verts        Samples       Length    Mean Time to generate  Result\n");
  printf("------------------------------------------------------------------\n");

  for(i=0;i<nSizes;i++) {

    localPASS = true;

    start = clock();

    for(j=0;j<nPolygons;j++) {

      test = polygen(Sizes[i]); 
      length = plc_arclength(test,NULL);
      
      if (fabs(length - 2.0)  > 1e-10) { localPASS = false; }
      plc_free(test);

    }

    end = clock();
    cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;
    cpu_time_used /= (double)(nPolygons);
    
    printf("%-10d   %-5d         %-5g     %-13.4g          ",Sizes[i],nPolygons,length,cpu_time_used);
    
    if (localPASS) {
      printf(" pass \n");
    } else {
      printf(" FAIL \n"); PASS = false;
    }
    
  }

  printf("------------------------------------------------------------------\n\n");
  
  return PASS;

}

bool chordlength_test(int nPolygons,int nEdges,int nSkips,int *Skips,plCurve *polygen(int nEdges),double prediction(int n, int k),char *prediction_string,char *polygon_type)

/* Given a generation function and a prediction for average chordlength in terms of number of edges and skip, do test. */

{
  int i;
  bool PASS = true;
  clock_t start,end;
  double cpu_time_used;

  printf("Mean Squared Chordlength test for %d %d-edge %s at %d skips.\n\t Skip list: ",
	 nPolygons,nEdges,polygon_type,nSkips);
  for(i=0;i<nSkips;i++) {
    printf("%5d ",Skips[i]);
    if (i % 10 == 0 && i > 0) { printf("\n"); }
  }
  
  printf("\nComputing data...");
  fflush(stdout);
 
  start = clock();

  double *allchords,*thispolychords;
  allchords = calloc(sizeof(double),nSkips);
  
  for(i=0;i<nPolygons;i++) {
    
    /* Generate polygon and compute chord data */
    plCurve *L;
    L = polygen(nEdges);
    thispolychords = plc_mean_squared_chordlengths(L,0,Skips,nSkips);
    
    /* Add to running total */
    int j;
    for(j=0;j<nSkips;j++) {
      allchords[j] += thispolychords[j];
    }
    
    /* Free memory */
    free(thispolychords);
    plc_free(L);

  }

  for(i=0;i<nSkips;i++) { allchords[i] /= (double)(nPolygons); }
  end = clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("done (%g sec).\n\n",cpu_time_used); 

  /* Now display results. */

  printf("According to CDS, the predicted value of squared chord length \nfor %s at skip k is %s\n\n",polygon_type,prediction_string);
  printf("Skip      Computed Average    Predicted Average   Difference  Result\n"
	 "---------------------------------------------------------------------------\n");

  for(i=0;i<nSkips;i++) {
    
    double predicted  = prediction(nEdges,Skips[i]);
    double percenterr = 100*(fabs(allchords[i] - predicted)/predicted);

    printf("%-5d     %-13.7g       %-13.7g       %3.3f %%",Skips[i],allchords[i],predicted,percenterr);

    if (percenterr > 1.0) { printf("      FAIL (> 1%%).\n"); PASS = false; }
    else { printf("      pass (< 1%%).\n"); }
    
  }

  printf("---------------------------------------------------------------------------\n");
  printf("\n\n");
  free(allchords);
  return PASS;

}

/* These prediction functions are for chordlength_test */

double space_pol_prediction(int n,int k) {
  return (double)((n-k)*6*k)/(double)((n-1)*n*(n+1));
}
char space_pol_predstring[256] = "((n-k)/n) (6k/(n-1)(n+1))";

double space_arm_prediction(int n,int k) {
  return (double)(6*k)/(double)(n*(n+0.5));
}
char space_arm_predstring[256] = "6k/(n(n+1/2))";

double plane_pol_prediction(int n,int k) {
  return (double)((n-k)*8*k)/(double)((n-1)*n*(n+2));
}
char plane_pol_predstring[256] = "((n-k)/n) (8k/(n(n+2)))";

double plane_arm_prediction(int n,int k) {
  return (double)(8*k)/(double)(n*(n+1));
}
char plane_arm_predstring[256] = "8k/(n(n+1))";

bool gyradius_test(int nPolygons,int nSizes,int *Sizes,plCurve *polygen(int nEdges),double prediction(int n),char *prediction_string,char *polygon_type)

  /* Given a generation function and a prediction for gyradius in terms of n, do test. */

{
  int i,j;
  bool PASS = true;
  clock_t start,end;
  double cpu_time_used;

  printf("(Squared) Radius of Gyration test for %d %s polygons at %d numbers of verts.\n\tSizelist: ",
	 nPolygons,polygon_type,nSizes);
  for(i=0;i<nSizes;i++) {
    printf("%5d ",Sizes[i]);
    if (i % 10 == 0 && i > 0) { printf("\n"); }
  }
  
  printf("\nComputing data...");
  fflush(stdout);
 
  start = clock();

  double *allgyradius;
  allgyradius = calloc(sizeof(double),nSizes);

  for(j=0;j<nSizes;j++) {
  
    for(i=0;i<nPolygons;i++) {
    
      /* Generate polygon and compute gyradius */
      plCurve *L;
      L = polygen(Sizes[j]);
      allgyradius[j] += plc_gyradius(L);
  
      /* Free memory */
      plc_free(L);
      
    }

    allgyradius[j] /= (double)(nPolygons);

  }
  end = clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("done (%g sec).\n\n",cpu_time_used); 

  /* Now display results. */

  printf("According to CDS, the expected value of gyradius for n edges is is %s\n\n",prediction_string);
  printf("n (num verts)    Mean Gyradius    Predicted Mean Gyradius   Difference  Result\n"
	 "------------------------------------------------------------------------------\n");

  for(i=0;i<nSizes;i++) {
    
    double predicted  = prediction(Sizes[i]);
    double percenterr = 100*(fabs(allgyradius[i] - predicted)/predicted);

    printf("%-5d            %-13.7g    %-13.7g             %7.3f %%",Sizes[i],allgyradius[i],predicted,percenterr);

    if (percenterr > 1.0) { printf("   FAIL (> 1%%).\n"); PASS = false; }
    else { printf("   pass (< 1%%).\n"); }
    
  }

  printf("\n\n");
  free(allgyradius);

  return PASS;

}

/* These prediction functions are for gyradius_test */

double space_pol_gyradius_prediction(int n) {
  return (double)(1.0)/(double)(2.0*n);
}
char space_pol_gyradius_predstring[256] = "(1/2) (1/n)";

double space_arm_gyradius_prediction(int n) {
  return (double)(n+2)/(double)((n+1)*(n+0.5));
}
char space_arm_gyradius_predstring[256] = "(n+2)/(n+1)(n+1/2)";

double plane_pol_gyradius_prediction(int n) {
  return (double)(2.0*(n+1))/(double)(3*n*(n+2));
}
char plane_pol_gyradius_predstring[256] = "(2/3) (n+1)/(n(n+2))";

double plane_arm_gyradius_prediction(int n) {
  return (double)(4.0*(n+2))/(double)(3*(n+1)*(n+1));
}
char plane_arm_gyradius_predstring[256] = "(4/3) (n+2)/(n+1)^2";


bool ftc_test(int nPolygons,int nSizes,int *Sizes,plCurve *polygen(int nEdges),double prediction(int n),char *prediction_string,char *polygon_type)
  
/* Tests to see if the given size of polygon has correct failure to close. */

{
  int i,j;
  bool PASS = true;
  clock_t start,end;
  double cpu_time_used;

  printf("Failure-to-close test for %d %s polygons at %d numbers of verts.\n\tSizelist: ",
	 nPolygons,polygon_type,nSizes);
  for(i=0;i<nSizes;i++) {
    printf("%5d ",Sizes[i]);
    if (i % 10 == 0 && i > 0) { printf("\n"); }
  }
  
  printf("\nComputing data...");
  fflush(stdout);
 
  start = clock();

  double *all_lastedgelength;
  all_lastedgelength = calloc(sizeof(double),nSizes);

  for(j=0;j<nSizes;j++) {
  
    for(i=0;i<nPolygons;i++) {
    
      /* Generate polygon and compute gyradius */
      plCurve *L;
      L = polygen(Sizes[j]);
      all_lastedgelength[j] += plc_M_sq_dist(L->cp[0].vt[0],L->cp[0].vt[L->cp[0].nv-1]);
  
      /* Free memory */
      plc_free(L);
      
    }

    all_lastedgelength[j] /= (double)(nPolygons);

  }
  end = clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("done (%g sec).\n\n",cpu_time_used); 

  /* Now display results. */

  printf("According to CDS, the expected value of squared last edgelength for n edge %s polygons is %s\n\n",polygon_type,prediction_string);
  printf("n (num verts)    Mean Sq FTC    Predicted Mean Sq FTC   Difference  Result\n"
	 "--------------------------------------------------------------------------\n");

  for(i=0;i<nSizes;i++) {
    
    double predicted  = prediction(Sizes[i]);
    double percenterr = 100*(fabs(all_lastedgelength[i] - predicted)/predicted);

    printf("%-5d            %-13.7g  %-13.7g           %3.3f %%",Sizes[i],all_lastedgelength[i],predicted,percenterr);

    if (percenterr > 1.0) { printf("     FAIL (> 1%%).\n"); PASS = false; }
    else { printf("     pass (< 1%%).\n"); }
    
  }

  printf("\n\n");

  free(all_lastedgelength);
  return PASS;

}
  
/* These prediction functions are for FTC_test */

double space_pol_ftc_prediction(int n) {
  return (double)(6.0)/(double)(n*(n+1));
}
char space_pol_ftc_predstring[256] = "6/n(n+1)";

double space_arm_ftc_prediction(int n) {
  return (double)(6.0*n)/(double)(n*(n+1));
}
char space_arm_ftc_predstring[256] = "6/(n+1)";

double plane_pol_ftc_prediction(int n) {
  return (double)(8.0)/(double)(n*(n+2));
}
char plane_pol_ftc_predstring[256] = "8/(n(n+2))";

double plane_arm_ftc_prediction(int n) {
   return (double)(8.0*n)/(double)(n*(n+2));
}
char plane_arm_ftc_predstring[256] = "8/(n+2)";
  
int main(int argc, char *argv[]) {

  bool PASS = {true};
  srand48(time(0));
  bool PAPERMODE = {false};

  if (argc > 1) { PAPERMODE = true; }

  printf("Random Polygon Generation Tests \n"
	 "------------------------------- \n"
	 "plCurve generates random polygons in four classes by direct sampling from the symmetric measure of Cantarella, Deguchi, Shonkwiler:\n"
	 "\n"
	 "\t Type                    Call                   \n"
	 "\t -------------------------------------------------------------------\n"
	 "\t closed space polygons   plc_random_closed_polygon(int nEdges) \n"
	 "\t open space polygons     plc_random_open_polygon(int nEdges) \n"
	 "\t closed plane polygons   plc_random_closed_plane_polygon(int nEdges) \n"
	 "\t open plane polygons     plc_random_open_plane_polygon(int nEdges) \n"
	 "\n"
	 "This program tests the polygons generated against theoretical calculations about the symmetric measure. All polygons generated\n"
	 "are guaranteed to have length 2. (Rescaling, if desired, requires you to pick a probability distribution function on polygon length and so is\n"
	 "left to the user).\n\n");

  test_gaussianarray();
  test_hdot();

  /* Timing and Length = 2.0 tests. */

  int Sizes[100] = {20,200,2000,20000};
  int nSizes = 4;
  int nPolygons = 40;

  if (!timing_and_length2_test(nPolygons,nSizes,Sizes,plc_random_closed_polygon_selfcheck,"closed space polygon")
      || !timing_and_length2_test(nPolygons,nSizes,Sizes,plc_random_open_polygon_selfcheck,"open space polygon")
      || !timing_and_length2_test(nPolygons,nSizes,Sizes,plc_random_closed_plane_polygon_selfcheck,"closed plane polygon")
      || !timing_and_length2_test(nPolygons,nSizes,Sizes,plc_random_open_plane_polygon_selfcheck,"open plane polygon")) {
    PASS = false;
  }

  /* Mean squared chordlength tests. */

  int Skips[100] = {50,100,150,250,300,350,450};
  int nSkips = 7;
  nPolygons = 200;
  int nEdges = 500;

  if (!chordlength_test(nPolygons,nEdges,nSkips,Skips,plc_random_closed_polygon,space_pol_prediction,space_pol_predstring,"closed space polygon")
      || !chordlength_test(nPolygons,nEdges,nSkips,Skips,plc_random_open_polygon,space_arm_prediction,space_arm_predstring,"open space polygon")
      || !chordlength_test(nPolygons,nEdges,nSkips,Skips,plc_random_closed_plane_polygon,plane_pol_prediction,plane_pol_predstring,"closed plane polygon")
      || !chordlength_test(nPolygons,nEdges,nSkips,Skips,plc_random_open_plane_polygon,plane_arm_prediction,plane_arm_predstring,"open plane polygon")) {
    PASS = false;
  }

  /* Gyradius tests. */

  int gySizes[100] = {100,200,300,400,500};
  int ngySizes = 5;
  nPolygons = 200;

  if (!gyradius_test(nPolygons,ngySizes,gySizes,plc_random_closed_polygon,space_pol_gyradius_prediction,space_pol_gyradius_predstring,"closed space polygon")
      || !gyradius_test(nPolygons,ngySizes,gySizes,plc_random_open_polygon,space_arm_gyradius_prediction,space_arm_gyradius_predstring,"open space polygon")
      || !gyradius_test(nPolygons,ngySizes,gySizes,plc_random_closed_plane_polygon,plane_pol_gyradius_prediction,plane_pol_gyradius_predstring,"closed plane polygon")
      || !gyradius_test(nPolygons,ngySizes,gySizes,plc_random_open_plane_polygon,plane_arm_gyradius_prediction,plane_arm_gyradius_predstring,"open plane polygon")) {
    PASS = false;
  };
  
   /* Failure to close tests. */

  int ftcSizes[100] = {100,200,300,400,500};
  int nftcSizes = 5;
  nPolygons = 40000;

  if (!ftc_test(nPolygons,nftcSizes,ftcSizes,plc_random_closed_polygon,space_pol_ftc_prediction,space_pol_ftc_predstring,"closed space polygon")
      || !ftc_test(nPolygons,nftcSizes,ftcSizes,plc_random_open_polygon,space_arm_ftc_prediction,space_arm_ftc_predstring,"open space polygon")
      || !ftc_test(nPolygons,nftcSizes,ftcSizes,plc_random_closed_plane_polygon,plane_pol_ftc_prediction,plane_pol_ftc_predstring,"closed plane polygon")
      || !ftc_test(nPolygons,nftcSizes,ftcSizes,plc_random_open_plane_polygon,plane_arm_ftc_prediction,plane_arm_ftc_predstring,"open plane polygon")) {
    PASS = false;
  };
 
   /*  outfile = fopen("space_arm_chord_table.tex","w"); */
/*     if (outfile == NULL) { printf("randompolygon_test: Couldn't open space_arm_chord_table.tex"); exit(1); } */
/*     fprintf(outfile,"%% Data for %d samples of %d edge OPEN SPACE polygons.\n",N,NEDGE); */
/*     fprintf(outfile,"%% Run generated:"); */

/*     time_t curtime; */
/*     struct tm *loctime; */
    
/*     /\* Get the current time. *\/ */
/*     curtime = time (NULL); */
    
/*     /\* Convert it to local time representation. *\/ */
/*     loctime = localtime (&curtime); */
    
/*     /\* Print out the date and time in the standard format. *\/ */
/*     fputs (asctime (loctime), outfile);  */

/*     fprintf(outfile,"%% Skip & Measured & Predicted & Percent Error\n"); */
  
/*   } */
 
  if (PASS) { exit(0); } else { exit(1); }

}
