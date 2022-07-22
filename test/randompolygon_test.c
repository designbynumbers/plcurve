/* 

   randompolygon_test.c : Test code for the random polygon generation functions in plCurve. 

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

#ifdef HAVE_ASSERT_H
#include<assert.h>
#endif

#ifdef HAVE_GSL_GSL_RNG_H
#include<gsl/gsl_rng.h>
#endif

int PD_VERBOSE=0;

gsl_rng *r;  /* The global random number generator */

bool PAPERMODE;
FILE *outfile;

plCurve *plc_random_closed_polygon_selfcheck(gsl_rng *r,int nEdges);
plCurve *plc_random_open_plane_polygon_selfcheck(gsl_rng *r,int nEdges);

plCurve *plc_random_open_polygon_selfcheck(gsl_rng *r,int nEdges);
plCurve *plc_random_closed_plane_polygon_selfcheck(gsl_rng *r,int nEdges);

plCurve *plc_random_polygon_selfcheck(int nEdges,bool selfcheck); /* A private version which turns on debugging code. */

double *gaussian_array(gsl_rng *r,int n);

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
    testvar = gaussian_array(r,N);
    end = clock();
    cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;
    
    printf("done. (%3.2g sec). ",cpu_time_used);
    
    double mean=0,moments[5]={0,0,0,0};
    int i;
    
    for(i=0;i<2*N;i++) { mean += testvar[i]; }
    mean /= 2.0*N;
    
    for(i=0;i<2*N;i++) {
      
      moments[2] += pow(testvar[i] - mean,2.0);
      moments[3] += pow(testvar[i] - mean,3.0);
      moments[4] += pow(testvar[i] - mean,4.0);
      
    }
    
    moments[2] /= 2.0*N; moments[3] /= 2.0*N; moments[4] /= 2.0*N;
    
    double skewness,kurtosis;
    
    skewness = moments[3]/pow(moments[2],1.5);
    kurtosis = moments[4]/pow(moments[2],2.0);
    
    double JB; /* Jarque-Bera statistic, from Wikipedia */
    
    JB = (1.0*N/6.0)*(pow(skewness,2.0) + (1.0/4.0)*pow(kurtosis - 3.0,2.0));
    
    printf("Moments = %3g,%3g,%3g,%3g S = %+7.5f. K = %+7.5f. JB statistic = %7.5g.\n",
	   mean,moments[2],moments[3],moments[4],skewness,kurtosis,JB);
    
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

bool timing_and_length2_test(int nPolygons,int nSizes,int *Sizes,
			     plCurve *polygen(gsl_rng *r,int nEdges),
			     char *typestring,char *shortstring,bool length2) {

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

  if (PAPERMODE) {
    fprintf(outfile,"TimingData%s ",shortstring);
  }

  for(i=0;i<nSizes;i++) {

    localPASS = true;

    start = clock();

    for(j=0;j<nPolygons;j++) {

      test = polygen(r,Sizes[i]); 
      length = plc_arclength(test,NULL);

      if (length2) {
	if (fabs(length - 2.0)  > 1e-10) { localPASS = false; }
      } else {
	if (fabs(length - Sizes[i]) > 1e-10) { localPASS = false; }
      }

      plc_free(test);

    }

    end = clock();
    cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;
    cpu_time_used /= (double)(nPolygons);
    
    printf("%-10d   %-5d         %-5g     %-13.4g          ",
	   Sizes[i],nPolygons,length,cpu_time_used);

    if (PAPERMODE) {
      fprintf(outfile,", %d , %g ",Sizes[i],cpu_time_used);
    }
    
    if (localPASS) {
      printf(" pass \n");
    } else {
      printf(" FAIL \n"); PASS = false;
    }
    
  }

  printf("------------------------------------------------------------------\n\n");

  if (PAPERMODE) {
    fprintf(outfile,"\n");
  }
  
  return PASS;

}

bool chordlength_test(int nPolygons,int nEdges,int nSkips,int *Skips,plCurve *polygen(gsl_rng *r,int nEdges),double prediction(int n, int k),
		      char *prediction_string,char *polygon_type,char *shortstring)

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
  
  if (PAPERMODE) {
    fprintf(outfile,"MeanSquareChordlengthData%s, \" %d samples\", \" %d gons\" ",shortstring,nPolygons,nEdges);
  }

  printf("\nComputing data...");
  fflush(stdout);
 
  start = clock();

  double *allchords,*allSE;
  allchords = calloc(sizeof(double),nSkips);
  allSE = calloc(sizeof(double),nSkips);
  
  double **datasets = calloc(sizeof(double *),nSkips);
  for(i=0;i<nSkips;i++) {
    datasets[i] = calloc(sizeof(double),nPolygons);
    assert(datasets[i] != NULL);
  }
  
  for(i=0;i<nPolygons;i++) {
    
    /* Generate polygon and compute chord data */   
    /* Add to running total and to data set*/
    
    int j;
    for(j=0;j<nSkips;j++) {
      plCurve *L;
      L = polygen(r,nEdges);
      double dist = plc_M_sq_dist(L->cp[0].vt[0],L->cp[0].vt[Skips[j]]);
      allchords[j] += dist;
      datasets[j][i] = dist;
      plc_free(L);
    }
    
  }

  for(i=0;i<nSkips;i++) {

    allchords[i] /= (double)(nPolygons); 

    /* Now we're going to compute the standard error for each data set. */
    double sampvar = 0;
    int j;
    
    for(j=0;j<nPolygons;j++) {

      sampvar += (datasets[i][j] - allchords[i])*(datasets[i][j] - allchords[i]);
      
    }
    
    sampvar /= (double)(nPolygons-1); /* Bessel's correction */
    double SE = sqrt(sampvar/(double)(nPolygons)); /* Standard error. */
    
    allSE[i] = SE;
    free(datasets[i]);
    
  }
  free(datasets);
  
  end = clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("done (%g sec).\n\n",cpu_time_used); 

  /* Now display results. */

  printf("According to CDS, the predicted value of squared chord length \nfor %s at skip k is %s\n\n",polygon_type,prediction_string);
  printf("Skip      Computed Average    Predicted Average   Difference     95%% Confidence     Result\n"
	 "----------------------------------------------------------------------------------------------------------------------\n");

  for(i=0;i<nSkips;i++) {
    
    double predicted  = prediction(nEdges,Skips[i]);
    double err = allchords[i] - predicted;
    double SE = allSE[i];

    printf("%-5d     %-13.7g       %-13.7g       %+-13.7g  %-13.7g",Skips[i],allchords[i],predicted,err,1.96*SE);

    if (PAPERMODE) {
      fprintf(outfile,", %d , %g ",Skips[i],allchords[i]);
    }

    if (fabs(err) > 1.96*SE) { printf("      WARN (> 95%% confidence interval).\n"); PASS = true; }
    else if (fabs(err) > 3.20*SE) { printf("      FAIL (> 99.9%% confidence interval).\n"); PASS = false; }
    else { printf("      pass (in 95%% confidence interval).\n"); }
    
  }

  printf("----------------------------------------------------------------------------------------------------------------------\n");
  printf("\n\n");

  if (PAPERMODE) {
    fprintf(outfile,"\n");
  }
    
  free(allchords);
  free(allSE);
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
char plane_pol_predstring[256] = "((n-k)/n) (8k/((n-1)(n+2)))";

double plane_arm_prediction(int n,int k) {
  return (double)(8*k)/(double)(n*(n+1));
}
char plane_arm_predstring[256] = "8k/(n(n+1))";

double eq_pol_prediction(int n,int k) {
  return (double)((n-k)*k)/(double)((n-1));
}
char eq_pol_predstring[256] = "(k(n-k)/(n-1))";

double eq_arm_prediction(int n,int k) {
  return (double)(4*k)/(double)(n*n);
};
char eq_arm_predstring[256] = "4k/n^2";

bool gyradius_test(int nPolygons,int nSizes,int *Sizes,plCurve *polygen(gsl_rng *r,int nEdges),double prediction(int n),char *prediction_string,char *polygon_type,char *shortstring)

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
  
  if (PAPERMODE) {
    fprintf(outfile,"Gyradius%s, \" %d samples \"",shortstring,nPolygons);
  }   

  printf("\nComputing data...");
  fflush(stdout);
 
  start = clock();

  double *allgyradius;
  allgyradius = calloc(sizeof(double),nSizes);
  double *allError;
  allError = calloc(sizeof(double),nSizes);

  for(j=0;j<nSizes;j++) {

    double *data = calloc(nPolygons,sizeof(double));
  
    for(i=0;i<nPolygons;i++) {
    
      /* Generate polygon and compute gyradius */
      plCurve *L;
      L = polygen(r,Sizes[j]);
      double gyr = plc_gyradius(L);
      
      allgyradius[j] += gyr;
      data[i] = gyr;
  
      /* Free memory */
      plc_free(L);
      
    }

    /* Now we need to compute the 95% confidence interval
       for this sample to see if it's in spec. To do that,
       we compute sample variance. */

    allgyradius[j] /= (double)(nPolygons);

    double sampvar = 0;
    for(i=0;i<nPolygons;i++) {

      sampvar += (data[i] - allgyradius[j])*(data[i] - allgyradius[j]);

    }

    sampvar /= (double)(nPolygons-1); /* Bessel's correction */
    double SE = sqrt(sampvar/(double)(nPolygons)); /* Standard error. */

    /* The 95% confidence interval is given by the sample mean +- 1.96 * SE. */

    allError[j] = 1.96*SE;
    free(data);
  }
  end = clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("done (%g sec).\n\n",cpu_time_used); 

  /* Now display results. */

  printf("According to CDS, the expected value of gyradius for n edges is is %s\n\n",prediction_string);
  printf("n (num verts)    Mean Gyradius    Predicted Mean Gyradius   Difference      95%% Error       Result\n"
	 "---------------------------------------------------------------------------------------------------------------\n");

  for(i=0;i<nSizes;i++) {
    
    double predicted  = prediction(Sizes[i]);
    double err = allgyradius[i] - predicted;
    double okerr = allError[i];

    printf("%-5d            %-13.7g    %-13.7g             %-+13.6g   %-13.6g",Sizes[i],allgyradius[i],predicted,err,okerr);

    if (PAPERMODE) {
      fprintf(outfile,", %d , %g",Sizes[i],allgyradius[i]);
    }

    if (fabs(err) > okerr) { printf("   WARN (outside 95 %% confidence).\n"); PASS = true; }
    else if (fabs(err) > (3.290/1.96)*okerr) { printf("   FAIL (outside 99.9 %% confidence).\n"); PASS = false; }
    else { printf("   pass (within 95 %%).\n"); }
    
  }

  printf("\n\n");
  free(allgyradius);

  if (PAPERMODE) { fprintf(outfile,"\n"); }

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

double eq_pol_gyradius_prediction(int n) {
  return (double)(n+1)/(double)(12.0);
}
char eq_pol_gyradius_predstring[256] = "(n+1)/12";

double eq_arm_gyradius_prediction(int n) {
  return (double)(2.0*(n+2))/(double)(3.0*n*(n+1));
}
char eq_arm_gyradius_predstring[256] = "(2/3) (n+2)/n(n+1)";


bool totalcurv_test(int nPolygons,int nSizes,int *Sizes,plCurve *polygen(gsl_rng *r,int nEdges),double prediction(int n),char *prediction_string,char *polygon_type,char *shortstring)

/* Given a generation function and a prediction for total curvature in terms of n, do test. */

{
  int i,j;
  bool PASS = true;
  clock_t start,end;
  double cpu_time_used;

  printf("Total curvature test for %d %s polygons at %d numbers of verts.\n\tSizelist: ",
	 nPolygons,polygon_type,nSizes);
  for(i=0;i<nSizes;i++) {
    printf("%5d ",Sizes[i]);
    if (i % 10 == 0 && i > 0) { printf("\n"); }
  }
  
  if (PAPERMODE) {
    fprintf(outfile,"Totalcurvature%s, \" %d samples \"",shortstring,nPolygons);
  }   

  printf("\nComputing data...");
  fflush(stdout);
 
  start = clock();

  double *alltc;
  alltc = calloc(sizeof(double),nSizes);
  double *allSE;
  allSE = calloc(sizeof(double),nSizes);

  for(j=0;j<nSizes;j++) {

    double *data = calloc(nPolygons,sizeof(double));
  
    for(i=0;i<nPolygons;i++) {
    
      /* Generate polygon and compute tc */
      plCurve *L;
      L = polygen(r,Sizes[j]);
      double tc = plc_totalcurvature(L,NULL);
      
      alltc[j] += tc;
      data[i] = tc;
  
      /* Free memory */
      plc_free(L);
      
    }

    /* Now we need to compute the 95% confidence interval
       for this sample to see if it's in spec. To do that,
       we compute sample variance. */

    alltc[j] /= (double)(nPolygons);

    double sampvar = 0;
    for(i=0;i<nPolygons;i++) {

      sampvar += (data[i] - alltc[j])*(data[i] - alltc[j]);

    }

    sampvar /= (double)(nPolygons-1); /* Bessel's correction */
    double SE = sqrt(sampvar/(double)(nPolygons)); /* Standard error. */

    /* The 95% confidence interval is given by the 
       sample mean +- 1.96 * SE. */

    allSE[j] = SE;
    free(data);
  }
  end = clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("done (%g sec).\n\n",cpu_time_used); 

  /* Now display results. */

  printf("According to CDS, the expected value of tc for n edges is is %s\n\n",prediction_string);
  printf("n (num verts)    Mean TotalK   Predicted Mean TotalK     Difference      95%% Error       Result\n"
	 "---------------------------------------------------------------------------------------------------------------\n");

  for(i=0;i<nSizes;i++) {
    
    double predicted  = prediction(Sizes[i]);
    double err = alltc[i] - predicted;
    double okerr = 1.96*allSE[i];

    printf("%-5d            %-13.7g %-13.7g             %-+13.6g   %-13.6g",Sizes[i],alltc[i],predicted,err,okerr);

    if (PAPERMODE) {
      fprintf(outfile,", %d , %g",Sizes[i],alltc[i]);
    }

    if (fabs(err) > 1.96*okerr) { printf("   WARN (outside 95 %% confidence).\n"); PASS = true; }
    else if (fabs(err) > 3.290*okerr) { printf("   FAIL (outside 99.9 %% confidence).\n"); PASS = false; }
    else { printf("   pass (within 95 %%).\n"); }
    
  }

  printf("\n\n");
  free(alltc);
  free(allSE);

  if (PAPERMODE) { fprintf(outfile,"\n"); }

  return PASS;

}

double PI = 3.1415926535897932385;

/* These prediction functions are for totalcurvature_test */

double space_pol_totalcurv_prediction(int n) {
return (double)(n*PI)/(double)(2.0) + (PI/4.0)*((double)(2*n)/(double)(2*n-3));
}
char space_pol_totalcurv_predstring[256] = "n pi/2 + (pi/4)(2n/(2n-3))";

double space_arm_totalcurv_prediction(int n) {
return (double)(n-1)*(PI/2.0);
}
char space_arm_totalcurv_predstring[256] = "(n-1)(pi/2)";

double plane_arm_totalcurv_prediction(int n) {
return (double)(n-1)*(PI/2.0);
}
char plane_arm_totalcurv_predstring[256] =  "(n-1)(pi/2)";

double eq_pol_totalcurv_prediction(int n)
{
if (n == 7) {
return 12.3689737;
} else if (n==8) {
return 13.9143058;
} else if (n==15) {
return 24.8244132;
} else if (n==16) {
return 26.3895713;
} else if (n==31) {
return 49.9121020;
} else if (n==32) {
return 51.4816285;
} else if (n==63) {
return 100.1572789;
} else if (n==64) {
return 101.7277732;
} else {
return 0.0;
}
}

char eq_pol_totalcurv_predstring[256] = "(values computed from integral formula)";

double eq_arm_totalcurv_prediction(int n) {
return (double)(n-1)*(PI/2.0);
}
char eq_arm_totalcurv_predstring[256] = "(n-1)(pi/2)";


bool ftc_test(int nPolygons,int nSizes,int *Sizes,plCurve *polygen(gsl_rng *r,int nEdges),double prediction(int n),char *prediction_string,char *polygon_type,char *shortstring)
  
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

  if (PAPERMODE) {
    fprintf(outfile,"FailureToClose%s, \"%d samples\"",shortstring, nPolygons);
  } 
  
  printf("\nComputing data...");
  fflush(stdout);
 
  start = clock();

  double *all_lastedgelength;
  all_lastedgelength = calloc(nSizes,sizeof(double));

  for(j=0;j<nSizes;j++) {

    all_lastedgelength[j] = 0;
  
    for(i=0;i<nPolygons;i++) {
    
      /* Generate polygon and compute ftc */
      plCurve *L;
      L = polygen(r,Sizes[j]);
      all_lastedgelength[j] += pow(plc_distance(L->cp[0].vt[0],L->cp[0].vt[L->cp[0].nv-1]),2.0);
  
      /* Free memory */
      plc_free(L);
      
    }

    all_lastedgelength[j] /= (double)(nPolygons);

  }
  end = clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("done (%g sec).\n\n",cpu_time_used); 

  /* Now display results. */

  printf("According to CDS, the expected value of failure-to-close for n edge %s polygons is %s\n\n",polygon_type,prediction_string);
  printf("n (num verts)    Mean Sq FTC    Predicted Mean Sq FTC   Difference   Result\n"
	 "---------------------------------------------------------------------------\n");

  for(i=0;i<nSizes;i++) {
    
    double predicted  = prediction(Sizes[i]);
    double percenterr = 100*(fabs(all_lastedgelength[i] - predicted)/predicted);

    printf("%-5d            %-13.7g  %-13.7g           %3.3f %%",Sizes[i],all_lastedgelength[i],predicted,percenterr);

    if (PAPERMODE) {
      fprintf(outfile,", %d, %g ",Sizes[i],all_lastedgelength[i]);
    }

    if (percenterr > 1.0) { printf("     FAIL (> 1%%).\n"); PASS = false; }
    else { printf("     pass (< 1%%).\n"); }
    
  }

  printf("Used %g seconds of cpu time (%g per polygon).\n",cpu_time_used,cpu_time_used/(double)(nSizes*nPolygons));
  printf("\n\n");

  if (PAPERMODE) { fprintf(outfile,"\n"); }

  free(all_lastedgelength);

  return PASS;

}
  
/* These prediction functions are for FTC_test */

double space_arm_ftc_prediction(int n) {
  return (double)(6.0)/(double)(n+0.5);
}
char space_arm_ftc_predstring[256] = "6/(n+1/2)";

double plane_arm_ftc_prediction(int n) {
   return (double)(8.0)/(double)(n+1);
}
char plane_arm_ftc_predstring[256] = "8/(n+1)";
  
int main(int argc, char *argv[]) {

  bool PASS = {true};

  const gsl_rng_type * T;
     
  gsl_rng_env_setup();
  
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  int seed = time(0);
  gsl_rng_set(r,seed);

  if (argc > 1) { PAPERMODE = true; }

  printf("Random Polygon Generation Tests \n"
	 "------------------------------- \n"
	 "plCurve generates fixed total length 2.0 random polygons in four\n"
	 "classes by direct sampling from the symmetric measure of Cantarella,\n"
	 "Deguchi, Shonkwiler [CPAM, 2012].\n"
	 "\n"
	 "We also generate equilateral length 2.0 random polygons in space\n"
	 "(closed and open) using the \"fast ergodic algorithm\" of Varela,\n"
	 "Hinson, Arsuaga, and Diao [J. Phys. A, 2009].\n"
	 "\n"
	 "This test suite code is using the random number generator %s \n"
	 "with seed %d.\n"
	 "\n"
	 "\t Type                     Call                   \n"
	 "\t -------------------------------------------------------------------\n"
	 "\t closed space polygons    plc_random_closed_polygon(gsl_rng *r,int nEdges) \n"
	 "\t open space polygons      plc_random_open_polygon(gsl_rng *r,int nEdges) \n"
	 "\n"
	 "\t closed plane polygons    plc_random_closed_plane_polygon(gsl_rng *r,int nEdges) \n"
	 "\t open plane polygons      plc_random_open_plane_polygon(gsl_rng *r,int nEdges) \n"
	 "\n"
	 "\t closed eq space polygons plc_random_equilateral_polygon(gsl_rng *r,int nEdges)\n"
	 "\t open eq space polygons   plc_random_equilateral_open_polygon(gsl_rng *r,int nEdges)\n"
	 "\n"
	 "This program tests the polygons generated against theoretical calculations\n"
	 "about the symmetric measure. All polygons generated are guaranteed to have\n"
	 "length 2. (Rescaling, if desired, requires you to pick a probability \n"
	 "distribution function on polygon length and so is left to the user).\n\n",gsl_rng_name(r),seed);

  test_gaussianarray();
  test_hdot();

  /* Timing and Length tests. */

  int Sizes[100] = {20,200,2000,20000};
  int nSizes = 3;
  int nPolygons = 40;
  int nEQSizes = 4;
  int eqSizes[100] = {10,100,500,1000};
  int oddeqSizes[100] = {9,99,499,999};

  if (PAPERMODE) { 
    outfile = fopen("timing_and_length2.csv","w");
    timestamp(outfile);
  }

  if (!timing_and_length2_test(nPolygons,nSizes,Sizes,plc_random_closed_polygon_selfcheck,"closed space ","CS",true)
      || !timing_and_length2_test(nPolygons,nSizes,Sizes,plc_random_open_polygon_selfcheck,"open space ","OS",true)
      || !timing_and_length2_test(nPolygons,nSizes,Sizes,plc_random_closed_plane_polygon_selfcheck,"closed plane ","CP",true)
      || !timing_and_length2_test(nPolygons,nSizes,Sizes,plc_random_open_plane_polygon_selfcheck,"open plane ","OP",true) 
      || !timing_and_length2_test(nPolygons,nEQSizes,eqSizes,plc_random_equilateral_closed_polygon,"closed equilateral space ","CES",false)
      || !timing_and_length2_test(nPolygons,nEQSizes,oddeqSizes,plc_random_equilateral_closed_polygon,"closed equilateral space ","CES",false)
      || !timing_and_length2_test(nPolygons,nEQSizes,eqSizes,plc_random_equilateral_open_polygon,"open equilateral space ","OES",true)) {
    PASS = false;
  }

  if (PAPERMODE) { fclose(outfile); }

/* Total curvature tests. */

  int tcSizes[100] = {150,150,200,250,300,350,400,450,500};
  int ntcSizes = 1;
  nPolygons = 60000;

  int nEQtcSizes = 4;
  int EQtcSizes[100] = {8,16,32,64};
  int oddEQtcSizes[100] = {7,15,31,63};
  int nEQPolygons = 60000;

  if (PAPERMODE) { 
    outfile = fopen("totalcurvature.csv","w");
    timestamp(outfile);
  }

  if (!totalcurv_test(nPolygons,ntcSizes,tcSizes,plc_random_closed_polygon,space_pol_totalcurv_prediction,space_pol_totalcurv_predstring,"closed space","CS")
      || !totalcurv_test(nPolygons,ntcSizes,tcSizes,plc_random_open_polygon,space_arm_totalcurv_prediction,space_arm_totalcurv_predstring,"open space","OS")
      || !totalcurv_test(nPolygons,ntcSizes,tcSizes,plc_random_open_plane_polygon,plane_arm_totalcurv_prediction,plane_arm_totalcurv_predstring,"open plane","OP")) {
    
    PASS = false;

  }
  
  if (!totalcurv_test(nEQPolygons,nEQtcSizes,EQtcSizes,plc_random_equilateral_closed_polygon,eq_pol_totalcurv_prediction,eq_pol_totalcurv_predstring,"closed equilateral space","CES")
      || !totalcurv_test(nEQPolygons,nEQtcSizes,oddEQtcSizes,plc_random_equilateral_closed_polygon,eq_pol_totalcurv_prediction,eq_pol_totalcurv_predstring,"closed equilateral space","CES")
      || !totalcurv_test(nEQPolygons,nEQtcSizes,EQtcSizes,plc_random_equilateral_open_polygon,eq_arm_totalcurv_prediction,eq_arm_totalcurv_predstring,"open equilateral space","OES")) 
    {
      PASS = false;
    }
  
  if (PAPERMODE) { fclose(outfile); }


  /* Mean squared chordlength tests. */

  int Skips[100] = {10,20,30,40,50,60,70,80,90};
  int nSkips = 9;
  nPolygons = 60000;
  int nEdges = 250;

  int nEQEdges = 250;
  int EQSkips[100] = {10,20,30,40,50,60,70,80,90};
  int noddEQEdges = 249;
  int nEQSkips = 9;
  nEQPolygons = 5000;
  
  if (PAPERMODE) { 
    outfile = fopen("mean_squared_chordlength.csv","w");
    timestamp(outfile);
  }

  if (!chordlength_test(nPolygons,nEdges,nSkips,Skips,plc_random_closed_polygon,space_pol_prediction,space_pol_predstring,"closed space polygon","CS")
      || !chordlength_test(nPolygons,nEdges,nSkips,Skips,plc_random_open_polygon,space_arm_prediction,space_arm_predstring,"open space polygon","OS")
      || !chordlength_test(nPolygons,nEdges,nSkips,Skips,plc_random_closed_plane_polygon,plane_pol_prediction,plane_pol_predstring,"closed plane polygon","CP")
      || !chordlength_test(nPolygons,nEdges,nSkips,Skips,plc_random_open_plane_polygon,plane_arm_prediction,plane_arm_predstring,"open plane polygon","OP"))
    {
      PASS = false;
      
    }

  if (!chordlength_test(nEQPolygons,nEQEdges,nEQSkips,EQSkips,plc_random_equilateral_closed_polygon,eq_pol_prediction,eq_pol_predstring,"closed equilateral space polygon","CES")
      || !chordlength_test(nEQPolygons,noddEQEdges,nEQSkips,EQSkips,plc_random_equilateral_closed_polygon,eq_pol_prediction,eq_pol_predstring,"closed equilateral space polygon","CES")
      || !chordlength_test(nEQPolygons,nEQEdges,nEQSkips,EQSkips,plc_random_equilateral_open_polygon,eq_arm_prediction,eq_arm_predstring,"open equilateral space polygon","OES")) 
    {
      PASS = false;
    }

  if (PAPERMODE) { fclose(outfile); }

  /* Gyradius tests. */

  int gySizes[100] = {150,150,200,250,300,350,400,450,500};
  int ngySizes = 1;
  nPolygons = 40000;
  int nEQgySizes = 5;
  int EQgySizes[100] = {50,100,150,200,250};
  int oddEQgySizes[100] = {49,99,149,199,249};
  nEQPolygons = 2000;

  if (PAPERMODE) { 
    outfile = fopen("gyradius.csv","w");
    timestamp(outfile);
  }

  if (!gyradius_test(nPolygons,ngySizes,gySizes,plc_random_closed_polygon,space_pol_gyradius_prediction,space_pol_gyradius_predstring,"closed space","CS")
      || !gyradius_test(nPolygons,ngySizes,gySizes,plc_random_open_polygon,space_arm_gyradius_prediction,space_arm_gyradius_predstring,"open space","OS")
      || !gyradius_test(nPolygons,ngySizes,gySizes,plc_random_closed_plane_polygon,plane_pol_gyradius_prediction,plane_pol_gyradius_predstring,"closed plane","CP")
      || !gyradius_test(nPolygons,ngySizes,gySizes,plc_random_open_plane_polygon,plane_arm_gyradius_prediction,plane_arm_gyradius_predstring,"open plane","OP")) {
    
    PASS = false;

  }
  
  if (!gyradius_test(nEQPolygons,nEQgySizes,EQgySizes,plc_random_equilateral_closed_polygon,eq_pol_gyradius_prediction,eq_pol_gyradius_predstring,"closed equilateral space","CES")
      || !gyradius_test(nEQPolygons,nEQgySizes,oddEQgySizes,plc_random_equilateral_closed_polygon,eq_pol_gyradius_prediction,eq_pol_gyradius_predstring,"closed equilateral space","CES")
      || !gyradius_test(nEQPolygons,nEQgySizes,EQgySizes,plc_random_equilateral_open_polygon,eq_arm_gyradius_prediction,eq_arm_gyradius_predstring,"open equilateral space","OES")) 
    {
      PASS = false;
    }
  
  if (PAPERMODE) { fclose(outfile); }
  
  /* Failure to close tests. */
  
  int ftcSizes[100] = {7,150,200,250,300,350,400,450,500};
  int nftcSizes = 1;
  nPolygons = 80000;
  
  if (PAPERMODE) { 
    outfile = fopen("failure_to_close.csv","w");
    timestamp(outfile);
  }  

  if (!ftc_test(nPolygons,nftcSizes,ftcSizes,plc_random_open_polygon,
		space_arm_ftc_prediction,space_arm_ftc_predstring,"open space","OS")
      || !ftc_test(nPolygons,nftcSizes,ftcSizes,plc_random_open_plane_polygon,
		   plane_arm_ftc_prediction,plane_arm_ftc_predstring,"open plane","OP")) {
    PASS = false;
  };
 
  if (PAPERMODE) { fclose(outfile); }

  gsl_rng_free(r);

  printf("==========================================\n");
 
  if (PASS) { 

    printf("Random Polygon Test Suite: PASS\n");
    exit(0); 

  } else { 

    printf("Random Polygon Test Suite: FAIL\n");
    exit(1); 

  }

}
