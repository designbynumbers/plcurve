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
   int i,j;
   
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

int main() {

  bool PASS = {true};
  srand48(time(0));

  printf("Random Polygon Generation Tests \n"
	 "------------------------------- \n"
	 "plCurve generates random closed polygons using the direct sampling from the symmetric measure algorithm\n"
	 "of Cantarella, Deguchi, and Shonkwiler. The call plc_random_closed_polygon(n) generates a polygon directly \n"
	 "sampled from this distribution. The polygon is guaranteed to be closed and have length 2.\n\n");

  plCurve *test;
  int i,j;
  int verts = 200;

  clock_t start,end;
  double cpu_time_used;

  double length;

  test_gaussianarray();
  test_hdot();

  printf("Verts        Length    Time to generate  Result\n");
  printf("-------------------------------------------------\n");

  for(i=0;i<4;i++,verts*=10) {

    for(j=0;j<5;j++) {

      start = clock();
      test = plc_random_closed_polygon_selfcheck(verts,true);  /* This is a private interface which turns out on some internal debugging code. */
      end = clock();
      cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

      length = plc_arclength(test,NULL);
      printf("%-10d   %-5g     %-13.4g    ",verts,length,cpu_time_used);
      if (fabs(length - 2.0) <1e-10) { 
	printf(" pass \n");
      } else {
	printf(" FAIL \n"); PASS = false;
      }

      plc_free(test);

    }

    printf("\n");

  }

  /* Distribution of last edgelength experiment */

  int N = 1000000;
  int NEDGE = 5;
  int J = 1;

  printf("Generating last edge data for %d %d-edge polygons to test closure...",N,NEDGE);

  double lastedge = 0;

  start = clock();

  for(i=0;i<N;i++) {

    test = plc_random_closed_polygon(NEDGE);
    lastedge += plc_distance(test->cp[0].vt[0],test->cp[0].vt[NEDGE-1]);

  }

  lastedge /= N;

  end = clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  if (fabs(lastedge - 2.0/NEDGE) > 1e-2) { 
    printf("FAIL (average length %g, should be 2/%d = %g).\n",lastedge,NEDGE,2.0/NEDGE);
    exit(1);
  } else {
    printf("done (%g sec)\n",cpu_time_used);
    printf("pass: Average length %g, close to predicted value of 2/%d = %g.\n\n",lastedge,NEDGE,2.0/NEDGE);
  }

  /* Random Polygon Chordlength Experiment */

  /* Generate chordlength data for N samples of NEDGE polygons at skips J, 2J, 3J, .... */

  N = 50000;

  printf("Generating chordlength data for %d %d-edge polygons at %d skips (%d, %d, %d ...) ... ",
	 N,NEDGE,(int)(NEDGE/J),J,2*J,3*J);
 
  start = clock();
  double *chords;
  chords = mean_squared_chordlengths(N,NEDGE,J);
  end = clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("done (%g sec).\n\n",cpu_time_used); 

  printf("According to CDS, the expected value of squared chord length at skip k for a closed n-gon is\n\n"
	 "  6 k (n-k) \n"
	 "------------- \n"
	 "(n-1) n (n+1) \n\n");

  printf("Skip      Computed Average    Predicted Average   Difference \n"
	 "-----------------------------------------------------------\n");

  for(i=J;i<NEDGE;i+=J) {
    
    double predicted = (double)(6*i*(NEDGE-i))/(double)(((NEDGE-1)*(NEDGE)*(NEDGE+1)));
    printf("%-5d     %-13.4g       %-13.4g       %3.2g%%\n",i,chords[i],predicted,100*(fabs(chords[i] - predicted)/predicted));
    
  }

  printf("\n\n"); 
 
  if (PASS) { exit(0); } else { exit(1); }

}
