/* 

   random_pe_polygon.c : Test code for the pseudoequilateral random polygon generation functions in plCurve. 

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

bool timing_and_length_test(int nPolygons,int nEdges,plCurve *polygen(int nEdges,double LOWER,double UPPER,int *nattempts),int nIntervals,double LOWER[],double UPPER[],char *typestring,char *shortstring) {

  clock_t start,end;
  double cpu_time_used;
  int i,j;
  plCurve *test;
  bool PASS = true;
  double length = 0;
  bool localPASS;

  printf("Timing and length test for pseudoequilateral %s %d-gons.\n",typestring,nEdges);

  printf("Edgelengths        Accepted/Generated    Ratio    Mean Time to Accepted  Second Moment of Edgelength  Quality Check\n");
  printf("-------------------------------------------------------------------------------------------------------------------\n");

  if (PAPERMODE) {
    fprintf(outfile," \"PE_timing_data%s\", \"%d edge polygons\" ",shortstring,nEdges);
  }

  for(i=0;i<nIntervals;i++) {

    int attempts = 0;
    int success = 0;
    double second_moment = 0;
    
    localPASS = true;

    start = clock();

    for(j=0;j<nPolygons;j++) {

      int n_attempts;

      test = polygen(nEdges,LOWER[i],UPPER[i],&n_attempts);
      attempts += n_attempts;
      
      if (test != NULL) {

	success++;
	length = plc_arclength(test,NULL);      
	if (fabs(length - 2.0)  > 1e-10) { localPASS = false; }
	
	double longedge,shortedge,meanedge,moment2edge;
	plc_edgelength_stats(test,&longedge,&shortedge,&meanedge,&moment2edge);
	if (fabs(meanedge - 2.0/(double)(nEdges)) > 1e-10) { localPASS = false; }
	if (longedge > UPPER[i] || shortedge < LOWER[i]) { localPASS = false; }

	second_moment += moment2edge;

	plc_free(test);

      } 

    }

    end = clock();
    cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;
    cpu_time_used /= (double)(success);
    second_moment /= (double)(success);

    printf("[%-5g,%-5g]      %5d/%-5d           %-7.5g     %-13.4g  (sec)   %7.5g ",LOWER[i],UPPER[i],success,attempts,(double)(success)/(double)(attempts),cpu_time_used,second_moment);

    if (PAPERMODE) {
      fprintf(outfile,", %g , %g , %g , %g, %g, ",LOWER[i],UPPER[i],(double)(nPolygons)/(double)(attempts),cpu_time_used,second_moment);
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

bool chordlength_test(int nPolygons,int nEdges,int nSkips,int *Skips,plCurve *polygen(int nEdges,double LOWER,double UPPER),
		      double LOWER,double UPPER,double prediction(int n, int k),
		      char *prediction_string,char *polygon_type,char *shortstring)

/* Given a generation function and a prediction for average chordlength in terms of number of edges and skip, do test. */

{
  int i;
  bool PASS = true;
  clock_t start,end;
  double cpu_time_used;

  printf("Mean Squared Chordlength test for %d %d-edge %s with edgelengths in [%g,%g] at %d skips.\n\t Skip list: ",
	 nPolygons,nEdges,polygon_type,LOWER,UPPER,nSkips);
  for(i=0;i<nSkips;i++) {
    printf("%5d ",Skips[i]);
    if (i % 10 == 0 && i > 0) { printf("\n"); }
  }
  
  if (PAPERMODE) {
    fprintf(outfile,"MeanSquareChordlengthData%s, \" %d samples\", \" %d gons\", \" edgelengths in [%g, %g] \" ",shortstring,nPolygons,nEdges,LOWER,UPPER);
  }

  printf("\nComputing data...");
  fflush(stdout);
 
  start = clock();

  double *allchords,*thispolychords;
  allchords = calloc(sizeof(double),nSkips);
  
  for(i=0;i<nPolygons;i++) {
    
    /* Generate polygon and compute chord data */
    plCurve *L;
    L = polygen(nEdges,LOWER,UPPER);
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

  printf("According to CDS, the predicted value of squared chord length \nfor equilateral %s at skip k is %s\n\n",polygon_type,prediction_string);
  printf("Skip      Computed Average    Predicted Average   Difference  \n"
	 "---------------------------------------------------------------------------\n");

  for(i=0;i<nSkips;i++) {
    
    double predicted  = prediction(nEdges,Skips[i]);
    double percenterr = 100*(fabs(allchords[i] - predicted)/predicted);

    printf("%-5d     %-13.7g       %-13.7g       %3.3f %% \n",Skips[i],allchords[i],predicted,percenterr);

    if (PAPERMODE) {
      fprintf(outfile,", %d , %g ",Skips[i],allchords[i]);
    }
    
  }

  printf("---------------------------------------------------------------------------\n");
  printf("\n\n");

  if (PAPERMODE) {
    fprintf(outfile,"\n");
  }
    
  free(allchords);
  return PASS;

}

/* These prediction functions are for chordlength_test */

double equilateral_pol_prediction(int n, int k) {
  double res,den; 
  res = ((double)(k*(n-k)*4.0));
  den = (double)(n)*(double)(n)*(double)(n-1);
  res /= den;
  return res;
}

char equilateral_pol_predstring[256] = "(k(n-k))/n 4/(n(n-1))";

double equilateral_arm_prediction(int n,int k) {
  return (double)(k*4.0)/(double)(n*n);
}

char equilateral_arm_predstring[256] = "4k/n^2";


bool gyradius_test(int nPolygons,int nSizes,int *Sizes,plCurve *polygen(int nEdges, double LOWER, double UPPER),double LOWER, double UPPERMULT,double prediction(int n),char *prediction_string,char *polygon_type,char *shortstring)

  /* Given a generation function and a prediction for gyradius in terms of n, do test. */

{
  int i,j;
  bool PASS = true;
  clock_t start,end;
  double cpu_time_used;

  printf("(Squared) Radius of Gyration test for %d %s polygons at %d numbers of verts with edgelengths in [%g,%g * 2.0/nEdges].\n\tSizelist: ",
	 nPolygons,polygon_type,nSizes,LOWER,UPPERMULT);
  for(i=0;i<nSizes;i++) {
    printf("%5d ",Sizes[i]);
    if (i % 10 == 0 && i > 0) { printf("\n"); }
  }
  
  if (PAPERMODE) {
    fprintf(outfile,"Gyradius%s, \" %d samples \", \" edgelengths in [%g,%g * 2.0/nEdges] \" ",shortstring,nPolygons,LOWER,UPPERMULT);
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
      L = polygen(Sizes[j],LOWER,UPPERMULT*2.0/(double)(Sizes[j]));
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

  printf("According to CDS, the expected value of gyradius for equilateral n edge polygons is %s\n\n",prediction_string);
  printf("n (num verts)    Mean Gyradius    Predicted Mean Gyradius   Difference  \n"
	 "------------------------------------------------------------------------------\n");

  for(i=0;i<nSizes;i++) {
    
    double predicted  = prediction(Sizes[i]);
    double percenterr = 100*(fabs(allgyradius[i] - predicted)/predicted);

    printf("%-5d            %-13.7g    %-13.7g             %7.3f %% \n",Sizes[i],allgyradius[i],predicted,percenterr);

    if (PAPERMODE) {
      fprintf(outfile,", %d , %g",Sizes[i],allgyradius[i]);
    }
    
  }

  printf("\n\n");
  free(allgyradius);

  if (PAPERMODE) { fprintf(outfile,"\n"); }

  return PASS;

}

/* These prediction functions are for gyradius_test */

double equilateral_pol_gyradius_prediction(int n) {
  return (double)(n+1)/(double)(3*n*n);
}

char equilateral_pol_gyradius_predstring[256] = "(1/3) (n+1)/n^2";

double equilateral_arm_gyradius_prediction(int n) {
  return (double)(2.0*(n+2))/(double)(3.0*n*(n+1));
}
char equilateral_arm_gyradius_predstring[256] = "(2/3) (n+2/(n(n+1))";
  
int main(int argc, char *argv[]) {

  bool PASS = {true};
  srand48(time(0));

  if (argc > 1) { PAPERMODE = true; }

  printf("Random Pseudoequilateral Polygon Generation Tests \n"
	 "------------------------------- \n"
	 "plCurve generates pseudoequilateral random polygons in four classes by rejection sampling from the symmetric measure of Cantarella, Deguchi, Shonkwiler."
	 "these polygons are sampled from the symmetric measure, but polygons with max edgelength > UPPER or min edgelength < LOWER are rejected. The mean edgelength"
	 "for an n-edge polygon is 2.0/n, so it is required that this be included in the interval [LOWER,UPPER].\n"
	 "\n"
	 "\t Type                    Call                   \n"
	 "\t -------------------------------------------------------------------\n"
	 "\t closed space polygons   plc_random_closed_polygon_PE(int nEdges,double LOWER,double UPPER) \n"
	 "\t open space polygons     plc_random_open_polygon_PE(int nEdges,double LOWER,double UPPER) \n"
	 "\t closed plane polygons   plc_random_closed_plane_polygon_PE(int nEdges,double LOWER,double UPPER) \n"
	 "\t open plane polygons     plc_random_open_plane_polygon_PE(int nEdges,double LOWER,double UPPER) \n"
	 "\n"
	 "This program tests the polygons generated against theoretical calculations about the symmetric measure. All polygons generated"
	 "are guaranteed to have length 2. (Rescaling, if desired, requires you to pick a probability distribution function on polygon length and so is"
	 "left to the user).\n\n");

  /* Timing and Length = 2.0 tests. */

  int nEdges = 200000;
  int nIntervals = 14;

  double LOWERS[100] = {0,0,0,0,0,0,0,0,0,0};
  double UPPERS[100] = {20*2.0/(double)(nEdges),10*2.0/(double)(nEdges),9*2.0/(double)(nEdges),8*2.0/(double)(nEdges),7*2.0/(double)(nEdges),
			6.5*2.0/(double)(nEdges),6*2.0/(double)(nEdges),5.5*2.0/(double)(nEdges),5*2.0/(double)(nEdges),4.5*2.0/(double)(nEdges),
			4.25*2.0/(double)(nEdges),4.125*2.0/(double)(nEdges),4.0*2.0/(double)(nEdges),3.95*2.0/(double)(nEdges)};

  int nPolygons = 200;
 
  if (PAPERMODE) { 
    outfile = fopen("pe_timing.csv","w");
    timestamp(outfile);
  }

  if (!timing_and_length_test(nPolygons,nEdges,plc_random_closed_polygon_PE_selfcheck,nIntervals,LOWERS,UPPERS,"closed space polygon","CS")
      || !timing_and_length_test(nPolygons,nEdges,plc_random_open_polygon_PE_selfcheck,nIntervals,LOWERS,UPPERS,"open space polygon","OS")) {
    PASS = false;
  };

  double NEWUPPERS[100] = {11*2.0/(double)(nEdges),10.5*2.0/(double)(nEdges),10*2.0/(double)(nEdges),9.5*2.0/(double)(nEdges),9*2.0/(double)(nEdges),
			   8.5*2.0/(double)(nEdges),8*2.0/(double)(nEdges),7.5*2.0/(double)(nEdges),7*2.0/(double)(nEdges),6.5*2.0/(double)(nEdges),
			   6.25*2.0/(double)(nEdges),6.125*2.0/(double)(nEdges),6.0*2.0/(double)(nEdges),5.95*2.0/(double)(nEdges)};
  
  if(!timing_and_length_test(nPolygons,nEdges,plc_random_closed_plane_polygon_PE_selfcheck,nIntervals,LOWERS,NEWUPPERS,"closed plane polygon","CP")
      || !timing_and_length_test(nPolygons,nEdges,plc_random_open_plane_polygon_PE_selfcheck,nIntervals,LOWERS,NEWUPPERS,"open plane polygon","OP")) {
    PASS = false;
  }

  if (PAPERMODE) { fclose(outfile); }

  exit(0);

  /* Mean squared chordlength tests. */

  int Skips[100] = {100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900};
  int nSkips = 19;
  nPolygons = 2000;
  nEdges = 2000;
  double EdgeMult = 5.0;
  
  if (PAPERMODE) { 
    outfile = fopen("pe_mean_squared_chordlength.csv","w");
    timestamp(outfile);
  }

  if (!chordlength_test(nPolygons,nEdges,nSkips,Skips,plc_random_closed_polygon_PE,0.0,EdgeMult*(2.0/nEdges),equilateral_pol_prediction,equilateral_pol_predstring,"closed space polygon","CS")
      || !chordlength_test(nPolygons,nEdges,nSkips,Skips,plc_random_open_polygon_PE,0.0,EdgeMult*(2.0/nEdges),equilateral_arm_prediction,equilateral_arm_predstring,"open space polygon","OS") ) {
    PASS = false;
  }

  EdgeMult = 7.0;

  if (!chordlength_test(nPolygons,nEdges,nSkips,Skips,plc_random_closed_plane_polygon_PE,0.0,EdgeMult*(2.0/nEdges),equilateral_pol_prediction,equilateral_pol_predstring,"closed plane polygon","CP")
      || !chordlength_test(nPolygons,nEdges,nSkips,Skips,plc_random_open_plane_polygon_PE,0.0,EdgeMult*(2.0/nEdges),equilateral_arm_prediction,equilateral_arm_predstring,"open plane polygon","OP")) {
    PASS = false;
  }

  if (PAPERMODE) { fclose(outfile); }

  /* Gyradius tests. */

  int gySizes[100] = {2000};
  int ngySizes = 1;
  nPolygons = 2000;

  if (PAPERMODE) { 
    outfile = fopen("pe_gyradius.csv","w");
    timestamp(outfile);
  }

  EdgeMult = 5.0;

  if (!gyradius_test(nPolygons,ngySizes,gySizes,plc_random_closed_polygon_PE,0.0,EdgeMult,
		     equilateral_pol_gyradius_prediction,equilateral_pol_gyradius_predstring,"closed space polygon","CS")
      || !gyradius_test(nPolygons,ngySizes,gySizes,plc_random_open_polygon_PE,0.0,EdgeMult,
			equilateral_arm_gyradius_prediction,equilateral_arm_gyradius_predstring,"open space polygon","OS")) {

    PASS = false;

  }

  EdgeMult = 7.0;

  if (!gyradius_test(nPolygons,ngySizes,gySizes,plc_random_closed_plane_polygon_PE,0.0,EdgeMult,
		     equilateral_pol_gyradius_prediction,equilateral_pol_gyradius_predstring,"closed plane polygon","CP")
      || !gyradius_test(nPolygons,ngySizes,gySizes,plc_random_open_plane_polygon_PE,0.0,EdgeMult,
			equilateral_arm_gyradius_prediction,equilateral_arm_gyradius_predstring,"open plane polygon","OP")) {
    PASS = false;
  };
  
  if (PAPERMODE) { fclose(outfile); }
 
  if (PASS) { exit(0); } else { exit(1); }

}
