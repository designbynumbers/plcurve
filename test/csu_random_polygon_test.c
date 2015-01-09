/* 

   csu_random_polygon_test.c : Test code for the CSU algorithm random polygon 
                               generation functions in plCurve. 

*/

#include<plCurve.h>
#include<tsmcmc.h>
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


plCurve *fantriangulation_action_angle(int n,double *theta, double *d);
plCurve *plc_random_equilateral_closed_polygon(gsl_rng *r,int nEdges);

bool test_fantriangulation_action_angle(int n,double *theta,double *d,bool verbose)
{

  if (verbose) {
  
    if (n <= 10) {

      int j;
      printf("testing the case   d[i] = (");
      for(j=0;j<n-1;j++) {

	printf("%3.2f%s",d[i],(j < n-2) ? ",":")\n");

      }

      printf("               theta[i] = (");
      for(j=0;j<n-3;j++) {

	printf("%3.2f%s",theta[i],(j < n-4) ? ",":")\n");

      }

    } else {

      printf("testing %d edge case d[i] = (");

      for(j=0;j<4;j++) {
	printf("%3.2f,",d[i]);
      }
      printf("...");
      for(j=n-4;j<n-1;j++) {
	printf("%3.2f%s",d[i],(j<n-2) ? ",":")\n");
      }

      printf("               theta[i] = (");
      for(j=0;j<3;j++) {
	printf("%3.2f,",theta[i]);
      }
      printf("...");
      for(j=n-6;j<n-3;j++) {
	printf("%3.2f%s",theta[i],(j<n-4) ? ",":")\n");
      }

    }

  }

  plCurve *L = fantriangulation_action_angle(n,theta,d);

  if (verbose) {

    printf("pass (didn't crash)\n");

    printf("\tchecking vt[0] at origin...");

  }
  
  if (fabs(plc_distance(L->cp[0].vt[0],plc_build_vect(0,0,0))) > 1e-8) {

    if (!verbose) {

      printf("FAIL in test of %d edge case\n",n);
      printf("\tchecking vt[0] at origin...");

    }
    
    printf("fail (|vt[0]| = %g != 0)\n",
	   plc_norm(L->cp[0].vt[0]));
    return false;

  }

  if (verbose) {
    
    printf("pass\n");
    printf("\tchecking edgelengths...");

  }

  int i;
  for(i=1;i<=n;i++) {

    if (fabs(plc_distance(L->cp[0].vt[i],L->cp[0].vt[i-1]) - 1.0) > 1e-8) {

      if (!verbose) {

	printf("FAIL in test of %d edge case\n",n);
	printf("\tchecking edgelengths...");

      }
      
      printf("fail (|vt[%d] -> vt[%d]| = %g != 1.0)\n",i-1,i,
	     plc_distance(L->cp[0].vt[i],L->cp[0].vt[i-1]));
      return false;

    }
  }


  if (verbose) {
    
    printf("pass\n");
    printf("\tchecking diagonals...");

  }
  
  for(i=1;i<n;i++) {
    
    if (fabs(plc_distance(L->cp[0].vt[i],L->cp[0].vt[0]) - d[i-1]) > 1e-8) {

      if (!verbose) {

	printf("FAIL in test of %d edge case\n",n);
	printf("\tchecking diagonals...");

      }
      
      printf("fail (|vt[%d] - vt[0]| = %g != d[%d] = %g)\n",
	     i,plc_distance(L->cp[0].vt[i],L->cp[0].vt[0]),
	     i-1,d[i-1]);
      return false;

    }

  }

  if (verbose) {

    printf("pass\n");
    printf("\tchecking dihedrals...");

  }

  plc_vector *normals = calloc((n-2)*sizeof(plc_vector));
  assert(normals != NULL);
  bool ok;

  for(i=2;i<n;i++) {

    normals[i-2] = plc_normalize(plc_cross_prod(L->cp[0].vt[i-1],L->cp[0].vt[i]),&ok);

    if (!ok) {

      if (!verbose) {

	printf("FAIL in test of %d edge case\n",n);
	printf("\tchecking dihedrals...");

      }
      
      printf("fail (couldn't generate a normal vector for vt[0] - vt[%d] - vt[%d] triangle)\n",
	     i-1,i);
      return false;

    }

  }

  double PI = 3.1415926535897932385;
  
  for(i=1;i<n-2;i++) {

    if (fabs(plc_angle(normals[i-1],normals[i])-theta[i]) > 1e-8 &&
	fabs(plc_angle(normals[i-1],normals[i])-theta[i]-PI) > 1e-8 &&
	fabs(plc_angle(normals[i-1],normals[i])-theta[i]+PI) > 1e-8) {

      if (!verbose) {

	printf("FAIL in test of %d edge case\n",n);
	printf("\tchecking dihedrals...");
	
      }
      
      printf("fail\n"
	     "angle between n[%d] and n[%d] is %g\n"
	     "theta[%d] = %g, theta[%d] + pi = %g, theta[%d] - pi = %g)\n",
	     i-1,i,plc_angle(normals[i-1],normals[i]),
	     i,theta[i],
	     i,theta[i] + PI,
	     i,theta[i] - PI);
      return false;

    }

  }

  if (verbose) {

    printf("pass\n");

  }

  /* Now we're actually going to build the polygon and compare vertex by vertex. */

  
  if (verbose) {

    printf("\tbuilding polygon with old code...");

  }

  /********* ADD NEW STUFF HERE! ************/
  
}

bool assembly_test()
{
  printf("----------------------------------------------\n"
	 "testing polygon assembly code                 \n"
	 "----------------------------------------------\n");

  printf("testing the case     d[i] = (1.0,1.2,1.8,0.9,1.0) \n"
	 "                 theta[i] =     (0.6,3.1,0.2)     \n"
	 "\tgenerated polygon....");

  double d[5] = {1.0,1.2,1.8,0.9,1.0};
  double theta[3] = {0.6,3.1,0.2};
  double e[6] = {1.0,1.0,1.0,1.0,1.0,1.0};

 
  

      
  

  

}


int main() {

  printf("test_csu_randompolygon generation (%s)\n",PACKAGE_STRING);
  printf("--------------------------------------------------------\n"
	 "Unit tests for direct sampler of equilateral random polygons. \n"
	 "========================================================\n");

  if (!assembly_test() || !randompoly_quality_test()) {

    printf("=======================================================\n");
    printf("test_csu_randompolygon generation:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================================\n");
    printf("test_csu_randompolygon generation:  PASS.\n");
    exit(0);

  }

  return 0;

}
