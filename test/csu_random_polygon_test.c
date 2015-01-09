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
    
    printf("testing the case     d[i] = (1.0,1.2,1.8,0.9,1.0) \n"
	   "                 theta[i] =     (0.6,3.1,0.2)     \n"
	   "\tgenerated polygon....");

  }

  plCurve *L = fantriangulation_action_angle(6,theta,d);

  printf("pass (didn't crash)\n");

  printf("\tchecking vt[0] at origin...");
  
  if (fabs(plc_distance(L->cp[0].vt[0],plc_build_vect(0,0,0))) > 1e-8) {

    printf("fail (|vt[0]| = %g != 0)\n",
	   plc_norm(L->cp[0].vt[0]));
    return false;

  }

  printf("pass\n");

  printf("\tchecking edgelengths...");

  int i;
  for(i=1;i<7;i++) {

    if (fabs(plc_distance(L->cp[0].vt[i],L->cp[0].vt[i-1]) - 1.0) > 1e-8) {

      printf("fail (|vt[%d] -> vt[%d]| = %g != 1.0)\n",i-1,i,
	     plc_distance(L->cp[0].vt[i],L->cp[0].vt[i-1]));
      return false;

    }
  }

  printf("pass\n");


  printf("\tchecking diagonals...");
  
  for(i=1;i<6;i++) {
    
    if (fabs(plc_distance(L->cp[0].vt[i],L->cp[0].vt[0]) - d[i-1]) > 1e-8) {

      printf("fail (|vt[%d] - vt[0]| = %g != d[%d] = %g)\n",
	     i,plc_distance(L->cp[0].vt[i],L->cp[0].vt[0]),
	     i-1,d[i-1]);
      return false;

    }

  }

  printf("\tchecking dihedrals...");

  plc_vector normals[4];
  bool ok;

  for(i=2;i<6;i++) {

    normals[i-2] = plc_normalize(plc_cross_prod(L->cp[0].vt[i-1],L->cp[0].vt[i]),&ok);

    if (!ok) {

      printf("fail (couldn't generate a normal vector for vt[0] - vt[%d] - vt[%d] triangle)\n",
	     i-1,i);
      return false;

    }

  }

  double PI = 3.1415926535897932385;
  
  for(i=1;i<4;i++) {

    if (fabs(plc_angle(normals[i-1],normals[i])-theta[i]) > 1e-8 &&
	fabs(plc_angle(normals[i-1],normals[i])-theta[i]-PI) > 1e-8 &&
	fabs(plc_angle(normals[i-1],normals[i])-theta[i]+PI) > 1e-8) {

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
