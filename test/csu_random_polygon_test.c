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

#ifdef HAVE_ASSERT_H
#include<assert.h>
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

/* We expect to be provided with ALL the diagonal lengths (starting with 
   |v0-v1| = d[0] = 1.0 and ending with |v0-v{n-1}| = d[n-1] = 1.0), and 
   the n-3 dihedral angles. */
{
  if (verbose) {
  
    if (n <= 10) {

      int j;
      printf("\ntesting the case   d[i] = (");
      for(j=0;j<n-1;j++) {

	printf("%3.2f%s",d[j],(j < n-2) ? ",":")\n");

      }

      printf("               theta[i] = (");
      for(j=0;j<n-3;j++) {

	printf("%3.2f%s",theta[j],(j < n-4) ? ",":")\n");

      }

    } else {
      int j;

      printf("\ntesting %d edge case d[i] = (",n);

      for(j=0;j<4;j++) {
	printf("%3.2f,",d[j]);
      }
      printf("...");
      for(j=n-4;j<n-1;j++) {
	printf("%3.2f%s",d[j],(j<n-2) ? ",":")\n");
      }

      printf("               theta[i] = (");
      for(j=0;j<3;j++) {
	printf("%3.2f,",theta[j]);
      }
      printf("...");
      for(j=n-6;j<n-3;j++) {
	printf("%3.2f%s",theta[j],(j<n-4) ? ",":")\n");
      }

    }

    printf("\tgenerating polygon...");

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

  plc_vector *normals = calloc((n-2),sizeof(plc_vector));
  assert(normals != NULL);
  bool ok;

  for(i=2;i<n;i++) {

    normals[i-2] = plc_normalize_vect(
				      plc_cross_prod(L->cp[0].vt[i-1],
						     L->cp[0].vt[i]),&ok);

    if (!ok) {

      if (!verbose) {

	printf("FAIL in test of %d edge case\n",n);
	printf("\tchecking dihedrals...");

      }
      
      printf("fail (couldn't generate a normal vector "
	     "for vt[0] - vt[%d] - vt[%d] triangle)\n",
	     i-1,i);
      return false;

    }

  }

  double PI = 3.1415926535897932385;
  
  for(i=1;i<n-2;i++) {

    bool ok1,ok2,ok3;
    
    if (fabs(plc_angle(normals[i-1],normals[i],&ok1)-theta[i-1]) > 1e-8 &&
	fabs(plc_angle(normals[i-1],normals[i],&ok2)-(2*PI - theta[i-1])) > 1e-8 &&
	fabs(plc_angle(normals[i-1],normals[i],&ok3)-(2*PI + theta[i-1])) > 1e-8) {

      if (!verbose) {

	printf("FAIL in test of %d edge case\n",n);
	printf("\tchecking dihedrals...");
	
      }
      
      printf("fail\n"
	     "angle between n[%d] and n[%d] is %g\n"
	     "theta[%d] = %g, theta[%d] + pi = %g, theta[%d] - pi = %g)\n",
	     i-1,i,plc_angle(normals[i-1],normals[i],&ok),
	     i,theta[i-1],
	     i,theta[i-1] + PI,
	     i,theta[i-1] - PI);
      return false;

    }

  }

  if (verbose) {

    printf("pass\n");

  }

  /* Now we're actually going to build the polygon and compare vertex by vertex. */

  
  /* if (verbose) { */

  /*   printf("\tbuilding polygon with old code..."); */

  /* } */

  /* tsmcmc_triangulation_t T; */
  /* T = tsmcmc_fan_triangulation(n); */
  /* double *edge_lengths = calloc(n,sizeof(double)); */
  /* for(i=0;i<n;i++) { edge_lengths[i] = 1.0; } */
  
  /* plCurve *Lcmp = tsmcmc_embed_polygon(T,edge_lengths,&(d[1]),theta); */

  /* if (verbose) { */

  /*   printf("pass (didn't crash)\n"); */

  /* } */

  /* if (verbose) { */

  /*   printf("\tcomparing vertex-by-vertex..."); */

  /* } */

  /* for(i=0;i<L->cp[0].nv;i++) { */

  /*   if (plc_distance(L->cp[0].vt[i],Lcmp->cp[0].vt[i]) > 1e-6) { */

  /*     if (verbose) { */

  /* 	printf("FAIL in test of %d edge case\n",n); */
  /* 	printf("\tchecking polygon vertex-by-vertex against old code..."); */
	
  /*     } */

  /*     printf("FAIL\n" */
  /* 	     "vertex %d via new code: (%g,%g,%g)\n" */
  /* 	     "vertex %d via old code: (%g,%g,%g)\n", */
  /* 	     i,plc_M_clist(L->cp[0].vt[i]), */
  /* 	     i,plc_M_clist(Lcmp->cp[0].vt[i])); */

  /*     return false; */

  /*   } */

  /* } */

  if (verbose) {

    printf("\thousekeeping...");

  }

  plc_free(L);
  /* plc_free(Lcmp); */
  /*free(edge_lengths);*/
  free(normals);
  /* tsmcmc_triangulation_free(T); */

  if (verbose) {

    printf("pass\n");
    printf("entire case test...pass\n");
  }

  return true;
  
}

bool assembly_test()
{
  printf("----------------------------------------------\n"
	 "testing polygon assembly code                 \n"
	 "----------------------------------------------\n");

  /* bool test_fantriangulation_action_angle(int n,double *theta,
                                             double *d,bool verbose) */
  
  /* We expect to be provided with ALL the diagonal lengths (starting with 
     |v0-v1| = d[0] = 1.0 and ending with |v0-v{n-1}| = d[n-1] = 1.0), and 
     the n-3 dihedral angles. */

  double equal_d[4] = {1.0,1.7,1.7,1.0};
  double flat_theta[2] = {0,0};

  if (!test_fantriangulation_action_angle(5,flat_theta,equal_d,true)) {

    return false;

  }

  double random_d[5] = {1.0,1.2,1.8,0.9,1.0};
  double random_theta[3] = {0.6,3.1,0.2};

  if (!test_fantriangulation_action_angle(6,random_theta,random_d,true)) {

    return false;

  }

  double long_d[22] = {1.0,1.4,2.1,2.3,1.7,2.3,2.9,3.5,3.2,4.0,4.7,
		       5.1,4.5,3.7,2.9,2.1,1.4,0.7,1.1,1.2,1.4,1.0};
  
  double long_theta[20] = {0.4,2.1,3.1,0.6,1.6,2.3,1.0,4.7,5.2,6.0,
			   0.4,2.1,3.1,0.6,1.6,2.3,1.0,4.7,5.2,6.0};

  if (!test_fantriangulation_action_angle(23,long_theta,long_d,true)) {

    return false;

  }

  return true;
  
}


int main() {

  printf("test_csu_randompolygon generation (%s)\n",PACKAGE_STRING);
  printf("--------------------------------------------------------\n"
	 "Unit tests for direct sampler of equilateral random polygons. \n"
	 "========================================================\n");

  if (!assembly_test() /*|| !randompoly_quality_test() */) {

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
