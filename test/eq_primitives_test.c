/* 

   eq_primitives_test.c : Run some sanity checks on the primitives for
   equilateral polygon generation.

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

#ifdef HAVE_STDIO_H
#include<stdio.h>
#endif

#ifdef HAVE_GSL_GSL_RNG_H
#include<gsl/gsl_rng.h>
#endif

#ifdef HAVE_GSL_GSL_RANDIST_H
#include<gsl/gsl_randist.h>
#endif

gsl_rng *r; /* The global random number generator */

// Process has done x out of n rounds,
// and we want a bar of width w and resolution r.

static inline void pbar(double current,int x, int n, int r, int w)
{
  int i;
  
  // Only update r times.
  if (n > 2*r) {  if ( x % (n/r) != 0 ) return; }
  
  // Calculuate the ratio of complete-to-incomplete.
  float ratio = x/(float)(n);
  int   c     = ratio * w;
  
  // Show the percentage complete.
  printf("%3d%% [", (int)(ratio*100) );
  
  // Show the load bar.
  for (i=0; i<c; i++)
    printf("=");
  
  for (i=c; i<w; i++)
    printf(" ");
  
  if (x==n) { printf("]"); }

  // ANSI Control codes to go back to the
  // previous line and clear it.
  // printf("]\n33[F33[J");
  
  printf("\r"); // Move to the first column
  fflush(stdout);
}

plc_vector plc_subsphere_sample(gsl_rng *r,plc_vector axis,double height);


bool test_height_and_sphereness(gsl_rng *r,plc_vector axis,double height,int nsamps) 

{
  bool ok;
  int i;
  plc_vector samp;

  axis = plc_normalize_vect(axis,&ok);

  printf("Generating samples to test sphereness, position.\n"
	 "Axis = (%g,%g,%g) \n"
	 "Height = %g \n"
	 ,plc_M_clist(axis),height);

  for(i=0;i<nsamps;i++) {

    double dist;

    samp = plc_subsphere_sample(r,axis,height);
    dist = plc_norm(samp);

    if (fabs(dist - 1.0) > 1e-8) { 

      printf("FAIL: Distance %g of sample %d from origin not close to 1.0\n",
	     dist,i);
      return false;

    }

    if (plc_dot_prod(axis,samp) < -1.0 || plc_dot_prod(axis,samp) > height) {

      printf("FAIL: Height %g of sample %d along axis not below target %g\n",
	     plc_dot_prod(axis,samp),i,height);
      return false;

    }
    
    pbar(0,i,nsamps,100,20);

  }

  pbar(0,nsamps,nsamps,100,20);


  printf("\n"
	 "pass: sphereness and height ok over %d samples.\n\n",nsamps);
  return true;

}

bool test_mean(gsl_rng *r,plc_vector axis,double height,int nsamps) 

{
  bool ok;
  int i;
  plc_vector samp;

  axis = plc_normalize_vect(axis,&ok);

  printf("Generating samples to test mean.\n"
	 "Axis = (%g,%g,%g) \n"
	 "Height = %g \n"
	 ,plc_M_clist(axis),height);

  double axis_height_mean = 0; 

  for(i=0;i<nsamps;i++) {

    samp = plc_subsphere_sample(r,axis,height);
    axis_height_mean += plc_dot_prod(axis,samp);
    pbar(0,i,nsamps,100,20);

  }

  pbar(0,nsamps,nsamps,100,20);
  axis_height_mean /= (double)(nsamps);

  if (fabs((height - 1.0)/2.0) < 1e-5) { /* Mean close to 0, use absolute error */

    if (fabs(axis_height_mean) > 1e-2) { 

      printf("FAIL: Height mean of %g not close to 0.\n",axis_height_mean);
      return false;

    } else {
      
      printf("pass: Height mean %g close to target of 0 over %d samps.\n",
	     axis_height_mean,nsamps);
      return true;

    }

  } else {

    double target_mean;
    target_mean = (height - 1.0)/2.0;
    double error; 

    error = fabs((axis_height_mean - target_mean)/target_mean);

    if (error > 0.01) {

      printf("FAIL: Height mean of %g not close to %g.\n",
	     axis_height_mean,target_mean);
      return false;

    } else {
      
      printf("pass: Height mean %g close to target of %g over %d samps.\n\n",
	     axis_height_mean,target_mean,nsamps);
      return true;

    }

  }

}

bool test_plc_subsphere_sample(gsl_rng *r,int nSamples) {

  plc_vector axis = {{1.0/sqrt(3.0),1.0/sqrt(3.0),1.0/sqrt(3.0)}};
  plc_vector z_axis = {{0,0,1.0}};
  plc_vector minusz = {{0,0,-1.0}};

  printf("Subsphere_sample TEST SUITE\n"
	 "---------------------------------------------\n");

  return (test_height_and_sphereness(r,axis,-0.5,nSamples) &&
	  test_height_and_sphereness(r,axis,0.5,nSamples) &&
	  test_height_and_sphereness(r,axis,-0.8,nSamples) &&
	  test_height_and_sphereness(r,axis,1.0,nSamples) 
	  &&
	  test_height_and_sphereness(r,z_axis,-0.5,nSamples) &&
	  test_height_and_sphereness(r,z_axis,0.5,nSamples) &&
	  test_height_and_sphereness(r,z_axis,-0.8,nSamples) &&
	  test_height_and_sphereness(r,z_axis,1.0,nSamples) 
	  &&
	  test_height_and_sphereness(r,minusz,-0.5,nSamples) &&
	  test_height_and_sphereness(r,minusz,0.5,nSamples) &&
	  test_height_and_sphereness(r,minusz,-0.8,nSamples) &&
	  test_height_and_sphereness(r,minusz,1.0,nSamples) 
	  &&
	  test_mean(r,axis,-0.25,nSamples) &&
	  test_mean(r,axis,-0.75,nSamples) &&
	  test_mean(r,axis,-0.625,nSamples) &&
	  test_mean(r,axis,0.33,nSamples)  
	  &&
	  test_mean(r,z_axis,-0.25,nSamples) &&
	  test_mean(r,z_axis,-0.75,nSamples) &&
	  test_mean(r,z_axis,-0.625,nSamples) &&
	  test_mean(r,z_axis,0.33,nSamples)  
	  &&
	  test_mean(r,minusz,-0.25,nSamples) &&
	  test_mean(r,minusz,-0.75,nSamples) &&
	  test_mean(r,minusz,-0.625,nSamples) &&
	  test_mean(r,minusz,0.33,nSamples)  
	  );
}


plc_vector plc_sphere_intersection_sample(gsl_rng *r,plc_vector A, plc_vector B);

bool test_center_dists(gsl_rng *r,plc_vector A,plc_vector B,int nsamps) {

  int i;
  plc_vector samp;

  printf("Generating samples to test distance from spheres.\n"
	 "Center A = (%g,%g,%g) \n"
	 "Center B = (%g,%g,%g) \n"
	 ,plc_M_clist(A),plc_M_clist(B));

  for(i=0;i<nsamps;i++) {

    double distA,distB;

    samp = plc_sphere_intersection_sample(r,A,B);

    distA = plc_distance(samp,A);
    distB = plc_distance(samp,B);

    if (fabs(distA - 1.0) > 1e-5 || fabs(distB - 1.0) > 1e-5) {

      printf("FAIL: sample %d distances from centers are %g, %g, not 1.0\n",
	     i,distA,distB);
      return false;

    }

    pbar(0,i,nsamps,100,20);

  }

  pbar(0,nsamps,nsamps,100,20);
       
  printf("pass: Sample distances from both centers close to 1 over %d samps.\n\n",
	 nsamps);
  return true;
  
}

bool test_angle_mean(gsl_rng *r,plc_vector A,plc_vector B,int nsamps) 

{
  int i;
  plc_vector samp,samp2;
  double pi = 3.14159265358979;
  double ABdist;

  printf("Generating samples to test mean angle.\n"
	 "Center A = (%g,%g,%g) \n"
	 "Center B = (%g,%g,%g) \n"
	 ,plc_M_clist(A),plc_M_clist(B));

  ABdist = plc_distance(A,B);

  if (ABdist > 1.99) { 

    printf("pass: Centers are too far apart to run this test.\n");

  }

  if (ABdist < 1e-2) {

    printf("pass: Centers are too close to run this test.\n");

  }

  printf("Center-center dist: %g\n",ABdist);

  bool ok;
  double angle_mean = 0; 
  plc_vector axis = plc_vect_diff(A,B);

  axis = plc_normalize_vect(axis,&ok);

  for(i=0;i<nsamps;i++) {

    bool ok;

    samp  = plc_vect_diff(plc_sphere_intersection_sample(r,A,B),B);
    samp2 = plc_vect_diff(plc_sphere_intersection_sample(r,A,B),B);

    angle_mean += 
      fabs(plc_angle(plc_cross_prod(axis,samp),plc_cross_prod(axis,samp2),&ok));
    
    if (!ok) {
      printf("FAIL: Couldn't measure angle.\n");
      exit(1);
    }
    
    pbar(0,i,nsamps,100,20);

  }

  pbar(0,nsamps,nsamps,100,20);
  angle_mean /= (double)(nsamps);

  if (fabs(angle_mean - pi/2.0) > 1e-2) { 

      printf("\n"
	     "FAIL: Angle mean of %g not close to pi/2.\n",angle_mean);
      return false;

  } else {
      
    printf("pass: Angle mean of %g close to target of pi/2 = %g over %d samps.\n\n",
	   angle_mean,pi/2.0,nsamps);
    return true;

  }

}


bool test_plc_sphere_intersection_sample(gsl_rng *r,int nSamples) {

  plc_vector A,B;

  gsl_ran_dir_3d(r,&(A.c[0]),&(A.c[1]),&(A.c[2]));
  gsl_ran_dir_3d(r,&(B.c[0]),&(B.c[1]),&(B.c[2]));
  B = plc_scale_vect(0.5,B);

  plc_vector AP = {{0,0,1.0}}, BP = {{0,0,3.0}};
  plc_vector CP = {{0,2.000001,3.0}};

  printf("Sphere_intersection_sample TEST SUITE\n"
	 "---------------------------------------------\n");

  return (
	  test_center_dists(r,A,B,nSamples) &&
	  test_center_dists(r,AP,BP,nSamples) &&
	  test_center_dists(r,AP,AP,nSamples) &&
	  test_center_dists(r,CP,BP,nSamples) &&
	  test_angle_mean(r,A,B,5*nSamples)
	  );
}

/* Double rotation tests. */

void plc_double_rotation(gsl_rng *r, 
			 plc_vector r1,plc_vector r2,plc_vector r3,
			 plc_vector *r1p, plc_vector *r2p, plc_vector *r3p);

bool test_double_rotation(gsl_rng *r,plc_vector r1,plc_vector r2,plc_vector r3,int nsamps) {
  
  int i;
  plc_vector r1p,r2p,r3p;

  printf("Generating double rotations to test geometry.\n"
	 "r1 = (%g,%g,%g) \n"
	 "r2 = (%g,%g,%g) \n"
	 "r2 = (%g,%g,%g) \n"
	 ,plc_M_clist(r1),plc_M_clist(r2),plc_M_clist(r3));

  plc_vector R;
  R = plc_vect_sum(r1,plc_vect_sum(r2,r3));

  for(i=0;i<nsamps;i++) {

    plc_vector Rp;

    plc_double_rotation(r,r1,r2,r3,&r1p,&r2p,&r3p);
    Rp = plc_vect_sum(r1p,plc_vect_sum(r2p,r3p));

    if (plc_distance(R,Rp) > 1e-5) {

      printf("FAIL: Double rotation didn't preserve sum.\n");
      return false;

    }

    if (fabs(plc_norm(r1p) - 1.0) > 1e-8 ||
	fabs(plc_norm(r2p) - 1.0) > 1e-8 ||
	fabs(plc_norm(r3p) - 1.0) > 1e-8) { 
     
      printf("FAIL: Double rotation didn't preserve norms.\n");
      return false;

    }
    
    pbar(0,i,nsamps,100,20);

  }

  pbar(0,nsamps,nsamps,100,20);
       
  printf("pass: Sum and vector norms preserved over %d samps.\n\n",
	 nsamps);
  return true;
  
}

bool test_plc_double_rotation(gsl_rng *r,int nSamples) {

  plc_vector A,B,C;

  gsl_ran_dir_3d(r,&(A.c[0]),&(A.c[1]),&(A.c[2]));
  gsl_ran_dir_3d(r,&(B.c[0]),&(B.c[1]),&(B.c[2]));
  gsl_ran_dir_3d(r,&(C.c[0]),&(C.c[1]),&(C.c[2]));

  plc_vector X = {{1.0,0,0}}, minusX = {{-1.0,0,0}};
  plc_vector Z = {{0,0,1.0}};

  printf("plc_double_rotation TEST SUITE\n"
	 "---------------------------------------------\n");

  return (
	  test_double_rotation(r,A,B,C,nSamples) &&
	  test_double_rotation(r,X,X,X,5) &&
	  test_double_rotation(r,X,minusX,X,5) &&
	  test_double_rotation(r,X,minusX,Z,nSamples) );
	 
}

/* Single rotation tests */

void plc_single_rotation(gsl_rng *r, 
			 plc_vector r1,plc_vector r2,
			 plc_vector *r1p, plc_vector *r2p);

bool test_single_rotation(gsl_rng *r,plc_vector r1,plc_vector r2,int nsamps) {
  
  int i;
  plc_vector r1p,r2p;

  printf("Generating single rotations to test geometry.\n"
	 "r1 = (%g,%g,%g) \n"
	 "r2 = (%g,%g,%g) \n"
	 ,plc_M_clist(r1),plc_M_clist(r2));

  plc_vector R;
  R = plc_vect_sum(r1,r2);

  for(i=0;i<nsamps;i++) {

    plc_vector Rp;

    plc_single_rotation(r,r1,r2,&r1p,&r2p);
    Rp = plc_vect_sum(r1p,r2p);

    if (plc_distance(R,Rp) > 1e-5) {

      printf("FAIL: Single rotation didn't preserve sum.\n");
      return false;

    }

    if (fabs(plc_norm(r1p) - 1.0) > 1e-8 ||
	fabs(plc_norm(r2p) - 1.0) > 1e-8 ) {
     
      printf("FAIL: Single rotation didn't preserve norms.\n");
      return false;

    }
    
    pbar(0,i,nsamps,100,20);

  }

  pbar(0,nsamps,nsamps,100,20);
       
  printf("pass: Sum and vector norms preserved over %d samps.\n\n",
	 nsamps);
  return true;
  
}

bool test_plc_single_rotation(gsl_rng *r,int nSamples) {

  plc_vector A,B;

  gsl_ran_dir_3d(r,&(A.c[0]),&(A.c[1]),&(A.c[2]));
  gsl_ran_dir_3d(r,&(B.c[0]),&(B.c[1]),&(B.c[2]));
 
  plc_vector X = {{1.0,0,0}}, minusX = {{-1.0,0,0}};
  plc_vector Z = {{0,0,1.0}};

  printf("plc_single_rotation TEST SUITE\n"
	 "---------------------------------------------\n");

  return (
	  test_single_rotation(r,A,B,nSamples) &&
	  test_single_rotation(r,X,X,5) &&
	  test_single_rotation(r,X,minusX,5) &&
	  test_single_rotation(r,X,Z,nSamples) 
	  );
	 
}



int main(int argc, char *argv[]) {

  bool PASS = {true};

  const gsl_rng_type * T;
     
  gsl_rng_env_setup();
  
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  int seed = time(0);
  gsl_rng_set(r,seed);

  printf("EQ_PRIMITIVES_TEST\n\n"
	 "Testing geometric primitives for equilateral polygon generation.\n");
  printf("with %s random number generator, seeded with %d.\n",gsl_rng_name(r),seed);
  printf("===================================================\n\n");

  PASS = test_plc_subsphere_sample(r,100000) &&
    test_plc_sphere_intersection_sample(r,100000) &&
    test_plc_double_rotation(r,100000) &&
    test_plc_single_rotation(r,100000);

  gsl_rng_free(r);

  printf("=================================\n");

  if (PASS) { 

    printf("EQ_PRIMITIVES_TEST: PASS\n");
    exit(0); 

  } else { 

    printf("EQ_PRIMITIVES_TEST: FAIL\n");
    exit(1); 

  }

}
