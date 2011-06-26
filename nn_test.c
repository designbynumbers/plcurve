/* 

   nn_test.c : Nearest neighbor test program. This program runs the plCurve nearest neighbor
   functions through a set of simple test cases to make sure that the code is working. 

*/

#include<config.h>
#include<plCurve.h>

#ifdef HAVE_MATH_H
  #include<math.h>
#endif
#ifdef HAVE_STDLIB_H
  #include<stdlib.h>
#endif
#ifdef HAVE_TIME_H
  #include<time.h>
#endif

#define VERBOSE 1

int gcd(int x, int y)
{
  int t;

  while (y) {
    t = x;
    x = y;
    y = t % y;
  }
  return(x);
}


plCurve *torusknot(int verts,int p,int q,double major_radius,double minor_radius) 

{

  int cp; 
  int vt;
  
  int *nv;
  bool *open;
  int *cc;

  plCurve *tk;

  cp = gcd(p,q);

  int i;

  nv = calloc(cp,sizeof(int));
  open = calloc(cp,sizeof(bool));
  cc = calloc(cp,sizeof(int));

  for (i=0;i<cp;i++) { 

    nv[i] = verts;
    open[i] = false;
    cc[i] = 0;

  }

  tk = plc_new(cp,nv,open,cc);

  /* We are now prepared to build the actual knot */

  double pofs;
  double theta; 
  double tstep;
  double pi = 3.1415926535897932384626433;
  double pangle,qangle;

  for (i=0;i<cp;i++) {

    for(vt=0,theta=0,pofs=i*(2*pi/cp),tstep = 2*pi/(double)(verts);
	vt<verts;
	vt++,theta+=tstep) {

      plc_vector loc;

      pangle = pofs + p*theta;
      qangle = q*theta;

      loc = plc_build_vect(
			   major_radius*cos(qangle)*(1+(minor_radius/major_radius)*cos(pangle)),
			   major_radius*sin(qangle)*(1+(minor_radius/major_radius)*cos(pangle)),
			   minor_radius*sin(pangle)
			   );

      tk->cp[i].vt[vt] = loc;

    }
    
  }

  plc_fix_wrap(tk);

  free(nv);
  free(open);
  free(cc);

  return tk;

}
  

plCurve *square() {

  int nv = 4;
  int cc = 0;
  bool open = {false};
  
  plCurve *build;

  build = plc_new(1,&nv,&open,&cc);

  build->cp[0].vt[0] = plc_build_vect(0,0,0);
  build->cp[0].vt[1] = plc_build_vect(1,0,0);
  build->cp[0].vt[2] = plc_build_vect(1,1,0);
  build->cp[0].vt[3] = plc_build_vect(0,1,0);

  return build;

}

plc_vector *sqbuf() {

  plc_vector *build;
  build = calloc(4,sizeof(plc_vector));

  build[0] = plc_build_vect(0,0,0);
  build[1] = plc_build_vect(1,0,0);
  build[2] = plc_build_vect(1,1,0);
  build[3] = plc_build_vect(0,1,0);

  return build;

}

int exhaustive_nearest_neighbor(plc_vector query,int n,plc_vector *buffer,int stride);

void square_buffer_tests(plc_vector pt,int target) {

  struct plc_nearest_neighbor_pc_data *pc_data = NULL;
  plc_vector *buffer = NULL;
  int plc_error = 0;

  buffer = sqbuf();

  /* Test 0: Exhaustive nearest neighbor on square*/
  
  if (exhaustive_nearest_neighbor(pt,4,buffer,1) != target) {

    exit(1);
    printf("nn_test: FAILED exhaustive test on (%g,%g,%g) and square buffer.\n",pt.c[0],pt.c[1],pt.c[2]);

  } else {

#ifdef VERBOSE
    printf("nn_test: passed exhaustive test on (%g,%g,%g) and square buffer.\n",pt.c[0],pt.c[1],pt.c[2]);
#endif
    
  }

  /* Test 1: Fancy nearest neighbor on square and 0.1, 0.1, 0.1 */

  if (plc_nearest_neighbor(pt,4,buffer,&pc_data,&plc_error) != target) {

    exit(1);
    printf("nn_test: FAILED plc_nearest_neighbor test on (%g,%g,%g) and square buffer.\n",pt.c[0],pt.c[1],pt.c[2]);

  } else {

#ifdef VERBOSE
    printf("nn_test: passed plc_nearest_neighbor test on (%g,%g,%g) and square buffer.\n",pt.c[0],pt.c[1],pt.c[2]);
#endif

  }

  /* Test 2: Fancy nearest neighbor on square and 0.1, 0.1, 0.1 using precomputed data. */

   if (plc_nearest_neighbor(pt,4,buffer,&pc_data,&plc_error) != target) {

    exit(1);
    printf("nn_test: FAILED plc_nearest_neighbor test on (%g,%g,%g) and square buffer using precomputed data.\n",pt.c[0],pt.c[1],pt.c[2]);

  } else {

#ifdef VERBOSE
    printf("nn_test: passed plc_nearest_neighbor test on (%g,%g,%g) and square buffer using precomputed data.\n",pt.c[0],pt.c[1],pt.c[2]);
#endif

  }

   /* Test 3: Fancy nearest neighbor on square and 0.1, 0.1, 0.1 using stale precomputed data. */

   pc_data->check_buffer = NULL; /* Force data to look stale */

   if (plc_nearest_neighbor(pt,4,buffer,&pc_data,&plc_error) != target) {

    exit(1);
    printf("nn_test: FAILED plc_nearest_neighbor test on (%g,%g,%g) and square buffer using STALE precomputed data.\n",pt.c[0],pt.c[1],pt.c[2]);

   } else {

#ifdef VERBOSE
    printf("nn_test: passed plc_nearest_neighbor test on (%g,%g,%g) and square buffer using STALE precomputed data.\n",pt.c[0],pt.c[1],pt.c[2]);
#endif

  }

   /* Test 4: Fancy nearest neighbor on square and 0.1, 0.1, 0.1 using precomputed data recovered from stale. */

   if (plc_nearest_neighbor(pt,4,buffer,&pc_data,&plc_error) != target) {

    exit(1);
    printf("nn_test: FAILED plc_nearest_neighbor test on (%g,%g,%g) and square buffer using precomputed data recovered from stale.\n",pt.c[0],pt.c[1],pt.c[2]);

  } else {

#ifdef VERBOSE
    printf("nn_test: passed plc_nearest_neighbor test on (%g,%g,%g) and square buffer using precomputed data recovered from stale.\n",pt.c[0],pt.c[1],pt.c[2]);
#endif

  }
  
   /* Cleanup */

   plc_nearest_neighbor_pc_data_free(&pc_data);
   free(buffer);

} 

void square_plCurve_tests(plc_vector pt,int target) {

  struct plc_nearest_vertex_pc_data *pc_data = NULL;
  plCurve *L = NULL, *nL = NULL;
  int plc_error = 0;
  int cp, vt;

  L = square();

  /* Test 0: Exhaustive nearest neighbor on square*/
  
  if (exhaustive_nearest_neighbor(pt,4,L->cp[0].vt,1) != target) {

    exit(1);
    printf("nn_test: FAILED exhaustive test on (%g,%g,%g) and square L.\n",pt.c[0],pt.c[1],pt.c[2]);

  } else {

#ifdef VERBOSE
    printf("nn_test: passed exhaustive test on (%g,%g,%g) and square L.\n",pt.c[0],pt.c[1],pt.c[2]);
#endif

  }

  /* Test 1: Fancy nearest neighbor on square and 0.1, 0.1, 0.1 */

  plc_nearest_vertex(pt,L,&cp,&vt,&pc_data,&plc_error);

  if (cp != 0 || vt != target) {

    exit(1);
    printf("nn_test: FAILED plc_nearest_vertex test on (%g,%g,%g) and square L.\n",pt.c[0],pt.c[1],pt.c[2]);

  } else {

#ifdef VERBOSE
    printf("nn_test: passed plc_nearest_vertex test on (%g,%g,%g) and square L.\n",pt.c[0],pt.c[1],pt.c[2]);
#endif

  }

  /* Test 2: Fancy nearest neighbor on square and 0.1, 0.1, 0.1 using precomputed data. */

   plc_nearest_vertex(pt,L,&cp,&vt,&pc_data,&plc_error);

  if (cp != 0 || vt != target) {

    exit(1);
    printf("nn_test: FAILED plc_nearest_vertex test on (%g,%g,%g) and square L using precomputed data.\n",pt.c[0],pt.c[1],pt.c[2]);

  } else {

#ifdef VERBOSE
    printf("nn_test: passed plc_nearest_vertex test on (%g,%g,%g) and square L using precomputed data.\n",pt.c[0],pt.c[1],pt.c[2]);
#endif

  }

   /* Test 3: Fancy nearest neighbor on square and 0.1, 0.1, 0.1 using stale precomputed data. */

  nL = plc_copy(L);
  pc_data->check_curve = nL; /* Force data to look stale (as if it applies to nL) */

  plc_nearest_vertex(pt,L,&cp,&vt,&pc_data,&plc_error);

  if (cp != 0 || vt != target) {

    exit(1);
    printf("nn_test: FAILED plc_nearest_vertex test on (%g,%g,%g) and square L using STALE precomputed data.\n",pt.c[0],pt.c[1],pt.c[2]);

  } else {

#ifdef VERBOSE
    printf("nn_test: passed plc_nearest_vertex test on (%g,%g,%g) and square L using STALE precomputed data.\n",pt.c[0],pt.c[1],pt.c[2]);
#endif

  }

   /* Test 4: Fancy nearest neighbor on square and 0.1, 0.1, 0.1 using precomputed data recovered from stale. */

   plc_nearest_vertex(pt,L,&cp,&vt,&pc_data,&plc_error);

  if (cp != 0 || vt != target) {

    exit(1);
    printf("nn_test: FAILED plc_nearest_vertex test on (%g,%g,%g) and square L using precomputed data recovered from stale.\n",pt.c[0],pt.c[1],pt.c[2]);

  } else {

#ifdef VERBOSE
    printf("nn_test: passed plc_nearest_vertex test on (%g,%g,%g) and square L using precomputed data recovered from stale.\n",pt.c[0],pt.c[1],pt.c[2]);
#endif

  }
  
   /* Cleanup */

   plc_nearest_vertex_pc_data_free(&pc_data);
   plc_free(L);
   plc_free(nL);

} 

void random_buffer_test(int npoints,int seed, int ntests) 

  /* Generates a random buffer of points and queries against this buffer, comparing results from exhaustion 
     with results from our function. Also does some rudimentary timing tests. */

{

  plc_vector *buffer;
  plc_vector pt;
  int i;

  srand(seed); /* We WANT this to be repeatable for debugging purposes. */

  #ifdef VERBOSE

  printf("Generating random buffer of %d points with seed %d.\n",npoints,seed);

  #endif

  buffer = calloc(npoints,sizeof(plc_vector));

  for(i=0;i<npoints;i++) {

    buffer[i] = plc_build_vect((1.0)*rand()/RAND_MAX,(1.0)*rand()/RAND_MAX,(1.0)*rand()/RAND_MAX);
    
  }

  pt = plc_build_vect((1.0)*rand()/RAND_MAX,(1.0)*rand()/RAND_MAX,(1.0)*rand()/RAND_MAX);


  printf("Testing random %d point buffer against exhaustive search.\n",npoints);

  clock_t start, end;
  double exhaustive_cpu_time_used,nn_cpu_time_used;
  int test;
  struct plc_nearest_neighbor_pc_data *pc_data = NULL;
  int plc_error;

  int exhaustive,nn;

  printf("Exhaustive (sec) | plCurve (sec) | Correct Results \n");
  printf("---------------------------------------------------\n");

  for(test=0;test<ntests;test++) {

    /* Compute using the exhaustive solution (presumed correct). */
     
    start = clock(); 
    exhaustive = exhaustive_nearest_neighbor(pt,npoints,buffer,1);
    end = clock();

    exhaustive_cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  
    /* Now compute with our method. */
  
    start = clock(); 
    nn = plc_nearest_neighbor(pt,npoints,buffer,&pc_data,&plc_error);
    end = clock();

    nn_cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    /* Print a report line. */

    if (nn == exhaustive) {

      printf("%17g  %14g   Pass\n",exhaustive_cpu_time_used,nn_cpu_time_used);

    } else {

      printf("Test FAILED with query point (%g,%g,%g) on %d point buffer generated with seed %d.\n",
	     pt.c[0],pt.c[1],pt.c[2],npoints,seed);
      exit(1);

    }

    /* Generate a new point for the next round. */

    pt = plc_build_vect((1.0)*rand()/RAND_MAX,(1.0)*rand()/RAND_MAX,(1.0)*rand()/RAND_MAX);

  }

  printf("Random round concluded: PASS\n\n");

  plc_nearest_neighbor_pc_data_free(&pc_data);
  free(buffer);

}  

void random_torusknot_test(int npoints,int seed, int p, int q, int ntests) 

  /* Generates a random buffer of points and queries against this buffer, comparing results from exhaustion 
     with results from our function. Also does some rudimentary timing tests. */

{

  plCurve *L;
  plc_vector pt;
  int i;

  srand(seed); /* We WANT this to be repeatable for debugging purposes. */

  #ifdef VERBOSE

  printf("Generating random (%d,%d) torus link of %d points with seed %d.\n",p,q,npoints,seed);

  #endif

  L = torusknot(npoints,p,q,2+((1.0)*rand())/RAND_MAX,(1.0*rand())/RAND_MAX);
  pt = plc_build_vect((6.0)*rand()/RAND_MAX - 3,(6.0)*rand()/RAND_MAX - 3,(6.0)*rand()/RAND_MAX - 3);

  if (gcd(p,q) == 1) {

    printf("Testing random %d point torus knot against exhaustive search.\n",npoints);

  } else {

    printf("Testing random %d point torus link for speed and stability.\n",npoints);

  }

  clock_t start, end;
  double exhaustive_cpu_time_used,nn_cpu_time_used;
  int test;
  struct plc_nearest_vertex_pc_data *pc_data = NULL;
  int plc_error;

  int exhaustive,nn;

  printf("Exhaustive (sec) | plCurve (sec) | Correct Results \n");
  printf("---------------------------------------------------\n");

  for(test=0;test<ntests;test++) {

    /* Compute using the exhaustive solution (presumed correct). */
   
    if (gcd(p,q) == 1) {

      start = clock(); 
      exhaustive = exhaustive_nearest_neighbor(pt,npoints,L->cp[0].vt,1);
      end = clock();
      
      exhaustive_cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    } else {

      exhaustive_cpu_time_used = 0;

    }
  
    /* Now compute with our method. */

    int cp, vt;
  
    start = clock(); 
    plc_nearest_vertex(pt,L,&cp,&vt,&pc_data,&plc_error);
    end = clock();

    nn_cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    /* Print a report line. */

    if (gcd(p,q) == 1) {

      if (vt == exhaustive) {
	
	printf("%17g  %14g   Pass\n",exhaustive_cpu_time_used,nn_cpu_time_used);
	
      } else {
	
	printf("Test FAILED with query point (%g,%g,%g) on %d point (%d,%d) torus knot generated with seed %d.\n",
	       pt.c[0],pt.c[1],pt.c[2],npoints,p,q,seed);
	exit(1);
	
      }

    } else {

     	printf("(not run)          %14g   Pass\n",nn_cpu_time_used);

    }

    /* Generate a new point for the next round. */

    pt = plc_build_vect((6.0)*rand()/RAND_MAX - 3,(6.0)*rand()/RAND_MAX - 3,(6.0)*rand()/RAND_MAX - 3);

  }

  printf("Random torusknot round concluded: PASS\n\n");

  plc_nearest_vertex_pc_data_free(&pc_data);
  plc_free(L);

}  

int main () {

  printf("nn_test: Nearest neighbor test suite for plCurve.\n");

  /* Square buffer tests. */

  square_buffer_tests(plc_build_vect(0.1,0.1,0.1),0);
  square_buffer_tests(plc_build_vect(-0.1,0.1,0.1),0);
  square_buffer_tests(plc_build_vect(0.1,0.1,-0.1),0);
  square_buffer_tests(plc_build_vect(-0.1,-0.1,-0.1),0);

  square_buffer_tests(plc_build_vect(1.1,0.1,0.1),1);
  square_buffer_tests(plc_build_vect(1.1,0.1,-0.1),1);
  square_buffer_tests(plc_build_vect(0.9,-0.1,0.1),1);
  square_buffer_tests(plc_build_vect(0.9,0.1,-0.1),1);

  square_buffer_tests(plc_build_vect(5,5,0.1),2);
  square_buffer_tests(plc_build_vect(1.1,1.1,1.1),2);
  square_buffer_tests(plc_build_vect(1.1,0.9,-0.1),2);
  square_buffer_tests(plc_build_vect(0.9,0.9,0.1),2);

  printf("Passed plc_nearest_neighbor tests on 4 vertex (square) buffer.\n");

  /* Square plCurve tests. */

  square_plCurve_tests(plc_build_vect(0.1,0.1,0.1),0);
  square_plCurve_tests(plc_build_vect(-0.1,0.1,0.1),0);
  square_plCurve_tests(plc_build_vect(0.1,0.1,-0.1),0);
  square_plCurve_tests(plc_build_vect(-0.1,-0.1,-0.1),0);

  square_plCurve_tests(plc_build_vect(1.1,0.1,0.1),1);
  square_plCurve_tests(plc_build_vect(1.1,0.1,-0.1),1);
  square_plCurve_tests(plc_build_vect(0.9,-0.1,0.1),1);
  square_plCurve_tests(plc_build_vect(0.9,0.1,-0.1),1);

  square_plCurve_tests(plc_build_vect(5,5,0.1),2);
  square_plCurve_tests(plc_build_vect(1.1,1.1,1.1),2);
  square_plCurve_tests(plc_build_vect(1.1,0.9,-0.1),2);
  square_plCurve_tests(plc_build_vect(0.9,0.9,0.1),2);

  printf("Passed plc_nearest_neighbor tests on 4 vertex (square) plCurve.\n");

  /* Random buffer tests. */

  random_buffer_test(100,117,10);
  random_buffer_test(500,131,10);
  random_buffer_test(1000,20397,10);
  random_buffer_test(10000,20397,10);
  random_buffer_test(30000,20397,10);
  random_buffer_test(50000,27397,10);

  printf("Passed all random buffer tests.\n");

  /* Random knot tests */

  random_torusknot_test(100,12898,3,2,10);
  random_torusknot_test(500,464,3,2,10);
  random_torusknot_test(1000,129387,3,2,10);
  random_torusknot_test(5000,872638,3,2,10);
  random_torusknot_test(10000,123897,3,2,10);

  random_torusknot_test(100,12898,6,14,10);
  random_torusknot_test(500,464,6,14,10);
  random_torusknot_test(1000,129387,6,14,10);
  random_torusknot_test(5000,872638,6,14,10);
  random_torusknot_test(10000,123897,6,14,10);

  printf("Passed all random torus knot tests.\n");


  printf("All nearest neighbor tests PASSED.\n");

  exit(0);

}
