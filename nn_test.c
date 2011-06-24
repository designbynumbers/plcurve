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

  struct plc_nearest_neighbor_pc_data *pc_data;
  plc_vector *buffer = NULL;
  int plc_error = 0;

  buffer = sqbuf();

  /* Test 0: Exhaustive nearest neighbor on square*/
  
  if (exhaustive_nearest_neighbor(pt,4,buffer,1) != target) {

    exit(1);
    printf("nn_test: FAILED exhaustive test on (%g,%g,%g) and square buffer.\n",pt.c[0],pt.c[1],pt.c[2]);

  } else {

    printf("nn_test: passed exhaustive test on (%g,%g,%g) and square buffer.\n",pt.c[0],pt.c[1],pt.c[2]);

  }

  /* Test 1: Fancy nearest neighbor on square and 0.1, 0.1, 0.1 */

  if (plc_nearest_neighbor(pt,4,buffer,&pc_data,&plc_error) != target) {

    exit(1);
    printf("nn_test: FAILED plc_nearest_neighbor test on (%g,%g,%g) and square buffer.\n",pt.c[0],pt.c[1],pt.c[2]);

  } else {

    printf("nn_test: passed plc_nearest_neighbor test on (%g,%g,%g) and square buffer.\n",pt.c[0],pt.c[1],pt.c[2]);

  }

  /* Test 2: Fancy nearest neighbor on square and 0.1, 0.1, 0.1 using precomputed data. */

   if (plc_nearest_neighbor(pt,4,buffer,&pc_data,&plc_error) != target) {

    exit(1);
    printf("nn_test: FAILED plc_nearest_neighbor test on (%g,%g,%g) and square buffer using precomputed data.\n",pt.c[0],pt.c[1],pt.c[2]);

  } else {

    printf("nn_test: passed plc_nearest_neighbor test on (%g,%g,%g) and square buffer using precomputed data.\n",pt.c[0],pt.c[1],pt.c[2]);

  }

   /* Test 3: Fancy nearest neighbor on square and 0.1, 0.1, 0.1 using stale precomputed data. */

   pc_data->check_buffer = NULL; /* Force data to look stale */

   if (plc_nearest_neighbor(pt,4,buffer,&pc_data,&plc_error) != target) {

    exit(1);
    printf("nn_test: FAILED plc_nearest_neighbor test on (%g,%g,%g) and square buffer using STALE precomputed data.\n",pt.c[0],pt.c[1],pt.c[2]);

  } else {

    printf("nn_test: passed plc_nearest_neighbor test on (%g,%g,%g) and square buffer using STALE precomputed data.\n",pt.c[0],pt.c[1],pt.c[2]);

  }

   /* Test 4: Fancy nearest neighbor on square and 0.1, 0.1, 0.1 using precomputed data recovered from stale. */

   if (plc_nearest_neighbor(pt,4,buffer,&pc_data,&plc_error) != target) {

    exit(1);
    printf("nn_test: FAILED plc_nearest_neighbor test on (%g,%g,%g) and square buffer using precomputed data recovered from stale.\n",pt.c[0],pt.c[1],pt.c[2]);

  } else {

    printf("nn_test: passed plc_nearest_neighbor test on (%g,%g,%g) and square buffer using precomputed data recovered from stale.\n",pt.c[0],pt.c[1],pt.c[2]);

  }
  
   /* Cleanup */

   plc_nearest_neighbor_pc_data_free(&pc_data);
   free(buffer);

} 

int main () {

  printf("nn_test: Nearest neighbor test suite for plCurve.\n");

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

  printf("All nearest neighbor tests PASSED.\n");

  exit(0);

}
