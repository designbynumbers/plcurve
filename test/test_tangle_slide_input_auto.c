
/*

   test_tangle_slide_input_auto.c : This code is auto-generated by 'assemble_tangleregentest.py' and shouldn't be edited by hand. It contains tests for the gatekeeper code in pd_tangle_slide.


*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
#endif

#ifdef HAVE_STDIO_H
   #include<stdio.h>
#endif

#ifdef HAVE_ASSERT_H
   #include<assert.h>
#endif

#ifdef HAVE_STRING_H
   #include<string.h>
#endif

#ifdef HAVE_STDINT_H
   #include<stdint.h>
#endif

#ifdef HAVE_STDLIB_H
   #include<stdlib.h>
#endif

#ifdef HAVE_STDBOOL_H
   #include<stdbool.h>
#endif

#include<plcTopology.h>

int PD_VERBOSE=15;

/* We need to include a prototype for the function we're testing, because
   it's not exposed in the header files. */

bool pdint_check_tslide_data_ok_and_find_te(pd_code_t *pd,pd_tangle_t *t,
					    pd_idx_t n,
					    pd_idx_t *overstrand_edges, 
					    pd_idx_t *border_faces,
					    pd_idx_t **tangle_slide_edges,
					    pd_idx_t **complementary_edges,
					    pd_boundary_or_t **complementary_or,
					    bool *overstrand_goes_OVER,
					    pd_or_t *overstrand_orientation);

int main() {

  printf("test_tangle_slide_input (%s)\n",PACKAGE_STRING);
  printf("--------------------------------------------------------\n"
	 "Unit tests for pdint_check_tslide_data_ok_and_find_te(). \n"
	 "========================================================\n");

  if () {

    printf("=======================================================\n");
    printf("test_tangle_slide_input:  PASS.\n");
    exit(0);

  } else {

    printf("=====================================================\n");
    printf("test_tangle_slide_input:  FAIL.\n");
    exit(1);

  }

  return 0;

}