/*

   test_homfly.c : Unit tests for the code in pdcode.c


*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
#endif

#ifdef HAVE_STDIO_H
   #include<stdio.h>
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

#ifdef HAVE_ARGTABLE2_H
  #include<argtable2.h>
#endif

#include<plcTopology.h>

#include<pd_multidx.h>
#include<pd_perm.h>
#include<pd_dihedral.h>

#include<pd_isomorphisms.h>
#include<pd_orientation.h>

int PD_VERBOSE=50;

char *pdcode_to_ccode(pd_code_t *pd);

bool trefoil_ccode_test() {

  printf("-------------------------------------------\n"
	 "testing diagrams based on trefoil knot\n"
	 "-------------------------------------------\n");

  pd_code_t *trefoil_pd = pd_build_torus_knot(2,3);

  /*

    According to the docs for build torus knot, this should produce the pd code

   +-----------------------------------------------+
   |                                               |
   +-- 2  ----\    /---0---\   /-- 4 --\   /---2---+
               \0 /         \1/         \2/
	        \            \           \
	       / \          / \         / \
   +-- 5   ---/   \---3  --/   \---1---/   \---5--+
   |                                              |
   +----------------------------------------------+

   which, using the convention from the documentation for lmpoly

       a
       |
       |
   b<--|-->d (meaning that this can go either way, depending on the orientation of cr)
       |
       V
       c

   So a crossing code representation of a plCurve is a char buffer
   containing lines of the form

   17+2b10c11c31a

   meaning that crossing 17 is a positive crossing

   connected in the a position to the b position of crossing 2,
   connected in the b position to the c position of crossing 10,
   connected in the c position to the c position of crossing 11 and
   connected in the d position to the a position of crossing 31.

   means we should generate the crossing code determined by

   +--------<------------------------<-------------+
   |                                               |
   +------>--a\    /d---->a\   /d---->a\   /d>-----+
               \1 /         \2/         \3/
	        \            \           \
	       / \          / \         / \
   +------>--b/   \c---->-b/   \c>----b/   \c>----+
   |                                              |
   +------<--------------<-------------------<----+

   or

  */

  char correct_ccode[2048] =
    "1+3d3c2b2a\n"
    "2+1d1c3b3a\n"
    "3+2d2c1b1a\n";

  printf("testing +trefoil pd code\n");
  printf("computing ccode...");

  char *ccode = pdcode_to_ccode(trefoil_pd);

  printf("done\n");
  printf("computed ccode:\n\n%s\n",ccode);

  printf("expected ccode:\n\n%s\n",correct_ccode);

  printf("comparing computed to expected ccode...");

  if (strcmp(ccode,correct_ccode)) {

    printf("FAIL\n\n");
    pd_printf("from PD %PD",trefoil_pd);

    return false;

  }

  printf("pass\n\n");
  free(ccode);

  printf("testing -trefoil pd code\n");
  pd_idx_t i;
  for(i=0;i<trefoil_pd->ncross;i++) {

    trefoil_pd->cross[i].sign = PD_NEG_ORIENTATION;

  }

  printf("computing ccode...");
  ccode = pdcode_to_ccode(trefoil_pd);
  printf("done\n");

  /*
    If we switch all of the crossings to negative,
    we get:

   +--------<------------------------<-------------+
   |                                               |
   +------>--d\    /c---->d\   /c---->d\   /c>-----+
               \1 /         \2/         \3/
	        /            /           /
	       / \          / \         / \
   +------>--a/   \b---->-a/   \b>----a/   \b>----+
   |                                              |
   +------<--------------<-------------------<----+

   or
  */

  char minustref_expected[4096] =
    "1-3b2a2d3c\n"
    "2-1b3a3d1c\n"
    "3-2b1a1d2c\n";

  printf("computed ccode:\n\n%s\n",ccode);

  printf("expected ccode:\n\n%s\n",minustref_expected);

  printf("comparing computed to expected ccode...");

  if (strcmp(ccode,minustref_expected)) {

    printf("FAIL\n\n");
    pd_printf("from PD %PD",trefoil_pd);

    return false;

  }

  printf("pass\n");
  free(ccode);

  printf("testing +-+ unknot pd code\n");
  for(i=0;i<trefoil_pd->ncross;i++) {

    trefoil_pd->cross[i].sign = PD_POS_ORIENTATION;

  }

  trefoil_pd->cross[1].sign = PD_NEG_ORIENTATION;

  printf("computing ccode...");
  ccode = pdcode_to_ccode(trefoil_pd);
  printf("done\n");

  /*
    If we switch ONLY THE MIDDLE CROSSING to negative,
    we get:

   +--------<------------------------<-------------+
   |                                               |
   +------>--a\    /d---->d\   /c---->a\   /d>-----+
               \1 /         \2/         \3/
	        \            /           \
	       / \          / \         / \
   +------>--b/   \c---->-a/   \b>----b/   \c>----+
   |                                              |
   +------<--------------<-------------------<----+

   or
  */

  char pmp_expected[4096] =
    "1+3d3c2a2d\n"
    "2-1c3b3a1d\n"
    "3+2c2b1b1a\n";

  printf("computed ccode:\n\n%s\n",ccode);

  printf("expected ccode:\n\n%s\n",pmp_expected);

  printf("comparing computed to expected ccode...");

  if (strcmp(ccode,pmp_expected)) {

    printf("FAIL\n\n");
    pd_printf("from PD %PD",trefoil_pd);

    return false;

  }

  printf("pass\n");
  free(ccode);

  pd_code_free(&trefoil_pd);
  printf("-----------------------------------------------\n"
	 "trefoil-based crossing code generation tests: PASS \n"
	 "-----------------------------------------------\n");

  return true;
}

bool test_all_signs(pd_code_t *pd) {

  pd_multidx_t *sign_iterator = pd_new_multidx(pd->ncross,NULL,orientation_ops);
  unsigned int max = pd_multidx_nvals(sign_iterator);
  unsigned int i;

  for(i=0;i<max;i++,pd_increment_multidx(sign_iterator)) {

    /* Actually set the crossings according to the information in the 
       iterator. */

    for(k=0;k<sign_iterator->nobj;k++) { 

      pd->cross[k].sign = (pd_orientation_t *)(sign_iterator->obj[k])->or;

    }

    pd_printf("\t testing crsigns %MULTIDX...",NULL,sign_iterator);

    if (!pd_ok(pd)) { 

      pd_printf("FAIL. pd %PD not ok after sign assignment.\n",pd);
      return false;

    }

    char *ccode = pdcode_to_ccode(pd);
    free(ccode); 

    printf("pass.\n");

  }
  
  pd_free_multidx(&sign_iterator);
  return true;

} 

bool unknot_generation_tests() {

  printf("-----------------------------------------------\n"
	 "testing diagrams based on unknots \n"
	 "-----------------------------------------------\n");

  printf("generating codes for unknots with 2-10 + crossings...\n"
	 "\t");

  pd_idx_t k;
  for(k=2;k<11;k++) {

    pd_code_t *pd = pd_build_unknot(k);
    char *ccode = pdcode_to_ccode(pd);
    pd_code_free(&pd);
    free(ccode);
    printf("%d ",k);

  }

  printf("testing all crossing signs for 7 crossing diagram...");

  pd_code_t *pd = pd_build_unknot(7);

  printf("testing all crossing signs for 7 crossing diagram...");
  if (test_all_signs(pd)) { 

    printf("pass (generated all ccodes w/o crashing)\n");


  } else {

    printf("FAIL.\n");
    return false;

  }
      
  pd_code_free(&pd);

  printf("testing all crossing signs for 6 crossing diagram...");
  pd = pd_build_unknot(6);

  if (test_all_signs(pd)) { 

    printf("pass (generated all ccodes w/o crashing)\n");

  } else {

    printf("FAIL.\n");
    return false;

  }
      
  pd_code_free(&pd);
  printf("pass (survived)\n");

  printf("-----------------------------------------------\n"
	 "unknot-based crossing code generation tests: PASS \n"
	 "-----------------------------------------------\n");

  return true;
}

bool test_ccode_conversion() {

  printf("tests for conversion of pd_code_t to Millett/Ewing crossing code\n");
  printf("tests function \"char *pdcode_to_ccode(pd_code_t *pd)\",\n"
	 "(not part of exposed API)\n\n");

  if (!trefoil_ccode_test()) { return false; }
  if (!unknot_generation_tests()) { return false; }
  return true;
}

int main() {

  printf("test_homfly (%s)\n",PACKAGE_STRING);
  printf("--------------------------------------------------------\n"
	 "Unit tests for computing HOMFLY polynomial from pdcode. \n"
	 "========================================================\n");

  if (!test_ccode_conversion()) {

    printf("=======================================================\n");
    printf("test_homfly:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================================\n");
    printf("test_homfly:  PASS.\n");
    exit(0);

  }

  return 0;

}
