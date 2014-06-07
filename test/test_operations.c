/*

   test_operations.c : Unit tests for the code in pd_operations.c. This is one of the longest 
   and most extensive test suites we've got, because this is where it's easiest to screw up.


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

#ifdef HAVE_ARGTABLE2_H
  #include<argtable2.h>
#endif

#include<plcTopology.h>

#include<pd_multidx.h>
#include<pd_perm.h>
#include<pd_dihedral.h>

#include<pd_isomorphisms.h>
#include<pd_orientation.h>
#include<polynomials.h>

int PD_VERBOSE=50;

bool trefoil_reverse_test() {
 
  printf("-------------------------------------------\n"
	 "testing diagrams based on torus knot/links\n"
	 "-------------------------------------------\n");

  pd_code_t *trefoil_pd = pd_build_torus_knot(2,3);

  /*
 
    According to the docs for build torus knot, this should produce
    the pd code

   +-----<-----------------------------------<-----+
   |                                               |
   +-- 2 >----\    />--0---\   /-> 4 --\   /-->2---+
               \0 /         \1/         \2/
	        \ +          \ +         \ +
	       / \          / \         / \
   +-- 5 >----/   \---3->--/   \->-1---/   \---5->+
   |                                              |
   +------------<-----------------------<---------+

   Reversing the orientation of component 0 should not change the sign
   of any of the crossings, however, we should renumber to get:

   +----->----------------------------------->-----+ 
   |                                               |
   +-- 4 <----\    /<--0---\   /-- 2<--\   /--<4---+
               \0 /         \1/         \2/
	        \ +          \ +         \ +
	       / \          / \         / \
   +-- 1 <----/   \---3-<--/   \-<-5---/   \---1-<+
   |                                              |
   +------------>----------1------------>---------+

   There are several tests we can do to check that this
   is the right thing:

  */  

  pd_printf("testing +trefoil pd code %PD",trefoil_pd);

  printf("reversing orientation...");
  pd_reorient_component(trefoil_pd,0,PD_NEG_ORIENTATION);
  if (pd_ok(trefoil_pd)) { 
    printf("pass (pd_ok after reorient)\n");
  } else {
    pd_printf("FAIL. pd %PD is not ok after reorient of component 0.\n",trefoil_pd);
    return false;
  }

  printf("checking crossing signs...");
  pd_idx_t cr;
  for(cr = 0; cr < 3;cr ++) { 
    if (trefoil_pd->cross[cr].sign != PD_POS_ORIENTATION) { 
      pd_printf("FAIL\n"
		"crossing %d is not in PD_POS_ORIENTATION in \n"
		"%PD",trefoil_pd,cr);
      return false;
    }
  }
  printf("pass (+++)\n");

  printf("checking for 3 bigon and 2 trigon faces...");
  pd_idx_t bigons[6],bigons_used = 0;
  pd_idx_t trigons[6],trigons_used = 0;
  pd_idx_t face;

  for(face=0;face<trefoil_pd->nfaces;face++) {

    if (trefoil_pd->face[face].nedges == 2) { 
      
      bigons[bigons_used++] = face;

    } else if (trefoil_pd->face[face].nedges == 3) {

      trigons[trigons_used++] = face;

    } else { 

      pd_printf("FAIL. Face %FACE is not a bigon or a trigon.",trefoil_pd,face);
      return false;

    }

  }

  if (bigons_used != 3 || trigons_used != 2) { 

    pd_printf("FAIL. Found %d bigons and %d trigons in pd %PD\n",trefoil_pd,bigons_used,trigons_used);
    return false;

  } 

  printf("pass (found 3 bigons, 3 trigons)\n");

  printf("checking that 0, 2, 4 occur positively in bigons...");

  pd_idx_t posface,posface_pos, negface, negface_pos;
  pd_idx_t elist[3] = {0,2,4};
  pd_idx_t edge;

  for(edge=0;edge<3;edge++) { 

    pd_face_and_pos(trefoil_pd,elist[edge],&posface,&posface_pos,&negface,&negface_pos);
    if (trefoil_pd->face[posface].nedges != 2) {

      pd_printf("FAIL. Edge %d occurs positively in non-bigon %FACE in pd %PD.\n",trefoil_pd,elist[edge],posface);
      return false;

    }

  }
  
  printf("pass\n");
  
  printf("checking that 0, 2, 4 occur negatively in trigons..."); 

  for(edge=0;edge<3;edge++) { 

    pd_face_and_pos(trefoil_pd,elist[edge],&posface,&posface_pos,&negface,&negface_pos);
    if (trefoil_pd->face[negface].nedges != 3) {

      pd_printf("FAIL. Edge %d occurs negatively in non-trigon %FACE in pd %PD.\n",trefoil_pd,elist[edge],negface);
      return false;

    }

  }
  
  printf("pass\n");

  printf("checking that 1, 3, 5 occur positively in trigons...");

  elist[0] = 1; elist[1] = 3; elist[2] = 5;

  for(edge=0;edge<3;edge++) { 

    pd_face_and_pos(trefoil_pd,elist[edge],&posface,&posface_pos,&negface,&negface_pos);
    if (trefoil_pd->face[posface].nedges != 3) {

      pd_printf("FAIL. Edge %d occurs positively in non-trigon %FACE in pd %PD.\n",trefoil_pd,elist[edge],posface);
      return false;

    }

  }
  
  printf("pass\n");
  
  printf("checking that 1, 3, 5 occur negatively in bigons..."); 

  for(edge=0;edge<3;edge++) { 

    pd_face_and_pos(trefoil_pd,elist[edge],&posface,&posface_pos,&negface,&negface_pos);
    if (trefoil_pd->face[negface].nedges != 2) {

      pd_printf("FAIL. Edge %d occurs negatively in non-bigon %FACE in pd %PD.\n",trefoil_pd,elist[edge],negface);
      return false;

    }

  }
  
  printf("pass\n");

  printf("+trefoil tests...pass.\n\n");

  return true;
}

bool test_all_component_reversals(pd_code_t *pd, int *tests_run) {

  pd_multidx_t *sign_iterator = pd_new_multidx(pd->ncomps,NULL,orientation_ops);
  unsigned int max = pd_multidx_nvals(sign_iterator);
  unsigned int i;

  printf("\t testing %d component reversals for %d component link...",max,pd->ncomps);

  for(i=0;i<max;i++,pd_increment_multidx(sign_iterator)) {

    /* Actually set the component orientations according to the information in the 
       iterator. */

    pd_code_t *working_copy = pd_copy(pd);

    pd_idx_t k;
    for(k=0;k<sign_iterator->nobj;k++) { 

      pd_reorient_component(working_copy,k,((pd_orientation_t *)(sign_iterator->obj[k]))->or);

    }

    if (!pd_ok(working_copy)) { 
      
      pd_printf("FAIL. pd %PD not ok \n"
		"after component orientation reversal %MULTIDX.\n",pd,sign_iterator);
      pd_printf("yielded pd %PD\n",working_copy);
      return false;
      
    }
    
    pd_code_free(&working_copy);
      
  }
    
  printf("pass.\n");
  pd_free_multidx(&sign_iterator);
  *tests_run = max;

  return true;

} 

/* We now write a test where we attempt to load a knot table and compute
   HOMFLYs for it. */

bool rolfsentabletest() 
{
 
  printf("\n--------------------------------------------------\n"
	 "Rolfsen/Thistlethwaite Link Table test\n"
	 "--------------------------------------------------\n"); 
  
  printf("Trying to determine srcdir from environment...");
  char *srcdir = getenv("srcdir");
  if (srcdir == NULL) { 

    printf("fail. (assuming data files local)\n");
    srcdir = calloc(4,sizeof(char));
    sprintf(srcdir,".");

  } else {

    printf("pass (srcdir = %s)\n",srcdir);

  }

  char *rolfsentable = calloc(4096,sizeof(char));

  sprintf(rolfsentable,"%s/rolfsentable.txt",srcdir);
  printf("Opening data file %s...",rolfsentable);

  FILE *infile;

  infile = fopen(rolfsentable,"r");
  
  if (infile != NULL) {

    printf("pass\n");

  } else { 

    printf("fail\n");
    return false;

  }

  printf("Loading Mathematica format pd codes...");

  int pd_codes_expected;
  if (!fscanf(infile,"pdcodes %d \n",&pd_codes_expected) == 1) { 

    printf("fail. (Couldn't read # of codes in file)");
    return false;

  } 

  pd_code_t **pdbuf;
  pdbuf = calloc(pd_codes_expected,sizeof(pd_code_t *));
  assert(pdbuf != NULL);

  int loaded;
  
  for(loaded = 0;loaded < pd_codes_expected && !feof(infile); loaded++) {

    pdbuf[loaded] = pd_read_KnotTheory(infile);

    if (pdbuf[loaded] == NULL) {

      printf("fail (on pd code %d)\n",loaded);
      return false;

    }

  }

  if (loaded != pd_codes_expected) { 

    printf("fail. (expected %d pd codes, got %d)\n",
	   pd_codes_expected,loaded);
    return false;

  }

  printf("pass (%d pd codes loaded, %d expected).\n",loaded,pd_codes_expected);
  fclose(infile);


  printf("trying all component reversals for %d pd codes...\n\n",loaded);
  clock_t start, end;
  double cpu_time_used;
  long int total_tests_run = 0;
  int tests_run;

  start = clock();
  int computed;

  for(computed = 0; computed < loaded; computed++) { 

    test_all_component_reversals(pdbuf[computed],&tests_run);
    total_tests_run += tests_run;
    pd_code_free(&(pdbuf[computed]));

  }

  end = clock();
  cpu_time_used =  ((double)(end - start))/CLOCKS_PER_SEC;

  printf("\n\n...pass (%ld component reversal tests computed in %2.2g sec).\n\n",total_tests_run,cpu_time_used);
  free(pdbuf);
  free(rolfsentable);

  char *thistlethwaitetable = calloc(4096,sizeof(char));

  sprintf(thistlethwaitetable,"%s/thistlethwaitetable.txt",srcdir);
  printf("Opening data file %s...",thistlethwaitetable);

  infile = fopen(thistlethwaitetable,"r");
  
  if (infile != NULL) {

    printf("pass\n");

  } else { 

    printf("fail\n");
    return false;

  }

  printf("Loading Mathematica format pd codes from file...");

  if (!fscanf(infile,"pdcodes %d \n",&pd_codes_expected) == 1) { 

    printf("fail. (Couldn't read # of codes in file)");
    return false;

  } 

  pdbuf = calloc(pd_codes_expected,sizeof(pd_code_t *));
  assert(pdbuf != NULL);

  for(loaded = 0;loaded < pd_codes_expected && !feof(infile); loaded++) {

    pdbuf[loaded] = pd_read_KnotTheory(infile);

    if (pdbuf[loaded] == NULL) {

      printf("fail (on pd code %d)\n",loaded);
      return false;

    }

  }

  if (loaded != pd_codes_expected) { 

    printf("fail. (expected %d pd codes, got %d)\n",
	   pd_codes_expected,loaded);
    return false;

  }

  printf("pass (%d pd codes loaded, %d expected).\n",loaded,pd_codes_expected);
  fclose(infile);

  printf("trying all component reversals for %d pd codes...\n\n",loaded);
  start = clock();
  total_tests_run = 0;

  for(computed = 0; computed < loaded; computed++) { 

    test_all_component_reversals(pdbuf[computed],&tests_run);
    total_tests_run += tests_run;
    pd_code_free(&(pdbuf[computed]));

  }

  end = clock();
  cpu_time_used =  ((double)(end - start))/CLOCKS_PER_SEC;

  printf("\n ... pass (%ld component reversal tests computed in %2.2g sec).\n\n",total_tests_run,cpu_time_used);

  free(pdbuf);
  free(thistlethwaitetable);
  
  printf("\n"
	 "--------------------------------------------------\n"
	 "Rolfsen/Thistlethwaite Table test: pass             \n"
	 "--------------------------------------------------\n\n"); 

  return true;

}
  

    

int main() {

  printf("test_operations (%s)\n",PACKAGE_STRING);
  printf("--------------------------------------------------------\n"
	 "Unit tests for pdcode operation primitives. \n"
	 "========================================================\n");

  if (!trefoil_reverse_test() || !rolfsentabletest()) {

    printf("=======================================================\n");
    printf("test_operations:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================================\n");
    printf("test_operations:  PASS.\n");
    exit(0);

  }

  return 0;

}