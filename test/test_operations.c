/*

   test_operations.c : Unit tests for the code in pd_operations.c. This is one of the longest 
   and most extensive test suites we've got, because this is where it's easiest to screw up.


*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
#endif

#include<stdio.h>
#include<assert.h>
#include<string.h>

#ifdef HAVE_STDINT_H
   #include<stdint.h>
#endif

#include<stdlib.h>

#ifdef HAVE_STDBOOL_H
   #include<stdbool.h>
#endif

#include<time.h>


#include<argtable2.h>
#include<plcTopology.h>

#include<pd_multidx.h>
#include<pd_perm.h>
#include<pd_dihedral.h>

#include<pd_isomorphisms.h>
#include<pd_orientation.h>
#include<polynomials.h>



bool isomorphism_matches_sign_iterator(pd_iso_t *iso,pd_code_t *pdA,pd_multidx_t *sign_iterator)
 
/* Runs a number of checks to make sure that the isomorphism matches the sign iterator given. 
   Of course, the most important check is that we actual HAVE an isomorphism: checking signs
   is secondary. 
*/
{
  /* Check 1: The component permutation is the identity. */
  
  if (!pd_perm_is_e(iso->compperm)) { 
    return false;
  }

  /* Check 2: The edgemap reverses or preserves all the orientations in each component. */

  pd_idx_t i,j;
  pd_or_t *signs_actual = calloc(pdA->ncomps,sizeof(pd_or_t));

  for(i=0;i<pdA->ncomps;i++) { 
    
    pd_or_t start_or = iso->edgemap->orient[ pdA->comp[i].edge[0]];
    for(j=1;j<pdA->comp[i].nedges;j++) { 

      if (iso->edgemap->orient[ pdA->comp[i].edge[j]] != start_or) { 

	free(signs_actual);
	return false;

      }

    }

    signs_actual[i] = start_or;

  }

  /* Now we check if the sign iterator matches. */

  for(i=0;i<sign_iterator->nobj;i++) { 
    
    if (signs_actual[i] != ((pd_orientation_t *)(sign_iterator->obj[i]))->orient) {

      free(signs_actual);
      return false;

    }

  }

  return true;
  
}

bool test_all_component_reversals(pd_code_t *pd, int *tests_run) {

  pd_multidx_t *sign_iterator = pd_new_multidx(pd->ncomps,NULL,orientation_ops);
  unsigned int max = pd_multidx_nvals(sign_iterator);
  unsigned int i;

  printf("\t testing %d component reversals for %d component link...\n",max,pd->ncomps);

  for(i=0;i<max;i++,pd_increment_multidx(sign_iterator)) {

    /* Actually set the component orientations according to the information in the 
       iterator. */

    pd_code_t *working_copy = pd_copy(pd);

    pd_idx_t k;
    for(k=0;k<sign_iterator->nobj;k++) { 

      pd_reorient_component(working_copy,k,((pd_orientation_t *)(sign_iterator->obj[k]))->orient);

    }

    if (!pd_ok(working_copy)) { 
      
      pd_printf("FAIL. pd %PD not ok \n"
		"after component orientation reversal %MULTIDX.\n",pd,sign_iterator);
      pd_printf("yielded pd %PD\n",working_copy);
      return false;
      
    }

    pd_printf("\t\t checking isos (all +) <-> %MULTIDX...",NULL,sign_iterator);

    pd_iso_t **isomorphisms;
    unsigned int nisos;
    isomorphisms = pd_build_isos(pd,working_copy,&nisos);

    if (nisos > 0) { 

      printf("done (found %d)\n",nisos);
    
    } else {
      
      pd_printf("FAIL. original pd %PD is not isomorphic to\n",pd);
      pd_printf("orientation-reversed pd %PD\n",working_copy);
      return false;
    
    }

    pd_printf("\t\t searching isos for sign matching %MULTIDX ...",NULL,sign_iterator);

    pd_idx_t i;
    bool found_ok_iso = false;
    pd_idx_t found_idx;
     
    for(i=0;i<nisos && !found_ok_iso;i++) { 
    
      if (isomorphism_matches_sign_iterator(isomorphisms[i],pd,sign_iterator)) {
      
	found_ok_iso = true;
	found_idx = i;

      }

    }
  
    if (found_ok_iso) { 

      printf("done (# %d matches)\n\n",found_idx);
    
    } else {

      printf("failed. (no isomorphism matches component signs)\n");
      return false;

    }

    pd_free_isos(&nisos,&isomorphisms);  
    pd_code_free(&working_copy);
      
  }
    
  printf("pass.\n");
  pd_free_multidx(&sign_iterator);
  *tests_run = max;

  return true;

} 


bool trefoil_reverse_test() {
 
  printf("-------------------------------------------\n"
	 "testing diagrams based on torus knot/links\n"
	 "-------------------------------------------\n");

  pd_code_t *trefoil_pd = pd_build_torus_knot(2,3);
  pd_printf("testing +trefoil pd code %PD",trefoil_pd);
  int tests_run;

  if (test_all_component_reversals(trefoil_pd,&tests_run)) {

    printf("\tdone (%d tests run)\n",tests_run);
    pd_code_free(&trefoil_pd);
    return true;

  } else {

    printf("\tfail (%d tests run)\n",tests_run);
    pd_code_free(&trefoil_pd);
    return false;

  }
    
  printf("+trefoil tests...pass.\n\n");

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
  if (!(fscanf(infile,"pdcodes %d \n",&pd_codes_expected) == 1)) { 

    printf("fail. (Couldn't read # of codes in file)");
    return false;

  } 

  pd_code_t **pdbuf;
  pdbuf = calloc(pd_codes_expected,sizeof(pd_code_t *));
  assert(pdbuf != NULL);

  int loaded;
  if (pd_codes_expected > 100) { pd_codes_expected = 100; }
  
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

  if (!(fscanf(infile,"pdcodes %d \n",&pd_codes_expected) == 1)) { 

    printf("fail. (Couldn't read # of codes in file)");
    return false;

  } 

  if (pd_codes_expected > 100) { pd_codes_expected = 100; }

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


bool connect_sum_tester(pd_code_t *pdA, pd_idx_t edgeA,
			pd_code_t *pdB, pd_idx_t edgeB,
			char *nameA,char *nameB,
			bool verbose)
{
  if (verbose) {
    
    printf("connect sum of %s and %s...",nameA,nameB);

  }
  
  pd_code_t *cs = pd_connect_sum(pdA,edgeA,pdB,edgeB);
 
  if (!pd_ok(cs)) {

    if (verbose) {

      printf("FAIL (connect sum didn't pass pd_ok)\n");

    }
    
    return false;

  }

  if (cs->ncross != pdA->ncross + pdB->ncross) {

    if (verbose) {

      printf("FAIL (connect sum has %d crossings != %d + %d)\n",
	     cs->ncross,pdA->ncross,pdB->ncross);

    }
    
    return false;

  }

  pd_idx_t i;

  for(i=0;i<cs->ncross;i++) {

    if (cs->cross[i].sign == PD_UNSET_ORIENTATION &&
	pdA->cross[0].sign != PD_UNSET_ORIENTATION &&
	pdB->cross[0].sign != PD_UNSET_ORIENTATION) {

      if (verbose) {

	printf("FAIL (connect sum lost crossing info)\n");

      }  

      return false;

    }

  }

  pd_code_free(&cs);

  if (verbose) {

    printf("pass\n");

  }

  return true;
  
}
  

bool connect_sum_tests()

{

  printf("\n--------------------------------------------------\n"
	 "Connect Sum tests\n"
	 "--------------------------------------------------\n"); 

  printf("generating example pdcodes (2,5) torus knot, 3-chain, 4-chain...");
  
  pd_code_t *torusKnot = pd_build_torus_knot(2,5);
  pd_code_t *simpleChainA = pd_build_simple_chain(3);
  pd_code_t *simpleChainB = pd_build_simple_chain(4);

  if (!pd_ok(torusKnot) || !pd_ok(simpleChainA) || !pd_ok(simpleChainB)) {

    printf("FAIL (didn't pass pd_ok)\n");
    return false;

  }

  printf("pass (pd_ok for all three)\n");

  if (!connect_sum_tester(torusKnot,3,simpleChainA,5,"(2,5) torus knot","3-chain",true)) {

    return false;

  }

  if (!connect_sum_tester(simpleChainB,9,simpleChainA,5,"4-chain","3-chain",true)) {

    return false;

  }

  printf("testing all connect sums of 3-chain and 4-chain...");

  pd_idx_t i,j;

  for(i=0;i<simpleChainA->nedges;i++) {

    for(j=0;j<simpleChainB->nedges;j++) {

      if (!connect_sum_tester(simpleChainA,i,simpleChainB,j,"3-chain","4-chain",false)) {

	printf("FAIL (on connect sum of 3-chain at %d and 4-chain at %d):\n",
	       i,j);

	connect_sum_tester(simpleChainA,i,simpleChainB,j,"3-chain","4-chain",true);

	return false;

      }

    }

  }
	
  printf("pass (all passed)\n");

  printf("testing all connect sums of torus knot and 4-chain...");

  for(i=0;i<torusKnot->nedges;i++) {

    for(j=0;j<simpleChainB->nedges;j++) {

      if (!connect_sum_tester(torusKnot,i,simpleChainB,j,"(2,5) torus knot","4-chain",false)) {

	printf("FAIL (on connect sum of (2,5) torus knot at %d and 4-chain at %d):\n",
	       i,j);

	connect_sum_tester(torusKnot,i,simpleChainB,j,"(2,5) torus knot","4-chain",true);

	return false;

      }

    }

  }
	
  printf("pass (all passed)\n");

  printf("housekeeping...");

  pd_code_free(&torusKnot);
  pd_code_free(&simpleChainA);
  pd_code_free(&simpleChainB);

  printf("done\n");
  
  return true;    
  
}
    

int main() {

  printf("test_operations (%s)\n",PACKAGE_STRING);
  printf("--------------------------------------------------------\n"
	 "Unit tests for pdcode operation primitives. \n"
	 "========================================================\n");

  if (!connect_sum_tests() || !trefoil_reverse_test() || !rolfsentabletest()) {

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
