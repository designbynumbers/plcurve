/*

   test_crossingops.c : Unit tests for the code in
   pd_crossingmoves.c. This is going to grow to be a relatively
   elaborate test suite because getting this stuff exactly right is
   one of the most important parts of the whole project.

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

int PD_VERBOSE=50;

void pd_compacting_copy(void *source, size_t obj_size, size_t nobj, 
			pd_idx_t ndeletions, pd_idx_t *deletions,
			void **target_var,
			pd_idx_t **target_idx,
			pd_idx_t *ntarget);  /* Returns the size of the target buffer*/

bool compacting_copy_int_test(size_t nobj, pd_idx_t ndeletions, pd_idx_t *deletions,
			      pd_idx_t nexpected, int *expected)

/* This runs a test on the compacting_copy API assuming that the objects are 
   all of type int, they are numbered 0..nobj-1, and that we give deletions, 
   and a buffer of expected results.  Just for our information... 

  void pd_compacting_copy(void *source, size_t obj_size, size_t nobj, 
			pd_idx_t ndeletions, pd_idx_t *deletions,
			void **target_var,
			pd_idx_t **target_idx,
			pd_idx_t *ntarget)  ( Returns the size of the target buffer.)
			
*/

{ 
  size_t obj_size = sizeof(int);
  void *source = calloc(nobj,obj_size);
  int  *intsource = (int *)(source);
  assert(source != NULL);
  pd_idx_t i;
  for(i=0;i<nobj;i++) { intsource[i] = i; }
 
  void *target_var;
  pd_idx_t ntarget;
  pd_idx_t *target_idx;
  
  for(i=0;i<10;i++) { intsource[i] = i; }

  printf("testing deleting ");
  for(i=0;i<ndeletions;i++) { printf("%d ",deletions[i]); }
  printf("from [0..%d]...\n\n",(int)(nobj));

  pd_compacting_copy(source,obj_size,nobj,
		     ndeletions,deletions,
		     &target_var,
		     &target_idx,&ntarget);

  printf("\t checking number of remaining elements == expected (%d)...",nexpected);
  
  if (ntarget != nexpected) { 

    printf("fail. (ntarget == %d != %d)\n",ntarget,nexpected);
    return false;
    
  } else { 

    printf("pass. (%d == %d)\n",ntarget,nexpected);

  }

  printf("\t checking compacted version == ( ");
  for(i=0;i<nexpected;i++) { printf("%d ",expected[i]); }
  printf(") ... ");

  int *int_targetvar = (int *)(target_var);

  for(i=0;i<nexpected;i++) { 

    if (int_targetvar[i] != expected[i]) { 

      printf("fail. (element %d == %d != %d)\n",i,int_targetvar[i],expected[i]);
      return false;

    }

  }

  printf("pass (compacted array ok)\n");
  printf("\t checking target idx ...");
  
  bool *used = calloc(nobj,sizeof(bool));
  for(i=0;i<nobj;i++) { used[i] = false; }
  
  for(i=0;i<nobj;i++) { 

    if ((target_idx[i] > ntarget) && target_idx[i] != PD_UNSET_IDX) { 
    
      printf("fail. (illegal index %d in map\n",target_idx[i]);
      int j;
      for(j=0;j<nobj;j++) { 
	if (target_idx[j] != PD_UNSET_IDX) { 
	  printf("\t\t%d -> %d \n",j,target_idx[j]); 
	} else {
	  printf("\t\t%d -> PD_UNSET_IDX \n",j);
	}
      }	      
      printf(")\n\n");
    
    }

    if (target_idx[i] != PD_UNSET_IDX) {

      if (used[target_idx[i]]) { 

	printf("fail. (two indices are assigned to target %d in map\n",target_idx[i]);
	int j;
	for(j=0;j<nobj;j++) { 
	  if (target_idx[j] != PD_UNSET_IDX) { 
	    printf("\t\t%d -> %d \n",j,target_idx[j]); 
	  } else {
	    printf("\t\t%d -> PD_UNSET_IDX \n",j);
	  }
	}	  
	printf(")\n\n");
	return false;
    
      }
  
      used[target_idx[i]] = true; 

    }

  }

  /* Now check that all the indices in the target array have been hit. */

  for(i=0;i<ntarget;i++) { 

    if (!used[i]) { 

      printf("fail. (missed target index %d in map",i);
      int j;
      for(j=0;j<nobj;j++) { 
	if (target_idx[j] != PD_UNSET_IDX) { 
	  printf("\t\t%d -> %d \n",j,target_idx[j]); 
	} else {
	  printf("\t\t%d -> PD_UNSET_IDX \n",j);
	}
      }	  
      printf(")\n\n");
      return false;

    }

  }

  printf("pass (target idx looks ok)\n");
  free(used);

  printf("\t housecleaning...");
  free(target_var);
  free(target_idx);
  free(source);
  printf("done\n\n");

  return true;

}
         
bool compacting_copy_tests() 
/* A sequence of tests for the pd_compacting_copy API designed to 
   stress most of the obvious cases (a lot of deletions, very few
   deletions, adjacent deletions, out-of-order deletions, repeated
   deletions and so forth). */
  
/* We're going to use the test harness:

   bool compacting_copy_int_test(size_t nobj, pd_idx_t ndeletions, pd_idx_t *deletions,
                                 pd_idx_t nexpected, int *expected_result) 
*/
{
  printf("------------------------------------------\n"
	 "Compacting-copy test suite \n"
	 "------------------------------------------\n");

  pd_idx_t deletions[15];
  int      expected_result[10];

  /* Single deletion */

  deletions[0] = 3;
  expected_result[0] = 0;
  expected_result[1] = 1;
  expected_result[2] = 2;
  expected_result[3] = 4;
  expected_result[4] = 5;
  expected_result[5] = 6;
  expected_result[6] = 7;
  expected_result[7] = 8;
  expected_result[8] = 9; 

  if (!compacting_copy_int_test(10,1,deletions,9,expected_result)) { return false; }

  /* Consecutive deletions. */

  deletions[0] = 3;
  deletions[1] = 4;
  deletions[2] = 5;

  expected_result[0] = 0;
  expected_result[1] = 1;
  expected_result[2] = 2;
  expected_result[3] = 6;
  expected_result[4] = 7;
  expected_result[5] = 8;
  expected_result[6] = 9;

  if (!compacting_copy_int_test(10,3,deletions,7,expected_result)) { return false; }

  /* Non consecutive deletions. */

  deletions[0] = 3;
  deletions[1] = 7;

  expected_result[0] = 0;
  expected_result[1] = 1;
  expected_result[2] = 2;
  expected_result[3] = 4;
  expected_result[4] = 5;
  expected_result[5] = 6;
  expected_result[6] = 8;
  expected_result[7] = 9;

  if (!compacting_copy_int_test(10,2,deletions,8,expected_result)) { return false; }

  /* Out of order deletions */

  deletions[0] = 7;
  deletions[1] = 3;

  expected_result[0] = 0;
  expected_result[1] = 1;
  expected_result[2] = 2;
  expected_result[3] = 4;
  expected_result[4] = 5;
  expected_result[5] = 6;
  expected_result[6] = 8;
  expected_result[7] = 9;

  if (!compacting_copy_int_test(10,2,deletions,8,expected_result)) { return false; }

  /* Repeated deletions */

  deletions[0] = 7;
  deletions[1] = 3;
  deletions[2] = 7;

  expected_result[0] = 0;
  expected_result[1] = 1;
  expected_result[2] = 2;
  expected_result[3] = 4;
  expected_result[4] = 5;
  expected_result[5] = 6;
  expected_result[6] = 8;
  expected_result[7] = 9;

  if (!compacting_copy_int_test(10,3,deletions,8,expected_result)) { return false; }

  /* Boundary index deletions */

  deletions[0] = 0;
  deletions[1] = 9;

  expected_result[0] = 1;
  expected_result[1] = 2;
  expected_result[2] = 3;
  expected_result[3] = 4;
  expected_result[4] = 5;
  expected_result[5] = 6;
  expected_result[6] = 7;
  expected_result[7] = 8;

  if (!compacting_copy_int_test(10,2,deletions,8,expected_result)) { return false; }

  /* No deletions at ALL. */

  deletions[0] = 0;

  expected_result[0] = 0;
  expected_result[1] = 1;
  expected_result[2] = 2;
  expected_result[3] = 3;
  expected_result[4] = 4;
  expected_result[5] = 5;
  expected_result[6] = 6;
  expected_result[7] = 7;
  expected_result[8] = 8;
  expected_result[9] = 9;

  if (!compacting_copy_int_test(10,0,NULL,10,expected_result)) { return false; } 

  printf("------------------------------------------\n"
	 "Compacting-copy test suite: PASS \n"
	 "------------------------------------------\n");

  return true;

}
  

bool find_simple_twist(pd_code_t *pd,pd_idx_t *cr) 
/* Finds a one-edge face (if it exists) and returns the corresponding crossing. */
{
  pd_idx_t i;

  for(i=0;i<pd->nfaces;i++) { 

    if (pd->face[i].nedges == 1) { 
      
      pd_idx_t edge = pd->face[i].edge[0];
      *cr = pd->edge[edge].head;
      return true;

    }

  }
  
  return false;
}

bool unravel_unknot_test(pd_idx_t n)
/* Attempts to do R1 loop deletions to unwind the "plectoneme" from the ends. */
{

  printf("------------------------------------------\n"
	 "Unravelling simple %2d-twist with R1 moves\n"
	 "------------------------------------------\n",
	 n);

  printf("building %2d-twist unknot...",n);
  pd_idx_t crossings_left;
  pd_code_t *pd = pd_build_unknot(n);
  
  if (pd_ok(pd)) { 

    printf("done.\n");

  } else { 

    printf("fail (couldn't generate twist).\n");
    return false;

  }

  printf("running successive untwists...\n\n");

  for(crossings_left=n;crossings_left > 2;crossings_left--) { 

    pd_idx_t target_cr;

    printf("\tlocating target crossing...");
    
    if (!find_simple_twist(pd,&target_cr)) { 

      pd_printf("fail (couldn't find R1 location in %PD)",pd);
      return false;

    } else {

      pd_printf("done. (will untwist at %CROSS).\n",pd,target_cr);
      
    }

    printf("\tuntwisting...");
    pd_code_t *untwisted = pd_R1_loopdeletion(pd,target_cr);
    
    if (!pd_ok(untwisted)) { 
      
      printf("fail. (resulting pd not ok)\n");
      return false;
      
    }

    if (untwisted->ncross != crossings_left-1) { 

      pd_printf("fail. (untwisted pd %PD has %d crossings, not %d)\n",pd,
		untwisted->ncross,crossings_left-1);
      return false;

    }

    pd_code_t *untwisted_compare = pd_build_unknot(crossings_left-1);
    
    if (!pd_isomorphic(untwisted,untwisted_compare)) { 

      printf("fail (not isomorphic to %d twist unknot)\n",crossings_left-1);
      return false;

    }

    pd_code_free(&untwisted_compare);
    pd_code_free(&pd);
    pd = untwisted;

    printf("\tpass (isomorphic to %d crossing twist)\n\n",crossings_left-1);

  }

  printf("------------------------------------------\n"
	 "Unravelling simple %2d-twist: pass         \n"
	 "------------------------------------------\n",
	 n);

  pd_code_free(&pd);

  return true;
    
}

bool r2_twist_tests(pd_idx_t n) {

  printf("-----------------------------------------\n"
	 "r2 bigon elimination tests on %2d-twist  \n"
	 "-----------------------------------------\n",n);

  printf("generating %2d twist...",n);
  pd_code_t *pd = pd_build_unknot(n);
  if (!pd_ok(pd)) { 
    printf("fail (couldn't generate unknot)\n");
    return false;
  } else {
    printf("pass\n");
  }

  printf("resetting crossing signs to +-+-+- along component 0...");
  pd_idx_t i;
  for(i=0;i<pd->ncross;i++) { 
    pd_idx_t cr = pd->edge[pd->comp[0].edge[i]].head;
    pd->cross[cr].sign = (i%2) ? PD_POS_ORIENTATION : PD_NEG_ORIENTATION;
  }
  if (!pd_ok(pd)) { 
    printf("fail (couldn't reset crossing signs)\n");
    return false;
  } else {
    printf("pass\n");
  }

  for(i=0;i<n-1;i++) { 

    pd_idx_t cr[2] = {pd->edge[pd->comp[0].edge[i]].head,
		      pd->edge[pd->comp[0].edge[i+1]].head};

    pd_code_t **outpd;
    pd_idx_t noutpd;
  
    pd_printf("\tperforming bigon elimination at %CROSS and %CROSS...",pd,cr[0],cr[1]);
    pd_R2_bigon_elimination(pd,cr,&noutpd,&outpd);
  
    if (noutpd != 1) { 
      printf("fail (constructed %d != 1 output pd codes)\n",noutpd);
      return false;
    } 

    if (!pd_ok(outpd[0])) { 
      pd_printf("fail (returned pd %PD does not pass pd_ok)\n",outpd[0]);
      return false;
    }

    printf("pass (returned ok pd code)\n");
    
    printf("\tchecking for isomorphism with %d twist knot...",n-2);
    pd_code_t *check_pd = pd_build_unknot(n-2);

    pd_idx_t j;
    
    for(j=0;j<(outpd[0])->ncross;j++) { 
      (outpd[0])->cross[j].sign = PD_POS_ORIENTATION; 
    }
    
    if (!pd_isomorphic((outpd[0]),check_pd)) { 
      printf("fail (not isomorphic)\n");
      return false;
    } 
    printf("pass\n");

    printf("\thousecleaning...");
  
    pd_code_free(&(outpd[0]));
    pd_code_free(&check_pd);
    free(outpd);
    
    printf("done\n\n");

  }

  pd_code_free(&pd);

  printf("----------------------------------------------\n"
	 "r2 bigon elimination tests on %2d-twist:PASS  \n"
	 "----------------------------------------------\n",n);

  return true;

}

bool is_0_crossing_diagram(pd_code_t *pd)
{

  if (pd->MAXVERTS == 0) { return false; }
  if (pd->MAXEDGES == 0) { return false; }
  if (pd->MAXFACES == 0) { return false; }
  if (pd->MAXCOMPONENTS == 0) { return false; }

  if (pd->ncross != 0) { return false; }
  if (pd->nedges != 1) { return false; }
  if (pd->nfaces != 2) { return false; }
  if (pd->ncomps != 1) { return false; }

  if (pd->cross == NULL) { return false; }
  
  if (pd->edge == NULL) { return false; }
  if (pd->edge[0].head != PD_UNSET_IDX) { return false; }
  if (pd->edge[0].headpos != PD_UNSET_POS) { return false; }
  if (pd->edge[0].tail != PD_UNSET_IDX) { return false; }
  if (pd->edge[0].tailpos != PD_UNSET_POS) { return false; }
  
  if (pd->face == NULL) { return false; }
  if (pd->face[0].nedges != 1) { return false; }
  if (pd->face[0].edge[0] != 0) { return false; }
  if (pd->face[0].or[0] != PD_POS_ORIENTATION) { return false; }
  if (pd->face[1].nedges != 1) { return false; }
  if (pd->face[1].edge[0] != 0) { return false; }
  if (pd->face[1].or[0] != PD_NEG_ORIENTATION) { return false; }
  
  if (pd->comp == NULL) { return false; }
  if (pd->comp[0].nedges != 1) { return false; }
  if (pd->comp[0].edge[0] != 0) { return false; }

  return true;

}

bool r2_tiny_tests() { 

  /* This tests the two cases where we pull apart linked rings or unfold a circle. */


  printf("--------------------------------------------\n"
	 "r2 bigon elimination tiny tests             \n"
	 "--------------------------------------------\n");

  printf("constructing 2-twist unknot...");
  pd_code_t *pd = pd_build_unknot(2);
  pd->cross[1].sign = PD_NEG_ORIENTATION;
  if (!pd_ok(pd)) { 
    printf("fail (returned pd not ok)\n");
    return false;
  } else {
    printf("done\n");
  }

  printf("applying r2 bigon elimination...");
  pd_idx_t noutpd;
  pd_code_t **outpd;
  pd_idx_t cr[2] = {0,1};

  pd_R2_bigon_elimination(pd,cr,&noutpd,&outpd);
  
  if (noutpd != 1) { 

    printf("fail (number of loops = %d != 1)\n",noutpd);
    return false;

  }

  if (!pd_ok(outpd[0])) { 

    pd_printf("fail (output pd %PD does not pass pd_ok)\n",outpd[0]);
    return false;

  }

  if (!is_0_crossing_diagram(outpd[0])) { 

    pd_printf("fail (output pd %PD doesn't match 0_crossing_diagram)\n",outpd[0]);
    return false;

  }

  if (outpd[0]->comp[0].tag != pd->comp[0].tag) { 

    printf("fail (tag of output %c != tag of input %c)\n",outpd[0]->comp[0].tag,pd->comp[0].tag);
    return false;

  }

  printf("pass (output matches 0 crossing loop)\n");
  pd_code_free(&outpd[0]);
  pd_code_free(&pd);
  free(outpd);

  printf("generating 2-link simple chain...");
  pd = pd_build_simple_chain(2);
  pd->cross[1].sign = PD_NEG_ORIENTATION;
  if (!pd_ok(pd)) { 
    printf("fail (returned pd not ok)\n");
    return false;
  } else {
    printf("done\n");
  }

  printf("applying r2 bigon elimination...");
  pd_R2_bigon_elimination(pd,cr,&noutpd,&outpd);
  
  if (noutpd != 2) { 

    printf("fail (number of loops = %d != 2)\n",noutpd);
    return false;

  }

  if (!pd_ok(outpd[0])) { 

    pd_printf("fail (output pd 0 %PD does not pass pd_ok)\n",outpd[0]);
    return false;

  }

  if (!pd_ok(outpd[1])) { 

    pd_printf("fail (output pd 1 %PD does not pass pd_ok)\n",outpd[0]);
    return false;

  }

  if (!is_0_crossing_diagram(outpd[0])) { 

    pd_printf("fail (output pd 0 %PD doesn't match 0_crossing_diagram)\n",outpd[0]);
    return false;

  }

  if (!is_0_crossing_diagram(outpd[1])) { 

    pd_printf("fail (output pd 1 %PD doesn't match 0_crossing_diagram)\n",outpd[0]);
    return false;

  }

  /* Try to figure out how the tags should be set. */

  pd_idx_t over_in,over_out;
  pd_idx_t over_comp,over_pos;

  pd_overstrand(pd,0,&over_in,&over_out);
  pd_component_and_pos(pd,over_in,&over_comp,&over_pos);

  pd_idx_t under_comp = (over_comp == 0) ? 1 : 0;

  if (outpd[0]->comp[0].tag != pd->comp[over_comp].tag) { 

    printf("fail (tag of output 0 %c != tag of input over component 0 %c)\n",outpd[0]->comp[0].tag,pd->comp[over_comp].tag);
    return false;

  }

  if (outpd[1]->comp[0].tag != pd->comp[under_comp].tag) { 

    printf("fail (tag of output 1 %c != tag of input under component 1 %c)\n",outpd[1]->comp[0].tag,pd->comp[under_comp].tag);
    return false;

  }

  printf("pass (output matches pair of 0 crossing loops)\n");
  pd_code_free(&outpd[0]);
  pd_code_free(&outpd[1]);
  pd_code_free(&pd);
  free(outpd);

  printf("--------------------------------------------\n"
	 "r2 bigon elimination tiny tests: PASS       \n"
	 "--------------------------------------------\n");
  
  return true;

} 
  
bool r2_tests() {

  if (!r2_twist_tests(4)) { return false; }
  if (!r2_twist_tests(5)) { return false; }
  if (!r2_tiny_tests()) { return false; }

  return true;
}

pd_code_t *pd_create_r1_testA_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(3);
assert(pd != NULL);
pd->ncross = 3;
pd->nedges = 6;
pd->ncomps = 1;
pd->nfaces = 5;
sprintf(pd->hash,"%s","AwYFBAMDAQEBBgAAAAAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 0;
pd->cross[0].edge[2] = 5;
pd->cross[0].edge[3] = 1;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 1;
pd->cross[1].edge[1] = 4;
pd->cross[1].edge[2] = 2;
pd->cross[1].edge[3] = 5;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 2;
pd->cross[2].edge[1] = 3;
pd->cross[2].edge[2] = 3;
pd->cross[2].edge[3] = 4;
pd->cross[2].sign = 1;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 1;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 1;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 3;

pd->edge[2].head = 2;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 1;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 2;
pd->edge[3].headpos = 1;
pd->edge[3].tail = 2;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 1;
pd->edge[4].headpos = 1;
pd->edge[4].tail = 2;
pd->edge[4].tailpos = 3;

pd->edge[5].head = 0;
pd->edge[5].headpos = 2;
pd->edge[5].tail = 1;
pd->edge[5].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 6;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;


/* Face data */

pd->face[0].nedges = 4;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 5;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 2;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 4;
pd->face[0].or[3] = 1;

pd->face[1].nedges = 3;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 0;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 1;
pd->face[1].or[1] = 1;

pd->face[1].edge[2] = 5;
pd->face[1].or[2] = 1;

pd->face[2].nedges = 3;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 2;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 4;
pd->face[2].or[1] = 0;

pd->face[2].edge[2] = 3;
pd->face[2].or[2] = 1;

pd->face[3].nedges = 1;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 0;
pd->face[3].or[0] = 1;

pd->face[4].nedges = 1;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 3;
pd->face[4].or[0] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}

pd_code_t *pd_create_r1_testA_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(2);
assert(pd != NULL);
pd->ncross = 2;
pd->nedges = 4;
pd->ncomps = 1;
pd->nfaces = 4;
sprintf(pd->hash,"%s","AgQEAwMBAQEEAAAAAAAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 1;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 2;
pd->cross[0].sign = 1;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 3;
pd->cross[1].edge[2] = 3;
pd->cross[1].edge[3] = 2;
pd->cross[1].sign = 1;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 1;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 0;
pd->edge[1].headpos = 1;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 1;
pd->edge[2].headpos = 3;
pd->edge[2].tail = 0;
pd->edge[2].tailpos = 3;

pd->edge[3].head = 1;
pd->edge[3].headpos = 2;
pd->edge[3].tail = 1;
pd->edge[3].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 4;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;


/* Face data */

pd->face[0].nedges = 3;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 2;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 1;
pd->face[0].or[2] = 1;

pd->face[1].nedges = 3;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 0;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 2;
pd->face[1].or[1] = 1;

pd->face[1].edge[2] = 3;
pd->face[1].or[2] = 0;

pd->face[2].nedges = 1;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 1;
pd->face[2].or[0] = 0;

pd->face[3].nedges = 1;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 3;
pd->face[3].or[0] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}

bool r1testA() 
{ 

/*                                                                            */
/*     _________________5_________<_____                                      */
/*    |                                /                                      */
/*    |                               /               (0)                     */
/*    |           ____0_____         /                                        */
/*    |    (1)    \        /        /                                         */
/*    |            \  (3) /        ^                                          */
/*    |             \    ^        /                                           */
/*    |              \  /        /                                            */
/*    |                /        /                                             */
/*    |______>________/ ____1__________>_______2___________                    */
/*                      0       1                         /                    */
/*                            /                          /                     */
/*                           /                          /                      */
/*                          /                          /                       */
/*                         /                (2)       /                        */
/*                        /                          /                         */
/*                       4   ____3_>___             /                          */
/*                      /    \        /            /                          */
/*                     /      \ (4)  /            /                           */
/*                    /        ^    /            /                            */
/*                   /          \  /            /                             */
/*                  /            \             /                              */
/*                 /______<______ \_____<_____/                               */
/*                                2                                           */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*R1 operation at vertex 0            					                      */
/*                                                                            */
/*       __________3___________                                               */
/*      |                     /                                               */
/*      |                    /                                                */
/*      |                   /                                                 */
/*      |       (2)        /                                                  */
/*      |                 ^                            (1)                    */
/*      |                /                                                    */
/*      |               /                                                     */
/*      |_________________>________0__________                                */
/*                   0                        /                               */
/*                    /                      /                                */
/*                   /          (0)         /                                 */
/*                  /                      /                                  */
/*                 /                      /                                   */
/*                / _____1____           /                                    */
/*               /  \        /          /                                     */
/*              2    \ (3)  /          /                                      */
/*             /      ^    /          /                                       */
/*            /        \  /          /                                        */
/*           /          \           /                                         */
/*          /______<____ \_________/                                          */
/*                       1                                                    */
/*                                                                            */

  printf("---------------------------------------\n"
	 "r1 test A\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r1_testA_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after configuration...");
  pd_code_t *after = pd_create_r1_testA_after_0();
  printf("done (passes pd_ok)\n");
  
  printf("performing r1 loop deletion at crossing 0...");
  pd_code_t *newpd = pd_R1_loopdeletion(before,0);
  if (!pd_ok(newpd)) { 
    printf("fail (output pd does not pass pd_ok)\n");
    return false;
  }
  printf("pass (output passes pd_ok)\n");

  printf("testing for isomorphism with after configuration...");
  if (!pd_isomorphic(newpd,after)) { 
    pd_printf("fail (output pd %PD ",newpd);
    pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after);
    return false;
  }

  printf("pass (output and after isomorphic)\n");
  
  printf("housecleaning...");
  pd_code_free(&before);
  pd_code_free(&after);
  pd_code_free(&newpd);
  printf("done\n");

  printf("---------------------------------------\n"
	 "r1 test A : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}

pd_code_t *pd_create_r1_testB_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(6);
assert(pd != NULL);
pd->ncross = 6;
pd->nedges = 12;
pd->ncomps = 2;
pd->nfaces = 8;
sprintf(pd->hash,"%s","BgwIBgQDAwMCAgECBgYAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 6;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 7;
pd->cross[0].sign = 1;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 7;
pd->cross[1].edge[2] = 5;
pd->cross[1].edge[3] = 8;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 3;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 2;
pd->cross[2].sign = 0;

pd->cross[3].edge[0] = 3;
pd->cross[3].edge[1] = 9;
pd->cross[3].edge[2] = 4;
pd->cross[3].edge[3] = 10;
pd->cross[3].sign = 0;

pd->cross[4].edge[0] = 4;
pd->cross[4].edge[1] = 11;
pd->cross[4].edge[2] = 5;
pd->cross[4].edge[3] = 10;
pd->cross[4].sign = 0;

pd->cross[5].edge[0] = 6;
pd->cross[5].edge[1] = 8;
pd->cross[5].edge[2] = 11;
pd->cross[5].edge[3] = 9;
pd->cross[5].sign = 0;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 1;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 2;
pd->edge[2].headpos = 3;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 3;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 2;
pd->edge[3].tailpos = 1;

pd->edge[4].head = 4;
pd->edge[4].headpos = 0;
pd->edge[4].tail = 3;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 1;
pd->edge[5].headpos = 2;
pd->edge[5].tail = 4;
pd->edge[5].tailpos = 2;

pd->edge[6].head = 0;
pd->edge[6].headpos = 1;
pd->edge[6].tail = 5;
pd->edge[6].tailpos = 0;

pd->edge[7].head = 1;
pd->edge[7].headpos = 1;
pd->edge[7].tail = 0;
pd->edge[7].tailpos = 3;

pd->edge[8].head = 5;
pd->edge[8].headpos = 1;
pd->edge[8].tail = 1;
pd->edge[8].tailpos = 3;

pd->edge[9].head = 3;
pd->edge[9].headpos = 1;
pd->edge[9].tail = 5;
pd->edge[9].tailpos = 3;

pd->edge[10].head = 4;
pd->edge[10].headpos = 3;
pd->edge[10].tail = 3;
pd->edge[10].tailpos = 3;

pd->edge[11].head = 5;
pd->edge[11].headpos = 2;
pd->edge[11].tail = 4;
pd->edge[11].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 6;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;

pd->comp[1].nedges = 6;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 6;
pd->comp[1].edge[1] = 7;
pd->comp[1].edge[2] = 8;
pd->comp[1].edge[3] = 9;
pd->comp[1].edge[4] = 10;
pd->comp[1].edge[5] = 11;


/* Face data */

pd->face[0].nedges = 6;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 2;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 3;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 10;
pd->face[0].or[3] = 1;

pd->face[0].edge[4] = 5;
pd->face[0].or[4] = 1;

pd->face[0].edge[5] = 7;
pd->face[0].or[5] = 0;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 6;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 9;
pd->face[1].or[2] = 1;

pd->face[1].edge[3] = 3;
pd->face[1].or[3] = 0;

pd->face[2].nedges = 3;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 8;
pd->face[2].or[1] = 1;

pd->face[2].edge[2] = 6;
pd->face[2].or[2] = 1;

pd->face[3].nedges = 3;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 4;
pd->face[3].or[0] = 0;

pd->face[3].edge[1] = 9;
pd->face[3].or[1] = 0;

pd->face[3].edge[2] = 11;
pd->face[3].or[2] = 0;

pd->face[4].nedges = 3;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 5;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 11;
pd->face[4].or[1] = 1;

pd->face[4].edge[2] = 8;
pd->face[4].or[2] = 0;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 0;
pd->face[5].or[0] = 1;

pd->face[5].edge[1] = 7;
pd->face[5].or[1] = 1;

pd->face[6].nedges = 2;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 4;
pd->face[6].or[0] = 1;

pd->face[6].edge[1] = 10;
pd->face[6].or[1] = 0;

pd->face[7].nedges = 1;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 2;
pd->face[7].or[0] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}

pd_code_t *pd_create_r1_testB_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(5);
assert(pd != NULL);
pd->ncross = 5;
pd->nedges = 10;
pd->ncomps = 2;
pd->nfaces = 7;
sprintf(pd->hash,"%s","BQoHBAMDAwMCAgIGBAAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 2;
pd->cross[0].edge[2] = 5;
pd->cross[0].edge[3] = 3;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 7;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 6;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 9;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 6;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 3;
pd->cross[3].edge[1] = 8;
pd->cross[3].edge[2] = 4;
pd->cross[3].edge[3] = 7;
pd->cross[3].sign = 1;

pd->cross[4].edge[0] = 4;
pd->cross[4].edge[1] = 8;
pd->cross[4].edge[2] = 5;
pd->cross[4].edge[3] = 9;
pd->cross[4].sign = 0;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 0;
pd->edge[2].headpos = 1;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 3;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 0;
pd->edge[3].tailpos = 3;

pd->edge[4].head = 4;
pd->edge[4].headpos = 0;
pd->edge[4].tail = 3;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 0;
pd->edge[5].headpos = 2;
pd->edge[5].tail = 4;
pd->edge[5].tailpos = 2;

pd->edge[6].head = 1;
pd->edge[6].headpos = 3;
pd->edge[6].tail = 2;
pd->edge[6].tailpos = 3;

pd->edge[7].head = 3;
pd->edge[7].headpos = 3;
pd->edge[7].tail = 1;
pd->edge[7].tailpos = 1;

pd->edge[8].head = 4;
pd->edge[8].headpos = 1;
pd->edge[8].tail = 3;
pd->edge[8].tailpos = 1;

pd->edge[9].head = 2;
pd->edge[9].headpos = 1;
pd->edge[9].tail = 4;
pd->edge[9].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 6;
pd->comp[0].tag = 'B';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;

pd->comp[1].nedges = 4;
pd->comp[1].tag = 'A';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 6;
pd->comp[1].edge[1] = 7;
pd->comp[1].edge[2] = 8;
pd->comp[1].edge[3] = 9;


/* Face data */

pd->face[0].nedges = 4;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 7;
pd->face[0].or[1] = 1;

pd->face[0].edge[2] = 4;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 9;
pd->face[0].or[3] = 1;

pd->face[1].nedges = 3;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 0;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 3;
pd->face[1].or[1] = 1;

pd->face[1].edge[2] = 7;
pd->face[1].or[2] = 0;

pd->face[2].nedges = 3;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 6;
pd->face[2].or[1] = 0;

pd->face[2].edge[2] = 2;
pd->face[2].or[2] = 1;

pd->face[3].nedges = 3;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 2;
pd->face[3].or[0] = 0;

pd->face[3].edge[1] = 9;
pd->face[3].or[1] = 0;

pd->face[3].edge[2] = 5;
pd->face[3].or[2] = 1;

pd->face[4].nedges = 3;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 3;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 5;
pd->face[4].or[1] = 0;

pd->face[4].edge[2] = 8;
pd->face[4].or[2] = 0;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 1;
pd->face[5].or[0] = 1;

pd->face[5].edge[1] = 6;
pd->face[5].or[1] = 1;

pd->face[6].nedges = 2;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 4;
pd->face[6].or[0] = 0;

pd->face[6].edge[1] = 8;
pd->face[6].or[1] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}

/*                                                                            */
/*                       ____2___                                             */
/*                       \      /                                             */
/*                        \(7) /                  (0)                         */
/*                         \  /                                               */
/*                           /                                                */
/*                          /                                                 */
/*                ____1____/ 2\_>__3____                                      */
/*               ^                     |                                      */
/*              0|      (1)   __9___>__________                               */
/*         __<__ | __6_      /           3     |                              */
/*        |      |     \    /          |       |                              */
/*        |      |      ^  /    (3)    |      10                              */
/*        |      0 (2)   \             4       |                              */
/*        7 (5)  |        \            | (6)   |                              */
/*        |      ^      / 5\           |       |                              */
/*        |      |     /    \__11__<__ |_______                               */
/*        |________>_8/                |4                                     */
/*              1                      |                                      */
/*               |       (4)          \/                                      */
/*               |                     |                                      */
/*               |____<____5___________|                                      */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*R1 operation at vertex 2	                                                  */  
/*                                                                            */
/*                                            (0)                             */
/*                _______7__>___________                                      */
/*               |                      |                                     */
/*               |      (1)    _3____>_________                                */
/*         __<__ | _0__      /           3     |                              */
/*        |     1|     \    /           |      |                              */
/*        |      |      ^  /            |      4                              */
/*        |      6        \0   (4)      8 (6)  |                              */
/*        1      | (2)     \            |      |                              */
/*        |      ^      /   \           |      |                              */
/*        | (5)  |     /     \___5__<__ |______|                              */
/*        |________> 2/                 |4                                    */
/*              2|                      |                                     */
/*               |        (3)          \/                                     */
/*               |                      |                                     */
/*               |__________<_____9_____|                                     */
/*                                                                            */
/*                                                                            */
bool r1testB() {
  
  printf("---------------------------------------\n"
	 "r1 test B\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r1_testB_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after configuration...");
  pd_code_t *after = pd_create_r1_testB_after_0();
  printf("done (passes pd_ok)\n");
  
  printf("performing r1 loop deletion at crossing 2...");
  pd_code_t *newpd = pd_R1_loopdeletion(before,2);
  if (!pd_ok(newpd)) { 
    printf("fail (output pd does not pass pd_ok)\n");
    return false;
  }
  printf("pass (output passes pd_ok)\n");

  printf("testing for isomorphism with after configuration...");
  if (!pd_isomorphic(newpd,after)) { 
    pd_printf("fail (output pd %PD ",newpd);
    pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after);
    return false;
  }

  printf("pass (output and after isomorphic)\n");
  
  printf("housecleaning...");
  pd_code_free(&before);
  pd_code_free(&after);
  pd_code_free(&newpd);
  printf("done\n");

  printf("---------------------------------------\n"
	 "r1 test B : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}

pd_code_t *pd_create_r1_testC_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(8);
assert(pd != NULL);
pd->ncross = 8;
pd->nedges = 16;
pd->ncomps = 2;
pd->nfaces = 10;
sprintf(pd->hash,"%s","CBAKBQQEBAMDAwMCAQIKBgAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 5;
pd->cross[0].edge[2] = 9;
pd->cross[0].edge[3] = 4;
pd->cross[0].sign = 1;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 10;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 15;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 7;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 6;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 13;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 14;
pd->cross[3].sign = 1;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 8;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 9;
pd->cross[4].sign = 1;

pd->cross[5].edge[0] = 5;
pd->cross[5].edge[1] = 15;
pd->cross[5].edge[2] = 6;
pd->cross[5].edge[3] = 14;
pd->cross[5].sign = 1;

pd->cross[6].edge[0] = 7;
pd->cross[6].edge[1] = 10;
pd->cross[6].edge[2] = 8;
pd->cross[6].edge[3] = 11;
pd->cross[6].sign = 0;

pd->cross[7].edge[0] = 11;
pd->cross[7].edge[1] = 12;
pd->cross[7].edge[2] = 12;
pd->cross[7].edge[3] = 13;
pd->cross[7].sign = 0;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 0;
pd->edge[4].headpos = 3;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 5;
pd->edge[5].headpos = 0;
pd->edge[5].tail = 0;
pd->edge[5].tailpos = 1;

pd->edge[6].head = 2;
pd->edge[6].headpos = 3;
pd->edge[6].tail = 5;
pd->edge[6].tailpos = 2;

pd->edge[7].head = 6;
pd->edge[7].headpos = 0;
pd->edge[7].tail = 2;
pd->edge[7].tailpos = 1;

pd->edge[8].head = 4;
pd->edge[8].headpos = 1;
pd->edge[8].tail = 6;
pd->edge[8].tailpos = 2;

pd->edge[9].head = 0;
pd->edge[9].headpos = 2;
pd->edge[9].tail = 4;
pd->edge[9].tailpos = 3;

pd->edge[10].head = 6;
pd->edge[10].headpos = 1;
pd->edge[10].tail = 1;
pd->edge[10].tailpos = 1;

pd->edge[11].head = 7;
pd->edge[11].headpos = 0;
pd->edge[11].tail = 6;
pd->edge[11].tailpos = 3;

pd->edge[12].head = 7;
pd->edge[12].headpos = 1;
pd->edge[12].tail = 7;
pd->edge[12].tailpos = 2;

pd->edge[13].head = 3;
pd->edge[13].headpos = 1;
pd->edge[13].tail = 7;
pd->edge[13].tailpos = 3;

pd->edge[14].head = 5;
pd->edge[14].headpos = 3;
pd->edge[14].tail = 3;
pd->edge[14].tailpos = 3;

pd->edge[15].head = 1;
pd->edge[15].headpos = 3;
pd->edge[15].tail = 5;
pd->edge[15].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 10;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;

pd->comp[1].nedges = 6;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 10;
pd->comp[1].edge[1] = 11;
pd->comp[1].edge[2] = 12;
pd->comp[1].edge[3] = 13;
pd->comp[1].edge[4] = 14;
pd->comp[1].edge[5] = 15;


/* Face data */

pd->face[0].nedges = 5;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 3;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 13;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 12;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 11;
pd->face[0].or[3] = 0;

pd->face[0].edge[4] = 8;
pd->face[0].or[4] = 1;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 0;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 4;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 8;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 10;
pd->face[1].or[3] = 0;

pd->face[2].nedges = 4;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 2;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 7;
pd->face[2].or[1] = 1;

pd->face[2].edge[2] = 11;
pd->face[2].or[2] = 1;

pd->face[2].edge[3] = 13;
pd->face[2].or[3] = 1;

pd->face[3].nedges = 4;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 3;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 9;
pd->face[3].or[1] = 1;

pd->face[3].edge[2] = 5;
pd->face[3].or[2] = 1;

pd->face[3].edge[3] = 14;
pd->face[3].or[3] = 0;

pd->face[4].nedges = 3;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 0;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 15;
pd->face[4].or[1] = 0;

pd->face[4].edge[2] = 5;
pd->face[4].or[2] = 0;

pd->face[5].nedges = 3;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 1;
pd->face[5].or[0] = 1;

pd->face[5].edge[1] = 6;
pd->face[5].or[1] = 0;

pd->face[5].edge[2] = 15;
pd->face[5].or[2] = 1;

pd->face[6].nedges = 3;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 1;
pd->face[6].or[0] = 0;

pd->face[6].edge[1] = 10;
pd->face[6].or[1] = 1;

pd->face[6].edge[2] = 7;
pd->face[6].or[2] = 0;

pd->face[7].nedges = 3;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 2;
pd->face[7].or[0] = 1;

pd->face[7].edge[1] = 14;
pd->face[7].or[1] = 1;

pd->face[7].edge[2] = 6;
pd->face[7].or[2] = 1;

pd->face[8].nedges = 2;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].or = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].or != NULL);

pd->face[8].edge[0] = 4;
pd->face[8].or[0] = 1;

pd->face[8].edge[1] = 9;
pd->face[8].or[1] = 0;

pd->face[9].nedges = 1;
pd->face[9].edge = calloc(pd->face[9].nedges,sizeof(pd_idx_t));
pd->face[9].or = calloc(pd->face[9].nedges,sizeof(pd_or_t));
assert(pd->face[9].edge != NULL);
assert(pd->face[9].or != NULL);

pd->face[9].edge[0] = 12;
pd->face[9].or[0] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}

pd_code_t *pd_create_r1_testC_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(7);
assert(pd != NULL);
pd->ncross = 7;
pd->nedges = 14;
pd->ncomps = 2;
pd->nfaces = 9;
sprintf(pd->hash,"%s","Bw4JBAQDAwMDAwMCAgoEAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 5;
pd->cross[0].edge[2] = 9;
pd->cross[0].edge[3] = 4;
pd->cross[0].sign = 1;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 10;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 13;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 7;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 6;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 11;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 12;
pd->cross[3].sign = 1;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 8;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 9;
pd->cross[4].sign = 1;

pd->cross[5].edge[0] = 5;
pd->cross[5].edge[1] = 13;
pd->cross[5].edge[2] = 6;
pd->cross[5].edge[3] = 12;
pd->cross[5].sign = 1;

pd->cross[6].edge[0] = 7;
pd->cross[6].edge[1] = 10;
pd->cross[6].edge[2] = 8;
pd->cross[6].edge[3] = 11;
pd->cross[6].sign = 0;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 0;
pd->edge[4].headpos = 3;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 5;
pd->edge[5].headpos = 0;
pd->edge[5].tail = 0;
pd->edge[5].tailpos = 1;

pd->edge[6].head = 2;
pd->edge[6].headpos = 3;
pd->edge[6].tail = 5;
pd->edge[6].tailpos = 2;

pd->edge[7].head = 6;
pd->edge[7].headpos = 0;
pd->edge[7].tail = 2;
pd->edge[7].tailpos = 1;

pd->edge[8].head = 4;
pd->edge[8].headpos = 1;
pd->edge[8].tail = 6;
pd->edge[8].tailpos = 2;

pd->edge[9].head = 0;
pd->edge[9].headpos = 2;
pd->edge[9].tail = 4;
pd->edge[9].tailpos = 3;

pd->edge[10].head = 6;
pd->edge[10].headpos = 1;
pd->edge[10].tail = 1;
pd->edge[10].tailpos = 1;

pd->edge[11].head = 3;
pd->edge[11].headpos = 1;
pd->edge[11].tail = 6;
pd->edge[11].tailpos = 3;

pd->edge[12].head = 5;
pd->edge[12].headpos = 3;
pd->edge[12].tail = 3;
pd->edge[12].tailpos = 3;

pd->edge[13].head = 1;
pd->edge[13].headpos = 3;
pd->edge[13].tail = 5;
pd->edge[13].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 10;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;

pd->comp[1].nedges = 4;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 10;
pd->comp[1].edge[1] = 11;
pd->comp[1].edge[2] = 12;
pd->comp[1].edge[3] = 13;


/* Face data */

pd->face[0].nedges = 4;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 4;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 8;
pd->face[0].or[2] = 0;

pd->face[0].edge[3] = 10;
pd->face[0].or[3] = 0;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 3;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 9;
pd->face[1].or[1] = 1;

pd->face[1].edge[2] = 5;
pd->face[1].or[2] = 1;

pd->face[1].edge[3] = 12;
pd->face[1].or[3] = 0;

pd->face[2].nedges = 3;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 13;
pd->face[2].or[1] = 0;

pd->face[2].edge[2] = 5;
pd->face[2].or[2] = 0;

pd->face[3].nedges = 3;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 1;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 6;
pd->face[3].or[1] = 0;

pd->face[3].edge[2] = 13;
pd->face[3].or[2] = 1;

pd->face[4].nedges = 3;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 1;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 10;
pd->face[4].or[1] = 1;

pd->face[4].edge[2] = 7;
pd->face[4].or[2] = 0;

pd->face[5].nedges = 3;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 2;
pd->face[5].or[0] = 0;

pd->face[5].edge[1] = 7;
pd->face[5].or[1] = 1;

pd->face[5].edge[2] = 11;
pd->face[5].or[2] = 1;

pd->face[6].nedges = 3;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 2;
pd->face[6].or[0] = 1;

pd->face[6].edge[1] = 12;
pd->face[6].or[1] = 1;

pd->face[6].edge[2] = 6;
pd->face[6].or[2] = 1;

pd->face[7].nedges = 3;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 3;
pd->face[7].or[0] = 0;

pd->face[7].edge[1] = 11;
pd->face[7].or[1] = 0;

pd->face[7].edge[2] = 8;
pd->face[7].or[2] = 1;

pd->face[8].nedges = 2;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].or = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].or != NULL);

pd->face[8].edge[0] = 4;
pd->face[8].or[0] = 1;

pd->face[8].edge[1] = 9;
pd->face[8].or[1] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}

/*                                                                            */
/*                          ________9___>__________                           */
/*                          \                      /         (3)               */
/*                           \         (8)        /                            */
/*                            \                  /                             */
/*                  ___________________4______  / __>______________             */
/*                  \         4 \              /0                 /             */
/*                   \ (0) __>_  8   (1)      /                  /              */
/*                    \   |(9)  | \          0                  /               */
/*                     \  |_12__|_11____10_ /___       (4)     /                */
/*                      \     7 |   6\     / 1 |              5                 */
/*                       3      |     \(6)1    15            /                  */
/*                        \     | (2)  7   (5) |            /                   */
/*                         \    13      \      ^           /                    */
/*                          \   |      /2\     |          /                     */
/*                           \_____2__/   \_<6_| ________/                      */
/*                              |3             |5                              */
/*                              |              |                               */
/*                              |      (7)     |                               */
/*                              |__>_14________|                               */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*R1 operation at vertex 7		                                              */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                          ________9___>__________                           */
/*                          \                     /        (1)                */
/*                           \        (8)        /                            */
/*                            \                 /                             */
/*                  ___________________4______ / __>_____________             */
/*                  \         4 \             /0                /             */
/*                   \ (7)       8    (0)    /                 /              */
/*                    \           \         0                 /               */
/*                     \       __<____10__ /___     (2)      /                */
/*                      \      |   6\ (4) /1  |             5                 */
/*                       3     |     \   1    13           /                  */
/*                        \   11 (5)  7 / (3) |           /                   */
/*                         \   |       \      ^          /                    */
/*                          \  |      /2\     |         /                     */
/*                           \_____2_/   _<6__|________/                      */
/*                             |3             |5                              */
/*                             |     (6)      |                               */
/*                             |              |                               */
/*                             |__>____12_____|                               */
/*                                                                            */
bool r1testC() {
  
  printf("---------------------------------------\n"
	 "r1 test C\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r1_testC_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after configuration...");
  pd_code_t *after = pd_create_r1_testC_after_0();
  printf("done (passes pd_ok)\n");
  
  printf("performing r1 loop deletion at crossing 7...");
  pd_code_t *newpd = pd_R1_loopdeletion(before,7);
  if (!pd_ok(newpd)) { 
    printf("fail (output pd does not pass pd_ok)\n");
    return false;
  }
  printf("pass (output passes pd_ok)\n");

  printf("testing for isomorphism with after configuration...");
  if (!pd_isomorphic(newpd,after)) { 
    pd_printf("fail (output pd %PD ",newpd);
    pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after);
    return false;
  }

  printf("pass (output and after isomorphic)\n");
  
  printf("housecleaning...");
  pd_code_free(&before);
  pd_code_free(&after);
  pd_code_free(&newpd);
  printf("done\n");

  printf("---------------------------------------\n"
	 "r1 test C : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}

pd_code_t *pd_create_r1_testD_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(15);
assert(pd != NULL);
pd->ncross = 15;
pd->nedges = 30;
pd->ncomps = 3;
pd->nfaces = 17;
sprintf(pd->hash,"%s","Dx4RBwcFBAQEBAMDAwMDAwICAgEDDAw");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 7;
pd->cross[0].edge[2] = 11;
pd->cross[0].edge[3] = 6;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 8;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 7;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 20;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 21;
pd->cross[2].sign = 0;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 9;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 10;
pd->cross[3].sign = 1;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 26;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 25;
pd->cross[4].sign = 0;

pd->cross[5].edge[0] = 4;
pd->cross[5].edge[1] = 18;
pd->cross[5].edge[2] = 5;
pd->cross[5].edge[3] = 17;
pd->cross[5].sign = 1;

pd->cross[6].edge[0] = 5;
pd->cross[6].edge[1] = 29;
pd->cross[6].edge[2] = 6;
pd->cross[6].edge[3] = 24;
pd->cross[6].sign = 1;

pd->cross[7].edge[0] = 8;
pd->cross[7].edge[1] = 19;
pd->cross[7].edge[2] = 9;
pd->cross[7].edge[3] = 20;
pd->cross[7].sign = 0;

pd->cross[8].edge[0] = 10;
pd->cross[8].edge[1] = 16;
pd->cross[8].edge[2] = 11;
pd->cross[8].edge[3] = 15;
pd->cross[8].sign = 0;

pd->cross[9].edge[0] = 12;
pd->cross[9].edge[1] = 15;
pd->cross[9].edge[2] = 23;
pd->cross[9].edge[3] = 14;
pd->cross[9].sign = 1;

pd->cross[10].edge[0] = 12;
pd->cross[10].edge[1] = 22;
pd->cross[10].edge[2] = 13;
pd->cross[10].edge[3] = 21;
pd->cross[10].sign = 1;

pd->cross[11].edge[0] = 13;
pd->cross[11].edge[1] = 22;
pd->cross[11].edge[2] = 14;
pd->cross[11].edge[3] = 23;
pd->cross[11].sign = 1;

pd->cross[12].edge[0] = 16;
pd->cross[12].edge[1] = 25;
pd->cross[12].edge[2] = 17;
pd->cross[12].edge[3] = 24;
pd->cross[12].sign = 0;

pd->cross[13].edge[0] = 18;
pd->cross[13].edge[1] = 26;
pd->cross[13].edge[2] = 19;
pd->cross[13].edge[3] = 27;
pd->cross[13].sign = 1;

pd->cross[14].edge[0] = 27;
pd->cross[14].edge[1] = 28;
pd->cross[14].edge[2] = 28;
pd->cross[14].edge[3] = 29;
pd->cross[14].sign = 1;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 5;
pd->edge[4].headpos = 0;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 6;
pd->edge[5].headpos = 0;
pd->edge[5].tail = 5;
pd->edge[5].tailpos = 2;

pd->edge[6].head = 0;
pd->edge[6].headpos = 3;
pd->edge[6].tail = 6;
pd->edge[6].tailpos = 2;

pd->edge[7].head = 1;
pd->edge[7].headpos = 3;
pd->edge[7].tail = 0;
pd->edge[7].tailpos = 1;

pd->edge[8].head = 7;
pd->edge[8].headpos = 0;
pd->edge[8].tail = 1;
pd->edge[8].tailpos = 1;

pd->edge[9].head = 3;
pd->edge[9].headpos = 1;
pd->edge[9].tail = 7;
pd->edge[9].tailpos = 2;

pd->edge[10].head = 8;
pd->edge[10].headpos = 0;
pd->edge[10].tail = 3;
pd->edge[10].tailpos = 3;

pd->edge[11].head = 0;
pd->edge[11].headpos = 2;
pd->edge[11].tail = 8;
pd->edge[11].tailpos = 2;

pd->edge[12].head = 10;
pd->edge[12].headpos = 0;
pd->edge[12].tail = 9;
pd->edge[12].tailpos = 0;

pd->edge[13].head = 11;
pd->edge[13].headpos = 0;
pd->edge[13].tail = 10;
pd->edge[13].tailpos = 2;

pd->edge[14].head = 9;
pd->edge[14].headpos = 3;
pd->edge[14].tail = 11;
pd->edge[14].tailpos = 2;

pd->edge[15].head = 8;
pd->edge[15].headpos = 3;
pd->edge[15].tail = 9;
pd->edge[15].tailpos = 1;

pd->edge[16].head = 12;
pd->edge[16].headpos = 0;
pd->edge[16].tail = 8;
pd->edge[16].tailpos = 1;

pd->edge[17].head = 5;
pd->edge[17].headpos = 3;
pd->edge[17].tail = 12;
pd->edge[17].tailpos = 2;

pd->edge[18].head = 13;
pd->edge[18].headpos = 0;
pd->edge[18].tail = 5;
pd->edge[18].tailpos = 1;

pd->edge[19].head = 7;
pd->edge[19].headpos = 1;
pd->edge[19].tail = 13;
pd->edge[19].tailpos = 2;

pd->edge[20].head = 2;
pd->edge[20].headpos = 1;
pd->edge[20].tail = 7;
pd->edge[20].tailpos = 3;

pd->edge[21].head = 10;
pd->edge[21].headpos = 3;
pd->edge[21].tail = 2;
pd->edge[21].tailpos = 3;

pd->edge[22].head = 11;
pd->edge[22].headpos = 1;
pd->edge[22].tail = 10;
pd->edge[22].tailpos = 1;

pd->edge[23].head = 9;
pd->edge[23].headpos = 2;
pd->edge[23].tail = 11;
pd->edge[23].tailpos = 3;

pd->edge[24].head = 12;
pd->edge[24].headpos = 3;
pd->edge[24].tail = 6;
pd->edge[24].tailpos = 3;

pd->edge[25].head = 4;
pd->edge[25].headpos = 3;
pd->edge[25].tail = 12;
pd->edge[25].tailpos = 1;

pd->edge[26].head = 13;
pd->edge[26].headpos = 1;
pd->edge[26].tail = 4;
pd->edge[26].tailpos = 1;

pd->edge[27].head = 14;
pd->edge[27].headpos = 0;
pd->edge[27].tail = 13;
pd->edge[27].tailpos = 3;

pd->edge[28].head = 14;
pd->edge[28].headpos = 1;
pd->edge[28].tail = 14;
pd->edge[28].tailpos = 2;

pd->edge[29].head = 6;
pd->edge[29].headpos = 1;
pd->edge[29].tail = 14;
pd->edge[29].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 12;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;
pd->comp[0].edge[10] = 10;
pd->comp[0].edge[11] = 11;

pd->comp[1].nedges = 12;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 12;
pd->comp[1].edge[1] = 13;
pd->comp[1].edge[2] = 14;
pd->comp[1].edge[3] = 15;
pd->comp[1].edge[4] = 16;
pd->comp[1].edge[5] = 17;
pd->comp[1].edge[6] = 18;
pd->comp[1].edge[7] = 19;
pd->comp[1].edge[8] = 20;
pd->comp[1].edge[9] = 21;
pd->comp[1].edge[10] = 22;
pd->comp[1].edge[11] = 23;

pd->comp[2].nedges = 6;
pd->comp[2].tag = 'C';

pd->comp[2].edge = calloc(pd->comp[2].nedges,sizeof(pd_idx_t));
assert(pd->comp[2].edge != NULL);

pd->comp[2].edge[0] = 24;
pd->comp[2].edge[1] = 25;
pd->comp[2].edge[2] = 26;
pd->comp[2].edge[3] = 27;
pd->comp[2].edge[4] = 28;
pd->comp[2].edge[5] = 29;


/* Face data */

pd->face[0].nedges = 7;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 6;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 29;
pd->face[0].or[2] = 0;

pd->face[0].edge[3] = 28;
pd->face[0].or[3] = 1;

pd->face[0].edge[4] = 27;
pd->face[0].or[4] = 0;

pd->face[0].edge[5] = 19;
pd->face[0].or[5] = 1;

pd->face[0].edge[6] = 8;
pd->face[0].or[6] = 0;

pd->face[1].nedges = 7;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 21;
pd->face[1].or[1] = 1;

pd->face[1].edge[2] = 13;
pd->face[1].or[2] = 1;

pd->face[1].edge[3] = 23;
pd->face[1].or[3] = 1;

pd->face[1].edge[4] = 15;
pd->face[1].or[4] = 1;

pd->face[1].edge[5] = 11;
pd->face[1].or[5] = 1;

pd->face[1].edge[6] = 7;
pd->face[1].or[6] = 1;

pd->face[2].nedges = 5;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 2;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 10;
pd->face[2].or[1] = 1;

pd->face[2].edge[2] = 15;
pd->face[2].or[2] = 0;

pd->face[2].edge[3] = 12;
pd->face[2].or[3] = 1;

pd->face[2].edge[4] = 21;
pd->face[2].or[4] = 0;

pd->face[3].nedges = 4;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 3;
pd->face[3].or[0] = 0;

pd->face[3].edge[1] = 9;
pd->face[3].or[1] = 0;

pd->face[3].edge[2] = 19;
pd->face[3].or[2] = 0;

pd->face[3].edge[3] = 26;
pd->face[3].or[3] = 0;

pd->face[4].nedges = 4;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 3;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 25;
pd->face[4].or[1] = 0;

pd->face[4].edge[2] = 16;
pd->face[4].or[2] = 0;

pd->face[4].edge[3] = 10;
pd->face[4].or[3] = 0;

pd->face[5].nedges = 4;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 5;
pd->face[5].or[0] = 0;

pd->face[5].edge[1] = 18;
pd->face[5].or[1] = 1;

pd->face[5].edge[2] = 27;
pd->face[5].or[2] = 1;

pd->face[5].edge[3] = 29;
pd->face[5].or[3] = 1;

pd->face[6].nedges = 4;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 6;
pd->face[6].or[0] = 1;

pd->face[6].edge[1] = 11;
pd->face[6].or[1] = 0;

pd->face[6].edge[2] = 16;
pd->face[6].or[2] = 1;

pd->face[6].edge[3] = 24;
pd->face[6].or[3] = 0;

pd->face[7].nedges = 3;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 1;
pd->face[7].or[0] = 0;

pd->face[7].edge[1] = 8;
pd->face[7].or[1] = 1;

pd->face[7].edge[2] = 20;
pd->face[7].or[2] = 1;

pd->face[8].nedges = 3;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].or = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].or != NULL);

pd->face[8].edge[0] = 2;
pd->face[8].or[0] = 0;

pd->face[8].edge[1] = 20;
pd->face[8].or[1] = 0;

pd->face[8].edge[2] = 9;
pd->face[8].or[2] = 1;

pd->face[9].nedges = 3;
pd->face[9].edge = calloc(pd->face[9].nedges,sizeof(pd_idx_t));
pd->face[9].or = calloc(pd->face[9].nedges,sizeof(pd_or_t));
assert(pd->face[9].edge != NULL);
assert(pd->face[9].or != NULL);

pd->face[9].edge[0] = 4;
pd->face[9].or[0] = 1;

pd->face[9].edge[1] = 17;
pd->face[9].or[1] = 0;

pd->face[9].edge[2] = 25;
pd->face[9].or[2] = 1;

pd->face[10].nedges = 3;
pd->face[10].edge = calloc(pd->face[10].nedges,sizeof(pd_idx_t));
pd->face[10].or = calloc(pd->face[10].nedges,sizeof(pd_or_t));
assert(pd->face[10].edge != NULL);
assert(pd->face[10].or != NULL);

pd->face[10].edge[0] = 4;
pd->face[10].or[0] = 0;

pd->face[10].edge[1] = 26;
pd->face[10].or[1] = 1;

pd->face[10].edge[2] = 18;
pd->face[10].or[2] = 0;

pd->face[11].nedges = 3;
pd->face[11].edge = calloc(pd->face[11].nedges,sizeof(pd_idx_t));
pd->face[11].or = calloc(pd->face[11].nedges,sizeof(pd_or_t));
assert(pd->face[11].edge != NULL);
assert(pd->face[11].or != NULL);

pd->face[11].edge[0] = 5;
pd->face[11].or[0] = 1;

pd->face[11].edge[1] = 24;
pd->face[11].or[1] = 1;

pd->face[11].edge[2] = 17;
pd->face[11].or[2] = 1;

pd->face[12].nedges = 3;
pd->face[12].edge = calloc(pd->face[12].nedges,sizeof(pd_idx_t));
pd->face[12].or = calloc(pd->face[12].nedges,sizeof(pd_or_t));
assert(pd->face[12].edge != NULL);
assert(pd->face[12].or != NULL);

pd->face[12].edge[0] = 12;
pd->face[12].or[0] = 0;

pd->face[12].edge[1] = 14;
pd->face[12].or[1] = 0;

pd->face[12].edge[2] = 22;
pd->face[12].or[2] = 0;

pd->face[13].nedges = 2;
pd->face[13].edge = calloc(pd->face[13].nedges,sizeof(pd_idx_t));
pd->face[13].or = calloc(pd->face[13].nedges,sizeof(pd_or_t));
assert(pd->face[13].edge != NULL);
assert(pd->face[13].or != NULL);

pd->face[13].edge[0] = 0;
pd->face[13].or[0] = 1;

pd->face[13].edge[1] = 7;
pd->face[13].or[1] = 0;

pd->face[14].nedges = 2;
pd->face[14].edge = calloc(pd->face[14].nedges,sizeof(pd_idx_t));
pd->face[14].or = calloc(pd->face[14].nedges,sizeof(pd_or_t));
assert(pd->face[14].edge != NULL);
assert(pd->face[14].or != NULL);

pd->face[14].edge[0] = 13;
pd->face[14].or[0] = 0;

pd->face[14].edge[1] = 22;
pd->face[14].or[1] = 1;

pd->face[15].nedges = 2;
pd->face[15].edge = calloc(pd->face[15].nedges,sizeof(pd_idx_t));
pd->face[15].or = calloc(pd->face[15].nedges,sizeof(pd_or_t));
assert(pd->face[15].edge != NULL);
assert(pd->face[15].or != NULL);

pd->face[15].edge[0] = 14;
pd->face[15].or[0] = 1;

pd->face[15].edge[1] = 23;
pd->face[15].or[1] = 0;

pd->face[16].nedges = 1;
pd->face[16].edge = calloc(pd->face[16].nedges,sizeof(pd_idx_t));
pd->face[16].or = calloc(pd->face[16].nedges,sizeof(pd_or_t));
assert(pd->face[16].edge != NULL);
assert(pd->face[16].or != NULL);

pd->face[16].edge[0] = 28;
pd->face[16].or[0] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}

pd_code_t *pd_create_r1_testD_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(14);
assert(pd != NULL);
pd->ncross = 14;
pd->nedges = 28;
pd->ncomps = 3;
pd->nfaces = 16;
sprintf(pd->hash,"%s","DhwQBwUFBAQEAwMDAwMDAwICAgMMDAQ");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 7;
pd->cross[0].edge[2] = 11;
pd->cross[0].edge[3] = 6;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 8;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 7;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 20;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 21;
pd->cross[2].sign = 0;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 9;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 10;
pd->cross[3].sign = 1;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 26;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 25;
pd->cross[4].sign = 0;

pd->cross[5].edge[0] = 4;
pd->cross[5].edge[1] = 18;
pd->cross[5].edge[2] = 5;
pd->cross[5].edge[3] = 17;
pd->cross[5].sign = 1;

pd->cross[6].edge[0] = 5;
pd->cross[6].edge[1] = 27;
pd->cross[6].edge[2] = 6;
pd->cross[6].edge[3] = 24;
pd->cross[6].sign = 1;

pd->cross[7].edge[0] = 8;
pd->cross[7].edge[1] = 19;
pd->cross[7].edge[2] = 9;
pd->cross[7].edge[3] = 20;
pd->cross[7].sign = 0;

pd->cross[8].edge[0] = 10;
pd->cross[8].edge[1] = 16;
pd->cross[8].edge[2] = 11;
pd->cross[8].edge[3] = 15;
pd->cross[8].sign = 0;

pd->cross[9].edge[0] = 12;
pd->cross[9].edge[1] = 15;
pd->cross[9].edge[2] = 23;
pd->cross[9].edge[3] = 14;
pd->cross[9].sign = 1;

pd->cross[10].edge[0] = 12;
pd->cross[10].edge[1] = 22;
pd->cross[10].edge[2] = 13;
pd->cross[10].edge[3] = 21;
pd->cross[10].sign = 1;

pd->cross[11].edge[0] = 13;
pd->cross[11].edge[1] = 22;
pd->cross[11].edge[2] = 14;
pd->cross[11].edge[3] = 23;
pd->cross[11].sign = 1;

pd->cross[12].edge[0] = 16;
pd->cross[12].edge[1] = 25;
pd->cross[12].edge[2] = 17;
pd->cross[12].edge[3] = 24;
pd->cross[12].sign = 1;

pd->cross[13].edge[0] = 18;
pd->cross[13].edge[1] = 26;
pd->cross[13].edge[2] = 19;
pd->cross[13].edge[3] = 27;
pd->cross[13].sign = 1;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 5;
pd->edge[4].headpos = 0;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 6;
pd->edge[5].headpos = 0;
pd->edge[5].tail = 5;
pd->edge[5].tailpos = 2;

pd->edge[6].head = 0;
pd->edge[6].headpos = 3;
pd->edge[6].tail = 6;
pd->edge[6].tailpos = 2;

pd->edge[7].head = 1;
pd->edge[7].headpos = 3;
pd->edge[7].tail = 0;
pd->edge[7].tailpos = 1;

pd->edge[8].head = 7;
pd->edge[8].headpos = 0;
pd->edge[8].tail = 1;
pd->edge[8].tailpos = 1;

pd->edge[9].head = 3;
pd->edge[9].headpos = 1;
pd->edge[9].tail = 7;
pd->edge[9].tailpos = 2;

pd->edge[10].head = 8;
pd->edge[10].headpos = 0;
pd->edge[10].tail = 3;
pd->edge[10].tailpos = 3;

pd->edge[11].head = 0;
pd->edge[11].headpos = 2;
pd->edge[11].tail = 8;
pd->edge[11].tailpos = 2;

pd->edge[12].head = 10;
pd->edge[12].headpos = 0;
pd->edge[12].tail = 9;
pd->edge[12].tailpos = 0;

pd->edge[13].head = 11;
pd->edge[13].headpos = 0;
pd->edge[13].tail = 10;
pd->edge[13].tailpos = 2;

pd->edge[14].head = 9;
pd->edge[14].headpos = 3;
pd->edge[14].tail = 11;
pd->edge[14].tailpos = 2;

pd->edge[15].head = 8;
pd->edge[15].headpos = 3;
pd->edge[15].tail = 9;
pd->edge[15].tailpos = 1;

pd->edge[16].head = 12;
pd->edge[16].headpos = 0;
pd->edge[16].tail = 8;
pd->edge[16].tailpos = 1;

pd->edge[17].head = 5;
pd->edge[17].headpos = 3;
pd->edge[17].tail = 12;
pd->edge[17].tailpos = 2;

pd->edge[18].head = 13;
pd->edge[18].headpos = 0;
pd->edge[18].tail = 5;
pd->edge[18].tailpos = 1;

pd->edge[19].head = 7;
pd->edge[19].headpos = 1;
pd->edge[19].tail = 13;
pd->edge[19].tailpos = 2;

pd->edge[20].head = 2;
pd->edge[20].headpos = 1;
pd->edge[20].tail = 7;
pd->edge[20].tailpos = 3;

pd->edge[21].head = 10;
pd->edge[21].headpos = 3;
pd->edge[21].tail = 2;
pd->edge[21].tailpos = 3;

pd->edge[22].head = 11;
pd->edge[22].headpos = 1;
pd->edge[22].tail = 10;
pd->edge[22].tailpos = 1;

pd->edge[23].head = 9;
pd->edge[23].headpos = 2;
pd->edge[23].tail = 11;
pd->edge[23].tailpos = 3;

pd->edge[24].head = 12;
pd->edge[24].headpos = 3;
pd->edge[24].tail = 6;
pd->edge[24].tailpos = 3;

pd->edge[25].head = 4;
pd->edge[25].headpos = 3;
pd->edge[25].tail = 12;
pd->edge[25].tailpos = 1;

pd->edge[26].head = 13;
pd->edge[26].headpos = 1;
pd->edge[26].tail = 4;
pd->edge[26].tailpos = 1;

pd->edge[27].head = 6;
pd->edge[27].headpos = 1;
pd->edge[27].tail = 13;
pd->edge[27].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 12;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;
pd->comp[0].edge[10] = 10;
pd->comp[0].edge[11] = 11;

pd->comp[1].nedges = 12;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 12;
pd->comp[1].edge[1] = 13;
pd->comp[1].edge[2] = 14;
pd->comp[1].edge[3] = 15;
pd->comp[1].edge[4] = 16;
pd->comp[1].edge[5] = 17;
pd->comp[1].edge[6] = 18;
pd->comp[1].edge[7] = 19;
pd->comp[1].edge[8] = 20;
pd->comp[1].edge[9] = 21;
pd->comp[1].edge[10] = 22;
pd->comp[1].edge[11] = 23;

pd->comp[2].nedges = 4;
pd->comp[2].tag = 'C';

pd->comp[2].edge = calloc(pd->comp[2].nedges,sizeof(pd_idx_t));
assert(pd->comp[2].edge != NULL);

pd->comp[2].edge[0] = 24;
pd->comp[2].edge[1] = 25;
pd->comp[2].edge[2] = 26;
pd->comp[2].edge[3] = 27;


/* Face data */

pd->face[0].nedges = 7;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 21;
pd->face[0].or[1] = 1;

pd->face[0].edge[2] = 13;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 23;
pd->face[0].or[3] = 1;

pd->face[0].edge[4] = 15;
pd->face[0].or[4] = 1;

pd->face[0].edge[5] = 11;
pd->face[0].or[5] = 1;

pd->face[0].edge[6] = 7;
pd->face[0].or[6] = 1;

pd->face[1].nedges = 5;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 0;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 6;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 27;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 19;
pd->face[1].or[3] = 1;

pd->face[1].edge[4] = 8;
pd->face[1].or[4] = 0;

pd->face[2].nedges = 5;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 2;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 10;
pd->face[2].or[1] = 1;

pd->face[2].edge[2] = 15;
pd->face[2].or[2] = 0;

pd->face[2].edge[3] = 12;
pd->face[2].or[3] = 1;

pd->face[2].edge[4] = 21;
pd->face[2].or[4] = 0;

pd->face[3].nedges = 4;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 3;
pd->face[3].or[0] = 0;

pd->face[3].edge[1] = 9;
pd->face[3].or[1] = 0;

pd->face[3].edge[2] = 19;
pd->face[3].or[2] = 0;

pd->face[3].edge[3] = 26;
pd->face[3].or[3] = 0;

pd->face[4].nedges = 4;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 3;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 25;
pd->face[4].or[1] = 0;

pd->face[4].edge[2] = 16;
pd->face[4].or[2] = 0;

pd->face[4].edge[3] = 10;
pd->face[4].or[3] = 0;

pd->face[5].nedges = 4;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 6;
pd->face[5].or[0] = 1;

pd->face[5].edge[1] = 11;
pd->face[5].or[1] = 0;

pd->face[5].edge[2] = 16;
pd->face[5].or[2] = 1;

pd->face[5].edge[3] = 24;
pd->face[5].or[3] = 0;

pd->face[6].nedges = 3;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 1;
pd->face[6].or[0] = 0;

pd->face[6].edge[1] = 8;
pd->face[6].or[1] = 1;

pd->face[6].edge[2] = 20;
pd->face[6].or[2] = 1;

pd->face[7].nedges = 3;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 2;
pd->face[7].or[0] = 0;

pd->face[7].edge[1] = 20;
pd->face[7].or[1] = 0;

pd->face[7].edge[2] = 9;
pd->face[7].or[2] = 1;

pd->face[8].nedges = 3;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].or = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].or != NULL);

pd->face[8].edge[0] = 4;
pd->face[8].or[0] = 1;

pd->face[8].edge[1] = 17;
pd->face[8].or[1] = 0;

pd->face[8].edge[2] = 25;
pd->face[8].or[2] = 1;

pd->face[9].nedges = 3;
pd->face[9].edge = calloc(pd->face[9].nedges,sizeof(pd_idx_t));
pd->face[9].or = calloc(pd->face[9].nedges,sizeof(pd_or_t));
assert(pd->face[9].edge != NULL);
assert(pd->face[9].or != NULL);

pd->face[9].edge[0] = 4;
pd->face[9].or[0] = 0;

pd->face[9].edge[1] = 26;
pd->face[9].or[1] = 1;

pd->face[9].edge[2] = 18;
pd->face[9].or[2] = 0;

pd->face[10].nedges = 3;
pd->face[10].edge = calloc(pd->face[10].nedges,sizeof(pd_idx_t));
pd->face[10].or = calloc(pd->face[10].nedges,sizeof(pd_or_t));
assert(pd->face[10].edge != NULL);
assert(pd->face[10].or != NULL);

pd->face[10].edge[0] = 5;
pd->face[10].or[0] = 0;

pd->face[10].edge[1] = 18;
pd->face[10].or[1] = 1;

pd->face[10].edge[2] = 27;
pd->face[10].or[2] = 1;

pd->face[11].nedges = 3;
pd->face[11].edge = calloc(pd->face[11].nedges,sizeof(pd_idx_t));
pd->face[11].or = calloc(pd->face[11].nedges,sizeof(pd_or_t));
assert(pd->face[11].edge != NULL);
assert(pd->face[11].or != NULL);

pd->face[11].edge[0] = 5;
pd->face[11].or[0] = 1;

pd->face[11].edge[1] = 24;
pd->face[11].or[1] = 1;

pd->face[11].edge[2] = 17;
pd->face[11].or[2] = 1;

pd->face[12].nedges = 3;
pd->face[12].edge = calloc(pd->face[12].nedges,sizeof(pd_idx_t));
pd->face[12].or = calloc(pd->face[12].nedges,sizeof(pd_or_t));
assert(pd->face[12].edge != NULL);
assert(pd->face[12].or != NULL);

pd->face[12].edge[0] = 12;
pd->face[12].or[0] = 0;

pd->face[12].edge[1] = 14;
pd->face[12].or[1] = 0;

pd->face[12].edge[2] = 22;
pd->face[12].or[2] = 0;

pd->face[13].nedges = 2;
pd->face[13].edge = calloc(pd->face[13].nedges,sizeof(pd_idx_t));
pd->face[13].or = calloc(pd->face[13].nedges,sizeof(pd_or_t));
assert(pd->face[13].edge != NULL);
assert(pd->face[13].or != NULL);

pd->face[13].edge[0] = 0;
pd->face[13].or[0] = 1;

pd->face[13].edge[1] = 7;
pd->face[13].or[1] = 0;

pd->face[14].nedges = 2;
pd->face[14].edge = calloc(pd->face[14].nedges,sizeof(pd_idx_t));
pd->face[14].or = calloc(pd->face[14].nedges,sizeof(pd_or_t));
assert(pd->face[14].edge != NULL);
assert(pd->face[14].or != NULL);

pd->face[14].edge[0] = 13;
pd->face[14].or[0] = 0;

pd->face[14].edge[1] = 22;
pd->face[14].or[1] = 1;

pd->face[15].nedges = 2;
pd->face[15].edge = calloc(pd->face[15].nedges,sizeof(pd_idx_t));
pd->face[15].or = calloc(pd->face[15].nedges,sizeof(pd_or_t));
assert(pd->face[15].edge != NULL);
assert(pd->face[15].or != NULL);

pd->face[15].edge[0] = 14;
pd->face[15].or[0] = 1;

pd->face[15].edge[1] = 23;
pd->face[15].or[1] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}

/*   __23___>_______                           _________>___11___________     */
/* ____(15)_14____ /_______15__>_____________ / ___>____16_________     |     */
/* |  \ 11        /9                        8/                    /     |     */
/* |   \  (12)   /                          /       (4)          /      |     */
/* |    \       /                          /                    /  (6)  |     */
/* |     ^    12            (2)           ^     _____25__<_________     |     */
/* |     22   /                          /     |              /12 |     |     */
/* |       \ /                         10      |     (9)     /   24     |    */
/* ^ (14)   \                          /       |           17 (11)|     |    */
/* |       / \                _____________3________4_>__  /__5______6_____>_ */
/* |_13___/ 10\                \      3\      4|         5/      6/    0|   | */
/*             \                \       \      |  (10)  18  (5) 29      |   | */
/*              \                \       \     |26_>___ / _27__ /       /   | */
/*               \                \  (8)  \            /13     \14     /    | */
/*                \                \       ^   (3)    /      /  \     /     | */
/*                 \                \       9        /      /(16)\   /      | */
/*                  \                \       \     19      /__28__\ /       | */
/*                   \                ^       \    /               /        | */
/*                    ^                2       \  /      (0)      0         | */
/*                     \                \      7 /               /          | */
/*                     21                \      / ____<__8_     /    (13)   | */
/*                       \                \    /           \   /            | */
/*                        \                \ 20             \ /             | */
/*                         \                 /     (7)       \              | */
/*          (1)             \               /               /1\             | */
/*                           \_____________/2 \            /   \            | */
/*                                             \          /     \           | */
/*                                              \___<__1_/       \____<__7__| */
/*                                                                            */
/*                                                                            */
/*R1 operation at vertex 14													  */

/*   __23___>_______                           _________>___11___________     */
/* ___(15)__14____ /_______15__>_____________ / ___>____16_________     |     */
/* |  \ 11        /9                        8/                    /     |     */
/* |   \  (12)   /                          /      (4)           /      |     */
/* |    \       /                          /                    /       |     */
/* |     ^    12          (2)             ^     _____25__<_________ (5) |     */
/* |     22   /                          /     |              /12 |     |     */
/* |       \ /                         10      |     (8)     /   24     |    */
/* ^  (14)  \                          /       |           17 (11)|     |    */
/* |       / \                _____________3________4_>__  /__5______6_____>_ */
/* |_13___/ 10\                \      3\      4|         5/      6|    0|   | */
/*             \                \       \      |  (9)   18  (10)  |     |   | */
/*              \                \       \     |26_>___ / ___27_>_|     /   | */
/*               \                \  (7)  \            /13             /    | */
/*                \                \       ^    (3)   /               /     | */
/*                 \                \       9        /               /      | */
/*                  \                \       \     19      (1)      /       | */
/*                   \                ^       \    /               /        | */
/*                    ^                2       \  /               0         | */
/*                     \                \      7 /               /          | */
/*                     21                \      / ____<__8_     /           | */
/*                       \                \    /           \   /     (13)   | */
/*                        \                \ 20             \ /             | */
/*       (0)               \                 /               \              | */
/*                          \               /        (6)    /1\             | */
/*                           \_____________/2 \            /   \            | */
/*                                             \          /     \           | */
/*                                              \___<__1_/       \____<__7__| */
/*                                                                            */
/*                                                                            */
bool r1testD() {
  
  printf("---------------------------------------\n"
	 "r1 test D\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r1_testD_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after configuration...");
  pd_code_t *after = pd_create_r1_testD_after_0();
  printf("done (passes pd_ok)\n");
  
  printf("performing r1 loop deletion at crossing 14...");
  pd_code_t *newpd = pd_R1_loopdeletion(before,14);
  if (!pd_ok(newpd)) { 
    printf("fail (output pd does not pass pd_ok)\n");
    return false;
  }
  printf("pass (output passes pd_ok)\n");

  printf("testing for isomorphism with after configuration...");
  if (!pd_isomorphic(newpd,after)) { 
    pd_printf("fail (output pd %PD ",newpd);
    pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after);
    return false;
  }

  printf("pass (output and after isomorphic)\n");
  
  printf("housecleaning...");
  pd_code_free(&before);
  pd_code_free(&after);
  pd_code_free(&newpd);
  printf("done\n");

  printf("---------------------------------------\n"
	 "r1 test D : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}

pd_code_t *pd_create_r1_testE_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(4);
assert(pd != NULL);
pd->ncross = 4;
pd->nedges = 8;
pd->ncomps = 1;
pd->nfaces = 6;
sprintf(pd->hash,"%s","BAgGBQMDAgIBAQgAAAAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 2;
pd->cross[0].edge[2] = 7;
pd->cross[0].edge[3] = 3;
pd->cross[0].sign = 1;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 5;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 6;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 7;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 6;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 3;
pd->cross[3].edge[1] = 4;
pd->cross[3].edge[2] = 4;
pd->cross[3].edge[3] = 5;
pd->cross[3].sign = 1;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 0;
pd->edge[2].headpos = 1;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 3;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 0;
pd->edge[3].tailpos = 3;

pd->edge[4].head = 3;
pd->edge[4].headpos = 1;
pd->edge[4].tail = 3;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 1;
pd->edge[5].headpos = 1;
pd->edge[5].tail = 3;
pd->edge[5].tailpos = 3;

pd->edge[6].head = 2;
pd->edge[6].headpos = 3;
pd->edge[6].tail = 1;
pd->edge[6].tailpos = 3;

pd->edge[7].head = 0;
pd->edge[7].headpos = 2;
pd->edge[7].tail = 2;
pd->edge[7].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 8;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;


/* Face data */

pd->face[0].nedges = 5;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 5;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 4;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 3;
pd->face[0].or[3] = 0;

pd->face[0].edge[4] = 7;
pd->face[0].or[4] = 0;

pd->face[1].nedges = 3;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 0;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 3;
pd->face[1].or[1] = 1;

pd->face[1].edge[2] = 5;
pd->face[1].or[2] = 1;

pd->face[2].nedges = 3;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 6;
pd->face[2].or[1] = 1;

pd->face[2].edge[2] = 2;
pd->face[2].or[2] = 1;

pd->face[3].nedges = 2;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 1;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 6;
pd->face[3].or[1] = 0;

pd->face[4].nedges = 2;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 2;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 7;
pd->face[4].or[1] = 1;

pd->face[5].nedges = 1;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 4;
pd->face[5].or[0] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}

pd_code_t *pd_create_r1_testE_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(3);
assert(pd != NULL);
pd->ncross = 3;
pd->nedges = 6;
pd->ncomps = 1;
pd->nfaces = 5;
sprintf(pd->hash,"%s","AwYFAwMCAgIBBgAAAAAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 2;
pd->cross[0].edge[2] = 5;
pd->cross[0].edge[3] = 3;
pd->cross[0].sign = 1;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 3;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 4;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 5;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 4;
pd->cross[2].sign = 1;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 0;
pd->edge[2].headpos = 1;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 1;
pd->edge[3].headpos = 1;
pd->edge[3].tail = 0;
pd->edge[3].tailpos = 3;

pd->edge[4].head = 2;
pd->edge[4].headpos = 3;
pd->edge[4].tail = 1;
pd->edge[4].tailpos = 3;

pd->edge[5].head = 0;
pd->edge[5].headpos = 2;
pd->edge[5].tail = 2;
pd->edge[5].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 6;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;


/* Face data */

pd->face[0].nedges = 3;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 4;
pd->face[0].or[1] = 1;

pd->face[0].edge[2] = 2;
pd->face[0].or[2] = 1;

pd->face[1].nedges = 3;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 3;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 5;
pd->face[1].or[2] = 0;

pd->face[2].nedges = 2;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 3;
pd->face[2].or[1] = 1;

pd->face[3].nedges = 2;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 1;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 4;
pd->face[3].or[1] = 0;

pd->face[4].nedges = 2;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 2;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 5;
pd->face[4].or[1] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}

/*                                                                  (2)       */
/*                          __________________2_____>________                 */
/*                         \                                 |                */
/*                          \                                |                */
/*                           ^              (4)              |                */
/*                            \                              |                */
/*                   ____________________>____7_____________ | ___>____       */
/*                  |          2                            0|        |       */
/*                  |            \                           |        |       */
/*                  |             \          (0)             |        |       */
/*                  |              \                         /        |       */
/*                  |               \                       <         |       */
/*                  |                \                     /          |       */
/*                  |                 \          __4__    /           |       */
/*                  ^                  1        |     |  3            |       */
/*                  |                   \       ^     | /             0       */
/*                  |                    \      | (5)  /              |       */
/*                  |                     ^     |_____/               |       */
/*                  6        (3)           \           3      (1)     |       */
/*                  |                       \       /                 |       */
/*                  |                        \     5                  |       */
/*                  |                         \   /                   |       */
/*                  |                          \ /                    |       */
/*                  |                           \                     |       */
/*                  |                            \                    |       */
/*                  |__________<______________/ 1 \_________<__________       */
/*                                                                            */
/*R1 operation at vertex 3												      */
/*                                                                  (0)       */
/*                          __________________2_____>________                 */
/*                         \                                 |                */
/*                          \                                |                */
/*                           ^              (4)              |                */
/*                            \                              |                */
/*                   ____________________>____5_____________ | ___>____       */
/*                  |          2                            0|        |       */
/*                  |            \                           |        |       */
/*                  |             \          (1)             |        |       */
/*                  |              \                         /        |       */
/*                  |               \                       /         |       */
/*                  |                \                     /          |       */
/*                  |                 \                   /           |       */
/*                  ^                  1                3             |       */
/*                  |                   \               /             0       */
/*                  |                    \             /              |       */
/*                  |                     ^           <               |       */
/*                  4        (3)           \         /        (2)     |       */
/*                  |                       \       /                 |       */
/*                  |                        \     /                  |       */
/*                  |                         \   /                   |       */
/*                  |                          \ /                    |       */
/*                  |                           \                     |       */
/*                  |                            \                    |       */
/*                  |__________<______________/ 1 \_________<__________       */
/*                                                                            */

bool r1testE() {
  
  printf("---------------------------------------\n"
	 "r1 test E\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r1_testE_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after configuration...");
  pd_code_t *after = pd_create_r1_testE_after_0();
  printf("done (passes pd_ok)\n");
  
  printf("performing r1 loop deletion at crossing 3...");
  pd_code_t *newpd = pd_R1_loopdeletion(before,3);
  if (!pd_ok(newpd)) { 
    printf("fail (output pd does not pass pd_ok)\n");
    return false;
  }
  printf("pass (output passes pd_ok)\n");

  printf("testing for isomorphism with after configuration...");
  if (!pd_isomorphic(newpd,after)) { 
    pd_printf("fail (output pd %PD ",newpd);
    pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after);
    return false;
  }

  printf("pass (output and after isomorphic)\n");
  
  printf("housecleaning...");
  pd_code_free(&before);
  pd_code_free(&after);
  pd_code_free(&newpd);
  printf("done\n");

  printf("---------------------------------------\n"
	 "r1 test E : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}

/* r2 test A assembled c code */

pd_code_t *pd_create_r2_testA_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(6);
assert(pd != NULL);
pd->ncross = 6;
pd->nedges = 12;
pd->ncomps = 2;
pd->nfaces = 8;
sprintf(pd->hash,"%s","BgwIBQUDAwICAgICCgIAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 6;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 5;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 7;
pd->cross[1].edge[2] = 9;
pd->cross[1].edge[3] = 6;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 9;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 8;
pd->cross[2].sign = 0;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 7;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 8;
pd->cross[3].sign = 1;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 11;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 10;
pd->cross[4].sign = 0;

pd->cross[5].edge[0] = 4;
pd->cross[5].edge[1] = 11;
pd->cross[5].edge[2] = 5;
pd->cross[5].edge[3] = 10;
pd->cross[5].sign = 1;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 1;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 5;
pd->edge[4].headpos = 0;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 0;
pd->edge[5].headpos = 3;
pd->edge[5].tail = 5;
pd->edge[5].tailpos = 2;

pd->edge[6].head = 1;
pd->edge[6].headpos = 3;
pd->edge[6].tail = 0;
pd->edge[6].tailpos = 1;

pd->edge[7].head = 3;
pd->edge[7].headpos = 1;
pd->edge[7].tail = 1;
pd->edge[7].tailpos = 1;

pd->edge[8].head = 2;
pd->edge[8].headpos = 3;
pd->edge[8].tail = 3;
pd->edge[8].tailpos = 3;

pd->edge[9].head = 1;
pd->edge[9].headpos = 2;
pd->edge[9].tail = 2;
pd->edge[9].tailpos = 1;

pd->edge[10].head = 4;
pd->edge[10].headpos = 3;
pd->edge[10].tail = 5;
pd->edge[10].tailpos = 3;

pd->edge[11].head = 5;
pd->edge[11].headpos = 1;
pd->edge[11].tail = 4;
pd->edge[11].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 10;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;

pd->comp[1].nedges = 2;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 10;
pd->comp[1].edge[1] = 11;


/* Face data */

pd->face[0].nedges = 5;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 5;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 11;
pd->face[0].or[2] = 0;

pd->face[0].edge[3] = 3;
pd->face[0].or[3] = 0;

pd->face[0].edge[4] = 7;
pd->face[0].or[4] = 0;

pd->face[1].nedges = 5;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 8;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 3;
pd->face[1].or[2] = 1;

pd->face[1].edge[3] = 10;
pd->face[1].or[3] = 0;

pd->face[1].edge[4] = 5;
pd->face[1].or[4] = 1;

pd->face[2].nedges = 3;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 1;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 6;
pd->face[2].or[1] = 1;

pd->face[2].edge[2] = 9;
pd->face[2].or[2] = 0;

pd->face[3].nedges = 3;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 2;
pd->face[3].or[0] = 0;

pd->face[3].edge[1] = 9;
pd->face[3].or[1] = 1;

pd->face[3].edge[2] = 7;
pd->face[3].or[2] = 1;

pd->face[4].nedges = 2;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 0;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 6;
pd->face[4].or[1] = 0;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 2;
pd->face[5].or[0] = 1;

pd->face[5].edge[1] = 8;
pd->face[5].or[1] = 1;

pd->face[6].nedges = 2;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 4;
pd->face[6].or[0] = 1;

pd->face[6].edge[1] = 10;
pd->face[6].or[1] = 1;

pd->face[7].nedges = 2;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 4;
pd->face[7].or[0] = 0;

pd->face[7].edge[1] = 11;
pd->face[7].or[1] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



pd_code_t *pd_create_r2_testA_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(4);
assert(pd != NULL);
pd->ncross = 4;
pd->nedges = 8;
pd->ncomps = 1;
pd->nfaces = 6;
sprintf(pd->hash,"%s","BAgGAwMDAwICAQgAAAAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 4;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 3;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 5;
pd->cross[1].edge[2] = 7;
pd->cross[1].edge[3] = 4;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 7;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 6;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 5;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 6;
pd->cross[3].sign = 1;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 1;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 0;
pd->edge[3].headpos = 3;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 1;
pd->edge[4].headpos = 3;
pd->edge[4].tail = 0;
pd->edge[4].tailpos = 1;

pd->edge[5].head = 3;
pd->edge[5].headpos = 1;
pd->edge[5].tail = 1;
pd->edge[5].tailpos = 1;

pd->edge[6].head = 2;
pd->edge[6].headpos = 3;
pd->edge[6].tail = 3;
pd->edge[6].tailpos = 3;

pd->edge[7].head = 1;
pd->edge[7].headpos = 2;
pd->edge[7].tail = 2;
pd->edge[7].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 8;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;


/* Face data */

pd->face[0].nedges = 3;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 3;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 5;
pd->face[0].or[2] = 0;

pd->face[1].nedges = 3;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 4;
pd->face[1].or[1] = 1;

pd->face[1].edge[2] = 7;
pd->face[1].or[2] = 0;

pd->face[2].nedges = 3;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 1;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 6;
pd->face[2].or[1] = 0;

pd->face[2].edge[2] = 3;
pd->face[2].or[2] = 1;

pd->face[3].nedges = 3;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 2;
pd->face[3].or[0] = 0;

pd->face[3].edge[1] = 7;
pd->face[3].or[1] = 1;

pd->face[3].edge[2] = 5;
pd->face[3].or[2] = 1;

pd->face[4].nedges = 2;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 0;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 4;
pd->face[4].or[1] = 0;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 2;
pd->face[5].or[0] = 1;

pd->face[5].edge[1] = 6;
pd->face[5].or[1] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}

pd_code_t *pd_create_r2_testA_after_1() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(2);
assert(pd != NULL);
pd->ncross = 0;
pd->nedges = 1;
pd->ncomps = 1;
pd->nfaces = 2;
sprintf(pd->hash,"%s","AAECAQEBAQAAAAAAAAAAAAAAAAAAAAA");

/* Crossing data. */


/* Edge data */

pd->edge[0].head = -1;
pd->edge[0].headpos = 4;
pd->edge[0].tail = -1;
pd->edge[0].tailpos = 4;


/* Component Data */

pd->comp[0].nedges = 1;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;


/* Face data */

pd->face[0].nedges = 1;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 1;

pd->face[1].nedges = 1;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 0;
pd->face[1].or[0] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



/*                                                                            */
/*                                                                            */
/*                                                                            */
/*       ----0->---|         __________<_______________    (0)                */
/*      1\ (4)     |         |                        |                       */
/*    ______<_6___ | __<_5________<_                  |                       */
/*   /             |0       5      |                  |                       */
/*  /      \  (2)  /         |     |                  |                       */
/* |        ^     1          |     |                  |                       */
/* |         9   <           |     |                  |                       */
/* |          \ /            |     |                  |                       */
/* |           \ 2          10     ^      (7)         |                       */
/* | (3)        \            | (6) |                  |                       */
/* |         /   \           |     4                  |                       */
/* |        /     \          |     |                  |                       */
/* |       /      /          |     |                  |                       */
/* |      2      /           |     |                  |                       */
/* |     /      8            |     /                  |                       */
/* |     | (5) /             |    /                   |                       */
/* \     |    ^    (1)       |   /                    |                       */
/*  7    \   /               |  /                     |                       */
/*   \    \ /                  /                      |                       */
/*    \    \                  /                       |                       */
/*     \  / \____3_>_________/                        |                       */
/*      \/ 3               4                          |                       */
/*                           |                        |                       */
/*                           |_________11___>_________|                       */
/*                                                                            */
/*                                                                            */
/* R2 on face 6 around crossings 4 and 5. Get 2 pds.                          */
/*                                                                            */
/*                                                                            */
/*                                (0)                                         */
/*       ----0->---|       				__________<_______________    */
/*       \ (4)     |         				|                        |    */
/*    ______<_4___ | ______________ 	                |                        |    */
/*   /   1         |0              |    		|                        |    */
/*  /      \  (1)  /               |    		|                        |    */
/* |        ^     1                |    		|                        |    */
/* |         7   <                 |    		|                        |    */
/* |          \ /                  |    		|                        |    */
/* |           \ 2                 ^    		|  (pd_unknot(0))        |    */
/* | (3)        \          (2)     |    		|                        |    */
/* |         /   \                 |    		|                        |    */
/* |        /     \                |    		|                        |    */
/* |       /      /                |    		|                        |    */
/* |      2      /                 |    		|                        |    */
/* |     /      6                  /    		|                        |    */
/* |     | (5) /                  /     		|                        |    */
/* \     |    ^                  /      		|                        |    */
/*  5    \   /                  /       		|                        |    */
/*   \    \ /                  /        		|                        |    */
/*    \    \                  /         		|                        |    */
/*     \  / \____3_>_________/          		|                        |    */
/*      \/ 3                            		|                        |    */
/*                                     			|                        |    */
/*                            				|______________>_________|    */
/*                                                                            */
/*                                                                            */
/* R2 on face 6 around crossings 4 and 5                                      */


bool r2testA() {
  
  printf("---------------------------------------\n"
	 "r2 test A\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r2_testA_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after-0 configuration...");
  pd_code_t *after0 = pd_create_r2_testA_after_0();
  printf("done (passes pd_ok)\n");

  printf("creating after-1 configuration...");
  pd_code_t *after1 = pd_create_r2_testA_after_1();
  printf("done (passes pd_ok)\n");

  pd_code_t *after[2] = {after0,after1};
    
  printf("performing r2 at crossings 4 and 5...\n");

  pd_idx_t  noutpd;
  pd_code_t **newpd;
  pd_idx_t  cr[2] = {4,5};

  pd_R2_bigon_elimination(before,cr,&noutpd,&newpd);

  if (noutpd != 2) { 
    printf("fail (didn't generate 2 output pd codes)\n");
    return false;
  }

  pd_idx_t ct;

  for(ct=0;ct<noutpd;ct++) { 

    if (!pd_ok(newpd[ct])) { 
      printf("fail (output pd %d does not pass pd_ok)\n",ct);
      return false;
    }
    printf("\t output pd %d passes pd_ok...pass\n",ct);
    printf("\t testing for isomorphism with after configuration %d...",ct);
    if (!pd_isomorphic(newpd[ct],after[ct])) { 
      pd_printf("fail (output pd %PD ",newpd[ct]);
      pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after[ct]);
      return false;
    }

    printf("pass (output %d and after %d isomorphic)\n",ct,ct);

  }

  printf("performing r2 at crossings 4 and 5...pass\n");
  printf("housecleaning...");

  pd_code_free(&before);
  pd_code_free(&after0);
  pd_code_free(&after1);
  pd_code_free(&newpd[0]);
  pd_code_free(&newpd[1]);

  free(newpd);

  printf("done\n");

  printf("---------------------------------------\n"
	 "r2 test A : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}

/* r2 test B assembled c code */

pd_code_t *pd_create_r2_testB_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(6);
assert(pd != NULL);
pd->ncross = 6;
pd->nedges = 12;
pd->ncomps = 1;
pd->nfaces = 8;
sprintf(pd->hash,"%s","BgwIBQQDAwMCAgIBDAAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 6;
pd->cross[0].edge[2] = 11;
pd->cross[0].edge[3] = 7;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 7;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 8;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 5;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 4;
pd->cross[2].sign = 0;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 10;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 9;
pd->cross[3].sign = 1;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 8;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 9;
pd->cross[4].sign = 0;

pd->cross[5].edge[0] = 5;
pd->cross[5].edge[1] = 11;
pd->cross[5].edge[2] = 6;
pd->cross[5].edge[3] = 10;
pd->cross[5].sign = 0;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 2;
pd->edge[4].headpos = 3;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 5;
pd->edge[5].headpos = 0;
pd->edge[5].tail = 2;
pd->edge[5].tailpos = 1;

pd->edge[6].head = 0;
pd->edge[6].headpos = 1;
pd->edge[6].tail = 5;
pd->edge[6].tailpos = 2;

pd->edge[7].head = 1;
pd->edge[7].headpos = 1;
pd->edge[7].tail = 0;
pd->edge[7].tailpos = 3;

pd->edge[8].head = 4;
pd->edge[8].headpos = 1;
pd->edge[8].tail = 1;
pd->edge[8].tailpos = 3;

pd->edge[9].head = 3;
pd->edge[9].headpos = 3;
pd->edge[9].tail = 4;
pd->edge[9].tailpos = 3;

pd->edge[10].head = 5;
pd->edge[10].headpos = 3;
pd->edge[10].tail = 3;
pd->edge[10].tailpos = 1;

pd->edge[11].head = 0;
pd->edge[11].headpos = 2;
pd->edge[11].tail = 5;
pd->edge[11].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 12;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;
pd->comp[0].edge[10] = 10;
pd->comp[0].edge[11] = 11;


/* Face data */

pd->face[0].nedges = 5;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 8;
pd->face[0].or[1] = 1;

pd->face[0].edge[2] = 3;
pd->face[0].or[2] = 0;

pd->face[0].edge[3] = 10;
pd->face[0].or[3] = 1;

pd->face[0].edge[4] = 6;
pd->face[0].or[4] = 1;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 7;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 11;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 5;
pd->face[1].or[3] = 0;

pd->face[2].nedges = 3;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 1;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 4;
pd->face[2].or[1] = 0;

pd->face[2].edge[2] = 8;
pd->face[2].or[2] = 0;

pd->face[3].nedges = 3;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 2;
pd->face[3].or[0] = 0;

pd->face[3].edge[1] = 5;
pd->face[3].or[1] = 1;

pd->face[3].edge[2] = 10;
pd->face[3].or[2] = 0;

pd->face[4].nedges = 3;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 2;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 9;
pd->face[4].or[1] = 0;

pd->face[4].edge[2] = 4;
pd->face[4].or[2] = 1;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 0;
pd->face[5].or[0] = 0;

pd->face[5].edge[1] = 7;
pd->face[5].or[1] = 1;

pd->face[6].nedges = 2;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 3;
pd->face[6].or[0] = 1;

pd->face[6].edge[1] = 9;
pd->face[6].or[1] = 1;

pd->face[7].nedges = 2;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 6;
pd->face[7].or[0] = 0;

pd->face[7].edge[1] = 11;
pd->face[7].or[1] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



pd_code_t *pd_create_r2_testB_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(4);
assert(pd != NULL);
pd->ncross = 4;
pd->nedges = 8;
pd->ncomps = 1;
pd->nfaces = 6;
sprintf(pd->hash,"%s","BAgGBAQDAgIBAQgAAAAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 4;
pd->cross[0].edge[2] = 7;
pd->cross[0].edge[3] = 5;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 5;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 6;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 3;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 2;
pd->cross[2].sign = 0;

pd->cross[3].edge[0] = 3;
pd->cross[3].edge[1] = 7;
pd->cross[3].edge[2] = 4;
pd->cross[3].edge[3] = 6;
pd->cross[3].sign = 0;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 2;
pd->edge[2].headpos = 3;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 3;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 2;
pd->edge[3].tailpos = 1;

pd->edge[4].head = 0;
pd->edge[4].headpos = 1;
pd->edge[4].tail = 3;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 1;
pd->edge[5].headpos = 1;
pd->edge[5].tail = 0;
pd->edge[5].tailpos = 3;

pd->edge[6].head = 3;
pd->edge[6].headpos = 3;
pd->edge[6].tail = 1;
pd->edge[6].tailpos = 3;

pd->edge[7].head = 0;
pd->edge[7].headpos = 2;
pd->edge[7].tail = 3;
pd->edge[7].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 8;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;


/* Face data */

pd->face[0].nedges = 4;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 2;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 3;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 6;
pd->face[0].or[3] = 0;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 5;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 7;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 3;
pd->face[1].or[3] = 0;

pd->face[2].nedges = 3;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 6;
pd->face[2].or[1] = 1;

pd->face[2].edge[2] = 4;
pd->face[2].or[2] = 1;

pd->face[3].nedges = 2;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 0;
pd->face[3].or[0] = 0;

pd->face[3].edge[1] = 5;
pd->face[3].or[1] = 1;

pd->face[4].nedges = 2;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 4;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 7;
pd->face[4].or[1] = 1;

pd->face[5].nedges = 1;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 2;
pd->face[5].or[0] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



/*                                                                            */
/*                   ______     ___0___                                       */
/*                  /      \   ^       \            (0)                       */
/*                 /        \0/         \                                     */
/*                6          /           \                                    */
/*               /    (7)   / \     (5)   \                                   */
/*               \         /   \          /                                   */
/*                \       /     7        /                                    */
/*                 ^     11      \      /                                     */
/*                  \   /         >    /                                      */
/*                   \ /           \  /                                       */
/*                    \5     (1)    \                                         */
/*                   / \             \                                        */
/*                  /   ^          /1 \                                       */
/*                 /     \        /    \                                      */
/*                /       5      1      \                                     */
/*               /         \    /        \                                    */
/*              ^           \2 <          >                                   */
/*             /              /            \                                  */
/*            /     (3)      /       (2)    \                                 */
/*           /              2  \             \                                */
/*          /              /(4) 4             \                               */
/*         /____10____<_______9________<___8___\                              */
/*                        3      4                                            */
/*                      /          \                                          */
/*                      \    (6)   /                                          */
/*                       \        /                                           */
/*                        \      /                                            */
/*                         \_3_>/                                             */
/*                                                                            */
/* R2 on face 6, edges 3 and 9, on vertices 3 and 4                           */
/*                                                                            */
/*                   ______     ___0_>_                                       */
/*                  /      \   /       \            (2)                       */
/*                 ^        \0/         \                                     */
/*                4          /           \                                    */
/*               /    (4)   / \     (3)   \                                   */
/*               \         /   \          /                                   */
/*                \       /     5        /                                    */
/*                 \     7       >      /                                     */
/*                  \   /         \    /                                      */
/*                   \ /           \  /                                       */
/*                    \3     (1)    \                                         */
/*                   / \             \                                        */
/*                  /   \          /1 \                                       */
/*                 /     ^        /    \                                      */
/*                /       3      1      \                                     */
/*               /         \    <        \                                    */
/*              /           \2 /          \                                   */
/*             /              /            \                                  */
/*            /              /              \                                 */
/*           /              /  \             \                                */
/*          /              /    \             \                               */
/*         /              /      \             \                              */
/*        |              /        \            |                              */
/*        |             |          \           |                              */
/*        |              \    (5)   |          |                              */
/*        |               \        /           |                              */
/*        |                \      /            |                              */
/*        |                 \_2_>/             |                              */
/*        |                                    |                              */
/*        |                 (0)                |                              */
/*        |                                    |                              */
/*        |______________<_____6_______________|                              */
/*                                                                            */
/*                                                                            */


    bool r2testB() {
  
  printf("---------------------------------------\n"
	 "r2 test B\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r2_testB_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after-0 configuration...");
  pd_code_t *after0 = pd_create_r2_testB_after_0();
  printf("done (passes pd_ok)\n");

  pd_code_t *after[1] = {after0};
    
  printf("performing r2 at crossings 3 and 4...\n");

  pd_idx_t  noutpd;
  pd_code_t **newpd;
  pd_idx_t  cr[2] = {3,4};

  pd_R2_bigon_elimination(before,cr,&noutpd,&newpd);

  if (noutpd != 1) { 
    printf("fail (did not generate 1 output pd codes)\n");
    return false;
  }

  pd_idx_t ct;

  for(ct=0;ct<noutpd;ct++) { 

    if (!pd_ok(newpd[ct])) { 
      printf("fail (output pd %d does not pass pd_ok)\n",ct);
      return false;
    }
    printf("\t output pd %d passes pd_ok...pass\n",ct);
    printf("\t testing for isomorphism with after configuration %d...",ct);
    if (!pd_isomorphic(newpd[ct],after[ct])) { 
      pd_printf("fail (output pd %PD ",newpd[ct]);
      pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after[ct]);
      return false;
    }

    printf("pass (output %d and after %d isomorphic)\n",ct,ct);

  }

  printf("performing r2 at crossings 3 and 4...pass\n");
  printf("housecleaning...");

  pd_code_free(&before);
  pd_code_free(&after0);
  pd_code_free(&newpd[0]);

  free(newpd);

  printf("done\n");

  printf("---------------------------------------\n"
	 "r2 test B : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}

/* r2 test C assembled c code */

pd_code_t *pd_create_r2_testC_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(6);
assert(pd != NULL);
pd->ncross = 6;
pd->nedges = 12;
pd->ncomps = 2;
pd->nfaces = 8;
sprintf(pd->hash,"%s","BgwIBgYCAgICAgICBgYAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 6;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 7;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 11;
pd->cross[1].edge[2] = 5;
pd->cross[1].edge[3] = 6;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 8;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 7;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 8;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 9;
pd->cross[3].sign = 1;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 10;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 9;
pd->cross[4].sign = 1;

pd->cross[5].edge[0] = 4;
pd->cross[5].edge[1] = 10;
pd->cross[5].edge[2] = 5;
pd->cross[5].edge[3] = 11;
pd->cross[5].sign = 1;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 1;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 5;
pd->edge[4].headpos = 0;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 1;
pd->edge[5].headpos = 2;
pd->edge[5].tail = 5;
pd->edge[5].tailpos = 2;

pd->edge[6].head = 0;
pd->edge[6].headpos = 1;
pd->edge[6].tail = 1;
pd->edge[6].tailpos = 3;

pd->edge[7].head = 2;
pd->edge[7].headpos = 3;
pd->edge[7].tail = 0;
pd->edge[7].tailpos = 3;

pd->edge[8].head = 3;
pd->edge[8].headpos = 1;
pd->edge[8].tail = 2;
pd->edge[8].tailpos = 1;

pd->edge[9].head = 4;
pd->edge[9].headpos = 3;
pd->edge[9].tail = 3;
pd->edge[9].tailpos = 3;

pd->edge[10].head = 5;
pd->edge[10].headpos = 1;
pd->edge[10].tail = 4;
pd->edge[10].tailpos = 1;

pd->edge[11].head = 1;
pd->edge[11].headpos = 1;
pd->edge[11].tail = 5;
pd->edge[11].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 6;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;

pd->comp[1].nedges = 6;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 6;
pd->comp[1].edge[1] = 7;
pd->comp[1].edge[2] = 8;
pd->comp[1].edge[3] = 9;
pd->comp[1].edge[4] = 10;
pd->comp[1].edge[5] = 11;


/* Face data */

pd->face[0].nedges = 6;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 7;
pd->face[0].or[1] = 1;

pd->face[0].edge[2] = 2;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 9;
pd->face[0].or[3] = 1;

pd->face[0].edge[4] = 4;
pd->face[0].or[4] = 1;

pd->face[0].edge[5] = 11;
pd->face[0].or[5] = 1;

pd->face[1].nedges = 6;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 6;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 5;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 10;
pd->face[1].or[3] = 0;

pd->face[1].edge[4] = 3;
pd->face[1].or[4] = 0;

pd->face[1].edge[5] = 8;
pd->face[1].or[5] = 0;

pd->face[2].nedges = 2;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 6;
pd->face[2].or[1] = 1;

pd->face[3].nedges = 2;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 1;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 7;
pd->face[3].or[1] = 0;

pd->face[4].nedges = 2;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 2;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 8;
pd->face[4].or[1] = 1;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 3;
pd->face[5].or[0] = 1;

pd->face[5].edge[1] = 9;
pd->face[5].or[1] = 0;

pd->face[6].nedges = 2;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 4;
pd->face[6].or[0] = 0;

pd->face[6].edge[1] = 10;
pd->face[6].or[1] = 1;

pd->face[7].nedges = 2;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 5;
pd->face[7].or[0] = 1;

pd->face[7].edge[1] = 11;
pd->face[7].or[1] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



pd_code_t *pd_create_r2_testC_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(4);
assert(pd != NULL);
pd->ncross = 4;
pd->nedges = 8;
pd->ncomps = 2;
pd->nfaces = 6;
sprintf(pd->hash,"%s","BAgGBAQCAgICAgQEAAAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 4;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 5;
pd->cross[0].sign = 1;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 7;
pd->cross[1].edge[2] = 3;
pd->cross[1].edge[3] = 4;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 6;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 5;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 6;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 7;
pd->cross[3].sign = 1;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 1;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 1;
pd->edge[3].headpos = 2;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 0;
pd->edge[4].headpos = 1;
pd->edge[4].tail = 1;
pd->edge[4].tailpos = 3;

pd->edge[5].head = 2;
pd->edge[5].headpos = 3;
pd->edge[5].tail = 0;
pd->edge[5].tailpos = 3;

pd->edge[6].head = 3;
pd->edge[6].headpos = 1;
pd->edge[6].tail = 2;
pd->edge[6].tailpos = 1;

pd->edge[7].head = 1;
pd->edge[7].headpos = 1;
pd->edge[7].tail = 3;
pd->edge[7].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 4;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;

pd->comp[1].nedges = 4;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 4;
pd->comp[1].edge[1] = 5;
pd->comp[1].edge[2] = 6;
pd->comp[1].edge[3] = 7;


/* Face data */

pd->face[0].nedges = 4;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 5;
pd->face[0].or[1] = 1;

pd->face[0].edge[2] = 2;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 7;
pd->face[0].or[3] = 1;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 4;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 3;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 6;
pd->face[1].or[3] = 0;

pd->face[2].nedges = 2;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 4;
pd->face[2].or[1] = 1;

pd->face[3].nedges = 2;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 1;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 5;
pd->face[3].or[1] = 0;

pd->face[4].nedges = 2;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 2;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 6;
pd->face[4].or[1] = 1;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 3;
pd->face[5].or[0] = 1;

pd->face[5].edge[1] = 7;
pd->face[5].or[1] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



/*                                                                            */
/*                                                   (0)                      */
/*                                                                            */
/*                                                                            */
/*      _____  --0->--            __7_>_                                      */
/*     /     \ 1      \          /     |                                      */
/*    ^       \        \        /      |                                      */
/*  11      /  \  (2)   \      /       |                                      */
/*  / (7)  ^    \        \    /        |                                      */
/*  \     5      \__6_>______/         |                                      */
/*   \   /                0\           |                                      */
/*    \ /                   \    (3)   |                                      */
/*    5/                     \         |                                      */
/*    /                       1        |                                      */
/*   /  \                      \_>_____|                                      */
/*  /    ^                            2|\                                     */
/* \     10         (1)                | \                                    */
/*  \      \                           |  \                                   */
/*   ^  (6) \                          |(4)|                                  */
/*    \      \                         |   2                                  */
/*     4     |                         8   |                                  */
/*      \____|___                      |   |                                  */
/*           |4  \                     |   /                                  */
/*           |    \                    |  /                                   */
/*           |     \                   | /                                    */
/*           |      \______<_3__________/                                     */
/*           \                          3                                     */
/*            \           (5)          |                                      */
/*             \                       /                                      */
/*              \_________<_9_________/                                       */
/*                                                                            */
/*  R2 on face 2                                                              */
/*                                                                            */
/*                                                                            */
/*      ________________7_>_____________                                      */
/*     /                               |                                      */
/*    /                                |                (0)                   */
/*   /                  (5)            |                                      */
/*  /                                  |                                      */
/*  \     ________3__>_____            |                                      */
/*   \   /                 \           |                                      */
/*    \ /                   \          |                                      */
/*    3/                     \         |                                      */
/*    /                       \        |                                      */
/*   /  \                      \_______|                                      */
/*  /    \                            1|\                                     */
/* \      ^                            | \                                    */
/*  \      6         (1)               |  \                                   */
/*   \  (4) \                          |   |                                  */
/*    ^      \                         4   0                                  */
/*     2     |                         |(2)|                                  */
/*      \____|___                      |   |                                  */
/*          2|   \                     |  <                                   */
/*           |    \                    | /                                    */
/*           |     \                    /                                     */
/*           |      \________<_1_______/                                      */
/*           \                          0                                     */
/*            \          (3)           |                                      */
/*             \                       /                                      */
/*              \_________<_5_________/                                       */



bool r2testC() {
  
  printf("---------------------------------------\n"
	 "r2 test C\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r2_testC_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after-0 configuration...");
  pd_code_t *after0 = pd_create_r2_testC_after_0();
  printf("done (passes pd_ok)\n");

  pd_code_t *after[1] = {after0};
    
  printf("performing r2 at crossings 1 and 0...\n");

  pd_idx_t  noutpd;
  pd_code_t **newpd;
  pd_idx_t  cr[2] = {1,0};

  pd_R2_bigon_elimination(before,cr,&noutpd,&newpd);

  if (noutpd != 1) { 
    printf("fail (did not generate 1 output pd codes)\n");
    return false;
  }

  pd_idx_t ct;

  for(ct=0;ct<noutpd;ct++) { 

    if (!pd_ok(newpd[ct])) { 
      printf("fail (output pd %d does not pass pd_ok)\n",ct);
      return false;
    }
    printf("\t output pd %d passes pd_ok...pass\n",ct);
    printf("\t testing for isomorphism with after configuration %d...",ct);
    if (!pd_isomorphic(newpd[ct],after[ct])) { 
      pd_printf("fail (output pd %PD ",newpd[ct]);
      pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after[ct]);
      return false;
    }

    printf("pass (output %d and after %d isomorphic)\n",ct,ct);

  }

  printf("performing r2 at crossings 1 and 0...pass\n");
  printf("housecleaning...");

  pd_code_free(&before);
  pd_code_free(&after0);
  pd_code_free(&newpd[0]);

  free(newpd);

  printf("done\n");

  printf("---------------------------------------\n"
	 "r2 test C : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}

/* r2 test D assembled c code */

pd_code_t *pd_create_r2_testD_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(12);
assert(pd != NULL);
pd->ncross = 12;
pd->nedges = 24;
pd->ncomps = 3;
pd->nfaces = 14;
sprintf(pd->hash,"%s","DBgOCAcGBAMDAwICAgICAgIDCgoEAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 4;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 3;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 20;
pd->cross[1].edge[2] = 9;
pd->cross[1].edge[3] = 23;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 4;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 5;
pd->cross[2].sign = 0;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 6;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 5;
pd->cross[3].sign = 0;

pd->cross[4].edge[0] = 6;
pd->cross[4].edge[1] = 16;
pd->cross[4].edge[2] = 7;
pd->cross[4].edge[3] = 17;
pd->cross[4].sign = 0;

pd->cross[5].edge[0] = 7;
pd->cross[5].edge[1] = 18;
pd->cross[5].edge[2] = 8;
pd->cross[5].edge[3] = 17;
pd->cross[5].sign = 0;

pd->cross[6].edge[0] = 8;
pd->cross[6].edge[1] = 21;
pd->cross[6].edge[2] = 9;
pd->cross[6].edge[3] = 20;
pd->cross[6].sign = 0;

pd->cross[7].edge[0] = 10;
pd->cross[7].edge[1] = 14;
pd->cross[7].edge[2] = 11;
pd->cross[7].edge[3] = 13;
pd->cross[7].sign = 0;

pd->cross[8].edge[0] = 10;
pd->cross[8].edge[1] = 23;
pd->cross[8].edge[2] = 19;
pd->cross[8].edge[3] = 22;
pd->cross[8].sign = 0;

pd->cross[9].edge[0] = 11;
pd->cross[9].edge[1] = 14;
pd->cross[9].edge[2] = 12;
pd->cross[9].edge[3] = 15;
pd->cross[9].sign = 0;

pd->cross[10].edge[0] = 12;
pd->cross[10].edge[1] = 16;
pd->cross[10].edge[2] = 13;
pd->cross[10].edge[3] = 15;
pd->cross[10].sign = 0;

pd->cross[11].edge[0] = 18;
pd->cross[11].edge[1] = 22;
pd->cross[11].edge[2] = 19;
pd->cross[11].edge[3] = 21;
pd->cross[11].sign = 0;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 1;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 0;
pd->edge[3].headpos = 3;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 2;
pd->edge[4].headpos = 1;
pd->edge[4].tail = 0;
pd->edge[4].tailpos = 1;

pd->edge[5].head = 3;
pd->edge[5].headpos = 3;
pd->edge[5].tail = 2;
pd->edge[5].tailpos = 3;

pd->edge[6].head = 4;
pd->edge[6].headpos = 0;
pd->edge[6].tail = 3;
pd->edge[6].tailpos = 1;

pd->edge[7].head = 5;
pd->edge[7].headpos = 0;
pd->edge[7].tail = 4;
pd->edge[7].tailpos = 2;

pd->edge[8].head = 6;
pd->edge[8].headpos = 0;
pd->edge[8].tail = 5;
pd->edge[8].tailpos = 2;

pd->edge[9].head = 1;
pd->edge[9].headpos = 2;
pd->edge[9].tail = 6;
pd->edge[9].tailpos = 2;

pd->edge[10].head = 7;
pd->edge[10].headpos = 0;
pd->edge[10].tail = 8;
pd->edge[10].tailpos = 0;

pd->edge[11].head = 9;
pd->edge[11].headpos = 0;
pd->edge[11].tail = 7;
pd->edge[11].tailpos = 2;

pd->edge[12].head = 10;
pd->edge[12].headpos = 0;
pd->edge[12].tail = 9;
pd->edge[12].tailpos = 2;

pd->edge[13].head = 7;
pd->edge[13].headpos = 3;
pd->edge[13].tail = 10;
pd->edge[13].tailpos = 2;

pd->edge[14].head = 9;
pd->edge[14].headpos = 1;
pd->edge[14].tail = 7;
pd->edge[14].tailpos = 1;

pd->edge[15].head = 10;
pd->edge[15].headpos = 3;
pd->edge[15].tail = 9;
pd->edge[15].tailpos = 3;

pd->edge[16].head = 4;
pd->edge[16].headpos = 1;
pd->edge[16].tail = 10;
pd->edge[16].tailpos = 1;

pd->edge[17].head = 5;
pd->edge[17].headpos = 3;
pd->edge[17].tail = 4;
pd->edge[17].tailpos = 3;

pd->edge[18].head = 11;
pd->edge[18].headpos = 0;
pd->edge[18].tail = 5;
pd->edge[18].tailpos = 1;

pd->edge[19].head = 8;
pd->edge[19].headpos = 2;
pd->edge[19].tail = 11;
pd->edge[19].tailpos = 2;

pd->edge[20].head = 6;
pd->edge[20].headpos = 3;
pd->edge[20].tail = 1;
pd->edge[20].tailpos = 1;

pd->edge[21].head = 11;
pd->edge[21].headpos = 3;
pd->edge[21].tail = 6;
pd->edge[21].tailpos = 1;

pd->edge[22].head = 8;
pd->edge[22].headpos = 3;
pd->edge[22].tail = 11;
pd->edge[22].tailpos = 1;

pd->edge[23].head = 1;
pd->edge[23].headpos = 3;
pd->edge[23].tail = 8;
pd->edge[23].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 10;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;

pd->comp[1].nedges = 10;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 10;
pd->comp[1].edge[1] = 11;
pd->comp[1].edge[2] = 12;
pd->comp[1].edge[3] = 13;
pd->comp[1].edge[4] = 14;
pd->comp[1].edge[5] = 15;
pd->comp[1].edge[6] = 16;
pd->comp[1].edge[7] = 17;
pd->comp[1].edge[8] = 18;
pd->comp[1].edge[9] = 19;

pd->comp[2].nedges = 4;
pd->comp[2].tag = 'C';

pd->comp[2].edge = calloc(pd->comp[2].nedges,sizeof(pd_idx_t));
assert(pd->comp[2].edge != NULL);

pd->comp[2].edge[0] = 20;
pd->comp[2].edge[1] = 21;
pd->comp[2].edge[2] = 22;
pd->comp[2].edge[3] = 23;


/* Face data */

pd->face[0].nedges = 8;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 23;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 10;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 13;
pd->face[0].or[3] = 0;

pd->face[0].edge[4] = 16;
pd->face[0].or[4] = 1;

pd->face[0].edge[5] = 6;
pd->face[0].or[5] = 0;

pd->face[0].edge[6] = 2;
pd->face[0].or[6] = 0;

pd->face[0].edge[7] = 4;
pd->face[0].or[7] = 0;

pd->face[1].nedges = 7;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 7;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 16;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 12;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 14;
pd->face[1].or[3] = 0;

pd->face[1].edge[4] = 10;
pd->face[1].or[4] = 0;

pd->face[1].edge[5] = 22;
pd->face[1].or[5] = 0;

pd->face[1].edge[6] = 18;
pd->face[1].or[6] = 0;

pd->face[2].nedges = 6;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 3;
pd->face[2].or[1] = 0;

pd->face[2].edge[2] = 6;
pd->face[2].or[2] = 1;

pd->face[2].edge[3] = 17;
pd->face[2].or[3] = 1;

pd->face[2].edge[4] = 8;
pd->face[2].or[4] = 1;

pd->face[2].edge[5] = 20;
pd->face[2].or[5] = 0;

pd->face[3].nedges = 4;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 9;
pd->face[3].or[0] = 0;

pd->face[3].edge[1] = 21;
pd->face[3].or[1] = 1;

pd->face[3].edge[2] = 19;
pd->face[3].or[2] = 1;

pd->face[3].edge[3] = 23;
pd->face[3].or[3] = 1;

pd->face[4].nedges = 3;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 1;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 5;
pd->face[4].or[1] = 1;

pd->face[4].edge[2] = 3;
pd->face[4].or[2] = 1;

pd->face[5].nedges = 3;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 8;
pd->face[5].or[0] = 0;

pd->face[5].edge[1] = 18;
pd->face[5].or[1] = 1;

pd->face[5].edge[2] = 21;
pd->face[5].or[2] = 0;

pd->face[6].nedges = 3;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 11;
pd->face[6].or[0] = 1;

pd->face[6].edge[1] = 15;
pd->face[6].or[1] = 1;

pd->face[6].edge[2] = 13;
pd->face[6].or[2] = 1;

pd->face[7].nedges = 2;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 1;
pd->face[7].or[0] = 0;

pd->face[7].edge[1] = 4;
pd->face[7].or[1] = 1;

pd->face[8].nedges = 2;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].or = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].or != NULL);

pd->face[8].edge[0] = 2;
pd->face[8].or[0] = 1;

pd->face[8].edge[1] = 5;
pd->face[8].or[1] = 0;

pd->face[9].nedges = 2;
pd->face[9].edge = calloc(pd->face[9].nedges,sizeof(pd_idx_t));
pd->face[9].or = calloc(pd->face[9].nedges,sizeof(pd_or_t));
assert(pd->face[9].edge != NULL);
assert(pd->face[9].or != NULL);

pd->face[9].edge[0] = 7;
pd->face[9].or[0] = 1;

pd->face[9].edge[1] = 17;
pd->face[9].or[1] = 0;

pd->face[10].nedges = 2;
pd->face[10].edge = calloc(pd->face[10].nedges,sizeof(pd_idx_t));
pd->face[10].or = calloc(pd->face[10].nedges,sizeof(pd_or_t));
assert(pd->face[10].edge != NULL);
assert(pd->face[10].or != NULL);

pd->face[10].edge[0] = 9;
pd->face[10].or[0] = 1;

pd->face[10].edge[1] = 20;
pd->face[10].or[1] = 1;

pd->face[11].nedges = 2;
pd->face[11].edge = calloc(pd->face[11].nedges,sizeof(pd_idx_t));
pd->face[11].or = calloc(pd->face[11].nedges,sizeof(pd_or_t));
assert(pd->face[11].edge != NULL);
assert(pd->face[11].or != NULL);

pd->face[11].edge[0] = 11;
pd->face[11].or[0] = 0;

pd->face[11].edge[1] = 14;
pd->face[11].or[1] = 1;

pd->face[12].nedges = 2;
pd->face[12].edge = calloc(pd->face[12].nedges,sizeof(pd_idx_t));
pd->face[12].or = calloc(pd->face[12].nedges,sizeof(pd_or_t));
assert(pd->face[12].edge != NULL);
assert(pd->face[12].or != NULL);

pd->face[12].edge[0] = 12;
pd->face[12].or[0] = 1;

pd->face[12].edge[1] = 15;
pd->face[12].or[1] = 0;

pd->face[13].nedges = 2;
pd->face[13].edge = calloc(pd->face[13].nedges,sizeof(pd_idx_t));
pd->face[13].or = calloc(pd->face[13].nedges,sizeof(pd_or_t));
assert(pd->face[13].edge != NULL);
assert(pd->face[13].or != NULL);

pd->face[13].edge[0] = 19;
pd->face[13].or[0] = 0;

pd->face[13].edge[1] = 22;
pd->face[13].or[1] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



pd_code_t *pd_create_r2_testD_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(10);
assert(pd != NULL);
pd->ncross = 10;
pd->nedges = 20;
pd->ncomps = 3;
pd->nfaces = 12;
sprintf(pd->hash,"%s","ChQMCQcEAwMCAgICAgICAwoIAgAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 4;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 3;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 19;
pd->cross[1].edge[2] = 9;
pd->cross[1].edge[3] = 18;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 4;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 5;
pd->cross[2].sign = 0;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 6;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 5;
pd->cross[3].sign = 0;

pd->cross[4].edge[0] = 6;
pd->cross[4].edge[1] = 17;
pd->cross[4].edge[2] = 7;
pd->cross[4].edge[3] = 16;
pd->cross[4].sign = 0;

pd->cross[5].edge[0] = 7;
pd->cross[5].edge[1] = 17;
pd->cross[5].edge[2] = 8;
pd->cross[5].edge[3] = 10;
pd->cross[5].sign = 0;

pd->cross[6].edge[0] = 8;
pd->cross[6].edge[1] = 18;
pd->cross[6].edge[2] = 9;
pd->cross[6].edge[3] = 19;
pd->cross[6].sign = 0;

pd->cross[7].edge[0] = 10;
pd->cross[7].edge[1] = 14;
pd->cross[7].edge[2] = 11;
pd->cross[7].edge[3] = 13;
pd->cross[7].sign = 0;

pd->cross[8].edge[0] = 11;
pd->cross[8].edge[1] = 14;
pd->cross[8].edge[2] = 12;
pd->cross[8].edge[3] = 15;
pd->cross[8].sign = 0;

pd->cross[9].edge[0] = 12;
pd->cross[9].edge[1] = 16;
pd->cross[9].edge[2] = 13;
pd->cross[9].edge[3] = 15;
pd->cross[9].sign = 0;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 1;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 0;
pd->edge[3].headpos = 3;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 2;
pd->edge[4].headpos = 1;
pd->edge[4].tail = 0;
pd->edge[4].tailpos = 1;

pd->edge[5].head = 3;
pd->edge[5].headpos = 3;
pd->edge[5].tail = 2;
pd->edge[5].tailpos = 3;

pd->edge[6].head = 4;
pd->edge[6].headpos = 0;
pd->edge[6].tail = 3;
pd->edge[6].tailpos = 1;

pd->edge[7].head = 5;
pd->edge[7].headpos = 0;
pd->edge[7].tail = 4;
pd->edge[7].tailpos = 2;

pd->edge[8].head = 6;
pd->edge[8].headpos = 0;
pd->edge[8].tail = 5;
pd->edge[8].tailpos = 2;

pd->edge[9].head = 1;
pd->edge[9].headpos = 2;
pd->edge[9].tail = 6;
pd->edge[9].tailpos = 2;

pd->edge[10].head = 7;
pd->edge[10].headpos = 0;
pd->edge[10].tail = 5;
pd->edge[10].tailpos = 3;

pd->edge[11].head = 8;
pd->edge[11].headpos = 0;
pd->edge[11].tail = 7;
pd->edge[11].tailpos = 2;

pd->edge[12].head = 9;
pd->edge[12].headpos = 0;
pd->edge[12].tail = 8;
pd->edge[12].tailpos = 2;

pd->edge[13].head = 7;
pd->edge[13].headpos = 3;
pd->edge[13].tail = 9;
pd->edge[13].tailpos = 2;

pd->edge[14].head = 8;
pd->edge[14].headpos = 1;
pd->edge[14].tail = 7;
pd->edge[14].tailpos = 1;

pd->edge[15].head = 9;
pd->edge[15].headpos = 3;
pd->edge[15].tail = 8;
pd->edge[15].tailpos = 3;

pd->edge[16].head = 4;
pd->edge[16].headpos = 3;
pd->edge[16].tail = 9;
pd->edge[16].tailpos = 1;

pd->edge[17].head = 5;
pd->edge[17].headpos = 1;
pd->edge[17].tail = 4;
pd->edge[17].tailpos = 1;

pd->edge[18].head = 1;
pd->edge[18].headpos = 3;
pd->edge[18].tail = 6;
pd->edge[18].tailpos = 1;

pd->edge[19].head = 6;
pd->edge[19].headpos = 3;
pd->edge[19].tail = 1;
pd->edge[19].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 10;
pd->comp[0].tag = 'B';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;

pd->comp[1].nedges = 8;
pd->comp[1].tag = 'A';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 10;
pd->comp[1].edge[1] = 11;
pd->comp[1].edge[2] = 12;
pd->comp[1].edge[3] = 13;
pd->comp[1].edge[4] = 14;
pd->comp[1].edge[5] = 15;
pd->comp[1].edge[6] = 16;
pd->comp[1].edge[7] = 17;

pd->comp[2].nedges = 2;
pd->comp[2].tag = 'C';

pd->comp[2].edge = calloc(pd->comp[2].nedges,sizeof(pd_idx_t));
assert(pd->comp[2].edge != NULL);

pd->comp[2].edge[0] = 18;
pd->comp[2].edge[1] = 19;


/* Face data */

pd->face[0].nedges = 9;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 3;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 6;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 16;
pd->face[0].or[3] = 0;

pd->face[0].edge[4] = 12;
pd->face[0].or[4] = 0;

pd->face[0].edge[5] = 14;
pd->face[0].or[5] = 0;

pd->face[0].edge[6] = 10;
pd->face[0].or[6] = 0;

pd->face[0].edge[7] = 8;
pd->face[0].or[7] = 1;

pd->face[0].edge[8] = 19;
pd->face[0].or[8] = 0;

pd->face[1].nedges = 7;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 0;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 18;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 8;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 17;
pd->face[1].or[3] = 0;

pd->face[1].edge[4] = 6;
pd->face[1].or[4] = 0;

pd->face[1].edge[5] = 2;
pd->face[1].or[5] = 0;

pd->face[1].edge[6] = 4;
pd->face[1].or[6] = 0;

pd->face[2].nedges = 4;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 7;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 10;
pd->face[2].or[1] = 1;

pd->face[2].edge[2] = 13;
pd->face[2].or[2] = 0;

pd->face[2].edge[3] = 16;
pd->face[2].or[3] = 1;

pd->face[3].nedges = 3;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 1;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 5;
pd->face[3].or[1] = 1;

pd->face[3].edge[2] = 3;
pd->face[3].or[2] = 1;

pd->face[4].nedges = 3;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 11;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 15;
pd->face[4].or[1] = 1;

pd->face[4].edge[2] = 13;
pd->face[4].or[2] = 1;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 1;
pd->face[5].or[0] = 0;

pd->face[5].edge[1] = 4;
pd->face[5].or[1] = 1;

pd->face[6].nedges = 2;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 2;
pd->face[6].or[0] = 1;

pd->face[6].edge[1] = 5;
pd->face[6].or[1] = 0;

pd->face[7].nedges = 2;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 7;
pd->face[7].or[0] = 0;

pd->face[7].edge[1] = 17;
pd->face[7].or[1] = 1;

pd->face[8].nedges = 2;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].or = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].or != NULL);

pd->face[8].edge[0] = 9;
pd->face[8].or[0] = 0;

pd->face[8].edge[1] = 18;
pd->face[8].or[1] = 1;

pd->face[9].nedges = 2;
pd->face[9].edge = calloc(pd->face[9].nedges,sizeof(pd_idx_t));
pd->face[9].or = calloc(pd->face[9].nedges,sizeof(pd_or_t));
assert(pd->face[9].edge != NULL);
assert(pd->face[9].or != NULL);

pd->face[9].edge[0] = 9;
pd->face[9].or[0] = 1;

pd->face[9].edge[1] = 19;
pd->face[9].or[1] = 1;

pd->face[10].nedges = 2;
pd->face[10].edge = calloc(pd->face[10].nedges,sizeof(pd_idx_t));
pd->face[10].or = calloc(pd->face[10].nedges,sizeof(pd_or_t));
assert(pd->face[10].edge != NULL);
assert(pd->face[10].or != NULL);

pd->face[10].edge[0] = 11;
pd->face[10].or[0] = 0;

pd->face[10].edge[1] = 14;
pd->face[10].or[1] = 1;

pd->face[11].nedges = 2;
pd->face[11].edge = calloc(pd->face[11].nedges,sizeof(pd_idx_t));
pd->face[11].or = calloc(pd->face[11].nedges,sizeof(pd_or_t));
assert(pd->face[11].edge != NULL);
assert(pd->face[11].or != NULL);

pd->face[11].edge[0] = 12;
pd->face[11].or[0] = 1;

pd->face[11].edge[1] = 15;
pd->face[11].or[1] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



/*                           ___13___<__                                      */
/*                          /  (6)     \                                     */
/*                         /        __>_ \ ____             (0)               */
/*                      _____11_   15    10\   \                              */
/*                     /   7    \ /        |    \                             */
/*                    /  /       /    (12) /     \                            */
/*                   /   |(11)  /9\       /       \                           */
/*                  /    \     /   \     ^         \                          */
/*                 /     14   ^     \12_/           \                         */
/*                ^        \_/                      16                        */
/*               10                                   \                       */
/*              /                                      \                      */
/*             /                                        \                     */
/*            /                    (1)                   \                    */
/*           /                                           /                    */
/*     ______________<_22_____                   _<_7__ / _____6__<______     */
/*    /     8                 \                 / (9)  /4                |    */
/*   /    /          (13)      \               /      /                  |    */
/*  /    /_______<_19_______________<_18_____ / _<_17/                   |    */
/* |                          11             /5            __<_3___      |    */
/* |                             \   (5)    /             /        \     |    */
/* |                             21        8    (2)      /    (4)   \    |    */
/* |                         ___________<_/             /        _>__\ __|    */
/* |                        |      6\               _______1    5     \3      */
/* |             (3)        9 (10)   ^              /   0    \  /      |       */
/* \                        |      20             /  /        /  (8)  /       */
/*  \                       |     /              /  |  (7)   / 2     /        */
/*   \                      |_>__________0_>____/   \       /  \_2_>/         */
/*    \                          1                   \_4_>_/                  */
/*     \_______23_>____________/                                              */
/*                           		                                      */
/*   r2 at crossings 1 and 6 	                                              */
/*                           		                                      */
/*                           ____<_3____                                      */
/*                          /           \                    (0)               */
/*                         /    (3) __>_ \ ____                               */
/*                      _____1_>   5     3\    \                              */
/*                     /   0    \ /        |    \                             */
/*                    /  /       /    (6)  /     \                            */
/*                   /   | (5)  /2\       /       \                           */
/*                  /    \     /   \     ^         6                          */
/*                 ^      4   ^     \_2_/           \                         */
/*                0        \_/                       \                        */
/*               /                                    \                       */
/*              /                                      \                      */
/*             /                  (1)                   \                     */
/*            /                                          \                    */
/*           /                                           /                    */
/*     ______________<_18_____                   __<___ / _____<_16______     */
/*    /     1                 \                17      /4                |    */
/*   /    /          (8)       \               / (7)  7                  |    */
/*  /    /__________<_9______________<_8_____ / ___<_/                   |    */
/* |                           6             /5            __<_13__      |    */
/* |                             \          /             /        \     |    */
/* |                              \        /      (2)    /   (4)    \    |    */
/* |               (9)             \      /             /        __>_\ __|    */
/* |                                \    |          _____11_   15     \9      */
/* |                                /    |         /   7    \  /      |       */
/* \                               /     |        /  /(10)    / (11)  /       */
/*  \                             /      |       /  |        /       ^        */
/*   \                           /       |_10_>_/   \       ^ 8\_12_/         */
/*    \                         /                    \_14__/                  */
/*     \_____19_>______________/                                              */
/*                           				                                  */
/*                           				                                  */


    bool r2testD() {
  
  printf("---------------------------------------\n"
	 "r2 test D\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r2_testD_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after-0 configuration...");
  pd_code_t *after0 = pd_create_r2_testD_after_0();
  printf("done (passes pd_ok)\n");

  pd_code_t *after[1] = {after0};
    
  printf("performing r2 at crossings 1 and 6...\n");

  pd_idx_t  noutpd;
  pd_code_t **newpd;
  pd_idx_t  cr[2] = {1,6};

  pd_R2_bigon_elimination(before,cr,&noutpd,&newpd);

  if (noutpd != 1) { 
    printf("fail (did not generate 1 output pd codes)\n");
    return false;
  }

  pd_idx_t ct;

  for(ct=0;ct<noutpd;ct++) { 

    if (!pd_ok(newpd[ct])) { 
      printf("fail (output pd %d does not pass pd_ok)\n",ct);
      return false;
    }
    printf("\t output pd %d passes pd_ok...pass\n",ct);
    printf("\t testing for isomorphism with after configuration %d...",ct);
    if (!pd_isomorphic(newpd[ct],after[ct])) { 
      pd_printf("fail (output pd %PD ",newpd[ct]);
      pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after[ct]);
      return false;
    }

    printf("pass (output %d and after %d isomorphic)\n",ct,ct);

  }

  printf("performing r2 at crossings 1 and 6...pass\n");
  printf("housecleaning...");

  pd_code_free(&before);
  pd_code_free(&after0);
  pd_code_free(&newpd[0]);

  free(newpd);

  printf("done\n");

  printf("---------------------------------------\n"
	 "r2 test D : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}

/* r2 test E assembled c code */

pd_code_t *pd_create_r2_testE_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(2);
assert(pd != NULL);
pd->ncross = 2;
pd->nedges = 4;
pd->ncomps = 1;
pd->nfaces = 4;
sprintf(pd->hash,"%s","AgQEBAIBAQEEAAAAAAAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 0;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 3;
pd->cross[0].sign = 1;

pd->cross[1].edge[0] = 1;
pd->cross[1].edge[1] = 2;
pd->cross[1].edge[2] = 2;
pd->cross[1].edge[3] = 3;
pd->cross[1].sign = 0;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 1;

pd->edge[1].head = 1;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 1;
pd->edge[2].headpos = 1;
pd->edge[2].tail = 1;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 0;
pd->edge[3].headpos = 3;
pd->edge[3].tail = 1;
pd->edge[3].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 4;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;


/* Face data */

pd->face[0].nedges = 4;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 3;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 2;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 1;
pd->face[0].or[3] = 0;

pd->face[1].nedges = 2;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 3;
pd->face[1].or[1] = 1;

pd->face[2].nedges = 1;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 0;

pd->face[3].nedges = 1;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 2;
pd->face[3].or[0] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



pd_code_t *pd_create_r2_testE_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(2);
assert(pd != NULL);
pd->ncross = 0;
pd->nedges = 1;
pd->ncomps = 1;
pd->nfaces = 2;
sprintf(pd->hash,"%s","AAECAQEBAQAAAAAAAAAAAAAAAAAAAAA");

/* Crossing data. */


/* Edge data */

pd->edge[0].head = -1;
pd->edge[0].headpos = 4;
pd->edge[0].tail = -1;
pd->edge[0].tailpos = 4;


/* Component Data */

pd->comp[0].nedges = 1;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;


/* Face data */

pd->face[0].nedges = 1;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 1;

pd->face[1].nedges = 1;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 0;
pd->face[1].or[0] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



/*                                                                            */
/*                                                                            */
/*                                                       (0)                  */
/*                                                                            */
/*                        ____________0_>_________                            */
/*                        |                      |                            */
/*                        |        (2)           |                            */
/*                        |                      |                            */
/*                        |                      |                            */
/*                         \       ______________/                            */
/*                          \     |                                           */
/*                           \_______________                                 */
/*                                |0        |                                 */
/*                                |         |                                 */
/*                                |  (1)    ^                                 */
/*                                1         3                                 */
/*                                |         |                                 */
/*                             _____________|                                 */
/*                            /    1                                          */
/*                           /    |                                           */
/*                          /     |                                           */
/*                         /      |_____2_>______                             */
/*                         |                     |                            */
/*                         |                     |                            */
/*                         |        (3)          |                            */
/*                         |                     |                            */
/*                         |_____________________|                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                      ______________________________                        */
/*                      |                            |                        */
/*                      |                            |                        */
/*                      |                            |                        */
/*                      |                            |                        */
/*                      |                            |                        */
/*                      |                            |                        */
/*                      |                            |                        */
/*                      |                            |                        */
/*                      |                            |                        */
/*                      |                            |                        */
/*                      |                            |                        */
/*                      |                            |                        */
/*                      |                            |                        */
/*                      |                            |                        */
/*                      |                            |                        */
/*                      |                            |                        */
/*                      |                            |                        */
/*                      |                            |                        */
/*                      |                            |                        */
/*                      |____________________________|                        */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */


    bool r2testE() {
  
  printf("---------------------------------------\n"
	 "r2 test E\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r2_testE_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after-0 configuration...");
  pd_code_t *after0 = pd_create_r2_testE_after_0();
  printf("done (passes pd_ok)\n");

  pd_code_t *after[1] = {after0};
    
  printf("performing r2 at crossings 0 and 1...\n");

  pd_idx_t  noutpd;
  pd_code_t **newpd;
  pd_idx_t  cr[2] = {0,1};

  pd_R2_bigon_elimination(before,cr,&noutpd,&newpd);

  if (noutpd != 1) { 
    printf("fail (did not generate 1 output pd codes)\n");
    return false;
  }

  pd_idx_t ct;

  for(ct=0;ct<noutpd;ct++) { 

    if (!pd_ok(newpd[ct])) { 
      printf("fail (output pd %d does not pass pd_ok)\n",ct);
      return false;
    }
    printf("\t output pd %d passes pd_ok...pass\n",ct);
    printf("\t testing for isomorphism with after configuration %d...",ct);
    if (!pd_isomorphic(newpd[ct],after[ct])) { 
      pd_printf("fail (output pd %PD ",newpd[ct]);
      pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after[ct]);
      return false;
    }

    printf("pass (output %d and after %d isomorphic)\n",ct,ct);

  }

  printf("performing r2 at crossings 0 and 1...pass\n");
  printf("housecleaning...");

  pd_code_free(&before);
  pd_code_free(&after0);
  pd_code_free(&newpd[0]);

  free(newpd);

  printf("done\n");

  printf("---------------------------------------\n"
	 "r2 test E : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}

/* r2 test F assembled c code */

pd_code_t *pd_create_r2_testF_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(3);
assert(pd != NULL);
pd->ncross = 3;
pd->nedges = 6;
pd->ncomps = 1;
pd->nfaces = 5;
sprintf(pd->hash,"%s","AwYFAwMCAgIBBgAAAAAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 3;
pd->cross[0].edge[2] = 5;
pd->cross[0].edge[3] = 2;
pd->cross[0].sign = 1;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 4;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 3;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 4;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 5;
pd->cross[2].sign = 0;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 0;
pd->edge[2].headpos = 3;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 1;
pd->edge[3].headpos = 3;
pd->edge[3].tail = 0;
pd->edge[3].tailpos = 1;

pd->edge[4].head = 2;
pd->edge[4].headpos = 1;
pd->edge[4].tail = 1;
pd->edge[4].tailpos = 1;

pd->edge[5].head = 0;
pd->edge[5].headpos = 2;
pd->edge[5].tail = 2;
pd->edge[5].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 6;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;


/* Face data */

pd->face[0].nedges = 3;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 2;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 4;
pd->face[0].or[2] = 0;

pd->face[1].nedges = 3;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 5;
pd->face[1].or[1] = 1;

pd->face[1].edge[2] = 3;
pd->face[1].or[2] = 1;

pd->face[2].nedges = 2;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 3;
pd->face[2].or[1] = 0;

pd->face[3].nedges = 2;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 1;
pd->face[3].or[0] = 0;

pd->face[3].edge[1] = 4;
pd->face[3].or[1] = 1;

pd->face[4].nedges = 2;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 2;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 5;
pd->face[4].or[1] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



pd_code_t *pd_create_r2_testF_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(2);
assert(pd != NULL);
pd->ncross = 1;
pd->nedges = 2;
pd->ncomps = 1;
pd->nfaces = 3;
sprintf(pd->hash,"%s","AQIDAgEBAQIAAAAAAAAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 0;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 1;
pd->cross[0].sign = 0;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 1;

pd->edge[1].head = 0;
pd->edge[1].headpos = 3;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 2;


/* Component Data */

pd->comp[0].nedges = 2;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;


/* Face data */

pd->face[0].nedges = 2;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 1;
pd->face[0].or[1] = 0;

pd->face[1].nedges = 1;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 0;
pd->face[1].or[0] = 0;

pd->face[2].nedges = 1;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 1;
pd->face[2].or[0] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



/*                                                                            */
/*                                                                            */
/*                          _______5_>________                                */
/*                         /                  \                (1)            */
/*                        /                    \                              */
/*                       /                      \                             */
/*                      /            (4)         \                            */
/*                      \                         \                           */
/*                       \                        /                           */
/*                        \                      /                            */
/*         _______________ \ ________2_>_______ / _______3_>______            */
/*         \               2\                  /0                /            */
/*          \                \                /                 /             */
/*           \                \     (0)      /                 /              */
/*            \                \            <                 /               */
/*             \                ^          0                 /                */
/*              \                4        /                 /                 */
/*               \       (3)      \      /        (2)      /                  */
/*                \                \    /                 /                   */
/*                 \                \  /                 /                    */
/*                  \                 /                 /                     */
/*                   \               /                 /                      */
/*                    \             / 1\              /                       */
/*                     \           /    \            /                        */
/*                      \         /      \          /                         */
/*                       \____<1_/        \________/                          */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            *//*                       r2 at 2 and 0                                                     */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                         ___________1_>_________                            */
/*                        /                       \            (2)            */
/*                       /                         \                          */
/*                      /                           \                         */
/*                     /          __0_>___           \                        */
/*                    /          /        \           \                       */
/*                   /          /          \           \                      */
/*                  /          /    (1)     \           \                     */
/*                 /           \            /            \                    */
/*                /             \          /              \                   */
/*               /               \        /     (0)        \                  */
/*               \                \      /                 /                  */
/*                \                \    /                 /                   */
/*                 \                \ 0/                 /                    */
/*                  \                 /                 /                     */
/*                   \               /                 /                      */
/*                    \             /  \              /                       */
/*                     \           /    \            /                        */
/*                      \         /      \          /                         */
/*                       \_______/        \________/                          */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */



    bool r2testF() {
  
  printf("---------------------------------------\n"
	 "r2 test F\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r2_testF_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after-0 configuration...");
  pd_code_t *after0 = pd_create_r2_testF_after_0();
  printf("done (passes pd_ok)\n");

  pd_code_t *after[1] = {after0};
    
  printf("performing r2 at crossings 0 and 2...\n");

  pd_idx_t  noutpd;
  pd_code_t **newpd;
  pd_idx_t  cr[2] = {0,2};

  pd_R2_bigon_elimination(before,cr,&noutpd,&newpd);

  if (noutpd != 1) { 
    printf("fail (did not generate 1 output pd codes)\n");
    return false;
  }

  pd_idx_t ct;

  for(ct=0;ct<noutpd;ct++) { 

    if (!pd_ok(newpd[ct])) { 
      printf("fail (output pd %d does not pass pd_ok)\n",ct);
      return false;
    }
    printf("\t output pd %d passes pd_ok...pass\n",ct);
    printf("\t testing for isomorphism with after configuration %d...",ct);
    if (!pd_isomorphic(newpd[ct],after[ct])) { 
      pd_printf("fail (output pd %PD ",newpd[ct]);
      pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after[ct]);
      return false;
    }

    printf("pass (output %d and after %d isomorphic)\n",ct,ct);

  }

  printf("performing r2 at crossings 0 and 2...pass\n");
  printf("housecleaning...");

  pd_code_free(&before);
  pd_code_free(&after0);
  pd_code_free(&newpd[0]);

  free(newpd);

  printf("done\n");

  printf("---------------------------------------\n"
	 "r2 test F : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}

/* r2 test G assembled c code */

pd_code_t *pd_create_r2_testG_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(10);
assert(pd != NULL);
pd->ncross = 10;
pd->nedges = 20;
pd->ncomps = 2;
pd->nfaces = 12;
sprintf(pd->hash,"%s","ChQMCgUEBAMCAgICAgICAgwIAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 4;
pd->cross[0].edge[2] = 11;
pd->cross[0].edge[3] = 5;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 5;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 6;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 9;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 8;
pd->cross[2].sign = 0;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 9;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 10;
pd->cross[3].sign = 0;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 11;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 10;
pd->cross[4].sign = 0;

pd->cross[5].edge[0] = 6;
pd->cross[5].edge[1] = 14;
pd->cross[5].edge[2] = 7;
pd->cross[5].edge[3] = 15;
pd->cross[5].sign = 1;

pd->cross[6].edge[0] = 7;
pd->cross[6].edge[1] = 14;
pd->cross[6].edge[2] = 8;
pd->cross[6].edge[3] = 13;
pd->cross[6].sign = 0;

pd->cross[7].edge[0] = 12;
pd->cross[7].edge[1] = 17;
pd->cross[7].edge[2] = 19;
pd->cross[7].edge[3] = 16;
pd->cross[7].sign = 0;

pd->cross[8].edge[0] = 12;
pd->cross[8].edge[1] = 18;
pd->cross[8].edge[2] = 13;
pd->cross[8].edge[3] = 17;
pd->cross[8].sign = 0;

pd->cross[9].edge[0] = 15;
pd->cross[9].edge[1] = 18;
pd->cross[9].edge[2] = 16;
pd->cross[9].edge[3] = 19;
pd->cross[9].sign = 0;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 0;
pd->edge[4].headpos = 1;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 1;
pd->edge[5].headpos = 1;
pd->edge[5].tail = 0;
pd->edge[5].tailpos = 3;

pd->edge[6].head = 5;
pd->edge[6].headpos = 0;
pd->edge[6].tail = 1;
pd->edge[6].tailpos = 3;

pd->edge[7].head = 6;
pd->edge[7].headpos = 0;
pd->edge[7].tail = 5;
pd->edge[7].tailpos = 2;

pd->edge[8].head = 2;
pd->edge[8].headpos = 3;
pd->edge[8].tail = 6;
pd->edge[8].tailpos = 2;

pd->edge[9].head = 3;
pd->edge[9].headpos = 1;
pd->edge[9].tail = 2;
pd->edge[9].tailpos = 1;

pd->edge[10].head = 4;
pd->edge[10].headpos = 3;
pd->edge[10].tail = 3;
pd->edge[10].tailpos = 3;

pd->edge[11].head = 0;
pd->edge[11].headpos = 2;
pd->edge[11].tail = 4;
pd->edge[11].tailpos = 1;

pd->edge[12].head = 8;
pd->edge[12].headpos = 0;
pd->edge[12].tail = 7;
pd->edge[12].tailpos = 0;

pd->edge[13].head = 6;
pd->edge[13].headpos = 3;
pd->edge[13].tail = 8;
pd->edge[13].tailpos = 2;

pd->edge[14].head = 5;
pd->edge[14].headpos = 1;
pd->edge[14].tail = 6;
pd->edge[14].tailpos = 1;

pd->edge[15].head = 9;
pd->edge[15].headpos = 0;
pd->edge[15].tail = 5;
pd->edge[15].tailpos = 3;

pd->edge[16].head = 7;
pd->edge[16].headpos = 3;
pd->edge[16].tail = 9;
pd->edge[16].tailpos = 2;

pd->edge[17].head = 8;
pd->edge[17].headpos = 3;
pd->edge[17].tail = 7;
pd->edge[17].tailpos = 1;

pd->edge[18].head = 9;
pd->edge[18].headpos = 1;
pd->edge[18].tail = 8;
pd->edge[18].tailpos = 1;

pd->edge[19].head = 7;
pd->edge[19].headpos = 2;
pd->edge[19].tail = 9;
pd->edge[19].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 12;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;
pd->comp[0].edge[10] = 10;
pd->comp[0].edge[11] = 11;

pd->comp[1].nedges = 8;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 12;
pd->comp[1].edge[1] = 13;
pd->comp[1].edge[2] = 14;
pd->comp[1].edge[3] = 15;
pd->comp[1].edge[4] = 16;
pd->comp[1].edge[5] = 17;
pd->comp[1].edge[6] = 18;
pd->comp[1].edge[7] = 19;


/* Face data */

pd->face[0].nedges = 10;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 6;
pd->face[0].or[1] = 1;

pd->face[0].edge[2] = 15;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 19;
pd->face[0].or[3] = 1;

pd->face[0].edge[4] = 17;
pd->face[0].or[4] = 1;

pd->face[0].edge[5] = 13;
pd->face[0].or[5] = 1;

pd->face[0].edge[6] = 8;
pd->face[0].or[6] = 1;

pd->face[0].edge[7] = 2;
pd->face[0].or[7] = 1;

pd->face[0].edge[8] = 10;
pd->face[0].or[8] = 1;

pd->face[0].edge[9] = 4;
pd->face[0].or[9] = 1;

pd->face[1].nedges = 5;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 5;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 11;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 3;
pd->face[1].or[3] = 0;

pd->face[1].edge[4] = 9;
pd->face[1].or[4] = 0;

pd->face[2].nedges = 4;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 1;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 8;
pd->face[2].or[1] = 0;

pd->face[2].edge[2] = 14;
pd->face[2].or[2] = 1;

pd->face[2].edge[3] = 6;
pd->face[2].or[3] = 0;

pd->face[3].nedges = 4;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 7;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 13;
pd->face[3].or[1] = 0;

pd->face[3].edge[2] = 18;
pd->face[3].or[2] = 1;

pd->face[3].edge[3] = 15;
pd->face[3].or[3] = 0;

pd->face[4].nedges = 3;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 12;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 16;
pd->face[4].or[1] = 0;

pd->face[4].edge[2] = 18;
pd->face[4].or[2] = 0;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 0;
pd->face[5].or[0] = 0;

pd->face[5].edge[1] = 5;
pd->face[5].or[1] = 1;

pd->face[6].nedges = 2;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 2;
pd->face[6].or[0] = 0;

pd->face[6].edge[1] = 9;
pd->face[6].or[1] = 1;

pd->face[7].nedges = 2;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 3;
pd->face[7].or[0] = 1;

pd->face[7].edge[1] = 10;
pd->face[7].or[1] = 0;

pd->face[8].nedges = 2;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].or = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].or != NULL);

pd->face[8].edge[0] = 4;
pd->face[8].or[0] = 0;

pd->face[8].edge[1] = 11;
pd->face[8].or[1] = 1;

pd->face[9].nedges = 2;
pd->face[9].edge = calloc(pd->face[9].nedges,sizeof(pd_idx_t));
pd->face[9].or = calloc(pd->face[9].nedges,sizeof(pd_or_t));
assert(pd->face[9].edge != NULL);
assert(pd->face[9].or != NULL);

pd->face[9].edge[0] = 7;
pd->face[9].or[0] = 0;

pd->face[9].edge[1] = 14;
pd->face[9].or[1] = 0;

pd->face[10].nedges = 2;
pd->face[10].edge = calloc(pd->face[10].nedges,sizeof(pd_idx_t));
pd->face[10].or = calloc(pd->face[10].nedges,sizeof(pd_or_t));
assert(pd->face[10].edge != NULL);
assert(pd->face[10].or != NULL);

pd->face[10].edge[0] = 12;
pd->face[10].or[0] = 1;

pd->face[10].edge[1] = 17;
pd->face[10].or[1] = 0;

pd->face[11].nedges = 2;
pd->face[11].edge = calloc(pd->face[11].nedges,sizeof(pd_idx_t));
pd->face[11].or = calloc(pd->face[11].nedges,sizeof(pd_or_t));
assert(pd->face[11].edge != NULL);
assert(pd->face[11].or != NULL);

pd->face[11].edge[0] = 16;
pd->face[11].or[0] = 1;

pd->face[11].edge[1] = 19;
pd->face[11].or[1] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



pd_code_t *pd_create_r2_testG_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(5);
assert(pd != NULL);
pd->ncross = 5;
pd->nedges = 10;
pd->ncomps = 1;
pd->nfaces = 7;
sprintf(pd->hash,"%s","BQoHBQUCAgICAgEKAAAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 4;
pd->cross[0].edge[2] = 9;
pd->cross[0].edge[3] = 5;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 5;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 6;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 7;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 6;
pd->cross[2].sign = 0;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 7;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 8;
pd->cross[3].sign = 0;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 9;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 8;
pd->cross[4].sign = 0;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 0;
pd->edge[4].headpos = 1;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 1;
pd->edge[5].headpos = 1;
pd->edge[5].tail = 0;
pd->edge[5].tailpos = 3;

pd->edge[6].head = 2;
pd->edge[6].headpos = 3;
pd->edge[6].tail = 1;
pd->edge[6].tailpos = 3;

pd->edge[7].head = 3;
pd->edge[7].headpos = 1;
pd->edge[7].tail = 2;
pd->edge[7].tailpos = 1;

pd->edge[8].head = 4;
pd->edge[8].headpos = 3;
pd->edge[8].tail = 3;
pd->edge[8].tailpos = 3;

pd->edge[9].head = 0;
pd->edge[9].headpos = 2;
pd->edge[9].tail = 4;
pd->edge[9].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 10;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;


/* Face data */

pd->face[0].nedges = 5;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 6;
pd->face[0].or[1] = 1;

pd->face[0].edge[2] = 2;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 8;
pd->face[0].or[3] = 1;

pd->face[0].edge[4] = 4;
pd->face[0].or[4] = 1;

pd->face[1].nedges = 5;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 5;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 9;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 3;
pd->face[1].or[3] = 0;

pd->face[1].edge[4] = 7;
pd->face[1].or[4] = 0;

pd->face[2].nedges = 2;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 5;
pd->face[2].or[1] = 1;

pd->face[3].nedges = 2;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 1;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 6;
pd->face[3].or[1] = 0;

pd->face[4].nedges = 2;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 2;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 7;
pd->face[4].or[1] = 1;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 3;
pd->face[5].or[0] = 1;

pd->face[5].edge[1] = 8;
pd->face[5].or[1] = 0;

pd->face[6].nedges = 2;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 4;
pd->face[6].or[0] = 0;

pd->face[6].edge[1] = 9;
pd->face[6].or[1] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}

pd_code_t *pd_create_r2_testG_after_1() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(3);
assert(pd != NULL);
pd->ncross = 3;
pd->nedges = 6;
pd->ncomps = 1;
pd->nfaces = 5;
sprintf(pd->hash,"%s","AwYFAwMCAgIBBgAAAAAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 3;
pd->cross[0].edge[2] = 5;
pd->cross[0].edge[3] = 2;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 4;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 3;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 4;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 5;
pd->cross[2].sign = 0;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 0;
pd->edge[2].headpos = 3;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 1;
pd->edge[3].headpos = 3;
pd->edge[3].tail = 0;
pd->edge[3].tailpos = 1;

pd->edge[4].head = 2;
pd->edge[4].headpos = 1;
pd->edge[4].tail = 1;
pd->edge[4].tailpos = 1;

pd->edge[5].head = 0;
pd->edge[5].headpos = 2;
pd->edge[5].tail = 2;
pd->edge[5].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 6;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;


/* Face data */

pd->face[0].nedges = 3;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 2;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 4;
pd->face[0].or[2] = 0;

pd->face[1].nedges = 3;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 5;
pd->face[1].or[1] = 1;

pd->face[1].edge[2] = 3;
pd->face[1].or[2] = 1;

pd->face[2].nedges = 2;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 3;
pd->face[2].or[1] = 0;

pd->face[3].nedges = 2;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 1;
pd->face[3].or[0] = 0;

pd->face[3].edge[1] = 4;
pd->face[3].or[1] = 1;

pd->face[4].nedges = 2;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 2;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 5;
pd->face[4].or[1] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



/*                  ___                               _19_>__                 */
/*                 /   \                             /       \                */
/*                ^     \           (0)             /         \               */
/*               0       \                         /   (11)    \              */
/*              /         \                       /             \             */
/*             /   (5)     \            __15_>____| _____16_>_______________  */
/*             |            \           |         |9               7       |  */
/*       _____ | ______5_>________6_>__________   |               |        |  */
/*      /      |0            1\         |5     |  |               |        |  */
/*     /       |               >  (2)   ^      |  |               |        |  */
/*    /        |                \      14 (9)  7  |               |        |  */
/*   /   (8)   |                 1      |      |  |       (4)     |       17  */
/*   ^         ^                  \     |      |  ^              12        |  */
/*    4       11                  |  8______<__|   \              /        |  */
/*     \      /                   | /    6          \            /         |  */
/*      \    /        (1)         |/    |           18          <   (10)   |  */
/*       \  /                    2|     |        (3)  \        /           |  */
/*       4\                      /|     |              \      /            |  */
/*         \                    / |     \               \    /             |  */
/*       /  \                  9  |      \               \  /              /  */
/*      /    ^                /   |       \                /              /   */
/*     /      3              /    |        \              /8\            /    */
/*    /        \   _____<___/     |         \            /   \          /     */
/*   /          \ /               |          \__13_<____/     \____<___/      */
/*  |   (7)      /                |                                           */
/*  |           /3\     (6)      /                                            */
/*  |          /   \            /                                             */
/*  |         /     \          /                                              */
/*  |__10_<__/       \___<_2__/                                               */
/*                                                                            */
/* R2 on face 9  at crossings 5 and 6                                         */

/*                  ___             (0/1)              ___5_>_                */
/*                 /   \                              /       \               */
/*                ^     \                            /         \              */
/*               0       \                          /  (4)      \             */
/*              /         \                        /             \            */
/*             /  (2)      \             __________| _____2_>________________ */
/*             |            \            |        2|                0       | */
/*       ______|______5_>________6_>__   \         |               |        | */
/*      /      |0            1\       |   \        |               |        | */
/*     /       |               \ (3)  |    |       |               |        | */
/*    /        |                1     |    |       |      (0)      |        | */
/*   ^  (6)    |                 >    |    |       |               |        | */
/*   4         |                  \   |    |       \               |        | */
/*    \        ^                  |  _|    |        ^              0        | */
/*     \      9      (1)          | /      |         4            <         | */
/*      \    /                    |/       |          \          /          | */
/*       \  /                    2|        |           \        /           | */
/*        \                      /|        |  (3)       \      /            | */
/*       4 \                    / |        |             \    /             | */
/*       /  \                  /  |        |              \  /     (2)      / */
/*      /    ^                /   |        |                /              /  */
/*     /      3              /    |         \              /1\            /   */
/*    /        \   ____<_7__/     |          \            /   \          /    */
/*   /            /               |           \___<_1____/     \____<_3_/     */
/*  |  (5)       /3               |                                          */
/*  |           /        (4)     /                                            */
/*  |          /   \            /                                             */
/*  |         /     \          /                                              */
/*  |___<_8__/       \___<_2__/                                               */
/*                                                                            */


    bool r2testG() {
  
  printf("---------------------------------------\n"
	 "r2 test G\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r2_testG_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after-0 configuration...");
  pd_code_t *after0 = pd_create_r2_testG_after_0();
  printf("done (passes pd_ok)\n");

  printf("creating after-1 configuration...");
  pd_code_t *after1 = pd_create_r2_testG_after_1();
  printf("done (passes pd_ok)\n");

  pd_code_t *after[2] = {after0,after1};
    
  printf("performing r2 at crossings 5 and 6...\n");

  pd_idx_t  noutpd;
  pd_code_t **newpd;
  pd_idx_t  cr[2] = { 5,6 };

  pd_R2_bigon_elimination(before,cr,&noutpd,&newpd);

  if (noutpd != 2) { 
    printf("fail (did not generate 2 output pd codes)\n");
    return false;
  }

  pd_idx_t ct;

  for(ct=0;ct<noutpd;ct++) { 

    if (!pd_ok(newpd[ct])) { 
      printf("fail (output pd %d does not pass pd_ok)\n",ct);
      return false;
    }
    printf("\t output pd %d passes pd_ok...pass\n",ct);
    printf("\t testing for isomorphism with after configuration %d...",ct);
    if (!pd_isomorphic(newpd[ct],after[ct])) { 
      pd_printf("fail (output pd %PD ",newpd[ct]);
      pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after[ct]);
      return false;
    }

    printf("pass (output %d and after %d isomorphic)\n",ct,ct);

  }

  printf("performing r2 at crossings 5 and 6...pass\n");
  printf("housecleaning...");

  pd_code_free(&before);
  pd_code_free(&after0);
  pd_code_free(&after1);
  pd_code_free(&newpd[0]);
  pd_code_free(&newpd[1]);

  free(newpd);

  printf("done\n");

  printf("---------------------------------------\n"
	 "r2 test G : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}

/* r2 test H assembled c code */

pd_code_t *pd_create_r2_testH_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(10);
assert(pd != NULL);
pd->ncross = 10;
pd->nedges = 20;
pd->ncomps = 3;
pd->nfaces = 12;
sprintf(pd->hash,"%s","ChQMBQUEBAQEBAICAgICAwwEBAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 16;
pd->cross[0].edge[2] = 11;
pd->cross[0].edge[3] = 19;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 17;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 16;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 15;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 12;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 13;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 12;
pd->cross[3].sign = 0;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 10;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 11;
pd->cross[4].sign = 0;

pd->cross[5].edge[0] = 4;
pd->cross[5].edge[1] = 10;
pd->cross[5].edge[2] = 5;
pd->cross[5].edge[3] = 9;
pd->cross[5].sign = 0;

pd->cross[6].edge[0] = 5;
pd->cross[6].edge[1] = 13;
pd->cross[6].edge[2] = 6;
pd->cross[6].edge[3] = 14;
pd->cross[6].sign = 0;

pd->cross[7].edge[0] = 6;
pd->cross[7].edge[1] = 15;
pd->cross[7].edge[2] = 7;
pd->cross[7].edge[3] = 14;
pd->cross[7].sign = 0;

pd->cross[8].edge[0] = 7;
pd->cross[8].edge[1] = 17;
pd->cross[8].edge[2] = 8;
pd->cross[8].edge[3] = 18;
pd->cross[8].sign = 0;

pd->cross[9].edge[0] = 8;
pd->cross[9].edge[1] = 19;
pd->cross[9].edge[2] = 9;
pd->cross[9].edge[3] = 18;
pd->cross[9].sign = 0;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 5;
pd->edge[4].headpos = 0;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 6;
pd->edge[5].headpos = 0;
pd->edge[5].tail = 5;
pd->edge[5].tailpos = 2;

pd->edge[6].head = 7;
pd->edge[6].headpos = 0;
pd->edge[6].tail = 6;
pd->edge[6].tailpos = 2;

pd->edge[7].head = 8;
pd->edge[7].headpos = 0;
pd->edge[7].tail = 7;
pd->edge[7].tailpos = 2;

pd->edge[8].head = 9;
pd->edge[8].headpos = 0;
pd->edge[8].tail = 8;
pd->edge[8].tailpos = 2;

pd->edge[9].head = 5;
pd->edge[9].headpos = 3;
pd->edge[9].tail = 9;
pd->edge[9].tailpos = 2;

pd->edge[10].head = 4;
pd->edge[10].headpos = 1;
pd->edge[10].tail = 5;
pd->edge[10].tailpos = 1;

pd->edge[11].head = 0;
pd->edge[11].headpos = 2;
pd->edge[11].tail = 4;
pd->edge[11].tailpos = 3;

pd->edge[12].head = 3;
pd->edge[12].headpos = 3;
pd->edge[12].tail = 2;
pd->edge[12].tailpos = 3;

pd->edge[13].head = 6;
pd->edge[13].headpos = 1;
pd->edge[13].tail = 3;
pd->edge[13].tailpos = 1;

pd->edge[14].head = 7;
pd->edge[14].headpos = 3;
pd->edge[14].tail = 6;
pd->edge[14].tailpos = 3;

pd->edge[15].head = 2;
pd->edge[15].headpos = 1;
pd->edge[15].tail = 7;
pd->edge[15].tailpos = 1;

pd->edge[16].head = 1;
pd->edge[16].headpos = 3;
pd->edge[16].tail = 0;
pd->edge[16].tailpos = 1;

pd->edge[17].head = 8;
pd->edge[17].headpos = 1;
pd->edge[17].tail = 1;
pd->edge[17].tailpos = 1;

pd->edge[18].head = 9;
pd->edge[18].headpos = 3;
pd->edge[18].tail = 8;
pd->edge[18].tailpos = 3;

pd->edge[19].head = 0;
pd->edge[19].headpos = 3;
pd->edge[19].tail = 9;
pd->edge[19].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 12;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;
pd->comp[0].edge[10] = 10;
pd->comp[0].edge[11] = 11;

pd->comp[1].nedges = 4;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 12;
pd->comp[1].edge[1] = 13;
pd->comp[1].edge[2] = 14;
pd->comp[1].edge[3] = 15;

pd->comp[2].nedges = 4;
pd->comp[2].tag = 'C';

pd->comp[2].edge = calloc(pd->comp[2].nedges,sizeof(pd_idx_t));
assert(pd->comp[2].edge != NULL);

pd->comp[2].edge[0] = 16;
pd->comp[2].edge[1] = 17;
pd->comp[2].edge[2] = 18;
pd->comp[2].edge[3] = 19;


/* Face data */

pd->face[0].nedges = 5;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 12;
pd->face[0].or[1] = 1;

pd->face[0].edge[2] = 3;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 11;
pd->face[0].or[3] = 1;

pd->face[0].edge[4] = 16;
pd->face[0].or[4] = 1;

pd->face[1].nedges = 5;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 5;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 14;
pd->face[1].or[1] = 1;

pd->face[1].edge[2] = 7;
pd->face[1].or[2] = 1;

pd->face[1].edge[3] = 18;
pd->face[1].or[3] = 1;

pd->face[1].edge[4] = 9;
pd->face[1].or[4] = 1;

pd->face[2].nedges = 4;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 19;
pd->face[2].or[1] = 0;

pd->face[2].edge[2] = 8;
pd->face[2].or[2] = 0;

pd->face[2].edge[3] = 17;
pd->face[2].or[3] = 0;

pd->face[3].nedges = 4;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 1;
pd->face[3].or[0] = 0;

pd->face[3].edge[1] = 17;
pd->face[3].or[1] = 1;

pd->face[3].edge[2] = 7;
pd->face[3].or[2] = 0;

pd->face[3].edge[3] = 15;
pd->face[3].or[3] = 1;

pd->face[4].nedges = 4;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 2;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 15;
pd->face[4].or[1] = 0;

pd->face[4].edge[2] = 6;
pd->face[4].or[2] = 0;

pd->face[4].edge[3] = 13;
pd->face[4].or[3] = 0;

pd->face[5].nedges = 4;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 3;
pd->face[5].or[0] = 0;

pd->face[5].edge[1] = 13;
pd->face[5].or[1] = 1;

pd->face[5].edge[2] = 5;
pd->face[5].or[2] = 0;

pd->face[5].edge[3] = 10;
pd->face[5].or[3] = 1;

pd->face[6].nedges = 4;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 4;
pd->face[6].or[0] = 1;

pd->face[6].edge[1] = 9;
pd->face[6].or[1] = 0;

pd->face[6].edge[2] = 19;
pd->face[6].or[2] = 1;

pd->face[6].edge[3] = 11;
pd->face[6].or[3] = 0;

pd->face[7].nedges = 2;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 0;
pd->face[7].or[0] = 1;

pd->face[7].edge[1] = 16;
pd->face[7].or[1] = 0;

pd->face[8].nedges = 2;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].or = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].or != NULL);

pd->face[8].edge[0] = 2;
pd->face[8].or[0] = 1;

pd->face[8].edge[1] = 12;
pd->face[8].or[1] = 0;

pd->face[9].nedges = 2;
pd->face[9].edge = calloc(pd->face[9].nedges,sizeof(pd_idx_t));
pd->face[9].or = calloc(pd->face[9].nedges,sizeof(pd_or_t));
assert(pd->face[9].edge != NULL);
assert(pd->face[9].or != NULL);

pd->face[9].edge[0] = 4;
pd->face[9].or[0] = 0;

pd->face[9].edge[1] = 10;
pd->face[9].or[1] = 0;

pd->face[10].nedges = 2;
pd->face[10].edge = calloc(pd->face[10].nedges,sizeof(pd_idx_t));
pd->face[10].or = calloc(pd->face[10].nedges,sizeof(pd_or_t));
assert(pd->face[10].edge != NULL);
assert(pd->face[10].or != NULL);

pd->face[10].edge[0] = 6;
pd->face[10].or[0] = 1;

pd->face[10].edge[1] = 14;
pd->face[10].or[1] = 0;

pd->face[11].nedges = 2;
pd->face[11].edge = calloc(pd->face[11].nedges,sizeof(pd_idx_t));
pd->face[11].or = calloc(pd->face[11].nedges,sizeof(pd_or_t));
assert(pd->face[11].edge != NULL);
assert(pd->face[11].or != NULL);

pd->face[11].edge[0] = 8;
pd->face[11].or[0] = 1;

pd->face[11].edge[1] = 18;
pd->face[11].or[1] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



pd_code_t *pd_create_r2_testH_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(8);
assert(pd != NULL);
pd->ncross = 8;
pd->nedges = 16;
pd->ncomps = 3;
pd->nfaces = 10;
sprintf(pd->hash,"%s","CBAKBgUEBAMCAgICAgMKBAIAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 10;
pd->cross[0].edge[2] = 9;
pd->cross[0].edge[3] = 13;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 11;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 10;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 8;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 9;
pd->cross[2].sign = 0;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 8;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 7;
pd->cross[3].sign = 0;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 14;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 15;
pd->cross[4].sign = 0;

pd->cross[5].edge[0] = 4;
pd->cross[5].edge[1] = 14;
pd->cross[5].edge[2] = 5;
pd->cross[5].edge[3] = 15;
pd->cross[5].sign = 0;

pd->cross[6].edge[0] = 5;
pd->cross[6].edge[1] = 11;
pd->cross[6].edge[2] = 6;
pd->cross[6].edge[3] = 12;
pd->cross[6].sign = 0;

pd->cross[7].edge[0] = 6;
pd->cross[7].edge[1] = 13;
pd->cross[7].edge[2] = 7;
pd->cross[7].edge[3] = 12;
pd->cross[7].sign = 0;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 5;
pd->edge[4].headpos = 0;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 6;
pd->edge[5].headpos = 0;
pd->edge[5].tail = 5;
pd->edge[5].tailpos = 2;

pd->edge[6].head = 7;
pd->edge[6].headpos = 0;
pd->edge[6].tail = 6;
pd->edge[6].tailpos = 2;

pd->edge[7].head = 3;
pd->edge[7].headpos = 3;
pd->edge[7].tail = 7;
pd->edge[7].tailpos = 2;

pd->edge[8].head = 2;
pd->edge[8].headpos = 1;
pd->edge[8].tail = 3;
pd->edge[8].tailpos = 1;

pd->edge[9].head = 0;
pd->edge[9].headpos = 2;
pd->edge[9].tail = 2;
pd->edge[9].tailpos = 3;

pd->edge[10].head = 1;
pd->edge[10].headpos = 3;
pd->edge[10].tail = 0;
pd->edge[10].tailpos = 1;

pd->edge[11].head = 6;
pd->edge[11].headpos = 1;
pd->edge[11].tail = 1;
pd->edge[11].tailpos = 1;

pd->edge[12].head = 7;
pd->edge[12].headpos = 3;
pd->edge[12].tail = 6;
pd->edge[12].tailpos = 3;

pd->edge[13].head = 0;
pd->edge[13].headpos = 3;
pd->edge[13].tail = 7;
pd->edge[13].tailpos = 1;

pd->edge[14].head = 4;
pd->edge[14].headpos = 1;
pd->edge[14].tail = 5;
pd->edge[14].tailpos = 1;

pd->edge[15].head = 5;
pd->edge[15].headpos = 3;
pd->edge[15].tail = 4;
pd->edge[15].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 10;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;

pd->comp[1].nedges = 4;
pd->comp[1].tag = 'C';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 10;
pd->comp[1].edge[1] = 11;
pd->comp[1].edge[2] = 12;
pd->comp[1].edge[3] = 13;

pd->comp[2].nedges = 2;
pd->comp[2].tag = 'B';

pd->comp[2].edge = calloc(pd->comp[2].nedges,sizeof(pd_idx_t));
assert(pd->comp[2].edge != NULL);

pd->comp[2].edge[0] = 14;
pd->comp[2].edge[1] = 15;


/* Face data */

pd->face[0].nedges = 6;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 11;
pd->face[0].or[1] = 1;

pd->face[0].edge[2] = 5;
pd->face[0].or[2] = 0;

pd->face[0].edge[3] = 14;
pd->face[0].or[3] = 1;

pd->face[0].edge[4] = 3;
pd->face[0].or[4] = 0;

pd->face[0].edge[5] = 8;
pd->face[0].or[5] = 1;

pd->face[1].nedges = 5;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 3;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 15;
pd->face[1].or[1] = 1;

pd->face[1].edge[2] = 5;
pd->face[1].or[2] = 1;

pd->face[1].edge[3] = 12;
pd->face[1].or[3] = 1;

pd->face[1].edge[4] = 7;
pd->face[1].or[4] = 1;

pd->face[2].nedges = 4;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 13;
pd->face[2].or[1] = 0;

pd->face[2].edge[2] = 6;
pd->face[2].or[2] = 0;

pd->face[2].edge[3] = 11;
pd->face[2].or[3] = 0;

pd->face[3].nedges = 4;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 2;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 7;
pd->face[3].or[1] = 0;

pd->face[3].edge[2] = 13;
pd->face[3].or[2] = 1;

pd->face[3].edge[3] = 9;
pd->face[3].or[3] = 0;

pd->face[4].nedges = 3;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 1;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 9;
pd->face[4].or[1] = 1;

pd->face[4].edge[2] = 10;
pd->face[4].or[2] = 1;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 0;
pd->face[5].or[0] = 1;

pd->face[5].edge[1] = 10;
pd->face[5].or[1] = 0;

pd->face[6].nedges = 2;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 2;
pd->face[6].or[0] = 0;

pd->face[6].edge[1] = 8;
pd->face[6].or[1] = 0;

pd->face[7].nedges = 2;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 4;
pd->face[7].or[0] = 0;

pd->face[7].edge[1] = 14;
pd->face[7].or[1] = 0;

pd->face[8].nedges = 2;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].or = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].or != NULL);

pd->face[8].edge[0] = 4;
pd->face[8].or[0] = 1;

pd->face[8].edge[1] = 15;
pd->face[8].or[1] = 0;

pd->face[9].nedges = 2;
pd->face[9].edge = calloc(pd->face[9].nedges,sizeof(pd_idx_t));
pd->face[9].or = calloc(pd->face[9].nedges,sizeof(pd_or_t));
assert(pd->face[9].edge != NULL);
assert(pd->face[9].or != NULL);

pd->face[9].edge[0] = 6;
pd->face[9].or[0] = 1;

pd->face[9].edge[1] = 12;
pd->face[9].or[1] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



/*                    ____>________________                                   */
/*                   /    ____<_1______    \                                  */
/*      (1)         /    /          ___\ _17______________                    */
/*                 /    |          ^   1^    \8           \                   */
/*                7     |   (0)   /     0 (2) \     (11)   \                  */
/*               /  (3) |       16  (7) |      8           18                 */
/*              /       \         \______<_19___\ _____<___/                  */
/*    ____>___ /___15_>_ \ _____       /0       9\                            */
/*   /        /7    (4)  2\ (8) \     ^           \                           */
/* 14  (10)  6             2    12   /             |                          */
/*   \____________<_13______\ __<  11              |                          */
/*          6|              3\     /               |                          */
/*           |                3   /                |                          */
/*           |                 > /                 9                          */
/*           ^                  /                  |                          */
/*           |      (5)        /4\       (6)       |                          */
/*           |                /   \                |                          */
/*           |               ^     \               |                          */
/*           |              /       \              |                          */
/*           5             /         4             |                          */
/*           \            /    (9)    \            /                          */
/*            \          10           /           /                           */
/*             \           \         /           /                            */
/*              \           \       /           /                             */
/*               \           \     <           /                              */
/*                \           \   /           /                               */
/*                 \           \ /           /                                */
/*                  \           /           /                                 */
/*                   \_________/5\____<____/                                  */
/*                                                                            */
/* R2 on face 8 at crossings 2 and 3                                          */
/*                    ___5_____>___________                                   */
/*                   /    _____<_1_____    \                                  */
/*                  /    /          >__\ _11______________                    */
/*       (1)       /    |          /   1|    \6           \                   */
/*                /     |         /     0 (2) \     (9)    \                  */
/*               /      |       10  (5) |      6           12                 */
/*              /       \         \________<_13_\ __<______/                  */
/*    _____>__ /_______  \    (4)      /0       7\                            */
/*   /        /5       \  \           /           \                           */
/* 15   (8)  4   (7)   14  \         ^             |                          */
/*   \_____________<___/    \       /              |                          */
/*           |4              \     9               |                          */
/*           |                \   /                |                          */
/*           ^                 \ /        (3)      |                          */
/*           |                  /                  |                          */
/*           |      (0)        /2\                 |                          */
/*           |                /   \                |                          */
/*           |               ^     \               |                          */
/*           3              /       \              |                          */
/*           |             /         2             7                          */
/*           \            /     (6)   \            /                          */
/*            \           8           /           /                           */
/*             \           \         /           /                            */
/*              \           \       /           /                             */
/*               \           \     /           /                              */
/*                \           \   <           /                               */
/*                 \           \ /           /                                */
/*                  \           /           /                                 */
/*                   \_________/3\___<_____/                                  */
/*                                                                            */



    bool r2testH() {
  
  printf("---------------------------------------\n"
	 "r2 test H\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r2_testH_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after-0 configuration...");
  pd_code_t *after0 = pd_create_r2_testH_after_0();
  printf("done (passes pd_ok)\n");

  pd_code_t *after[1] = {after0};
    
  printf("performing r2 at crossings 2 and 3...\n");

  pd_idx_t  noutpd;
  pd_code_t **newpd;
  pd_idx_t  cr[2] = {2,3};

  pd_R2_bigon_elimination(before,cr,&noutpd,&newpd);

  if (noutpd != 1) { 
    printf("fail (did not generate 1 output pd codes)\n");
    return false;
  }

  pd_idx_t ct;

  for(ct=0;ct<noutpd;ct++) { 

    if (!pd_ok(newpd[ct])) { 
      printf("fail (output pd %d does not pass pd_ok)\n",ct);
      return false;
    }
    printf("\t output pd %d passes pd_ok...pass\n",ct);
    printf("\t testing for isomorphism with after configuration %d...",ct);
    if (!pd_isomorphic(newpd[ct],after[ct])) { 
      pd_printf("fail (output pd %PD ",newpd[ct]);
      pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after[ct]);
      return false;
    }

    printf("pass (output %d and after %d isomorphic)\n",ct,ct);

  }

  printf("performing r2 at crossings 2 and 3...pass\n");
  printf("housecleaning...");

  pd_code_free(&before);
  pd_code_free(&after0);
  pd_code_free(&newpd[0]);

  free(newpd);

  printf("done\n");

  printf("---------------------------------------\n"
	 "r2 test H : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}

/* r2 test I assembled c code */

pd_code_t *pd_create_r2_testI_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(10);
assert(pd != NULL);
pd->ncross = 10;
pd->nedges = 20;
pd->ncomps = 6;
pd->nfaces = 12;
sprintf(pd->hash,"%s","ChQMCgQEBAQCAgICAgICBgQEBAQCAgA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 16;
pd->cross[0].edge[2] = 3;
pd->cross[0].edge[3] = 17;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 17;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 16;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 4;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 5;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 4;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 7;
pd->cross[3].sign = 1;

pd->cross[4].edge[0] = 5;
pd->cross[4].edge[1] = 8;
pd->cross[4].edge[2] = 6;
pd->cross[4].edge[3] = 9;
pd->cross[4].sign = 1;

pd->cross[5].edge[0] = 6;
pd->cross[5].edge[1] = 8;
pd->cross[5].edge[2] = 7;
pd->cross[5].edge[3] = 11;
pd->cross[5].sign = 0;

pd->cross[6].edge[0] = 9;
pd->cross[6].edge[1] = 12;
pd->cross[6].edge[2] = 10;
pd->cross[6].edge[3] = 13;
pd->cross[6].sign = 1;

pd->cross[7].edge[0] = 10;
pd->cross[7].edge[1] = 12;
pd->cross[7].edge[2] = 11;
pd->cross[7].edge[3] = 15;
pd->cross[7].sign = 1;

pd->cross[8].edge[0] = 13;
pd->cross[8].edge[1] = 18;
pd->cross[8].edge[2] = 14;
pd->cross[8].edge[3] = 19;
pd->cross[8].sign = 1;

pd->cross[9].edge[0] = 14;
pd->cross[9].edge[1] = 18;
pd->cross[9].edge[2] = 15;
pd->cross[9].edge[3] = 19;
pd->cross[9].sign = 1;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 0;
pd->edge[3].headpos = 2;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 2;
pd->edge[4].headpos = 1;
pd->edge[4].tail = 3;
pd->edge[4].tailpos = 1;

pd->edge[5].head = 4;
pd->edge[5].headpos = 0;
pd->edge[5].tail = 2;
pd->edge[5].tailpos = 3;

pd->edge[6].head = 5;
pd->edge[6].headpos = 0;
pd->edge[6].tail = 4;
pd->edge[6].tailpos = 2;

pd->edge[7].head = 3;
pd->edge[7].headpos = 3;
pd->edge[7].tail = 5;
pd->edge[7].tailpos = 2;

pd->edge[8].head = 4;
pd->edge[8].headpos = 1;
pd->edge[8].tail = 5;
pd->edge[8].tailpos = 1;

pd->edge[9].head = 6;
pd->edge[9].headpos = 0;
pd->edge[9].tail = 4;
pd->edge[9].tailpos = 3;

pd->edge[10].head = 7;
pd->edge[10].headpos = 0;
pd->edge[10].tail = 6;
pd->edge[10].tailpos = 2;

pd->edge[11].head = 5;
pd->edge[11].headpos = 3;
pd->edge[11].tail = 7;
pd->edge[11].tailpos = 2;

pd->edge[12].head = 6;
pd->edge[12].headpos = 1;
pd->edge[12].tail = 7;
pd->edge[12].tailpos = 1;

pd->edge[13].head = 8;
pd->edge[13].headpos = 0;
pd->edge[13].tail = 6;
pd->edge[13].tailpos = 3;

pd->edge[14].head = 9;
pd->edge[14].headpos = 0;
pd->edge[14].tail = 8;
pd->edge[14].tailpos = 2;

pd->edge[15].head = 7;
pd->edge[15].headpos = 3;
pd->edge[15].tail = 9;
pd->edge[15].tailpos = 2;

pd->edge[16].head = 0;
pd->edge[16].headpos = 1;
pd->edge[16].tail = 1;
pd->edge[16].tailpos = 3;

pd->edge[17].head = 1;
pd->edge[17].headpos = 1;
pd->edge[17].tail = 0;
pd->edge[17].tailpos = 3;

pd->edge[18].head = 8;
pd->edge[18].headpos = 1;
pd->edge[18].tail = 9;
pd->edge[18].tailpos = 1;

pd->edge[19].head = 9;
pd->edge[19].headpos = 3;
pd->edge[19].tail = 8;
pd->edge[19].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 4;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;

pd->comp[1].nedges = 4;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 4;
pd->comp[1].edge[1] = 5;
pd->comp[1].edge[2] = 6;
pd->comp[1].edge[3] = 7;

pd->comp[2].nedges = 4;
pd->comp[2].tag = 'C';

pd->comp[2].edge = calloc(pd->comp[2].nedges,sizeof(pd_idx_t));
assert(pd->comp[2].edge != NULL);

pd->comp[2].edge[0] = 8;
pd->comp[2].edge[1] = 9;
pd->comp[2].edge[2] = 10;
pd->comp[2].edge[3] = 11;

pd->comp[3].nedges = 4;
pd->comp[3].tag = 'D';

pd->comp[3].edge = calloc(pd->comp[3].nedges,sizeof(pd_idx_t));
assert(pd->comp[3].edge != NULL);

pd->comp[3].edge[0] = 12;
pd->comp[3].edge[1] = 13;
pd->comp[3].edge[2] = 14;
pd->comp[3].edge[3] = 15;

pd->comp[4].nedges = 2;
pd->comp[4].tag = 'E';

pd->comp[4].edge = calloc(pd->comp[4].nedges,sizeof(pd_idx_t));
assert(pd->comp[4].edge != NULL);

pd->comp[4].edge[0] = 16;
pd->comp[4].edge[1] = 17;

pd->comp[5].nedges = 2;
pd->comp[5].tag = 'F';

pd->comp[5].edge = calloc(pd->comp[5].nedges,sizeof(pd_idx_t));
assert(pd->comp[5].edge != NULL);

pd->comp[5].edge[0] = 18;
pd->comp[5].edge[1] = 19;


/* Face data */

pd->face[0].nedges = 10;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 5;
pd->face[0].or[1] = 1;

pd->face[0].edge[2] = 9;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 13;
pd->face[0].or[3] = 1;

pd->face[0].edge[4] = 19;
pd->face[0].or[4] = 1;

pd->face[0].edge[5] = 15;
pd->face[0].or[5] = 1;

pd->face[0].edge[6] = 11;
pd->face[0].or[6] = 1;

pd->face[0].edge[7] = 7;
pd->face[0].or[7] = 1;

pd->face[0].edge[8] = 3;
pd->face[0].or[8] = 1;

pd->face[0].edge[9] = 16;
pd->face[0].or[9] = 0;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 17;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 3;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 4;
pd->face[1].or[3] = 1;

pd->face[2].nedges = 4;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 2;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 7;
pd->face[2].or[1] = 0;

pd->face[2].edge[2] = 8;
pd->face[2].or[2] = 1;

pd->face[2].edge[3] = 5;
pd->face[2].or[3] = 0;

pd->face[3].nedges = 4;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 6;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 11;
pd->face[3].or[1] = 0;

pd->face[3].edge[2] = 12;
pd->face[3].or[2] = 1;

pd->face[3].edge[3] = 9;
pd->face[3].or[3] = 0;

pd->face[4].nedges = 4;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 10;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 15;
pd->face[4].or[1] = 0;

pd->face[4].edge[2] = 18;
pd->face[4].or[2] = 1;

pd->face[4].edge[3] = 13;
pd->face[4].or[3] = 0;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 0;
pd->face[5].or[0] = 1;

pd->face[5].edge[1] = 16;
pd->face[5].or[1] = 1;

pd->face[6].nedges = 2;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 0;
pd->face[6].or[0] = 0;

pd->face[6].edge[1] = 17;
pd->face[6].or[1] = 1;

pd->face[7].nedges = 2;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 2;
pd->face[7].or[0] = 0;

pd->face[7].edge[1] = 4;
pd->face[7].or[1] = 0;

pd->face[8].nedges = 2;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].or = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].or != NULL);

pd->face[8].edge[0] = 6;
pd->face[8].or[0] = 0;

pd->face[8].edge[1] = 8;
pd->face[8].or[1] = 0;

pd->face[9].nedges = 2;
pd->face[9].edge = calloc(pd->face[9].nedges,sizeof(pd_idx_t));
pd->face[9].or = calloc(pd->face[9].nedges,sizeof(pd_or_t));
assert(pd->face[9].edge != NULL);
assert(pd->face[9].or != NULL);

pd->face[9].edge[0] = 10;
pd->face[9].or[0] = 0;

pd->face[9].edge[1] = 12;
pd->face[9].or[1] = 0;

pd->face[10].nedges = 2;
pd->face[10].edge = calloc(pd->face[10].nedges,sizeof(pd_idx_t));
pd->face[10].or = calloc(pd->face[10].nedges,sizeof(pd_or_t));
assert(pd->face[10].edge != NULL);
assert(pd->face[10].or != NULL);

pd->face[10].edge[0] = 14;
pd->face[10].or[0] = 0;

pd->face[10].edge[1] = 18;
pd->face[10].or[1] = 0;

pd->face[11].nedges = 2;
pd->face[11].edge = calloc(pd->face[11].nedges,sizeof(pd_idx_t));
pd->face[11].or = calloc(pd->face[11].nedges,sizeof(pd_or_t));
assert(pd->face[11].edge != NULL);
assert(pd->face[11].or != NULL);

pd->face[11].edge[0] = 14;
pd->face[11].or[0] = 1;

pd->face[11].edge[1] = 19;
pd->face[11].or[1] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



pd_code_t *pd_create_r2_testI_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(4);
assert(pd != NULL);
pd->ncross = 4;
pd->nedges = 8;
pd->ncomps = 3;
pd->nfaces = 6;
sprintf(pd->hash,"%s","BAgGBAQCAgICAwQCAgAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 4;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 5;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 5;
pd->cross[1].edge[2] = 3;
pd->cross[1].edge[3] = 4;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 6;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 7;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 6;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 7;
pd->cross[3].sign = 1;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 1;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 1;
pd->edge[3].headpos = 2;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 0;
pd->edge[4].headpos = 1;
pd->edge[4].tail = 1;
pd->edge[4].tailpos = 3;

pd->edge[5].head = 1;
pd->edge[5].headpos = 1;
pd->edge[5].tail = 0;
pd->edge[5].tailpos = 3;

pd->edge[6].head = 2;
pd->edge[6].headpos = 1;
pd->edge[6].tail = 3;
pd->edge[6].tailpos = 1;

pd->edge[7].head = 3;
pd->edge[7].headpos = 3;
pd->edge[7].tail = 2;
pd->edge[7].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 4;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;

pd->comp[1].nedges = 2;
pd->comp[1].tag = 'E';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 4;
pd->comp[1].edge[1] = 5;

pd->comp[2].nedges = 2;
pd->comp[2].tag = 'B';

pd->comp[2].edge = calloc(pd->comp[2].nedges,sizeof(pd_idx_t));
assert(pd->comp[2].edge != NULL);

pd->comp[2].edge[0] = 6;
pd->comp[2].edge[1] = 7;


/* Face data */

pd->face[0].nedges = 4;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 4;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 3;
pd->face[0].or[2] = 0;

pd->face[0].edge[3] = 6;
pd->face[0].or[3] = 1;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 7;
pd->face[1].or[1] = 1;

pd->face[1].edge[2] = 3;
pd->face[1].or[2] = 1;

pd->face[1].edge[3] = 5;
pd->face[1].or[3] = 0;

pd->face[2].nedges = 2;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 4;
pd->face[2].or[1] = 1;

pd->face[3].nedges = 2;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 0;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 5;
pd->face[3].or[1] = 1;

pd->face[4].nedges = 2;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 2;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 6;
pd->face[4].or[1] = 0;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 2;
pd->face[5].or[0] = 1;

pd->face[5].edge[1] = 7;
pd->face[5].or[1] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}

pd_code_t *pd_create_r2_testI_after_1() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(4);
assert(pd != NULL);
pd->ncross = 4;
pd->nedges = 8;
pd->ncomps = 3;
pd->nfaces = 6;
sprintf(pd->hash,"%s","BAgGBAQCAgICAwQCAgAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 4;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 5;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 5;
pd->cross[1].edge[2] = 3;
pd->cross[1].edge[3] = 4;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 6;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 7;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 6;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 7;
pd->cross[3].sign = 1;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 1;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 1;
pd->edge[3].headpos = 2;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 0;
pd->edge[4].headpos = 1;
pd->edge[4].tail = 1;
pd->edge[4].tailpos = 3;

pd->edge[5].head = 1;
pd->edge[5].headpos = 1;
pd->edge[5].tail = 0;
pd->edge[5].tailpos = 3;

pd->edge[6].head = 2;
pd->edge[6].headpos = 1;
pd->edge[6].tail = 3;
pd->edge[6].tailpos = 1;

pd->edge[7].head = 3;
pd->edge[7].headpos = 3;
pd->edge[7].tail = 2;
pd->edge[7].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 4;
pd->comp[0].tag = 'D';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;

pd->comp[1].nedges = 2;
pd->comp[1].tag = 'C';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 4;
pd->comp[1].edge[1] = 5;

pd->comp[2].nedges = 2;
pd->comp[2].tag = 'F';

pd->comp[2].edge = calloc(pd->comp[2].nedges,sizeof(pd_idx_t));
assert(pd->comp[2].edge != NULL);

pd->comp[2].edge[0] = 6;
pd->comp[2].edge[1] = 7;


/* Face data */

pd->face[0].nedges = 4;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 4;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 3;
pd->face[0].or[2] = 0;

pd->face[0].edge[3] = 6;
pd->face[0].or[3] = 1;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 7;
pd->face[1].or[1] = 1;

pd->face[1].edge[2] = 3;
pd->face[1].or[2] = 1;

pd->face[1].edge[3] = 5;
pd->face[1].or[3] = 0;

pd->face[2].nedges = 2;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 4;
pd->face[2].or[1] = 1;

pd->face[3].nedges = 2;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 0;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 5;
pd->face[3].or[1] = 1;

pd->face[4].nedges = 2;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 2;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 6;
pd->face[4].or[1] = 0;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 2;
pd->face[5].or[0] = 1;

pd->face[5].edge[1] = 7;
pd->face[5].or[1] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



/*                                                                            */
/*                                                                            */
/*                                                           (0)              */
/*     __<__   _1_>_   _5_>_   _9_>_   _13>_   _____                          */
/*    /     \1/     \2/     \4/     \6/     \8/     \                         */
/*   /       \  (1)  \  (2)  \  (3)  \  (4)  \       \                        */
/*  /       / \     / \     / \     / \     / \       \                       */
/* /       /   \   /   \   /   \   ^   \   /   \       \                      */
/* |      |     | |     | |     | |     | ^     |       |                     */
/* |      |     | |     | ^     | |     | |     |       |                     */
/* |      ^     ^ ^     | |     | |     | |     |       |                     */
/* |      |     | |     | |     | 12    | |    14       |                     */
/* | (5)  | (6) | 4 (7) | | (8) | | (9) | |(10) | (11)  |                     */
/* |      |     | |     | 8     6 |    10 |     |       |                     */
/*16      |    17 |     | |     | |     | |     |      19                     */
/* |      0     | |     2 |     | |     | |     |       |                     */
/* |      |     | |     | |     | |     | 18    |       |                     */
/* |      |     | |     | |     | |     | |     |       |                     */
/* |      |     | |     | |     | |     | |     |       |                     */
/* |      |     | |     | |     | |     | |     |       |                     */
/* |      |     | |     | |     | |     | |     |       |                     */
/* |      |     | |     | |     | |     | |     |       |                     */
/* \       \   /   \   <   \   <   \   <   \   <       /                      */
/*  \       \ /     \ /     \ /     \ /     \ /       /                       */
/*   \       \       \       /       \       \       /                        */
/*    \__>__/0\_<_3_/3\_<_7_/5\_<11_/7\_<15_/9\__<__/                         */
/*                                                                            */
/*                                                                            */
/*  R2 on face 8 at crossings 4 and 5                                         */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*     _____   __1_>   _____      (1)                                         */
/*    /     \0/     \2/     \                                                 */
/*   /       \  (0)  \       \                                                */
/*  /       / \     / \       \                                               */
/* /       /   \   /   \       \                                              */
/* |      ^     ^ ^     |       |                                             */
/* |      |     | |     |       |                                             */
/* |      |     | |     |       |                                             */
/* |      |     | |     |       |                                             */
/* |      |     | |     2       |                                             */
/* 5  (3) | (2) | 6 (4) |  (5)  7                    x2                       */
/* |      0     | |     |       |                                             */
/* |      |     4 |     |       |                                             */
/* |      |     | |     |       |                                             */
/* |      |     | |     |       |                                             */
/* |      |     | |     |       |                                             */
/* |      |     | |     |       |                                             */
/* |      |     | |     |       |                                             */
/* |      |     | |     |       |                                             */
/* \       \   /   \   <       /                                              */
/*  \       \ /     \ /       /                                               */
/*   \       \       \       /                                                */
/*    \___>_/1\_<_3_/3\__<__/                                                 */
/*                                                                            */
/*                                                                            */



    bool r2testI() {
  
  printf("---------------------------------------\n"
	 "r2 test I\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r2_testI_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after-0 configuration...");
  pd_code_t *after0 = pd_create_r2_testI_after_0();
  printf("done (passes pd_ok)\n");

  printf("creating after-1 configuration...");
  pd_code_t *after1 = pd_create_r2_testI_after_1();
  printf("done (passes pd_ok)\n");

  pd_code_t *after[2] = {after0,after1};
    
  printf("performing r2 at crossings 4 and 5...\n");

  pd_idx_t  noutpd;
  pd_code_t **newpd;
  pd_idx_t  cr[2] = { 4,5 };

  pd_R2_bigon_elimination(before,cr,&noutpd,&newpd);

  if (noutpd != 2) { 
    printf("fail (did not generate 2 output pd codes)\n");
    return false;
  }

  pd_idx_t ct;

  for(ct=0;ct<noutpd;ct++) { 

    if (!pd_ok(newpd[ct])) { 
      printf("fail (output pd %d does not pass pd_ok)\n",ct);
      return false;
    }
    printf("\t output pd %d passes pd_ok...pass\n",ct);
    printf("\t testing for isomorphism with after configuration %d...",ct);
    if (!pd_isomorphic(newpd[ct],after[ct])) { 
      pd_printf("fail (output pd %PD ",newpd[ct]);
      pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after[ct]);
      return false;
    }

    printf("pass (output %d and after %d isomorphic)\n",ct,ct);

  }

  printf("performing r2 at crossings 4 and 5...pass\n");
  printf("housecleaning...");

  pd_code_free(&before);
  pd_code_free(&after0);
  pd_code_free(&after1);
  pd_code_free(&newpd[0]);
  pd_code_free(&newpd[1]);

  free(newpd);

  printf("done\n");

  printf("---------------------------------------\n"
	 "r2 test I : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}

/* r2 test J assembled c code */

pd_code_t *pd_create_r2_testJ_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(11);
assert(pd != NULL);
pd->ncross = 11;
pd->nedges = 22;
pd->ncomps = 4;
pd->nfaces = 13;
sprintf(pd->hash,"%s","CxYNBQUEBAQEAwMDAwICAgQKBAQEAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 5;
pd->cross[0].edge[2] = 9;
pd->cross[0].edge[3] = 4;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 21;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 18;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 16;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 15;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 16;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 17;
pd->cross[3].sign = 1;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 21;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 20;
pd->cross[4].sign = 0;

pd->cross[5].edge[0] = 5;
pd->cross[5].edge[1] = 18;
pd->cross[5].edge[2] = 6;
pd->cross[5].edge[3] = 19;
pd->cross[5].sign = 1;

pd->cross[6].edge[0] = 6;
pd->cross[6].edge[1] = 11;
pd->cross[6].edge[2] = 7;
pd->cross[6].edge[3] = 12;
pd->cross[6].sign = 1;

pd->cross[7].edge[0] = 7;
pd->cross[7].edge[1] = 13;
pd->cross[7].edge[2] = 8;
pd->cross[7].edge[3] = 12;
pd->cross[7].sign = 1;

pd->cross[8].edge[0] = 8;
pd->cross[8].edge[1] = 20;
pd->cross[8].edge[2] = 9;
pd->cross[8].edge[3] = 19;
pd->cross[8].sign = 1;

pd->cross[9].edge[0] = 10;
pd->cross[9].edge[1] = 14;
pd->cross[9].edge[2] = 11;
pd->cross[9].edge[3] = 15;
pd->cross[9].sign = 0;

pd->cross[10].edge[0] = 10;
pd->cross[10].edge[1] = 17;
pd->cross[10].edge[2] = 13;
pd->cross[10].edge[3] = 14;
pd->cross[10].sign = 1;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 0;
pd->edge[4].headpos = 3;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 5;
pd->edge[5].headpos = 0;
pd->edge[5].tail = 0;
pd->edge[5].tailpos = 1;

pd->edge[6].head = 6;
pd->edge[6].headpos = 0;
pd->edge[6].tail = 5;
pd->edge[6].tailpos = 2;

pd->edge[7].head = 7;
pd->edge[7].headpos = 0;
pd->edge[7].tail = 6;
pd->edge[7].tailpos = 2;

pd->edge[8].head = 8;
pd->edge[8].headpos = 0;
pd->edge[8].tail = 7;
pd->edge[8].tailpos = 2;

pd->edge[9].head = 0;
pd->edge[9].headpos = 2;
pd->edge[9].tail = 8;
pd->edge[9].tailpos = 2;

pd->edge[10].head = 9;
pd->edge[10].headpos = 0;
pd->edge[10].tail = 10;
pd->edge[10].tailpos = 0;

pd->edge[11].head = 6;
pd->edge[11].headpos = 1;
pd->edge[11].tail = 9;
pd->edge[11].tailpos = 2;

pd->edge[12].head = 7;
pd->edge[12].headpos = 3;
pd->edge[12].tail = 6;
pd->edge[12].tailpos = 3;

pd->edge[13].head = 10;
pd->edge[13].headpos = 2;
pd->edge[13].tail = 7;
pd->edge[13].tailpos = 1;

pd->edge[14].head = 9;
pd->edge[14].headpos = 1;
pd->edge[14].tail = 10;
pd->edge[14].tailpos = 3;

pd->edge[15].head = 2;
pd->edge[15].headpos = 3;
pd->edge[15].tail = 9;
pd->edge[15].tailpos = 3;

pd->edge[16].head = 3;
pd->edge[16].headpos = 1;
pd->edge[16].tail = 2;
pd->edge[16].tailpos = 1;

pd->edge[17].head = 10;
pd->edge[17].headpos = 1;
pd->edge[17].tail = 3;
pd->edge[17].tailpos = 3;

pd->edge[18].head = 5;
pd->edge[18].headpos = 1;
pd->edge[18].tail = 1;
pd->edge[18].tailpos = 3;

pd->edge[19].head = 8;
pd->edge[19].headpos = 3;
pd->edge[19].tail = 5;
pd->edge[19].tailpos = 3;

pd->edge[20].head = 4;
pd->edge[20].headpos = 3;
pd->edge[20].tail = 8;
pd->edge[20].tailpos = 1;

pd->edge[21].head = 1;
pd->edge[21].headpos = 1;
pd->edge[21].tail = 4;
pd->edge[21].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 10;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;

pd->comp[1].nedges = 4;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 10;
pd->comp[1].edge[1] = 11;
pd->comp[1].edge[2] = 12;
pd->comp[1].edge[3] = 13;

pd->comp[2].nedges = 4;
pd->comp[2].tag = 'C';

pd->comp[2].edge = calloc(pd->comp[2].nedges,sizeof(pd_idx_t));
assert(pd->comp[2].edge != NULL);

pd->comp[2].edge[0] = 14;
pd->comp[2].edge[1] = 15;
pd->comp[2].edge[2] = 16;
pd->comp[2].edge[3] = 17;

pd->comp[3].nedges = 4;
pd->comp[3].tag = 'D';

pd->comp[3].edge = calloc(pd->comp[3].nedges,sizeof(pd_idx_t));
assert(pd->comp[3].edge != NULL);

pd->comp[3].edge[0] = 18;
pd->comp[3].edge[1] = 19;
pd->comp[3].edge[2] = 20;
pd->comp[3].edge[3] = 21;


/* Face data */

pd->face[0].nedges = 5;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 15;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 11;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 6;
pd->face[0].or[3] = 0;

pd->face[0].edge[4] = 18;
pd->face[0].or[4] = 0;

pd->face[1].nedges = 5;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 3;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 20;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 8;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 13;
pd->face[1].or[3] = 1;

pd->face[1].edge[4] = 17;
pd->face[1].or[4] = 0;

pd->face[2].nedges = 4;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 1;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 21;
pd->face[2].or[1] = 0;

pd->face[2].edge[2] = 3;
pd->face[2].or[2] = 0;

pd->face[2].edge[3] = 16;
pd->face[2].or[3] = 0;

pd->face[3].nedges = 4;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 2;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 17;
pd->face[3].or[1] = 1;

pd->face[3].edge[2] = 10;
pd->face[3].or[2] = 1;

pd->face[3].edge[3] = 15;
pd->face[3].or[3] = 1;

pd->face[4].nedges = 4;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 6;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 12;
pd->face[4].or[1] = 1;

pd->face[4].edge[2] = 8;
pd->face[4].or[2] = 1;

pd->face[4].edge[3] = 19;
pd->face[4].or[3] = 0;

pd->face[5].nedges = 4;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 7;
pd->face[5].or[0] = 0;

pd->face[5].edge[1] = 11;
pd->face[5].or[1] = 0;

pd->face[5].edge[2] = 14;
pd->face[5].or[2] = 0;

pd->face[5].edge[3] = 13;
pd->face[5].or[3] = 0;

pd->face[6].nedges = 3;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 0;
pd->face[6].or[0] = 0;

pd->face[6].edge[1] = 4;
pd->face[6].or[1] = 0;

pd->face[6].edge[2] = 21;
pd->face[6].or[2] = 1;

pd->face[7].nedges = 3;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 0;
pd->face[7].or[0] = 1;

pd->face[7].edge[1] = 18;
pd->face[7].or[1] = 1;

pd->face[7].edge[2] = 5;
pd->face[7].or[2] = 0;

pd->face[8].nedges = 3;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].or = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].or != NULL);

pd->face[8].edge[0] = 4;
pd->face[8].or[0] = 1;

pd->face[8].edge[1] = 9;
pd->face[8].or[1] = 0;

pd->face[8].edge[2] = 20;
pd->face[8].or[2] = 1;

pd->face[9].nedges = 3;
pd->face[9].edge = calloc(pd->face[9].nedges,sizeof(pd_idx_t));
pd->face[9].or = calloc(pd->face[9].nedges,sizeof(pd_or_t));
assert(pd->face[9].edge != NULL);
assert(pd->face[9].or != NULL);

pd->face[9].edge[0] = 5;
pd->face[9].or[0] = 1;

pd->face[9].edge[1] = 19;
pd->face[9].or[1] = 1;

pd->face[9].edge[2] = 9;
pd->face[9].or[2] = 1;

pd->face[10].nedges = 2;
pd->face[10].edge = calloc(pd->face[10].nedges,sizeof(pd_idx_t));
pd->face[10].or = calloc(pd->face[10].nedges,sizeof(pd_or_t));
assert(pd->face[10].edge != NULL);
assert(pd->face[10].or != NULL);

pd->face[10].edge[0] = 2;
pd->face[10].or[0] = 0;

pd->face[10].edge[1] = 16;
pd->face[10].or[1] = 1;

pd->face[11].nedges = 2;
pd->face[11].edge = calloc(pd->face[11].nedges,sizeof(pd_idx_t));
pd->face[11].or = calloc(pd->face[11].nedges,sizeof(pd_or_t));
assert(pd->face[11].edge != NULL);
assert(pd->face[11].or != NULL);

pd->face[11].edge[0] = 7;
pd->face[11].or[0] = 1;

pd->face[11].edge[1] = 12;
pd->face[11].or[1] = 0;

pd->face[12].nedges = 2;
pd->face[12].edge = calloc(pd->face[12].nedges,sizeof(pd_idx_t));
pd->face[12].or = calloc(pd->face[12].nedges,sizeof(pd_or_t));
assert(pd->face[12].edge != NULL);
assert(pd->face[12].or != NULL);

pd->face[12].edge[0] = 10;
pd->face[12].or[0] = 0;

pd->face[12].edge[1] = 14;
pd->face[12].or[1] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



pd_code_t *pd_create_r2_testJ_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(9);
assert(pd != NULL);
pd->ncross = 9;
pd->nedges = 18;
pd->ncomps = 4;
pd->nfaces = 11;
sprintf(pd->hash,"%s","CRILCAQEAwMDAwICAgIECgQCAgAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 5;
pd->cross[0].edge[2] = 9;
pd->cross[0].edge[3] = 4;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 13;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 10;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 16;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 17;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 16;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 17;
pd->cross[3].sign = 1;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 13;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 12;
pd->cross[4].sign = 0;

pd->cross[5].edge[0] = 5;
pd->cross[5].edge[1] = 10;
pd->cross[5].edge[2] = 6;
pd->cross[5].edge[3] = 11;
pd->cross[5].sign = 1;

pd->cross[6].edge[0] = 6;
pd->cross[6].edge[1] = 15;
pd->cross[6].edge[2] = 7;
pd->cross[6].edge[3] = 14;
pd->cross[6].sign = 1;

pd->cross[7].edge[0] = 7;
pd->cross[7].edge[1] = 15;
pd->cross[7].edge[2] = 8;
pd->cross[7].edge[3] = 14;
pd->cross[7].sign = 0;

pd->cross[8].edge[0] = 8;
pd->cross[8].edge[1] = 12;
pd->cross[8].edge[2] = 9;
pd->cross[8].edge[3] = 11;
pd->cross[8].sign = 1;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 0;
pd->edge[4].headpos = 3;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 5;
pd->edge[5].headpos = 0;
pd->edge[5].tail = 0;
pd->edge[5].tailpos = 1;

pd->edge[6].head = 6;
pd->edge[6].headpos = 0;
pd->edge[6].tail = 5;
pd->edge[6].tailpos = 2;

pd->edge[7].head = 7;
pd->edge[7].headpos = 0;
pd->edge[7].tail = 6;
pd->edge[7].tailpos = 2;

pd->edge[8].head = 8;
pd->edge[8].headpos = 0;
pd->edge[8].tail = 7;
pd->edge[8].tailpos = 2;

pd->edge[9].head = 0;
pd->edge[9].headpos = 2;
pd->edge[9].tail = 8;
pd->edge[9].tailpos = 2;

pd->edge[10].head = 5;
pd->edge[10].headpos = 1;
pd->edge[10].tail = 1;
pd->edge[10].tailpos = 3;

pd->edge[11].head = 8;
pd->edge[11].headpos = 3;
pd->edge[11].tail = 5;
pd->edge[11].tailpos = 3;

pd->edge[12].head = 4;
pd->edge[12].headpos = 3;
pd->edge[12].tail = 8;
pd->edge[12].tailpos = 1;

pd->edge[13].head = 1;
pd->edge[13].headpos = 1;
pd->edge[13].tail = 4;
pd->edge[13].tailpos = 1;

pd->edge[14].head = 7;
pd->edge[14].headpos = 3;
pd->edge[14].tail = 6;
pd->edge[14].tailpos = 3;

pd->edge[15].head = 6;
pd->edge[15].headpos = 1;
pd->edge[15].tail = 7;
pd->edge[15].tailpos = 1;

pd->edge[16].head = 3;
pd->edge[16].headpos = 1;
pd->edge[16].tail = 2;
pd->edge[16].tailpos = 1;

pd->edge[17].head = 2;
pd->edge[17].headpos = 3;
pd->edge[17].tail = 3;
pd->edge[17].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 10;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;

pd->comp[1].nedges = 4;
pd->comp[1].tag = 'D';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 10;
pd->comp[1].edge[1] = 11;
pd->comp[1].edge[2] = 12;
pd->comp[1].edge[3] = 13;

pd->comp[2].nedges = 2;
pd->comp[2].tag = 'B';

pd->comp[2].edge = calloc(pd->comp[2].nedges,sizeof(pd_idx_t));
assert(pd->comp[2].edge != NULL);

pd->comp[2].edge[0] = 14;
pd->comp[2].edge[1] = 15;

pd->comp[3].nedges = 2;
pd->comp[3].tag = 'C';

pd->comp[3].edge = calloc(pd->comp[3].nedges,sizeof(pd_idx_t));
assert(pd->comp[3].edge != NULL);

pd->comp[3].edge[0] = 16;
pd->comp[3].edge[1] = 17;


/* Face data */

pd->face[0].nedges = 8;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 17;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 3;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 12;
pd->face[0].or[3] = 0;

pd->face[0].edge[4] = 8;
pd->face[0].or[4] = 0;

pd->face[0].edge[5] = 15;
pd->face[0].or[5] = 1;

pd->face[0].edge[6] = 6;
pd->face[0].or[6] = 0;

pd->face[0].edge[7] = 10;
pd->face[0].or[7] = 0;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 13;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 3;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 16;
pd->face[1].or[3] = 0;

pd->face[2].nedges = 4;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 6;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 14;
pd->face[2].or[1] = 1;

pd->face[2].edge[2] = 8;
pd->face[2].or[2] = 1;

pd->face[2].edge[3] = 11;
pd->face[2].or[3] = 0;

pd->face[3].nedges = 3;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 0;
pd->face[3].or[0] = 0;

pd->face[3].edge[1] = 4;
pd->face[3].or[1] = 0;

pd->face[3].edge[2] = 13;
pd->face[3].or[2] = 1;

pd->face[4].nedges = 3;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 0;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 10;
pd->face[4].or[1] = 1;

pd->face[4].edge[2] = 5;
pd->face[4].or[2] = 0;

pd->face[5].nedges = 3;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 4;
pd->face[5].or[0] = 1;

pd->face[5].edge[1] = 9;
pd->face[5].or[1] = 0;

pd->face[5].edge[2] = 12;
pd->face[5].or[2] = 1;

pd->face[6].nedges = 3;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 5;
pd->face[6].or[0] = 1;

pd->face[6].edge[1] = 11;
pd->face[6].or[1] = 1;

pd->face[6].edge[2] = 9;
pd->face[6].or[2] = 1;

pd->face[7].nedges = 2;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 2;
pd->face[7].or[0] = 0;

pd->face[7].edge[1] = 16;
pd->face[7].or[1] = 1;

pd->face[8].nedges = 2;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].or = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].or != NULL);

pd->face[8].edge[0] = 2;
pd->face[8].or[0] = 1;

pd->face[8].edge[1] = 17;
pd->face[8].or[1] = 1;

pd->face[9].nedges = 2;
pd->face[9].edge = calloc(pd->face[9].nedges,sizeof(pd_idx_t));
pd->face[9].or = calloc(pd->face[9].nedges,sizeof(pd_or_t));
assert(pd->face[9].edge != NULL);
assert(pd->face[9].or != NULL);

pd->face[9].edge[0] = 7;
pd->face[9].or[0] = 1;

pd->face[9].edge[1] = 14;
pd->face[9].or[1] = 0;

pd->face[10].nedges = 2;
pd->face[10].edge = calloc(pd->face[10].nedges,sizeof(pd_idx_t));
pd->face[10].or = calloc(pd->face[10].nedges,sizeof(pd_or_t));
assert(pd->face[10].edge != NULL);
assert(pd->face[10].or != NULL);

pd->face[10].edge[0] = 7;
pd->face[10].or[0] = 0;

pd->face[10].edge[1] = 15;
pd->face[10].or[1] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



/*                                                                            */
/*     ____13______>_______________    ________________________<_17_          */
/*    /                            \ /                              \         */
/*   /                              /10                              \        */
/*  /              (5)             / \               (3)              \       */
/* /                              /  10                                \      */
/* \                             14(12)\                               /      */
/*  \   _____<____7____            >   <        ____2___>____         /       */
/*   \ /        (11)   \            \ /        / (10)        \       /        */
/*    \                 \            \__15_>_________16_>_____\ ____/         */
/*   /7\____<____12_____ \  __<_11__/9       /2               3\              */
/*  /                    6\                 ^                   \             */
/*  \                      ^      (0)      1                    /             */
/*   \        (4)           6             /           (2)      /              */
/*    \             ________ \ ___<__18_______<__             /               */
/*     \          19         5\         /1       \           3     (1)        */
/*      8         /            ^  (7)  ^          \         /                 */
/*       \       /              \     0            \       /                  */
/*        \     /                5   /             21     /                   */
/*         >   <       (9)        \ /                \   <                    */
/*          \ /                    \        (6)       \ /                     */
/*           /                    /0\                  /                      */
/*          /8\                  /   \                /4\                     */
/*         /   \                /     ^              /   \                    */
/*        /     \____9___>_____/      \_________4___/     \                   */
/*       /                                                 \                  */
/*      /                       (8)                         \                 */
/*      \                                                   /                 */
/*       \                                                 /                  */
/*        \_____________________20_______>________________/                   */
/*                                                                            */
/* R2 on face 12 at crossings 9 and 10                                        */
/*     ________15__________>_              ________<________17______          */
/*    /                      \            /                         \         */
/*   /                        \          /                           \        */
/*  /                          \        /             (8)             \       */
/* /           (10)             \      /                               \      */
/* \                             \     \                               /      */
/*  \   _______<__7____          /      \       ______2__>___         /       */
/*   \ /        (9)    \        /        \     / (7)         \       /        */
/*    /                 \      /          \_________16__>_____\ ____/         */
/*   /7\_____<__14______ \  __/             2/                3\              */
/*  /                    6\                 ^                   \             */
/*  \                      ^               1          (1)       /             */
/*   \        (2)           6             /                    /              */
/*    \             ________ \ ___<_10_______<___             3      (0)      */
/*     8           /         5\         /1       \           /                */
/*      \        11            ^  (4)  ^          \         /                 */
/*       \       /              \     /           13       /                  */
/*        >     /                5   0              \     <                   */
/*         \   <         (6)      \ /        (3)     \   /                    */
/*          \ /                    \                  \ /                     */
/*           /                    /0\                  /                      */
/*          /8\                  /   \                /4\                     */
/*         /   \                /     \              /   \                    */
/*        /     \______9__>____/       \____<__4____/     \                   */
/*       /                                                 \                  */
/*      /                                                   \                 */
/*      \                          (5)                      /                 */
/*       \                                                 /                  */
/*        \______________________12____>__________________/                   */



    bool r2testJ() {
  
  printf("---------------------------------------\n"
	 "r2 test J\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r2_testJ_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after-0 configuration...");
  pd_code_t *after0 = pd_create_r2_testJ_after_0();
  printf("done (passes pd_ok)\n");

  pd_code_t *after[1] = {after0};
    
  printf("performing r2 at crossings 9 and 10...\n");

  pd_idx_t  noutpd;
  pd_code_t **newpd;
  pd_idx_t  cr[2] = {9,10};

  pd_R2_bigon_elimination(before,cr,&noutpd,&newpd);

  if (noutpd != 1) { 
    printf("fail (did not generate 1 output pd codes)\n");
    return false;
  }

  pd_idx_t ct;

  for(ct=0;ct<noutpd;ct++) { 

    if (!pd_ok(newpd[ct])) { 
      printf("fail (output pd %d does not pass pd_ok)\n",ct);
      return false;
    }
    printf("\t output pd %d passes pd_ok...pass\n",ct);
    printf("\t testing for isomorphism with after configuration %d...",ct);
    if (!pd_isomorphic(newpd[ct],after[ct])) { 
      pd_printf("fail (output pd %PD ",newpd[ct]);
      pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after[ct]);
      return false;
    }

    printf("pass (output %d and after %d isomorphic)\n",ct,ct);

  }

  printf("performing r2 at crossings 9 and 10...pass\n");
  printf("housecleaning...");

  pd_code_free(&before);
  pd_code_free(&after0);
  pd_code_free(&newpd[0]);

  free(newpd);

  printf("done\n");

  printf("---------------------------------------\n"
	 "r2 test J : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}

/* r2 test K assembled c code */

pd_code_t *pd_create_r2_testK_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(6);
assert(pd != NULL);
pd->ncross = 6;
pd->nedges = 12;
pd->ncomps = 2;
pd->nfaces = 8;
sprintf(pd->hash,"%s","BgwIBAQEAwMCAgICCAQAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 8;
pd->cross[0].edge[2] = 7;
pd->cross[0].edge[3] = 11;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 9;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 8;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 6;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 7;
pd->cross[2].sign = 0;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 6;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 5;
pd->cross[3].sign = 0;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 9;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 10;
pd->cross[4].sign = 1;

pd->cross[5].edge[0] = 4;
pd->cross[5].edge[1] = 11;
pd->cross[5].edge[2] = 5;
pd->cross[5].edge[3] = 10;
pd->cross[5].sign = 1;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 5;
pd->edge[4].headpos = 0;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 3;
pd->edge[5].headpos = 3;
pd->edge[5].tail = 5;
pd->edge[5].tailpos = 2;

pd->edge[6].head = 2;
pd->edge[6].headpos = 1;
pd->edge[6].tail = 3;
pd->edge[6].tailpos = 1;

pd->edge[7].head = 0;
pd->edge[7].headpos = 2;
pd->edge[7].tail = 2;
pd->edge[7].tailpos = 3;

pd->edge[8].head = 1;
pd->edge[8].headpos = 3;
pd->edge[8].tail = 0;
pd->edge[8].tailpos = 1;

pd->edge[9].head = 4;
pd->edge[9].headpos = 1;
pd->edge[9].tail = 1;
pd->edge[9].tailpos = 1;

pd->edge[10].head = 5;
pd->edge[10].headpos = 3;
pd->edge[10].tail = 4;
pd->edge[10].tailpos = 3;

pd->edge[11].head = 0;
pd->edge[11].headpos = 3;
pd->edge[11].tail = 5;
pd->edge[11].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 8;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;

pd->comp[1].nedges = 4;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 8;
pd->comp[1].edge[1] = 9;
pd->comp[1].edge[2] = 10;
pd->comp[1].edge[3] = 11;


/* Face data */

pd->face[0].nedges = 4;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 11;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 4;
pd->face[0].or[2] = 0;

pd->face[0].edge[3] = 9;
pd->face[0].or[3] = 0;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 9;
pd->face[1].or[1] = 1;

pd->face[1].edge[2] = 3;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 6;
pd->face[1].or[3] = 1;

pd->face[2].nedges = 4;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 2;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 5;
pd->face[2].or[1] = 0;

pd->face[2].edge[2] = 11;
pd->face[2].or[2] = 1;

pd->face[2].edge[3] = 7;
pd->face[2].or[3] = 0;

pd->face[3].nedges = 3;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 1;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 7;
pd->face[3].or[1] = 1;

pd->face[3].edge[2] = 8;
pd->face[3].or[2] = 1;

pd->face[4].nedges = 3;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 3;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 10;
pd->face[4].or[1] = 1;

pd->face[4].edge[2] = 5;
pd->face[4].or[2] = 1;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 0;
pd->face[5].or[0] = 1;

pd->face[5].edge[1] = 8;
pd->face[5].or[1] = 0;

pd->face[6].nedges = 2;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 2;
pd->face[6].or[0] = 0;

pd->face[6].edge[1] = 6;
pd->face[6].or[1] = 0;

pd->face[7].nedges = 2;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 4;
pd->face[7].or[0] = 1;

pd->face[7].edge[1] = 10;
pd->face[7].or[1] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



pd_code_t *pd_create_r2_testK_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(4);
assert(pd != NULL);
pd->ncross = 4;
pd->nedges = 8;
pd->ncomps = 2;
pd->nfaces = 6;
sprintf(pd->hash,"%s","BAgGBgMCAgIBAgYCAAAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 0;
pd->cross[0].edge[2] = 5;
pd->cross[0].edge[3] = 1;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 1;
pd->cross[1].edge[1] = 5;
pd->cross[1].edge[2] = 2;
pd->cross[1].edge[3] = 4;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 2;
pd->cross[2].edge[1] = 7;
pd->cross[2].edge[2] = 3;
pd->cross[2].edge[3] = 6;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 3;
pd->cross[3].edge[1] = 7;
pd->cross[3].edge[2] = 4;
pd->cross[3].edge[3] = 6;
pd->cross[3].sign = 1;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 1;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 1;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 3;

pd->edge[2].head = 2;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 1;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 3;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 2;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 1;
pd->edge[4].headpos = 3;
pd->edge[4].tail = 3;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 0;
pd->edge[5].headpos = 2;
pd->edge[5].tail = 1;
pd->edge[5].tailpos = 1;

pd->edge[6].head = 3;
pd->edge[6].headpos = 3;
pd->edge[6].tail = 2;
pd->edge[6].tailpos = 3;

pd->edge[7].head = 2;
pd->edge[7].headpos = 1;
pd->edge[7].tail = 3;
pd->edge[7].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 6;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;

pd->comp[1].nedges = 2;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 6;
pd->comp[1].edge[1] = 7;


/* Face data */

pd->face[0].nedges = 6;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 1;
pd->face[0].or[1] = 1;

pd->face[0].edge[2] = 4;
pd->face[0].or[2] = 0;

pd->face[0].edge[3] = 7;
pd->face[0].or[3] = 1;

pd->face[0].edge[4] = 2;
pd->face[0].or[4] = 0;

pd->face[0].edge[5] = 5;
pd->face[0].or[5] = 1;

pd->face[1].nedges = 3;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 2;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 6;
pd->face[1].or[1] = 1;

pd->face[1].edge[2] = 4;
pd->face[1].or[2] = 1;

pd->face[2].nedges = 2;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 1;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 5;
pd->face[2].or[1] = 0;

pd->face[3].nedges = 2;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 3;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 6;
pd->face[3].or[1] = 0;

pd->face[4].nedges = 2;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 3;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 7;
pd->face[4].or[1] = 0;

pd->face[5].nedges = 1;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 0;
pd->face[5].or[0] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



/*                                                                            */
/*                      _________1___________                                 */
/*                     /                     \              (1)               */
/*                    /                       \                               */
/*                   /        (3)              \                              */
/*                  <                           \                             */
/*                 /               __8_____>___________                       */
/*                /               /     (5)    /1      \                      */
/*                \               \           ^         \                     */
/*                 \     7_____>__ \ ____0___/           \                    */
/*                  \   /          0\                     \                   */
/*                   \2/             \                     \                  */
/*                    /               ^          (0)        9                 */
/*                   / \               \                    /                 */
/*                  /   \              11                  <                  */
/*                 ^     \               \                /                   */
/*                /       \      (2)      \   _<____4______                   */
/*               /         2               \ /    (7)   /4 \                  */
/*              /           \               \          /    \                 */
/*              |            \             /5\_<__10__/      ^                */
/*              |             \           /                   \               */
/*              |              \         /                    /               */
/*              |       (6)     \       5                    /                */
/*              6                \     /                    /                 */
/*              |                 >   <         (4)        /                  */
/*              \                  \ /                    /                   */
/*               \                  \                    /                    */
/*                \                /3\__________________3                     */
/*                 \______________/                                           */
/* R2 on face 5  at crossings 0 and 1                                         */                                             
/*                      _____________________                                 */
/*                     /                     \                                */
/*                    /                       \                               */
/*                   <                         \                              */
/*                  /                           0                             */
/*                 /            (5)             /  __>_                       */
/*                /                            /  /    \                      */
/*                \                           /  /      \                     */
/*                 \     ____________________/  /        \                    */
/*                  \   /                      /          \                   */
/*                   \0/                      /            \                  */
/*                    /                      7    (4)       \                 */
/*                   / \        (0)         /               /                 */
/*                  /   \                  /               /                  */
/*                 ^     \                /               /                   */
/*                /       \               \   _<____3______                   */
/*               /         1               \ /   (3)    /2 \                  */
/*              /           \              3\          /    \                 */
/*              |            \             / \__<__6__/      ^                */
/*              |             \           /                   \               */
/*              |              \         4                    /               */
/*              |    (2)        \       /                    /                */
/*              5                >     /        (1)         /                 */
/*              |                 \   <                    2                  */
/*              \                  \ /                    /                   */
/*               \                  \                    /                    */
/*                \                /1\__________________/                     */
/*                 \______________/                                           */
/*                                                                            */



    bool r2testK() {
  
  printf("---------------------------------------\n"
	 "r2 test K\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r2_testK_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after-0 configuration...");
  pd_code_t *after0 = pd_create_r2_testK_after_0();
  printf("done (passes pd_ok)\n");

  pd_code_t *after[1] = {after0};
    
  printf("performing r2 at crossings 0 and 1...\n");

  pd_idx_t  noutpd;
  pd_code_t **newpd;
  pd_idx_t  cr[2] = {0,1};

  pd_R2_bigon_elimination(before,cr,&noutpd,&newpd);

  if (noutpd != 1) { 
    printf("fail (did not generate 1 output pd codes)\n");
    return false;
  }

  pd_idx_t ct;

  for(ct=0;ct<noutpd;ct++) { 

    if (!pd_ok(newpd[ct])) { 
      printf("fail (output pd %d does not pass pd_ok)\n",ct);
      return false;
    }
    printf("\t output pd %d passes pd_ok...pass\n",ct);
    printf("\t testing for isomorphism with after configuration %d...",ct);
    if (!pd_isomorphic(newpd[ct],after[ct])) { 
      pd_printf("fail (output pd %PD ",newpd[ct]);
      pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after[ct]);
      return false;
    }

    printf("pass (output %d and after %d isomorphic)\n",ct,ct);

  }

  printf("performing r2 at crossings 0 and 1...pass\n");
  printf("housecleaning...");

  pd_code_free(&before);
  pd_code_free(&after0);
  pd_code_free(&newpd[0]);

  free(newpd);

  printf("done\n");

  printf("---------------------------------------\n"
	 "r2 test K : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}

/* r2 test L assembled c code */

pd_code_t *pd_create_r2_testL_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(7);
assert(pd != NULL);
pd->ncross = 7;
pd->nedges = 14;
pd->ncomps = 2;
pd->nfaces = 9;
sprintf(pd->hash,"%s","Bw4JBQUEAwMCAgICAgoEAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 10;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 11;
pd->cross[0].sign = 1;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 13;
pd->cross[1].edge[2] = 9;
pd->cross[1].edge[3] = 10;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 7;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 6;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 7;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 8;
pd->cross[3].sign = 1;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 9;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 8;
pd->cross[4].sign = 1;

pd->cross[5].edge[0] = 4;
pd->cross[5].edge[1] = 13;
pd->cross[5].edge[2] = 5;
pd->cross[5].edge[3] = 12;
pd->cross[5].sign = 1;

pd->cross[6].edge[0] = 5;
pd->cross[6].edge[1] = 11;
pd->cross[6].edge[2] = 6;
pd->cross[6].edge[3] = 12;
pd->cross[6].sign = 1;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 1;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 5;
pd->edge[4].headpos = 0;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 6;
pd->edge[5].headpos = 0;
pd->edge[5].tail = 5;
pd->edge[5].tailpos = 2;

pd->edge[6].head = 2;
pd->edge[6].headpos = 3;
pd->edge[6].tail = 6;
pd->edge[6].tailpos = 2;

pd->edge[7].head = 3;
pd->edge[7].headpos = 1;
pd->edge[7].tail = 2;
pd->edge[7].tailpos = 1;

pd->edge[8].head = 4;
pd->edge[8].headpos = 3;
pd->edge[8].tail = 3;
pd->edge[8].tailpos = 3;

pd->edge[9].head = 1;
pd->edge[9].headpos = 2;
pd->edge[9].tail = 4;
pd->edge[9].tailpos = 1;

pd->edge[10].head = 0;
pd->edge[10].headpos = 1;
pd->edge[10].tail = 1;
pd->edge[10].tailpos = 3;

pd->edge[11].head = 6;
pd->edge[11].headpos = 1;
pd->edge[11].tail = 0;
pd->edge[11].tailpos = 3;

pd->edge[12].head = 5;
pd->edge[12].headpos = 3;
pd->edge[12].tail = 6;
pd->edge[12].tailpos = 3;

pd->edge[13].head = 1;
pd->edge[13].headpos = 1;
pd->edge[13].tail = 5;
pd->edge[13].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 10;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;

pd->comp[1].nedges = 4;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 10;
pd->comp[1].edge[1] = 11;
pd->comp[1].edge[2] = 12;
pd->comp[1].edge[3] = 13;


/* Face data */

pd->face[0].nedges = 5;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 10;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 9;
pd->face[0].or[2] = 0;

pd->face[0].edge[3] = 3;
pd->face[0].or[3] = 0;

pd->face[0].edge[4] = 7;
pd->face[0].or[4] = 0;

pd->face[1].nedges = 5;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 2;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 8;
pd->face[1].or[1] = 1;

pd->face[1].edge[2] = 4;
pd->face[1].or[2] = 1;

pd->face[1].edge[3] = 12;
pd->face[1].or[3] = 0;

pd->face[1].edge[4] = 6;
pd->face[1].or[4] = 1;

pd->face[2].nedges = 4;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 11;
pd->face[2].or[1] = 1;

pd->face[2].edge[2] = 5;
pd->face[2].or[2] = 0;

pd->face[2].edge[3] = 13;
pd->face[2].or[3] = 1;

pd->face[3].nedges = 3;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 1;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 6;
pd->face[3].or[1] = 0;

pd->face[3].edge[2] = 11;
pd->face[3].or[2] = 0;

pd->face[4].nedges = 3;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 4;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 9;
pd->face[4].or[1] = 1;

pd->face[4].edge[2] = 13;
pd->face[4].or[2] = 0;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 0;
pd->face[5].or[0] = 0;

pd->face[5].edge[1] = 10;
pd->face[5].or[1] = 1;

pd->face[6].nedges = 2;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 2;
pd->face[6].or[0] = 0;

pd->face[6].edge[1] = 7;
pd->face[6].or[1] = 1;

pd->face[7].nedges = 2;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 3;
pd->face[7].or[0] = 1;

pd->face[7].edge[1] = 8;
pd->face[7].or[1] = 0;

pd->face[8].nedges = 2;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].or = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].or != NULL);

pd->face[8].edge[0] = 5;
pd->face[8].or[0] = 1;

pd->face[8].edge[1] = 12;
pd->face[8].or[1] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



pd_code_t *pd_create_r2_testL_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(5);
assert(pd != NULL);
pd->ncross = 5;
pd->nedges = 10;
pd->ncomps = 2;
pd->nfaces = 7;
sprintf(pd->hash,"%s","BQoHBQQDAgICAgIIAgAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 2;
pd->cross[0].edge[2] = 7;
pd->cross[0].edge[3] = 3;
pd->cross[0].sign = 1;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 3;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 4;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 7;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 6;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 4;
pd->cross[3].edge[1] = 8;
pd->cross[3].edge[2] = 5;
pd->cross[3].edge[3] = 9;
pd->cross[3].sign = 1;

pd->cross[4].edge[0] = 5;
pd->cross[4].edge[1] = 8;
pd->cross[4].edge[2] = 6;
pd->cross[4].edge[3] = 9;
pd->cross[4].sign = 1;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 0;
pd->edge[2].headpos = 1;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 1;
pd->edge[3].headpos = 1;
pd->edge[3].tail = 0;
pd->edge[3].tailpos = 3;

pd->edge[4].head = 3;
pd->edge[4].headpos = 0;
pd->edge[4].tail = 1;
pd->edge[4].tailpos = 3;

pd->edge[5].head = 4;
pd->edge[5].headpos = 0;
pd->edge[5].tail = 3;
pd->edge[5].tailpos = 2;

pd->edge[6].head = 2;
pd->edge[6].headpos = 3;
pd->edge[6].tail = 4;
pd->edge[6].tailpos = 2;

pd->edge[7].head = 0;
pd->edge[7].headpos = 2;
pd->edge[7].tail = 2;
pd->edge[7].tailpos = 1;

pd->edge[8].head = 4;
pd->edge[8].headpos = 1;
pd->edge[8].tail = 3;
pd->edge[8].tailpos = 1;

pd->edge[9].head = 3;
pd->edge[9].headpos = 3;
pd->edge[9].tail = 4;
pd->edge[9].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 8;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;

pd->comp[1].nedges = 2;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 8;
pd->comp[1].edge[1] = 9;


/* Face data */

pd->face[0].nedges = 5;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 4;
pd->face[0].or[1] = 1;

pd->face[0].edge[2] = 9;
pd->face[0].or[2] = 0;

pd->face[0].edge[3] = 6;
pd->face[0].or[3] = 1;

pd->face[0].edge[4] = 2;
pd->face[0].or[4] = 1;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 6;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 8;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 4;
pd->face[1].or[3] = 0;

pd->face[2].nedges = 3;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 1;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 3;
pd->face[2].or[1] = 0;

pd->face[2].edge[2] = 7;
pd->face[2].or[2] = 0;

pd->face[3].nedges = 2;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 0;
pd->face[3].or[0] = 0;

pd->face[3].edge[1] = 3;
pd->face[3].or[1] = 1;

pd->face[4].nedges = 2;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 2;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 7;
pd->face[4].or[1] = 1;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 5;
pd->face[5].or[0] = 0;

pd->face[5].edge[1] = 8;
pd->face[5].or[1] = 1;

pd->face[6].nedges = 2;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 5;
pd->face[6].or[0] = 1;

pd->face[6].edge[1] = 9;
pd->face[6].or[1] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



/*                                                                            */
/*                              _____8___________                             */
/*                 ______>     /                 \                            */
/*                /       \   /                   >         (1)               */
/*               2         \ /         (7)         \                          */
/*              /    (6)    \                       \                         */
/*             /           /3\                      |                         */
/*     ____________7___>__/   \________3_____>_____ |____                     */
/*    /     2/                                     4|    \                    */
/*   /      /                                       |     \                   */
/*  /       \               (0)                    /       \                  */
/*  ^        ^                                    9         \                 */
/*  |         \             __<___10___          /           \                */
/*  |          \           /           \        /             \               */
/*  |           1         /    (5)      \      <               \              */
/*  |            \       /               \    /                 \             */
/*  |             \     /                 \  /                   \            */
/*  6              \___________<____0_______/                     \           */
/*   \                 |0                 1|         (4)          /           */
/*    \      (3)       |                   ^                     /            */
/*     \               |                   |                    4             */
/*      \             11        (2)        |                   /              */
/*       \             |                  13                  /               */
/*        \            |                   /                 /                */
/*         \           |                  /                 /                 */
/*          \          |                 /                 /                  */
/*           \______________<_____5_____/_____<___________/                   */
/*                    6\               /5                                     */
/*                      \      (8)    /                                       */
/*                       \___12____>_/                                        */
/* R2 on face 5 at crossings 0 and 1                                          */
/*                              _____0___________                             */
/*                 _______     /                 \                            */
/*                /       >   /                   \                           */
/*               /         \ /         (3)         >                          */
/*              2   (4)     \                       \         (0)             */
/*             /           /0\                      |                         */
/*     _____________7__>__/   \______3_________>___ |____                     */
/*    /     2/                                     1|    \                    */
/*   ^      /                                       |     \                   */
/*  /       \                                      /       \                  */
/*  |        \                                    /         \                 */
/*  |         ^                                  /           \                */
/*  |          \              (2)               /             \               */
/*  |           \                              /               \              */
/*  |            \                            /                 \             */
/*  |             \                          /                   \            */
/*  6              \______1_________________/                     \           */
/*   \                                                            4           */
/*    \                                      (1)                 /            */
/*     \                  ___<____8__                           /             */
/*      \                /           \                         /              */
/*       \              /             \                       /               */
/*        \            /               \                     /                */
/*         \           |      (5)      |                    /                 */
/*          \          |               |                   /                  */
/*           \____________<_______5___ | _____<___________/                   */
/*                    4\      (6)      /3                                     */
/*                      \             /                                       */
/*                       \_9_____>___/                                        */



    bool r2testL() {
  
  printf("---------------------------------------\n"
	 "r2 test L\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r2_testL_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after-0 configuration...");
  pd_code_t *after0 = pd_create_r2_testL_after_0();
  printf("done (passes pd_ok)\n");

  pd_code_t *after[1] = {after0};
    
  printf("performing r2 at crossings 0 and 1...\n");

  pd_idx_t  noutpd;
  pd_code_t **newpd;
  pd_idx_t  cr[2] = {0,1};

  pd_R2_bigon_elimination(before,cr,&noutpd,&newpd);

  if (noutpd != 1) { 
    printf("fail (did not generate 1 output pd codes)\n");
    return false;
  }

  pd_idx_t ct;

  for(ct=0;ct<noutpd;ct++) { 

    if (!pd_ok(newpd[ct])) { 
      printf("fail (output pd %d does not pass pd_ok)\n",ct);
      return false;
    }
    printf("\t output pd %d passes pd_ok...pass\n",ct);
    printf("\t testing for isomorphism with after configuration %d...",ct);
    if (!pd_isomorphic(newpd[ct],after[ct])) { 
      pd_printf("fail (output pd %PD ",newpd[ct]);
      pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after[ct]);
      return false;
    }

    printf("pass (output %d and after %d isomorphic)\n",ct,ct);

  }

  printf("performing r2 at crossings 0 and 1...pass\n");
  printf("housecleaning...");

  pd_code_free(&before);
  pd_code_free(&after0);
  pd_code_free(&newpd[0]);

  free(newpd);

  printf("done\n");

  printf("---------------------------------------\n"
	 "r2 test L : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}

/* r2 test N assembled c code */

pd_code_t *pd_create_r2_testN_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(7);
assert(pd != NULL);
pd->ncross = 7;
pd->nedges = 14;
pd->ncomps = 1;
pd->nfaces = 9;
sprintf(pd->hash,"%s","Bw4JBgQEAwMCAgICAQ4AAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 7;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 8;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 8;
pd->cross[1].edge[2] = 13;
pd->cross[1].edge[3] = 9;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 7;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 6;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 9;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 10;
pd->cross[3].sign = 1;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 13;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 12;
pd->cross[4].sign = 1;

pd->cross[5].edge[0] = 4;
pd->cross[5].edge[1] = 11;
pd->cross[5].edge[2] = 5;
pd->cross[5].edge[3] = 12;
pd->cross[5].sign = 1;

pd->cross[6].edge[0] = 5;
pd->cross[6].edge[1] = 11;
pd->cross[6].edge[2] = 6;
pd->cross[6].edge[3] = 10;
pd->cross[6].sign = 1;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 1;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 5;
pd->edge[4].headpos = 0;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 6;
pd->edge[5].headpos = 0;
pd->edge[5].tail = 5;
pd->edge[5].tailpos = 2;

pd->edge[6].head = 2;
pd->edge[6].headpos = 3;
pd->edge[6].tail = 6;
pd->edge[6].tailpos = 2;

pd->edge[7].head = 0;
pd->edge[7].headpos = 1;
pd->edge[7].tail = 2;
pd->edge[7].tailpos = 1;

pd->edge[8].head = 1;
pd->edge[8].headpos = 1;
pd->edge[8].tail = 0;
pd->edge[8].tailpos = 3;

pd->edge[9].head = 3;
pd->edge[9].headpos = 1;
pd->edge[9].tail = 1;
pd->edge[9].tailpos = 3;

pd->edge[10].head = 6;
pd->edge[10].headpos = 3;
pd->edge[10].tail = 3;
pd->edge[10].tailpos = 3;

pd->edge[11].head = 5;
pd->edge[11].headpos = 1;
pd->edge[11].tail = 6;
pd->edge[11].tailpos = 1;

pd->edge[12].head = 4;
pd->edge[12].headpos = 3;
pd->edge[12].tail = 5;
pd->edge[12].tailpos = 3;

pd->edge[13].head = 1;
pd->edge[13].headpos = 2;
pd->edge[13].tail = 4;
pd->edge[13].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 14;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;
pd->comp[0].edge[10] = 10;
pd->comp[0].edge[11] = 11;
pd->comp[0].edge[12] = 12;
pd->comp[0].edge[13] = 13;


/* Face data */

pd->face[0].nedges = 6;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 6;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 11;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 4;
pd->face[0].or[3] = 0;

pd->face[0].edge[4] = 13;
pd->face[0].or[4] = 1;

pd->face[0].edge[5] = 8;
pd->face[0].or[5] = 0;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 0;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 9;
pd->face[1].or[1] = 1;

pd->face[1].edge[2] = 2;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 7;
pd->face[1].or[3] = 1;

pd->face[2].nedges = 4;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 3;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 12;
pd->face[2].or[1] = 0;

pd->face[2].edge[2] = 5;
pd->face[2].or[2] = 1;

pd->face[2].edge[3] = 10;
pd->face[2].or[3] = 0;

pd->face[3].nedges = 3;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 2;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 10;
pd->face[3].or[1] = 1;

pd->face[3].edge[2] = 6;
pd->face[3].or[2] = 1;

pd->face[4].nedges = 3;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 3;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 9;
pd->face[4].or[1] = 0;

pd->face[4].edge[2] = 13;
pd->face[4].or[2] = 0;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 0;
pd->face[5].or[0] = 1;

pd->face[5].edge[1] = 8;
pd->face[5].or[1] = 1;

pd->face[6].nedges = 2;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 1;
pd->face[6].or[0] = 0;

pd->face[6].edge[1] = 7;
pd->face[6].or[1] = 0;

pd->face[7].nedges = 2;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 4;
pd->face[7].or[0] = 1;

pd->face[7].edge[1] = 12;
pd->face[7].or[1] = 1;

pd->face[8].nedges = 2;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].or = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].or != NULL);

pd->face[8].edge[0] = 5;
pd->face[8].or[0] = 0;

pd->face[8].edge[1] = 11;
pd->face[8].or[1] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



pd_code_t *pd_create_r2_testN_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(5);
assert(pd != NULL);
pd->ncross = 5;
pd->nedges = 10;
pd->ncomps = 1;
pd->nfaces = 7;
sprintf(pd->hash,"%s","BQoHBAQDAwICAgEKAAAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 5;
pd->cross[0].edge[2] = 9;
pd->cross[0].edge[3] = 4;
pd->cross[0].sign = 1;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 6;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 5;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 8;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 9;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 8;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 7;
pd->cross[3].sign = 1;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 6;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 7;
pd->cross[4].sign = 1;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 0;
pd->edge[4].headpos = 3;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 1;
pd->edge[5].headpos = 3;
pd->edge[5].tail = 0;
pd->edge[5].tailpos = 1;

pd->edge[6].head = 4;
pd->edge[6].headpos = 1;
pd->edge[6].tail = 1;
pd->edge[6].tailpos = 1;

pd->edge[7].head = 3;
pd->edge[7].headpos = 3;
pd->edge[7].tail = 4;
pd->edge[7].tailpos = 3;

pd->edge[8].head = 2;
pd->edge[8].headpos = 1;
pd->edge[8].tail = 3;
pd->edge[8].tailpos = 1;

pd->edge[9].head = 0;
pd->edge[9].headpos = 2;
pd->edge[9].tail = 2;
pd->edge[9].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 10;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;


/* Face data */

pd->face[0].nedges = 4;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 6;
pd->face[0].or[1] = 1;

pd->face[0].edge[2] = 3;
pd->face[0].or[2] = 0;

pd->face[0].edge[3] = 8;
pd->face[0].or[3] = 1;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 2;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 7;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 4;
pd->face[1].or[2] = 1;

pd->face[1].edge[3] = 9;
pd->face[1].or[3] = 0;

pd->face[2].nedges = 3;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 4;
pd->face[2].or[1] = 0;

pd->face[2].edge[2] = 6;
pd->face[2].or[2] = 0;

pd->face[3].nedges = 3;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 1;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 9;
pd->face[3].or[1] = 1;

pd->face[3].edge[2] = 5;
pd->face[3].or[2] = 1;

pd->face[4].nedges = 2;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 0;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 5;
pd->face[4].or[1] = 0;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 2;
pd->face[5].or[0] = 0;

pd->face[5].edge[1] = 8;
pd->face[5].or[1] = 0;

pd->face[6].nedges = 2;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 3;
pd->face[6].or[0] = 1;

pd->face[6].edge[1] = 7;
pd->face[6].or[1] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



/*                  ________   ________                                       */
/*                 /        \0/        1                                      */
/*                /   (5)    \   (6)    \     (0)                             */
/*                8         / \         |                                     */
/*                |        ^   \        |                                     */
/*                 \      /     ^       |                                     */
/*                  >    0       \     /                                      */
/*                   \  /         7   <                                       */
/*                  1 \            \ /                                        */
/*                     \    (1)     \2                                        */
/*                   /  \          / \                                        */
/*                  /    9        2   \                                       */
/*                 ^      \      /     ^                                      */
/*                /        >    /       \                                     */
/*               |    (4)   \  <   (3)   |                                    */
/*               13          _/          6                                    */
/*                \         /3          /                                     */
/*                 \       /   \       /                                      */
/*                  \     3    10     /                                       */
/*                   \   <       >   /                                        */
/*                    \ /         \ /                                         */
/*                    4\    (2)    \6                                         */
/*                    / \         / \                                         */
/*                   /   ^       ^   \                                        */
/*                  /     \     /     \                                       */
/*                 /      12   5       \                                      */
/*                /   (7)   \ /  (8)    \                                     */
/*               /           /          11                                    */
/*               \_4____>___/5\___<_____/                                     */
/*  R2 on face 5 at crossings 0 and 1                                         */
/*                     ________________                                       */
/*                    /                \                                      */
/*                   /  ____>_          \              (0)                    */
/*                  /  /      \         |                                     */
/*                 /  /        \        |                                     */
/*                /  /          \       |                                     */
/*               /   |           \     /                                      */
/*              /    5     (4)    \   /                                       */
/*             /      \            \ /                                        */
/*            /        \          1 \                                         */
/*           /          \          / \                                        */
/*           1           \        ^   \                                       */
/*            \   (3)     \      /     \                                      */
/*             \           \    0       \                                     */
/*              \           \  /         6                                    */
/*               \         0 _/   (2)    |                                    */
/*                \         /           /                                     */
/*                 \       ^   \       /                                      */
/*                  >     /     ^     <                                       */
/*                   \   9       4   /                                        */
/*                    \ /         \ /                                         */
/*                   2 \     (1)   \4                                         */
/*                    / \         / \                                         */
/*                   /   2       7   \                                        */
/*                  ^     \     /     ^                                       */
/*                 /       >   <       \                                      */
/*                8   (5)   \ /   (6)   \                                     */
/*               /           /          3                                     */
/*               \__________/3\_________/                                     */
/*                                                                            */

bool r2testN() {
  
  printf("---------------------------------------\n"
	 "r2 test N\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r2_testN_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after-0 configuration...");
  pd_code_t *after0 = pd_create_r2_testN_after_0();
  printf("done (passes pd_ok)\n");

  pd_code_t *after[1] = {after0};
    
  printf("performing r2 at crossings 0 and 1...\n");

  pd_idx_t  noutpd;
  pd_code_t **newpd;
  pd_idx_t  cr[2] = {0,1};

  pd_R2_bigon_elimination(before,cr,&noutpd,&newpd);

  if (noutpd != 1) { 
    printf("fail (did not generate 1 output pd codes)\n");
    return false;
  }

  pd_idx_t ct;

  for(ct=0;ct<noutpd;ct++) { 

    if (!pd_ok(newpd[ct])) { 
      printf("fail (output pd %d does not pass pd_ok)\n",ct);
      return false;
    }
    printf("\t output pd %d passes pd_ok...pass\n",ct);
    printf("\t testing for isomorphism with after configuration %d...",ct);
    if (!pd_isomorphic(newpd[ct],after[ct])) { 
      pd_printf("fail (output pd %PD ",newpd[ct]);
      pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after[ct]);
      return false;
    }

    printf("pass (output %d and after %d isomorphic)\n",ct,ct);

  }

  printf("performing r2 at crossings 0 and 1...pass\n");
  printf("housecleaning...");

  pd_code_free(&before);
  pd_code_free(&after0);
  pd_code_free(&newpd[0]);

  free(newpd);

  printf("done\n");

  printf("---------------------------------------\n"
	 "r2 test N : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}

/* r2 test O assembled c code */

pd_code_t *pd_create_r2_testO_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(9);
assert(pd != NULL);
pd->ncross = 9;
pd->nedges = 18;
pd->ncomps = 2;
pd->nfaces = 11;
sprintf(pd->hash,"%s","CRILBQQEAwMDAwMDAwICDAYAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 4;
pd->cross[0].edge[2] = 11;
pd->cross[0].edge[3] = 5;
pd->cross[0].sign = 1;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 13;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 14;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 6;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 7;
pd->cross[2].sign = 0;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 10;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 9;
pd->cross[3].sign = 0;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 17;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 16;
pd->cross[4].sign = 1;

pd->cross[5].edge[0] = 5;
pd->cross[5].edge[1] = 12;
pd->cross[5].edge[2] = 6;
pd->cross[5].edge[3] = 13;
pd->cross[5].sign = 1;

pd->cross[6].edge[0] = 7;
pd->cross[6].edge[1] = 15;
pd->cross[6].edge[2] = 8;
pd->cross[6].edge[3] = 14;
pd->cross[6].sign = 1;

pd->cross[7].edge[0] = 8;
pd->cross[7].edge[1] = 15;
pd->cross[7].edge[2] = 9;
pd->cross[7].edge[3] = 16;
pd->cross[7].sign = 0;

pd->cross[8].edge[0] = 10;
pd->cross[8].edge[1] = 12;
pd->cross[8].edge[2] = 11;
pd->cross[8].edge[3] = 17;
pd->cross[8].sign = 1;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 0;
pd->edge[4].headpos = 1;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 5;
pd->edge[5].headpos = 0;
pd->edge[5].tail = 0;
pd->edge[5].tailpos = 3;

pd->edge[6].head = 2;
pd->edge[6].headpos = 1;
pd->edge[6].tail = 5;
pd->edge[6].tailpos = 2;

pd->edge[7].head = 6;
pd->edge[7].headpos = 0;
pd->edge[7].tail = 2;
pd->edge[7].tailpos = 3;

pd->edge[8].head = 7;
pd->edge[8].headpos = 0;
pd->edge[8].tail = 6;
pd->edge[8].tailpos = 2;

pd->edge[9].head = 3;
pd->edge[9].headpos = 3;
pd->edge[9].tail = 7;
pd->edge[9].tailpos = 2;

pd->edge[10].head = 8;
pd->edge[10].headpos = 0;
pd->edge[10].tail = 3;
pd->edge[10].tailpos = 1;

pd->edge[11].head = 0;
pd->edge[11].headpos = 2;
pd->edge[11].tail = 8;
pd->edge[11].tailpos = 2;

pd->edge[12].head = 5;
pd->edge[12].headpos = 1;
pd->edge[12].tail = 8;
pd->edge[12].tailpos = 1;

pd->edge[13].head = 1;
pd->edge[13].headpos = 1;
pd->edge[13].tail = 5;
pd->edge[13].tailpos = 3;

pd->edge[14].head = 6;
pd->edge[14].headpos = 3;
pd->edge[14].tail = 1;
pd->edge[14].tailpos = 3;

pd->edge[15].head = 7;
pd->edge[15].headpos = 1;
pd->edge[15].tail = 6;
pd->edge[15].tailpos = 1;

pd->edge[16].head = 4;
pd->edge[16].headpos = 3;
pd->edge[16].tail = 7;
pd->edge[16].tailpos = 3;

pd->edge[17].head = 8;
pd->edge[17].headpos = 3;
pd->edge[17].tail = 4;
pd->edge[17].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 12;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;
pd->comp[0].edge[10] = 10;
pd->comp[0].edge[11] = 11;

pd->comp[1].nedges = 6;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 12;
pd->comp[1].edge[1] = 13;
pd->comp[1].edge[2] = 14;
pd->comp[1].edge[3] = 15;
pd->comp[1].edge[4] = 16;
pd->comp[1].edge[5] = 17;


/* Face data */

pd->face[0].nedges = 5;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 14;
pd->face[0].or[1] = 1;

pd->face[0].edge[2] = 8;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 16;
pd->face[0].or[3] = 1;

pd->face[0].edge[4] = 4;
pd->face[0].or[4] = 1;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 2;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 6;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 12;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 10;
pd->face[1].or[3] = 0;

pd->face[2].nedges = 4;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 2;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 9;
pd->face[2].or[1] = 0;

pd->face[2].edge[2] = 15;
pd->face[2].or[2] = 0;

pd->face[2].edge[3] = 7;
pd->face[2].or[3] = 0;

pd->face[3].nedges = 3;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 0;
pd->face[3].or[0] = 0;

pd->face[3].edge[1] = 5;
pd->face[3].or[1] = 1;

pd->face[3].edge[2] = 13;
pd->face[3].or[2] = 1;

pd->face[4].nedges = 3;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 1;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 7;
pd->face[4].or[1] = 1;

pd->face[4].edge[2] = 14;
pd->face[4].or[2] = 0;

pd->face[5].nedges = 3;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 1;
pd->face[5].or[0] = 0;

pd->face[5].edge[1] = 13;
pd->face[5].or[1] = 0;

pd->face[5].edge[2] = 6;
pd->face[5].or[2] = 1;

pd->face[6].nedges = 3;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 3;
pd->face[6].or[0] = 0;

pd->face[6].edge[1] = 10;
pd->face[6].or[1] = 1;

pd->face[6].edge[2] = 17;
pd->face[6].or[2] = 0;

pd->face[7].nedges = 3;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 3;
pd->face[7].or[0] = 1;

pd->face[7].edge[1] = 16;
pd->face[7].or[1] = 0;

pd->face[7].edge[2] = 9;
pd->face[7].or[2] = 1;

pd->face[8].nedges = 3;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].or = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].or != NULL);

pd->face[8].edge[0] = 4;
pd->face[8].or[0] = 0;

pd->face[8].edge[1] = 17;
pd->face[8].or[1] = 1;

pd->face[8].edge[2] = 11;
pd->face[8].or[2] = 1;

pd->face[9].nedges = 3;
pd->face[9].edge = calloc(pd->face[9].nedges,sizeof(pd_idx_t));
pd->face[9].or = calloc(pd->face[9].nedges,sizeof(pd_or_t));
assert(pd->face[9].edge != NULL);
assert(pd->face[9].or != NULL);

pd->face[9].edge[0] = 5;
pd->face[9].or[0] = 0;

pd->face[9].edge[1] = 11;
pd->face[9].or[1] = 0;

pd->face[9].edge[2] = 12;
pd->face[9].or[2] = 1;

pd->face[10].nedges = 2;
pd->face[10].edge = calloc(pd->face[10].nedges,sizeof(pd_idx_t));
pd->face[10].or = calloc(pd->face[10].nedges,sizeof(pd_or_t));
assert(pd->face[10].edge != NULL);
assert(pd->face[10].or != NULL);

pd->face[10].edge[0] = 8;
pd->face[10].or[0] = 0;

pd->face[10].edge[1] = 15;
pd->face[10].or[1] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



pd_code_t *pd_create_r2_testO_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(7);
assert(pd != NULL);
pd->ncross = 7;
pd->nedges = 14;
pd->ncomps = 2;
pd->nfaces = 9;
sprintf(pd->hash,"%s","Bw4JBAQDAwMDAwMCAgoEAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 4;
pd->cross[0].edge[2] = 9;
pd->cross[0].edge[3] = 5;
pd->cross[0].sign = 1;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 11;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 12;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 6;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 7;
pd->cross[2].sign = 0;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 8;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 7;
pd->cross[3].sign = 0;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 13;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 12;
pd->cross[4].sign = 1;

pd->cross[5].edge[0] = 5;
pd->cross[5].edge[1] = 10;
pd->cross[5].edge[2] = 6;
pd->cross[5].edge[3] = 11;
pd->cross[5].sign = 1;

pd->cross[6].edge[0] = 8;
pd->cross[6].edge[1] = 10;
pd->cross[6].edge[2] = 9;
pd->cross[6].edge[3] = 13;
pd->cross[6].sign = 1;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 0;
pd->edge[4].headpos = 1;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 5;
pd->edge[5].headpos = 0;
pd->edge[5].tail = 0;
pd->edge[5].tailpos = 3;

pd->edge[6].head = 2;
pd->edge[6].headpos = 1;
pd->edge[6].tail = 5;
pd->edge[6].tailpos = 2;

pd->edge[7].head = 3;
pd->edge[7].headpos = 3;
pd->edge[7].tail = 2;
pd->edge[7].tailpos = 3;

pd->edge[8].head = 6;
pd->edge[8].headpos = 0;
pd->edge[8].tail = 3;
pd->edge[8].tailpos = 1;

pd->edge[9].head = 0;
pd->edge[9].headpos = 2;
pd->edge[9].tail = 6;
pd->edge[9].tailpos = 2;

pd->edge[10].head = 5;
pd->edge[10].headpos = 1;
pd->edge[10].tail = 6;
pd->edge[10].tailpos = 1;

pd->edge[11].head = 1;
pd->edge[11].headpos = 1;
pd->edge[11].tail = 5;
pd->edge[11].tailpos = 3;

pd->edge[12].head = 4;
pd->edge[12].headpos = 3;
pd->edge[12].tail = 1;
pd->edge[12].tailpos = 3;

pd->edge[13].head = 6;
pd->edge[13].headpos = 3;
pd->edge[13].tail = 4;
pd->edge[13].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 10;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;

pd->comp[1].nedges = 4;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 10;
pd->comp[1].edge[1] = 11;
pd->comp[1].edge[2] = 12;
pd->comp[1].edge[3] = 13;


/* Face data */

pd->face[0].nedges = 4;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 7;
pd->face[0].or[1] = 1;

pd->face[0].edge[2] = 3;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 12;
pd->face[0].or[3] = 0;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 2;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 6;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 10;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 8;
pd->face[1].or[3] = 0;

pd->face[2].nedges = 3;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 5;
pd->face[2].or[1] = 1;

pd->face[2].edge[2] = 11;
pd->face[2].or[2] = 1;

pd->face[3].nedges = 3;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 0;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 12;
pd->face[3].or[1] = 1;

pd->face[3].edge[2] = 4;
pd->face[3].or[2] = 1;

pd->face[4].nedges = 3;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 1;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 11;
pd->face[4].or[1] = 0;

pd->face[4].edge[2] = 6;
pd->face[4].or[2] = 1;

pd->face[5].nedges = 3;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 3;
pd->face[5].or[0] = 0;

pd->face[5].edge[1] = 8;
pd->face[5].or[1] = 1;

pd->face[5].edge[2] = 13;
pd->face[5].or[2] = 0;

pd->face[6].nedges = 3;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 4;
pd->face[6].or[0] = 0;

pd->face[6].edge[1] = 13;
pd->face[6].or[1] = 1;

pd->face[6].edge[2] = 9;
pd->face[6].or[2] = 1;

pd->face[7].nedges = 3;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 5;
pd->face[7].or[0] = 0;

pd->face[7].edge[1] = 9;
pd->face[7].or[1] = 0;

pd->face[7].edge[2] = 10;
pd->face[7].or[2] = 1;

pd->face[8].nedges = 2;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].or = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].or != NULL);

pd->face[8].edge[0] = 2;
pd->face[8].or[0] = 1;

pd->face[8].edge[1] = 7;
pd->face[8].or[1] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



/*              _______>__   __________                                       */
/*             /          \0/          0                                      */
/*            4            \            \                                     */
/*           /    (8)     / \     (3)    \                                    */
/*          /            ^   5            \                                   */
/*          |           11(9)  >           <                                   */
/*          ____17__>____12_>__\ _13_>   /                                    */
/*         / 4\  (6)  /8       5\ (5) \ /      (0)                            */
/*        /    ^    10    (1)    6     /1                                     */
/*       ^      3______<____2_____|__1/ \                                     */
/*      /          3|            2|      \                                    */
/*     /            |             |       \                                   */
/*    /             |             |        \                                  */
/*   /              |             |         \                                 */
/*  |      (7)      ^             |          \                                */
/*  |               |             |           \                               */
/*  |               |     (2)     |           14                               */
/*  |               |             7    (4)    |                               */
/*  |               |             |           |                               */
/*  |               |             |           |                               */
/*  |               9            \/           |                               */
/* 16               |             |           |                               */
/*  \               |             |          /                                */
/*   \              |             |         /                                 */
/*    \_________________<____15_______<____/                                  */
/*                 7\             /6                                          */
/*                   \   (10)    /                                            */
/*                    \__<___8__/                                             */
/*                                                                            */
/*                                                                            */
/*              ____>_____   __________                                       */
/*             /          \0/          \                                      */
/*            /            \            0                                     */
/*           4    (6)     / \    (2)     \                                    */
/*          /            9   5            \                                   */
/*          |           / (7) >           |                                   */
/*         4___13___>____10_>__\ _11_>   <                                    */
/*         /  \  (5)  /6       5\ (4) \ /                                     */
/*        /    ^     8   (1)     6     /1          (3)                        */
/*       /      3_______<____2____|__1/ \                                     */
/*      /          3|            2|      \                                    */
/*     /            |             |       \                                   */
/*    /             |    (8)      |        \                                  */
/*   /              |             |         \                                 */
/*  |               ^             /          \                                */
/*  |                \           /            \                               */
/*  |                 \_______7_/             |                               */
/*  |                                         |                               */
/*  |                                         |                               */
/*  |                          (0)            |                               */
/*  |                                         |                               */
/*  |                                         |                               */
/*  \                                        /                                */
/*   \                                      /                                 */
/*    \______<_____________12_____________ /                                  */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/* r2 at crossings 7 and 6                                                    */



bool r2testO() {
  
  printf("---------------------------------------\n"
	 "r2 test O\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r2_testO_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after-0 configuration...");
  pd_code_t *after0 = pd_create_r2_testO_after_0();
  printf("done (passes pd_ok)\n");

  pd_code_t *after[1] = {after0};
    
  printf("performing r2 at crossings 7 and 6...\n");

  pd_idx_t  noutpd;
  pd_code_t **newpd;
  pd_idx_t  cr[2] = {7,6};

  pd_R2_bigon_elimination(before,cr,&noutpd,&newpd);

  if (noutpd != 1) { 
    printf("fail (did not generate 1 output pd codes)\n");
    return false;
  }

  pd_idx_t ct;

  for(ct=0;ct<noutpd;ct++) { 

    if (!pd_ok(newpd[ct])) { 
      printf("fail (output pd %d does not pass pd_ok)\n",ct);
      return false;
    }
    printf("\t output pd %d passes pd_ok...pass\n",ct);
    printf("\t testing for isomorphism with after configuration %d...",ct);
    if (!pd_isomorphic(newpd[ct],after[ct])) { 
      pd_printf("fail (output pd %PD ",newpd[ct]);
      pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after[ct]);
      return false;
    }

    printf("pass (output %d and after %d isomorphic)\n",ct,ct);

  }

  printf("performing r2 at crossings 7 and 6...pass\n");
  printf("housecleaning...");

  pd_code_free(&before);
  pd_code_free(&after0);
  pd_code_free(&newpd[0]);

  free(newpd);

  printf("done\n");

  printf("---------------------------------------\n"
	 "r2 test O : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}

/* r2 test P assembled c code */

pd_code_t *pd_create_r2_testP_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(8);
assert(pd != NULL);
pd->ncross = 8;
pd->nedges = 16;
pd->ncomps = 3;
pd->nfaces = 10;
sprintf(pd->hash,"%s","CBAKBgQEBAQDAgICAQMIBAQAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 4;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 3;
pd->cross[0].sign = 1;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 8;
pd->cross[1].edge[2] = 7;
pd->cross[1].edge[3] = 11;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 15;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 12;
pd->cross[2].sign = 0;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 13;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 12;
pd->cross[3].sign = 0;

pd->cross[4].edge[0] = 4;
pd->cross[4].edge[1] = 6;
pd->cross[4].edge[2] = 5;
pd->cross[4].edge[3] = 5;
pd->cross[4].sign = 1;

pd->cross[5].edge[0] = 6;
pd->cross[5].edge[1] = 11;
pd->cross[5].edge[2] = 7;
pd->cross[5].edge[3] = 10;
pd->cross[5].sign = 1;

pd->cross[6].edge[0] = 8;
pd->cross[6].edge[1] = 13;
pd->cross[6].edge[2] = 9;
pd->cross[6].edge[3] = 14;
pd->cross[6].sign = 0;

pd->cross[7].edge[0] = 9;
pd->cross[7].edge[1] = 15;
pd->cross[7].edge[2] = 10;
pd->cross[7].edge[3] = 14;
pd->cross[7].sign = 0;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 1;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 0;
pd->edge[3].headpos = 3;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 4;
pd->edge[4].headpos = 0;
pd->edge[4].tail = 0;
pd->edge[4].tailpos = 1;

pd->edge[5].head = 4;
pd->edge[5].headpos = 3;
pd->edge[5].tail = 4;
pd->edge[5].tailpos = 2;

pd->edge[6].head = 5;
pd->edge[6].headpos = 0;
pd->edge[6].tail = 4;
pd->edge[6].tailpos = 1;

pd->edge[7].head = 1;
pd->edge[7].headpos = 2;
pd->edge[7].tail = 5;
pd->edge[7].tailpos = 2;

pd->edge[8].head = 6;
pd->edge[8].headpos = 0;
pd->edge[8].tail = 1;
pd->edge[8].tailpos = 1;

pd->edge[9].head = 7;
pd->edge[9].headpos = 0;
pd->edge[9].tail = 6;
pd->edge[9].tailpos = 2;

pd->edge[10].head = 5;
pd->edge[10].headpos = 3;
pd->edge[10].tail = 7;
pd->edge[10].tailpos = 2;

pd->edge[11].head = 1;
pd->edge[11].headpos = 3;
pd->edge[11].tail = 5;
pd->edge[11].tailpos = 1;

pd->edge[12].head = 3;
pd->edge[12].headpos = 3;
pd->edge[12].tail = 2;
pd->edge[12].tailpos = 3;

pd->edge[13].head = 6;
pd->edge[13].headpos = 1;
pd->edge[13].tail = 3;
pd->edge[13].tailpos = 1;

pd->edge[14].head = 7;
pd->edge[14].headpos = 3;
pd->edge[14].tail = 6;
pd->edge[14].tailpos = 3;

pd->edge[15].head = 2;
pd->edge[15].headpos = 1;
pd->edge[15].tail = 7;
pd->edge[15].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 8;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;

pd->comp[1].nedges = 4;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 8;
pd->comp[1].edge[1] = 9;
pd->comp[1].edge[2] = 10;
pd->comp[1].edge[3] = 11;

pd->comp[2].nedges = 4;
pd->comp[2].tag = 'C';

pd->comp[2].edge = calloc(pd->comp[2].nedges,sizeof(pd_idx_t));
assert(pd->comp[2].edge != NULL);

pd->comp[2].edge[0] = 12;
pd->comp[2].edge[1] = 13;
pd->comp[2].edge[2] = 14;
pd->comp[2].edge[3] = 15;


/* Face data */

pd->face[0].nedges = 6;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 4;
pd->face[0].or[1] = 1;

pd->face[0].edge[2] = 5;
pd->face[0].or[2] = 0;

pd->face[0].edge[3] = 6;
pd->face[0].or[3] = 1;

pd->face[0].edge[4] = 10;
pd->face[0].or[4] = 0;

pd->face[0].edge[5] = 15;
pd->face[0].or[5] = 1;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 0;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 3;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 13;
pd->face[1].or[2] = 1;

pd->face[1].edge[3] = 8;
pd->face[1].or[3] = 0;

pd->face[2].nedges = 4;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 11;
pd->face[2].or[1] = 0;

pd->face[2].edge[2] = 6;
pd->face[2].or[2] = 0;

pd->face[2].edge[3] = 4;
pd->face[2].or[3] = 0;

pd->face[3].nedges = 4;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 2;
pd->face[3].or[0] = 0;

pd->face[3].edge[1] = 15;
pd->face[3].or[1] = 0;

pd->face[3].edge[2] = 9;
pd->face[3].or[2] = 0;

pd->face[3].edge[3] = 13;
pd->face[3].or[3] = 0;

pd->face[4].nedges = 4;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 7;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 8;
pd->face[4].or[1] = 1;

pd->face[4].edge[2] = 14;
pd->face[4].or[2] = 1;

pd->face[4].edge[3] = 10;
pd->face[4].or[3] = 1;

pd->face[5].nedges = 3;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 1;
pd->face[5].or[0] = 1;

pd->face[5].edge[1] = 12;
pd->face[5].or[1] = 1;

pd->face[5].edge[2] = 3;
pd->face[5].or[2] = 1;

pd->face[6].nedges = 2;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 2;
pd->face[6].or[0] = 1;

pd->face[6].edge[1] = 12;
pd->face[6].or[1] = 0;

pd->face[7].nedges = 2;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 7;
pd->face[7].or[0] = 0;

pd->face[7].edge[1] = 11;
pd->face[7].or[1] = 1;

pd->face[8].nedges = 2;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].or = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].or != NULL);

pd->face[8].edge[0] = 9;
pd->face[8].or[0] = 1;

pd->face[8].edge[1] = 14;
pd->face[8].or[1] = 0;

pd->face[9].nedges = 1;
pd->face[9].edge = calloc(pd->face[9].nedges,sizeof(pd_idx_t));
pd->face[9].or = calloc(pd->face[9].nedges,sizeof(pd_or_t));
assert(pd->face[9].edge != NULL);
assert(pd->face[9].or != NULL);

pd->face[9].edge[0] = 5;
pd->face[9].or[0] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



pd_code_t *pd_create_r2_testP_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(6);
assert(pd != NULL);
pd->ncross = 6;
pd->nedges = 12;
pd->ncomps = 3;
pd->nfaces = 8;
sprintf(pd->hash,"%s","BgwICAQDAgICAgEDBgQCAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 4;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 3;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 5;
pd->cross[1].edge[2] = 5;
pd->cross[1].edge[3] = 4;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 9;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 6;
pd->cross[2].sign = 0;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 7;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 6;
pd->cross[3].sign = 0;

pd->cross[4].edge[0] = 7;
pd->cross[4].edge[1] = 11;
pd->cross[4].edge[2] = 8;
pd->cross[4].edge[3] = 10;
pd->cross[4].sign = 0;

pd->cross[5].edge[0] = 8;
pd->cross[5].edge[1] = 11;
pd->cross[5].edge[2] = 9;
pd->cross[5].edge[3] = 10;
pd->cross[5].sign = 0;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 1;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 0;
pd->edge[3].headpos = 3;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 1;
pd->edge[4].headpos = 3;
pd->edge[4].tail = 0;
pd->edge[4].tailpos = 1;

pd->edge[5].head = 1;
pd->edge[5].headpos = 2;
pd->edge[5].tail = 1;
pd->edge[5].tailpos = 1;

pd->edge[6].head = 3;
pd->edge[6].headpos = 3;
pd->edge[6].tail = 2;
pd->edge[6].tailpos = 3;

pd->edge[7].head = 4;
pd->edge[7].headpos = 0;
pd->edge[7].tail = 3;
pd->edge[7].tailpos = 1;

pd->edge[8].head = 5;
pd->edge[8].headpos = 0;
pd->edge[8].tail = 4;
pd->edge[8].tailpos = 2;

pd->edge[9].head = 2;
pd->edge[9].headpos = 1;
pd->edge[9].tail = 5;
pd->edge[9].tailpos = 2;

pd->edge[10].head = 4;
pd->edge[10].headpos = 3;
pd->edge[10].tail = 5;
pd->edge[10].tailpos = 3;

pd->edge[11].head = 5;
pd->edge[11].headpos = 1;
pd->edge[11].tail = 4;
pd->edge[11].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 6;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;

pd->comp[1].nedges = 4;
pd->comp[1].tag = 'C';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 6;
pd->comp[1].edge[1] = 7;
pd->comp[1].edge[2] = 8;
pd->comp[1].edge[3] = 9;

pd->comp[2].nedges = 2;
pd->comp[2].tag = 'B';

pd->comp[2].edge = calloc(pd->comp[2].nedges,sizeof(pd_idx_t));
assert(pd->comp[2].edge != NULL);

pd->comp[2].edge[0] = 10;
pd->comp[2].edge[1] = 11;


/* Face data */

pd->face[0].nedges = 8;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 3;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 7;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 10;
pd->face[0].or[3] = 0;

pd->face[0].edge[4] = 9;
pd->face[0].or[4] = 1;

pd->face[0].edge[5] = 1;
pd->face[0].or[5] = 0;

pd->face[0].edge[6] = 4;
pd->face[0].or[6] = 1;

pd->face[0].edge[7] = 5;
pd->face[0].or[7] = 0;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 2;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 9;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 11;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 7;
pd->face[1].or[3] = 0;

pd->face[2].nedges = 3;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 1;
pd->face[2].or[0] = 1;

pd->face[2].edge[1] = 6;
pd->face[2].or[1] = 1;

pd->face[2].edge[2] = 3;
pd->face[2].or[2] = 1;

pd->face[3].nedges = 2;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 0;
pd->face[3].or[0] = 0;

pd->face[3].edge[1] = 4;
pd->face[3].or[1] = 0;

pd->face[4].nedges = 2;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 2;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 6;
pd->face[4].or[1] = 0;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 8;
pd->face[5].or[0] = 1;

pd->face[5].edge[1] = 10;
pd->face[5].or[1] = 1;

pd->face[6].nedges = 2;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 8;
pd->face[6].or[0] = 0;

pd->face[6].edge[1] = 11;
pd->face[6].or[1] = 1;

pd->face[7].nedges = 1;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 5;
pd->face[7].or[0] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}



/*                                                                            */
/*            __<______9___                                                   */
/*           /             \                                                  */
/*          /        ______________________<__________                        */
/*         /  (8)  14      6                          \                       */
/*        /        /          \                        \                      */
/*        \       /            ^                        \                     */
/*         \     /              \                        \         (3)        */
/*          \   <     (4)        8              (1)       \                   */
/*           \ /              ____\ _________              \                  */
/*            \              /    1\         0             13                 */
/*           /7\            ^      /          \              \                */
/*          /   \          /  (7) /            \              \               */
/*         /    10        7      ^              \              \              */
/*         |      \      /      /                \              |             */
/*         |       >__________11                  \             | ____        */
/*         |         5                 (2)         \            |/    ^       */
/*         |          /                             \           /3     \      */
/*         |          ^                              >         /|       \     */
/*         |           6                              \       / |        2    */
/*         |            \ _______<__________4__________\___<_3  ^  (6)   /    */
/*         |             \                             0\       |       /     */
/*         |            /4\                              \ (5) 12      /      */
/*          \          /   \                              \     |     /       */
/*           \        /(9)  ^                              1_>_ /____/        */
/*            \      /      /      (0)                         /2             */
/*             \     \__5__/                                  /               */
/*              \                                            /                */
/*               \__________________15___________>__________/                 */
/*                                                                            */
/*                                                                            */
/*            ___11________                                                   */
/*           /             \                                                  */
/*          /        _________________<_____________7_                        */
/*         /  (6)   8      4                          \                       */
/*        /        /          \                        \         (1)          */
/*        \       /            \                        \                     */
/*         >     /              \                        \                    */
/*          \   <       (5)     ^                         \                   */
/*           \ /                /                          \                  */
/*           5\                /                            \                 */
/*           / \              /                              \                */
/*          /   \_________10_/                                \               */
/*         /                                                   \              */
/*         |             ___________0_________>___              |             */
/*         |            /                         \             | __<_        */
/*         |           /                           \            |/    \       */
/*         |          /              (3)            \           /      \      */
/*         |    (0)   \                              \         /|3      \     */
/*         |           \                              \       3 ^        \    */
/*         |            \ _______<________4____________\__<__/  |  (4)   2    */
/*         |            1\                             0\       |       /     */
/*         |            / \                              \ (2)  6      /      */
/*          \          /   ^                              1     |     /       */
/*           \        /     \                              \>__ /____/        */
/*            \      / (7)  /                                 2/              */
/*             \     \__5__/                                  /               */
/*              \                                            /                */
/*               \_________________9______________>_________/                 */
/*                                                                            */
/* r2 at crossings 1 and 5 */


bool r2testP() {
  
  printf("---------------------------------------\n"
	 "r2 test P\n"
	 "---------------------------------------\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_r2_testP_before_0();
  printf("done (passes pd_ok)\n");

  printf("creating after-0 configuration...");
  pd_code_t *after0 = pd_create_r2_testP_after_0();
  printf("done (passes pd_ok)\n");

  pd_code_t *after[1] = {after0};
    
  printf("performing r2 at crossings 1 and 5...\n");

  pd_idx_t  noutpd;
  pd_code_t **newpd;
  pd_idx_t  cr[2] = {1,5};

  pd_R2_bigon_elimination(before,cr,&noutpd,&newpd);

  if (noutpd != 1) { 
    printf("fail (did not generate 1 output pd codes)\n");
    return false;
  }

  pd_idx_t ct;

  for(ct=0;ct<noutpd;ct++) { 

    if (!pd_ok(newpd[ct])) { 
      printf("fail (output pd %d does not pass pd_ok)\n",ct);
      return false;
    }
    printf("\t output pd %d passes pd_ok...pass\n",ct);
    printf("\t testing for isomorphism with after configuration %d...",ct);
    if (!pd_isomorphic(newpd[ct],after[ct])) { 
      pd_printf("fail (output pd %PD ",newpd[ct]);
      pd_printf("is not isomorphism to expected \"after\" configuration %PD\n",after[ct]);
      return false;
    }

    printf("pass (output %d and after %d isomorphic)\n",ct,ct);

  }

  printf("performing r2 at crossings 1 and 5...pass\n");
  printf("housecleaning...");

  pd_code_free(&before);
  pd_code_free(&after0);
  pd_code_free(&newpd[0]);

  free(newpd);

  printf("done\n");

  printf("---------------------------------------\n"
	 "r2 test P : PASS\n"
	 "---------------------------------------\n");
  
  return true;

}



int main() {

  printf("test_crossingops (%s)\n",PACKAGE_STRING);
  printf("--------------------------------------------------------\n"
	 "Unit tests for pdcode Reidemeister move primitives. \n"
	 "========================================================\n");

  if (!r2testA() || !r2testB() || !r2testC() || !r2testD() || !r2testE() || !r2testF() || !r2testG() || !r2testH() || !r2testI() || !r2testJ() || !r2testL() || !r2testN() || !r2testO() || !r2testP() || !r1testA() || !r1testB() || !r1testC() || !r1testD() || !r1testE() || !r2_tests() || !compacting_copy_tests() || !unravel_unknot_test(10)) {

    printf("\n=======================================================\n");
    printf("test_crossingops:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================================\n");
    printf("test_crossingops:  PASS.\n");
    exit(0);

  }

  return 0;

}
    
