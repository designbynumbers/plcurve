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

  for(crossings_left=n;crossings_left > 1;crossings_left--) { 

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

  return true;
    
}

int main() {

  printf("test_crossingops (%s)\n",PACKAGE_STRING);
  printf("--------------------------------------------------------\n"
	 "Unit tests for pdcode Reidemeister move primitives. \n"
	 "========================================================\n");

  if (!compacting_copy_tests() || !unravel_unknot_test(10)) {

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
    
