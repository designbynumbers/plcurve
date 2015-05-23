/* 

   test_invariants.c : Unit tests for the code in pd_invariantss.c.


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

#ifdef HAVE_GSL_GSL_PERMUTATION_H
   #include<gsl/gsl_permutation.h>
#endif

#ifdef HAVE_STDBOOL_H
   #include<stdbool.h>
#endif

#ifdef HAVE_ASSERT_H
   #include<assert.h>
#endif

#ifdef HAVE_ARGTABLE2_H
   #include<argtable2.h> 
#endif

#include<plcTopology.h>

#include<pd_multidx.h>
#include<pd_perm.h>
#include<pd_dihedral.h>

#include<pd_invariants.h>

int PD_VERBOSE=50;

bool unknot_interlacement_test(pd_idx_t ncross)
/*

  Since the unknot we're generating looks like
  
  /-\/-\/-\...../-\
  \-/\-/\-/     \-/

  we should have no interlacements at all. 

*/

{
  pd_code_t *pdA;
  
  printf("%d-crossing unknot interlacement test  \n",ncross);
  printf("---------------------------------------\n");

  pdA = pd_build_unknot(ncross);

  printf("built %d crossing unknot...done\n",ncross);
  printf("checking number of interlaced crossing pairs...");
  unsigned int *interlacements = pd_interlaced_crossings(pdA);

  if (interlacements == NULL) {

    printf("FAIL (didn't return a buffer of interlacements)\n");
    return false;

  }

  if (interlacements[0] != 0) {

    printf("FAIL (expected 0 interlacements, got %u)\n",interlacements[0]);
    return false;

  }

  printf("pass (expected 0, got 0)\n");

  pd_code_free(&pdA);
  free(interlacements);
  
  printf("-----------------------------------------------\n");
  printf("%d-crossing unknot interlacement test ... PASS \n"
	 "-----------------------------------------------\n\n",ncross);
 
  return true;

}  

bool twist_interlacement_test(pd_idx_t ntwist)

/* The twist knot we're generating looks like this:
			      	              
    +----<--------+  +---+   +----+    +--<--+	  +--<-+
    |	      +---|------|->------|------->--|-------+ |
    |	      |	  +--+	 +---+	  +----+     +----+  | |
    |  +----+ |	  1  2 	 3     	       	       	  n  | |
    +--|->----+	       	       	       	       	     | |
     A |    | B	       	       	      		     | |
       |    +----------<----------------<------------- |       	       	     
       +-->------------------->------------------------+ 

   Now none of the crossing pairs <digit>-<digit> are interlaced,
   because we visit them in increasing order going "out" (left-right)
   and decreasing order coming "back" (right-left). 

   Whether A and B are interlaced with each other depends on 
   whether the number of twist crossings is even (NO) or odd (YES).

   On the other hand, any pair in the form <A or B>-<digit>
   ARE interlaced since we visit 

   LETTERS -> NUMBERS -> LETTERS -> NUMBERS

   in order around the curve. So I claim that there should be 2n
   interlacements in an n-twist knot with n even, and 2n+1
   interlacements in a n-twist knot with n odd.  

   (And n+2 crossings overall.)

*/		      
       	       	       	       	       	       	       	       	       	      
{							   
  pd_code_t *pdA;
  
  printf("%d-twist knot interlacement test  \n",ntwist);
  printf("-------------------------------------------\n");

  if (ntwist < 3) {

    printf("FAIL. n-twist knot test not implemented for n < 3.\n");
    return false;

  }

  pdA = pd_build_twist_knot(ntwist);

  printf("built %d twist knot...done\n",ntwist);

  printf("checking number of crossings...");
  if (pdA->ncross != ntwist+2) {

    printf("FAIL (expected %d+2 == %d crossings, got %d)\n",
	   ntwist+2,ntwist,pdA->ncross);
    return false;
    
  }
  printf("pass (expected %d, got %d)\n",
	 ntwist+2,pdA->ncross);
  
  printf("checking number of interlaced crossing pairs...");
  unsigned int *interlacements = pd_interlaced_crossings(pdA);

  if (interlacements == NULL) {

    printf("FAIL (didn't return a buffer of interlacements)\n");
    return false;

  }

  if (ntwist % 2 == 0) { 

    if (interlacements[0] != 2*ntwist) {

      printf("FAIL (expected %d interlacements, got %u)\n",
	     2*ntwist,interlacements[0]);
      return false;

    }

    printf("pass (expected %d, got %u)\n",2*ntwist,interlacements[0]);

  } else {

    if (interlacements[0] != 2*ntwist+1) {

      printf("FAIL (expected %d interlacements, got %u)\n",
	     2*ntwist+1,interlacements[0]);
      return false;

    }

    printf("pass (expected %d, got %u)\n",2*ntwist+1,interlacements[0]);

  }

  pd_code_free(&pdA);
  free(interlacements);
  
  printf("------------------------------------------\n");
  printf("%d-twist knot interlacement test ... PASS \n"
	 "------------------------------------------\n\n",ntwist);

  return true;

}  


bool torusknot_interlacement_test(pd_idx_t q)

/* 
   In the 2-q torus knot, ALL pairs of crossings are interlaced.
   So the answer should be q(q+1)/2.

   In the 2-q torus link, NO pairs of crossings are interlaced.
   So the answer should be (0,0).
*/		      
       	       	       	       	       	       	       	       	       	      
{							   
  pd_code_t *pdA;
  
  printf("-------------------------------------------\n"
	 "(2,%d)-torus knot/link interlacement test  \n",q);
  printf("-------------------------------------------\n");

  pdA = pd_build_torus_knot(2,q);
  printf("built (2,%d) twist knot...done\n",q);
  
  printf("checking number of interlaced crossing pairs...");
  unsigned int *interlacements = pd_interlaced_crossings(pdA);

  if (interlacements == NULL) {

    printf("FAIL (didn't return a buffer of interlacements)\n");
    return false;
 
  }

  if (q % 2 == 0) {

    int i; 

    for(i=0;i<2;i++) { 

      if (interlacements[i] != 0) {

	printf("FAIL (expected 0 interlacements on component %d, got %u)\n",
	       i,interlacements[0]);
	return false;

      }

    }

    printf("pass (expected {0,0}, got {%u,%u})\n",
	   interlacements[0],interlacements[1]);

  } else {

    if (interlacements[0] != (q*(q-1))/2) {
      
      printf("FAIL (expected %d interlacements, got %u)\n",
	     (q*(q-1))/2,interlacements[0]);
      return false;

    }

    printf("pass (expected %d, got %u)\n",(q*(q-1))/2,interlacements[0]);

  }

  pd_code_free(&pdA);
  free(interlacements);
  
  printf("--------------------------------------------\n");
  printf("(2,%d)-torus knot interlacement test...PASS \n"
	 "--------------------------------------------\n\n",q);

  return true;

}  

bool test_easy_examples() {

  int store_verb = PD_VERBOSE;
  PD_VERBOSE = 0;

  printf("----------------------------------\n"
	 "pd_interlaced_crossings test suite\n"
	 "----------------------------------\n");

  int i;
  for(i=2;i<16;i++) { 

    if (!unknot_interlacement_test(i)) { return false; }

  }

  for(i=3;i<16;i++) {

    if (!twist_interlacement_test(i)) { return false; }

  }

  for(i=2;i<17;i++) {

    if (!torusknot_interlacement_test(i)) { return false; }

  }
  
  printf("----------------------------------------\n"
	 "pd_interlaced_crossings test suite: PASS\n"
	 "----------------------------------------\n");

  PD_VERBOSE = store_verb;

  return true;

}

int main() {

  printf("test_invariants (%s)\n",PACKAGE_STRING);
  printf("---------------------------------------\n"
	 "Unit tests for pd_invariants.c\n"
	 "=======================================\n");

  if (!test_easy_examples()) {

    printf("=====================================\n");
    printf("test_invariants:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================\n");
    printf("test_invariants:  PASS.\n");
    exit(0);

  }

  return 0;

}
  
  

  

  

  
