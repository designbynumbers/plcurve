/* 

   test_invariants.c : Unit tests for the code in pd_invariantss.c.


*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
#endif

#include<stdio.h>
#include<string.h>


#ifdef HAVE_STDINT_H
   #include<stdint.h>
#endif

#include<stdlib.h>


#ifdef HAVE_GSL_GSL_PERMUTATION_H
   #include<gsl/gsl_permutation.h>
#endif

#ifdef HAVE_STDBOOL_H
   #include<stdbool.h>
#endif

#include<assert.h>

#ifdef HAVE_ARGTABLE2_H
   #include<argtable2.h> 
#endif

#include<plcTopology.h>

#include<pd_multidx.h>
#include<pd_perm.h>
#include<pd_dihedral.h>

#include<pd_invariants.h>



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
  int *interlacements = pd_interlaced_crossings(pdA);

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

   (And n+2 crossings overall.) Now we need to think about the signs 
   of these interlacements. The numbered crossing signs alternate-- 
   assuming the straight section happens first and is the "e1" tangent
   vector, the "e2" tangent vector is going to be going up or down 
   depending on which crossing number we're on. 

   The crossings A and B have opposite signs as well, for the same
   reason.

   So we have the following interlacement signs:

   AB (if interlaced), is always a -1 interlacement
   A(numbers) and B(numbers) must cancel, since A and B have opposite signs.

   So the overall score should be -1 (if AB are interlaced),
   0 (otherwise).
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
  
  printf("checking signed number of interlaced crossing pairs...");
  int *interlacements = pd_interlaced_crossings(pdA);

  if (interlacements == NULL) {

    printf("FAIL (didn't return a buffer of interlacements)\n");
    return false;

  }

  if (ntwist % 2 == 0) { 

    if (interlacements[0] != 0) {

      printf("FAIL (expected %d signed interlacements, got %d)\n",
	     0,interlacements[0]);
      return false;

    }

    printf("pass (expected %d, got %u)\n",0,interlacements[0]);

  } else {

    if (interlacements[0] != -1) {

      printf("FAIL (expected %d interlacements, got %d)\n",
	     -1,interlacements[0]);
      return false;

    }

    printf("pass (expected %d, got %d)\n",-1,interlacements[0]);

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
   In the 2-q torus knot, where q = 2n + 1, ALL pairs of crossings are
   interlaced.  So the answer should be q(q-1)/2. However, you have to
   think about the sign convention.

   Looking at the diagram, you see that the crossings alternate
   -+-+-
   and so forth. Now there are an odd number of crossings, so 
   generally you have (n+1) - signs and n (+) signs. This means
   that we have to split up the (2n+1)(2n)/2 crossings in a 
   careful way to make sure that we understand how the sign 
   convention is going to play out.

   I think that we should get n(n+1) interlacements which are -+,
   and so have sign - overall. Which means that we should have 
   
   ((2n+1)n - (n)(n+1)) - (n(n+1))

   as the overall signed sum. Using Mathematica, this simplies to -n.
   Putting this in terms of q, we see n = (q-1)/2, so -n = -(q-1)/2,
   and that should be the answer. It's comforting that this agrees with 
   our hand calculation of -1 for the trefoil diagram.

   In the 2-q torus link, NO pairs of crossings are interlaced.
   So the answer should be (0,0).
*/		      
       	       	       	       	       	       	       	       	       	      
{							   
  pd_code_t *pdA;
  
  printf("-------------------------------------------\n"
	 "(2,%d)-torus knot/link interlacement test  \n",q);
  printf("-------------------------------------------\n");

  pdA = pd_build_torus_knot(2,q);
  printf("built (2,%d) torus knot...done\n",q);
  
  printf("checking signed number of interlaced crossing pairs...");
  int *interlacements = pd_interlaced_crossings(pdA);

  if (interlacements == NULL) {

    printf("FAIL (didn't return a buffer of interlacements)\n");
    return false;
 
  }

  if (q % 2 == 0) {

    int i; 

    for(i=0;i<2;i++) { 

      if (interlacements[i] != 0) {

	printf("FAIL (expected 0 interlacements on component %d, got %d)\n",
	       i,interlacements[0]);
	return false;

      }

    }

    printf("pass (expected {0,0}, got {%d,%d})\n",
	   interlacements[0],interlacements[1]);

  } else {

    if (interlacements[0] != -((int)(q)-1)/2) {
      
      printf("FAIL (expected %d interlacements, got %d)\n",
	     -((int)(q)-1)/2,interlacements[0]);
      return false;

    }

    printf("pass (expected %d, got %d)\n",-((int)(q)-1)/2,interlacements[0]);

  }

  pd_code_free(&pdA);
  free(interlacements);
  
  printf("--------------------------------------------\n");
  printf("(2,%d)-torus knot interlacement test...PASS \n"
	 "--------------------------------------------\n\n",q);

  return true;

}  

bool test_easy_examples() {

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
  
  

  

  

  
