/* 

   test_2crossing.c : Unit tests for the code in pdcode.c


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

#include <plCurve.h>
#include<plcTopology.h>

#include<pd_multidx.h>
#include<pd_perm.h>
#include<pd_dihedral.h>

#include<pd_isomorphisms.h>

int PD_VERBOSE=50;

bool two_crossing() {

  printf("----------------------------------\n"
	 "two-crossing regeneration tests\n"
	 "----------------------------------\n");
  
  PD_VERBOSE = 50;

  printf("building 2 crossing diagram...");

  pd_code_t *pdA = pd_code_new(2);
  pdA->ncross = 2;
  pdA->cross[0].edge[0] = 0;
  pdA->cross[0].edge[1] = 2;
  pdA->cross[0].edge[2] = 1;
  pdA->cross[0].edge[3] = 3;
  pdA->cross[0].sign = PD_NEG_ORIENTATION;

  pdA->cross[1].edge[0] = 0;
  pdA->cross[1].edge[1] = 3;
  pdA->cross[1].edge[2] = 1;
  pdA->cross[1].edge[3] = 2;
  pdA->cross[1].sign = PD_NEG_ORIENTATION;
  pd_regenerate(pdA);

  if (pd_ok(pdA)) { 

    printf("pass (passes pd_ok)\n");

  } else { 

    pd_printf("fail (pd %PD doesn't pass pd_ok)\n",pdA);
    return false;

  }

  printf("testing whether 2-crossing pd is isomorphic to self...");

  if (pd_isomorphic(pdA, pdA)) { 

    printf("pass (pd isomorphic to self)\n");
    
  } else { 

    printf("fail (pd is not isomorphic to itself)\n");
    return false;

  }

  return true;

}

int main() {

  printf("test_2crossing (%s)\n",PACKAGE_STRING);
  printf("---------------------------------------\n"
	 "Unit tests to see how well two-crossing diagrams are handled.\n"
	 "=======================================\n");

  if (!two_crossing()) {

    printf("=====================================\n");
    printf("test_pdcode:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================\n");
    printf("test_pdcode:  PASS.\n");
    exit(0);

  }

  return 0;

}
