/*

   test_readkt.c : Unit tests for the code in pdcode.c which reads KnotTheory-style PD codes.

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

extern int PD_VERBOSE;

void create_inputfile() {

  printf("opening 6a5and6n1.kt for writing...");

  FILE *outfile;
  outfile = fopen("6a5and6n1.kt","w");

  if (outfile == NULL) { 

    printf("fail. (Couldn't open file for writing.)\n");
    exit(1);

  }

  printf("done (opened ok).\n");

  printf("writing data and closing file...");

  fprintf(outfile,
	  "PD[X[6, 1, 7, 2], X[10, 3, 11, 4], X[12, 7, 9, 8], X[8, 11, 5, 12], X[2, 5, 3, 6], X[4, 9, 1, 10]]\n"
	  "PD[X[6, 1, 7, 2], X[12, 8, 9, 7], X[4, 12, 1, 11], X[5, 11, 6, 10], X[3, 8, 4, 5], X[9, 3, 10, 2]]\n"
	  );

  fclose(outfile);

  printf("done\n");

}
  
bool test_6a5_and_6n1() { 

   printf("-------------------------------------------\n"
	  "testing 6a5 and 6n1 KnotTheory pd codes\n"
	  "-------------------------------------------\n");
  
   create_inputfile();

   FILE *f;
   pd_code_t *l6a5, *l6n1;
   
   printf("reading KnotTheory-style pd codes from 6a5and6n1.kt...");
    
   f = fopen("6a5and6n1.kt", "r");
   l6a5 = pd_read_KnotTheory(f);
   l6n1 = pd_read_KnotTheory(f);
   fclose(f);

   printf("done (read pd codes and passed pd_ok).\n");

   printf("checking that 6a5 is alternating...");

   if (!pd_is_alternating(l6a5)) { 

     printf("fail (not alternating)\n");
     return false;

   } 

   printf("pass (alternating)\n");

   printf("checking that 6n5 is not alternating...");

   if (pd_is_alternating(l6n1)) { 

     printf("fail (alternating)\n");
     return false;

   }

   printf("pass\n");

   printf("-------------------------------------------\n"
	  "6a5 and 6n1 KnotTheory pd code test: PASS\n"
	  "-------------------------------------------\n");
  
   return true;
}

int main() {

  PD_VERBOSE = 50;

  printf("test_readkt (%s)\n",PACKAGE_STRING);
  printf("--------------------------------------------------------\n"
	 "Unit tests for reading KnotTheory-style pd codes. \n"
	 "========================================================\n");

  if (!test_6a5_and_6n1()) {

    printf("=======================================================\n");
    printf("test_readkt:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================================\n");
    printf("test_readkt:  PASS.\n");
    exit(0);

  }

  return 0;

}
