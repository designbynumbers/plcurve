/* 

   test_bad_homfly.c: A test for a failure of the pd_homfly code.


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

#ifdef HAVE_ASSERT_H
  #include<assert.h>
#endif

#include<plcTopology.h>

#include<pd_multidx.h>
#include<pd_perm.h>
#include<pd_dihedral.h>

#include<pd_isomorphisms.h>

int PD_VERBOSE=50;

pd_code_t *pd_create_8_7_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(8);
assert(pd != NULL);
pd->ncross = 8;
pd->nedges = 16;
pd->ncomps = 1;
pd->nfaces = 10;
sprintf(pd->hash,"%s","CBAKBgUEAwMDAgICAgEQAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 3;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 4;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 11;
pd->cross[1].edge[2] = 15;
pd->cross[1].edge[3] = 10;
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

pd->cross[4].edge[0] = 4;
pd->cross[4].edge[1] = 12;
pd->cross[4].edge[2] = 5;
pd->cross[4].edge[3] = 11;
pd->cross[4].sign = 1;

pd->cross[5].edge[0] = 5;
pd->cross[5].edge[1] = 12;
pd->cross[5].edge[2] = 6;
pd->cross[5].edge[3] = 13;
pd->cross[5].sign = 1;

pd->cross[6].edge[0] = 6;
pd->cross[6].edge[1] = 14;
pd->cross[6].edge[2] = 7;
pd->cross[6].edge[3] = 13;
pd->cross[6].sign = 1;

pd->cross[7].edge[0] = 7;
pd->cross[7].edge[1] = 14;
pd->cross[7].edge[2] = 8;
pd->cross[7].edge[3] = 15;
pd->cross[7].sign = 1;


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
pd->edge[3].headpos = 1;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 4;
pd->edge[4].headpos = 0;
pd->edge[4].tail = 0;
pd->edge[4].tailpos = 3;

pd->edge[5].head = 5;
pd->edge[5].headpos = 0;
pd->edge[5].tail = 4;
pd->edge[5].tailpos = 2;

pd->edge[6].head = 6;
pd->edge[6].headpos = 0;
pd->edge[6].tail = 5;
pd->edge[6].tailpos = 2;

pd->edge[7].head = 7;
pd->edge[7].headpos = 0;
pd->edge[7].tail = 6;
pd->edge[7].tailpos = 2;

pd->edge[8].head = 2;
pd->edge[8].headpos = 3;
pd->edge[8].tail = 7;
pd->edge[8].tailpos = 2;

pd->edge[9].head = 3;
pd->edge[9].headpos = 1;
pd->edge[9].tail = 2;
pd->edge[9].tailpos = 1;

pd->edge[10].head = 1;
pd->edge[10].headpos = 3;
pd->edge[10].tail = 3;
pd->edge[10].tailpos = 3;

pd->edge[11].head = 4;
pd->edge[11].headpos = 3;
pd->edge[11].tail = 1;
pd->edge[11].tailpos = 1;

pd->edge[12].head = 5;
pd->edge[12].headpos = 1;
pd->edge[12].tail = 4;
pd->edge[12].tailpos = 1;

pd->edge[13].head = 6;
pd->edge[13].headpos = 3;
pd->edge[13].tail = 5;
pd->edge[13].tailpos = 3;

pd->edge[14].head = 7;
pd->edge[14].headpos = 1;
pd->edge[14].tail = 6;
pd->edge[14].tailpos = 1;

pd->edge[15].head = 1;
pd->edge[15].headpos = 2;
pd->edge[15].tail = 7;
pd->edge[15].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 16;
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
pd->comp[0].edge[14] = 14;
pd->comp[0].edge[15] = 15;


/* Face data */

pd->face[0].nedges = 6;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 8;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 14;
pd->face[0].or[2] = 0;

pd->face[0].edge[3] = 6;
pd->face[0].or[3] = 0;

pd->face[0].edge[4] = 12;
pd->face[0].or[4] = 0;

pd->face[0].edge[5] = 4;
pd->face[0].or[5] = 0;

pd->face[1].nedges = 5;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 5;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 13;
pd->face[1].or[1] = 1;

pd->face[1].edge[2] = 7;
pd->face[1].or[2] = 1;

pd->face[1].edge[3] = 15;
pd->face[1].or[3] = 1;

pd->face[1].edge[4] = 11;
pd->face[1].or[4] = 1;

pd->face[2].nedges = 4;
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

pd->face[2].edge[3] = 8;
pd->face[2].or[3] = 1;

pd->face[3].nedges = 3;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 0;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 4;
pd->face[3].or[1] = 1;

pd->face[3].edge[2] = 11;
pd->face[3].or[2] = 0;

pd->face[4].nedges = 3;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 0;
pd->face[4].or[0] = 0;

pd->face[4].edge[1] = 10;
pd->face[4].or[1] = 0;

pd->face[4].edge[2] = 3;
pd->face[4].or[2] = 1;

pd->face[5].nedges = 3;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 1;
pd->face[5].or[0] = 0;

pd->face[5].edge[1] = 3;
pd->face[5].or[1] = 0;

pd->face[5].edge[2] = 9;
pd->face[5].or[2] = 0;

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

pd->face[7].edge[0] = 5;
pd->face[7].or[0] = 0;

pd->face[7].edge[1] = 12;
pd->face[7].or[1] = 1;

pd->face[8].nedges = 2;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].or = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].or != NULL);

pd->face[8].edge[0] = 6;
pd->face[8].or[0] = 1;

pd->face[8].edge[1] = 13;
pd->face[8].or[1] = 0;

pd->face[9].nedges = 2;
pd->face[9].edge = calloc(pd->face[9].nedges,sizeof(pd_idx_t));
pd->face[9].or = calloc(pd->face[9].nedges,sizeof(pd_or_t));
assert(pd->face[9].edge != NULL);
assert(pd->face[9].or != NULL);

pd->face[9].edge[0] = 7;
pd->face[9].or[0] = 0;

pd->face[9].edge[1] = 14;
pd->face[9].or[1] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}

bool test_bad_homfly() {

  printf("creating (assumed) 8_7 pd code...");
  pd_code_t *pd = pd_create_8_7_0();
   
  if (!pd_ok(pd)) {

    printf("fail (doesn't pass pd_ok).\n");
    return false;

  }

  printf("pass (pd_ok is true)\n");

  printf("computing HOMFLY...");
  char *homfly = pd_homfly(pd);
  
  if (!strcmp(homfly,"1")) { 

    printf("fail (HOMFLY == %s == 1)\n",homfly);
    return false;

  } 

  printf("tentative pass (HOMFLY == %s != 1)\n",homfly);

  free(homfly);
  pd_code_free(&pd);

  return true;

}

int main() {

  printf("test_bad_homfly (%s)\n",PACKAGE_STRING);
  printf("---------------------------------------\n"
	 "Unit test for pd_homfly \n"
	 "=======================================\n");

  if (!test_bad_homfly()) {

    printf("=====================================\n");
    printf("test_bad_homfly:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================\n");
    printf("test_bad_homfly:  PASS.\n");
    exit(0);

  }

  return 0;

}
