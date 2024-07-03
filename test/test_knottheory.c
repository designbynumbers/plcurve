/* 

   test_knottheory.c : Unit tests for reading and writing pd codes in KnotTheory format.

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

#ifdef HAVE_STDBOOL_H
   #include<stdbool.h>
#endif

#ifdef HAVE_ARGTABLE2_H
  #include<argtable2.h>
#endif

#include<assert.h>

#include<plcTopology.h>

#include<pd_multidx.h>
#include<pd_perm.h>
#include<pd_dihedral.h>

#include<pd_isomorphisms.h>



pd_code_t *pd_create_badguy_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(9);
assert(pd != NULL);
pd->ncross = 9;
pd->nedges = 18;
pd->ncomps = 1;
pd->nfaces = 11;

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 5;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 6;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 12;
pd->cross[1].edge[2] = 17;
pd->cross[1].edge[3] = 13;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 10;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 11;
pd->cross[2].sign = 0;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 15;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 16;
pd->cross[3].sign = 1;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 9;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 8;
pd->cross[4].sign = 0;

pd->cross[5].edge[0] = 4;
pd->cross[5].edge[1] = 14;
pd->cross[5].edge[2] = 5;
pd->cross[5].edge[3] = 13;
pd->cross[5].sign = 1;

pd->cross[6].edge[0] = 6;
pd->cross[6].edge[1] = 11;
pd->cross[6].edge[2] = 7;
pd->cross[6].edge[3] = 12;
pd->cross[6].sign = 0;

pd->cross[7].edge[0] = 7;
pd->cross[7].edge[1] = 16;
pd->cross[7].edge[2] = 8;
pd->cross[7].edge[3] = 17;
pd->cross[7].sign = 1;

pd->cross[8].edge[0] = 9;
pd->cross[8].edge[1] = 15;
pd->cross[8].edge[2] = 10;
pd->cross[8].edge[3] = 14;
pd->cross[8].sign = 1;


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
pd->edge[5].headpos = 1;
pd->edge[5].tail = 5;
pd->edge[5].tailpos = 2;

pd->edge[6].head = 6;
pd->edge[6].headpos = 0;
pd->edge[6].tail = 0;
pd->edge[6].tailpos = 3;

pd->edge[7].head = 7;
pd->edge[7].headpos = 0;
pd->edge[7].tail = 6;
pd->edge[7].tailpos = 2;

pd->edge[8].head = 4;
pd->edge[8].headpos = 3;
pd->edge[8].tail = 7;
pd->edge[8].tailpos = 2;

pd->edge[9].head = 8;
pd->edge[9].headpos = 0;
pd->edge[9].tail = 4;
pd->edge[9].tailpos = 1;

pd->edge[10].head = 2;
pd->edge[10].headpos = 1;
pd->edge[10].tail = 8;
pd->edge[10].tailpos = 2;

pd->edge[11].head = 6;
pd->edge[11].headpos = 1;
pd->edge[11].tail = 2;
pd->edge[11].tailpos = 3;

pd->edge[12].head = 1;
pd->edge[12].headpos = 1;
pd->edge[12].tail = 6;
pd->edge[12].tailpos = 3;

pd->edge[13].head = 5;
pd->edge[13].headpos = 3;
pd->edge[13].tail = 1;
pd->edge[13].tailpos = 3;

pd->edge[14].head = 8;
pd->edge[14].headpos = 3;
pd->edge[14].tail = 5;
pd->edge[14].tailpos = 1;

pd->edge[15].head = 3;
pd->edge[15].headpos = 1;
pd->edge[15].tail = 8;
pd->edge[15].tailpos = 1;

pd->edge[16].head = 7;
pd->edge[16].headpos = 1;
pd->edge[16].tail = 3;
pd->edge[16].tailpos = 3;

pd->edge[17].head = 1;
pd->edge[17].headpos = 2;
pd->edge[17].tail = 7;
pd->edge[17].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 18;
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
pd->comp[0].edge[16] = 16;
pd->comp[0].edge[17] = 17;


/* Face data */

pd->face[0].nedges = 4;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].orient = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].orient != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].orient[0] = 0;

pd->face[0].edge[1] = 5;
pd->face[0].orient[1] = 0;

pd->face[0].edge[2] = 14;
pd->face[0].orient[2] = 1;

pd->face[0].edge[3] = 10;
pd->face[0].orient[3] = 1;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].orient = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].orient != NULL);

pd->face[1].edge[0] = 2;
pd->face[1].orient[0] = 1;

pd->face[1].edge[1] = 16;
pd->face[1].orient[1] = 1;

pd->face[1].edge[2] = 7;
pd->face[1].orient[2] = 0;

pd->face[1].edge[3] = 11;
pd->face[1].orient[3] = 0;

pd->face[2].nedges = 4;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].orient = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].orient != NULL);

pd->face[2].edge[0] = 4;
pd->face[2].orient[0] = 1;

pd->face[2].edge[1] = 13;
pd->face[2].orient[1] = 0;

pd->face[2].edge[2] = 17;
pd->face[2].orient[2] = 0;

pd->face[2].edge[3] = 8;
pd->face[2].orient[3] = 1;

pd->face[3].nedges = 3;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].orient = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].orient != NULL);

pd->face[3].edge[0] = 0;
pd->face[3].orient[0] = 1;

pd->face[3].edge[1] = 6;
pd->face[3].orient[1] = 1;

pd->face[3].edge[2] = 12;
pd->face[3].orient[2] = 1;

pd->face[4].nedges = 3;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].orient = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].orient != NULL);

pd->face[4].edge[0] = 0;
pd->face[4].orient[0] = 0;

pd->face[4].edge[1] = 13;
pd->face[4].orient[1] = 1;

pd->face[4].edge[2] = 5;
pd->face[4].orient[2] = 1;

pd->face[5].nedges = 3;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].orient = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].orient != NULL);

pd->face[5].edge[0] = 1;
pd->face[5].orient[0] = 1;

pd->face[5].edge[1] = 11;
pd->face[5].orient[1] = 1;

pd->face[5].edge[2] = 6;
pd->face[5].orient[2] = 0;

pd->face[6].nedges = 3;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].orient = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].orient != NULL);

pd->face[6].edge[0] = 2;
pd->face[6].orient[0] = 0;

pd->face[6].edge[1] = 10;
pd->face[6].orient[1] = 0;

pd->face[6].edge[2] = 15;
pd->face[6].orient[2] = 1;

pd->face[7].nedges = 3;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].orient = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].orient != NULL);

pd->face[7].edge[0] = 3;
pd->face[7].orient[0] = 1;

pd->face[7].edge[1] = 8;
pd->face[7].orient[1] = 0;

pd->face[7].edge[2] = 16;
pd->face[7].orient[2] = 0;

pd->face[8].nedges = 3;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].orient = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].orient != NULL);

pd->face[8].edge[0] = 3;
pd->face[8].orient[0] = 0;

pd->face[8].edge[1] = 15;
pd->face[8].orient[1] = 0;

pd->face[8].edge[2] = 9;
pd->face[8].orient[2] = 0;

pd->face[9].nedges = 3;
pd->face[9].edge = calloc(pd->face[9].nedges,sizeof(pd_idx_t));
pd->face[9].orient = calloc(pd->face[9].nedges,sizeof(pd_or_t));
assert(pd->face[9].edge != NULL);
assert(pd->face[9].orient != NULL);

pd->face[9].edge[0] = 4;
pd->face[9].orient[0] = 0;

pd->face[9].edge[1] = 9;
pd->face[9].orient[1] = 1;

pd->face[9].edge[2] = 14;
pd->face[9].orient[2] = 0;

pd->face[10].nedges = 3;
pd->face[10].edge = calloc(pd->face[10].nedges,sizeof(pd_idx_t));
pd->face[10].orient = calloc(pd->face[10].nedges,sizeof(pd_or_t));
assert(pd->face[10].edge != NULL);
assert(pd->face[10].orient != NULL);

pd->face[10].edge[0] = 7;
pd->face[10].orient[0] = 1;

pd->face[10].edge[1] = 17;
pd->face[10].orient[1] = 1;

pd->face[10].edge[2] = 12;
pd->face[10].orient[2] = 0;


/* End of data. */

pd_regenerate_hash(pd);
assert(pd_ok(pd));
return pd;

}

bool checkpd(pd_code_t *pd,char *name)
{
  char tmpName[256] = "/tmp/ktfileXXXXXX";
  FILE *tmpfile;
  int tmpfd;
  
  printf("writing pd code %s to knottheory format...",name);
  tmpfd = mkstemp(tmpName);
  if (tmpfd == -1) {
    printf("FAIL (couldn't create tempfile)\n");
    return false;
  }
  tmpfile = fdopen(tmpfd,"w");
  if (tmpfile == NULL) {
    printf("FAIL (got fd for tempfile, but couldn't open for writing)\n");
    return false;
  }

  pd_write_KnotTheory(tmpfile,pd);
  fclose(tmpfile);

  printf("done\n");

  printf("reopening temp file %s...",tmpName);
  tmpfile = fopen(tmpName,"r");
  if (tmpfile == NULL) {
    printf("FAIL (couldn't open file)\n");
    return false;
  }

  printf("reading KnotTheory format pdcode...");
  pd_code_t *newpd = pd_read_KnotTheory(tmpfile);
  if (newpd == NULL) {
    printf("FAIL (couldn't parse file)\n");
    return false;
  }
  if (!pd_ok(newpd)) {
    pd_printf("FAIL (returned pd %PD is not pd_ok)\n",pd);
    return false;
  }
  printf("pass\n");

  printf("comparing old and new pdcodes...");
  if (!pd_diagram_isotopic(pd,newpd)) {
    printf("FAIL (old and new pd are not diagram isotopic)\n\n");
    pd_printf("OLD pd\n===================\n%PD",pd);
    pd_printf("NEW pd\n===================\n%PD",newpd);
    return false;
  }
  printf("pass (pd codes are diagram-isotopic)\n");

  pd_code_free(&newpd);
  fclose(tmpfile);

  return true;
}

bool test_inouts() {

  printf("test writing out and reading in pd codes\n"
	 "----------------------------------------\n");

  pd_code_t *pd = pd_create_badguy_0();
  if (!checkpd(pd,"ambig 9 test case")) {
    printf("--------------------------------------\n"
	   "test writing out and reading in: FAIL\n");
    return false;
  }

  printf("-------------------------------------\n"
	 "test writing out and reading in: pass\n\n");
  return true;
}

int main() {

  printf("test_knottheory (%s)\n",PACKAGE_STRING);
  printf("---------------------------------------\n"
	 "Test read and write from KnotTheory.   \n"
	 "=======================================\n");

  if (!test_inouts()) {

    printf("=====================================\n");
    printf("test_knottheory:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================\n");
    printf("test_knottheory:  PASS.\n");
    exit(0);

  }

  return 0;

}
