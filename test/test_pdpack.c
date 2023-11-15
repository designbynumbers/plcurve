/* 

   test_pdpack.c : Unit tests for the code in pd_pack.c


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

#ifdef HAVE_ASSERT_H
   #include<assert.h>
#endif

#include<ordie.h>
#include<plcTopology.h>

#include<pd_multidx.h>
#include<pd_perm.h>
#include<pd_dihedral.h>

#include<pd_isomorphisms.h>
#include<pd_storage.h>
#include<pd_pack.h>

#include"argtable2.h" /* We use a local copy of argtable */



bool pd_identical(pd_code_t *pdA,pd_code_t *pdB)
/* This is function that should have basically no use outside debugging,
   but it compares, number for number, the information in the pd codes */
{
  int i,j;

  if (strcmp(pdA->hash,pdB->hash) != 0) {

    return false;

  }

  if (pdA->ncross != pdB->ncross) {

    return false;

  }

  for(i=0;i<pdA->ncross;i++) {

    for(j=0;j<4;j++) {
      
      if (pdA->cross[i].edge[j] != pdB->cross[i].edge[j]) {

	return false;

      }

    }

  }

  if (pdA->nedges != pdB->nedges) {

    return false;

  }

  for(i=0;i<pdA->nedges;i++) {

    if (pdA->edge[i].head != pdB->edge[i].head) {

      return false;

    }

    if (pdA->edge[i].headpos != pdB->edge[i].headpos) {

      return false;

    }

    if (pdA->edge[i].tail != pdB->edge[i].tail) {

      return false;

    }

    if (pdA->edge[i].tailpos != pdB->edge[i].tailpos) {

      return false;

    }

  }

  if (pdA->nfaces != pdB->nfaces) {

    return false;

  }

  for(i=0;i<pdA->nfaces;i++) {

    for(j=0;j<pdA->face[i].nedges;j++) {

      if (pdA->face[i].edge[j] != pdB->face[i].edge[j]) {

	return false;

      }

      if (pdA->face[i].orient[j] != pdB->face[i].orient[j]) {

	return false;

      }

    }

  }

  if (pdA->ncomps != pdB->ncomps) {

    return false;

  }

  for(i=0;i<pdA->ncomps;i++) {

    for(j=0;j<pdA->comp[i].nedges;j++) {

      if (pdA->comp[i].edge[j] != pdB->comp[i].edge[j]) {

	return false;

      }

    }

  }

  return true;

}
  

bool pack_torusknots(int nknots)

/* pack and unpack torus knots */

{
  int i;

  printf("packing and unpacking %d (distinct) torus knots...\n",
	 nknots);

  for(i=0;i<nknots;i++) {

    printf("\tbuilding (%d,%d) torus knot...",2,i+3);
    pd_code_t *pd;
    pd = pd_build_torus_knot(2,i+3);
    printf("done\n");

    unsigned int *packform;
    unsigned int packlength;
    pd_code_t *newpd;

    packform = pd_pack(pd,&packlength);

    printf("\tpacked %d byte pdcode to %d byte bitstring (%g%% compression)...done\n",
	   (int)(pd_size(pd)),(int)(sizeof(unsigned int))*packlength,
	   100.0 - 100.0*(double)(sizeof(unsigned int)*packlength)/(double)(pd_size(pd)));

    printf("\tunpacking to new pd code...");
    newpd = pd_unpack(packform);
    printf("done\n");

    printf("\tchecking that new and old pd are isomorphic...");
    if (!pd_isomorphic(newpd,pd)) {
      printf("FAIL (not isomorphic)\n");
      return false;
    }
    printf("pass\n");

    printf("\tchecking that new and old pd are diagram-isotopic...");
    if (!pd_diagram_isotopic(newpd,pd)) {
      printf("FAIL (not diagram-isotopic)\n");
      return false;
    }
    printf("pass\n");


    printf("\tchecking that new and old pd are identical...");
    if (!pd_identical(newpd,pd)) {
      printf("FAIL (not identical)\n");
      return false;
    }
    printf("pass\n");

    printf("\thousekeeping...");
    pd_code_free(&pd);
    pd_code_free(&newpd); free(packform);
    printf("done\n\n");
    
  }
    
  printf("pass.\n"); /* If we didn't crash, we pass. */

  return true;
  
}

bool pack_torusknots_unsigned(int nknots)

/* pack and unpack torus knots */

{
  int i;

  printf("packing and unpacking %d (distinct) torus knot SHADOWS...\n",
	 nknots);

  for(i=0;i<nknots;i++) {

    printf("\tbuilding (%d,%d) torus knot...",2,i+3);
    pd_code_t *pd;
    pd = pd_build_torus_knot(2,i+3);
    printf("done\n");

    int i;
    for(i=0;i<pd->ncross;i++) { pd->cross[i].sign = PD_UNSET_ORIENTATION; }

    unsigned int *packform;
    unsigned int packlength;
    pd_code_t *newpd;

    packform = pd_pack(pd,&packlength);

    printf("\tpacked %d byte pdcode to %d byte bitstring (%g%% compression)...done\n",
	   (int)(pd_size(pd)),(int)(sizeof(unsigned int))*packlength,
	   100.0 - 100.0*(double)(sizeof(unsigned int)*packlength)/(double)(pd_size(pd)));

    printf("\tunpacking to new pd code...");
    newpd = pd_unpack(packform);
    printf("done\n");

    printf("\tchecking that new and old pd are isomorphic...");
    if (!pd_isomorphic(newpd,pd)) {
      printf("FAIL (not isomorphic)\n");
      return false;
    }
    printf("pass\n");

    printf("\tchecking that new and old pd are identical...");
    if (!pd_identical(newpd,pd)) {
      printf("FAIL (not identical)\n");
      return false;
    }
    printf("pass\n");

    printf("\thousekeeping...");
    pd_code_free(&pd);
    pd_code_free(&newpd);
    free(packform);
    printf("done\n\n");
    
  }
    
  printf("pass.\n"); /* If we didn't crash, we pass. */

  return true;
  
}


struct sum_triple {

    pd_idx_t triple[3];

};

int cmp_pdx(const void *a,const void *b)
{
  pd_idx_t *A = (pd_idx_t *)(a);
  pd_idx_t *B = (pd_idx_t *)(b);

  return *A - *B;
}

int cmp_sum_triple(const void *a,const void *b)
{
  struct sum_triple *A = (struct sum_triple *)(a);
  struct sum_triple *B = (struct sum_triple *)(b);

  pd_idx_t i;

  for(i=0;i<3;i++) {

    if (A->triple[i] != B->triple[i]) { return A->triple[i] - B->triple[i]; }

  }

  return 0;

}


unsigned int count_unknotwye_isomorphism_types(int n)

/* # isomorphism types == # of (unordered) triples which sum to ntwists */

/* 
   This is best handled as a direct calculation. 

   Observe that given any such triple, we can write it in
   the form a + b + c = n.

   Regrouping, c = n - (a + b) is completely determined by
   a and b, so the problem is finding pairs which sum to
   each number between 0 and n (inclusive).

   Now if a + b = m (<= n), then a is between 0 and m and 
   b = m - a is completely determined by a. We can thus
   generate all such triples as follows:

*/
{
  struct sum_triple *stbuf;

  stbuf = calloc((n+1)*(n+1),sizeof(struct sum_triple));
  assert(stbuf != NULL);

  pd_idx_t a,m,st;

  for(st=0,m=0;m<=n;m++) {

    for(a=0;a<=m;a++,st++) { 

      stbuf[st].triple[0] = a;
      stbuf[st].triple[1] = m-a;
      stbuf[st].triple[2] = n-m;
      
    }
 
  }
   
  pd_idx_t this_st;
  
  for(this_st=0;this_st < st;this_st++) {

    qsort(stbuf[this_st].triple,3,sizeof(pd_idx_t),cmp_pdx);
    
  }

  qsort(stbuf,st,sizeof(struct sum_triple),cmp_sum_triple);
   
  /* Now we can count the distinct elements in the (sorted) buffer. */

  pd_idx_t unique_sts;

  for(this_st=1,unique_sts=1;
      this_st < st;
      this_st++) {

    if (stbuf[this_st].triple[0] != stbuf[this_st-1].triple[0] ||
	stbuf[this_st].triple[1] != stbuf[this_st-1].triple[1] ||
	stbuf[this_st].triple[2] != stbuf[this_st-1].triple[2]) {

      unique_sts++;

    }

  }
   
  /* Finally, free the buffer and return. */

  free(stbuf);
  return unique_sts;

}

bool pack_unknotwyes(int n)

/* Pack and unpack unknotted wyes with n twists. */

{
  pd_idx_t a,m;

  printf("packing/unpacking all unknotted wyes with %d twists...\n",
	 n);

  for(m=0;m<=n;m++) {

    for(a=0;a<=m;a++) {

      printf("\tbuilding %d-%d-%d unknotted wye...",a,m-a,n-m);
      pd_code_t *pd;
      pd = pd_build_unknot_wye(a,m-a,n-m);
      printf("done\n");
      
      unsigned int *packform;
      unsigned int packlength;
      pd_code_t *newpd;
      
      packform = pd_pack(pd,&packlength);
      
      printf("\tpacked %d byte pdcode to %d byte bitstring (%g%% compression)...done\n",
	     (int)(pd_size(pd)),(int)(sizeof(unsigned int))*packlength,
	     100.0 - 100.0*(double)(sizeof(unsigned int)*packlength)/(double)(pd_size(pd)));
      
      printf("\tunpacking to new pd code...");
      newpd = pd_unpack(packform);
      printf("done\n");
      
      printf("\tchecking that new and old pd are isomorphic...");
      if (!pd_isomorphic(newpd,pd)) {
	printf("FAIL (not isomorphic)\n");
	return false;
      }
      printf("pass\n");

      printf("\tchecking that new and old pd are diagram-isotopic...");
      if (!pd_diagram_isotopic(newpd,pd)) {
	printf("FAIL (not diagram-isotopic)\n");
	return false;
      }
      printf("pass\n");

      
      printf("\tchecking that new and old pd are identical...");
      if (!pd_identical(newpd,pd)) {
	printf("FAIL (not identical)\n");
	return false;
      }
      printf("pass\n");
      
      printf("\thousekeeping...");
      pd_code_free(&pd);
      pd_code_free(&newpd);
      free(packform);
      printf("done\n\n");
      
    }

  }

  printf("pass.\n"); /* If we didn't crash, we pass. */
  return true;
  
}

bool pack_unknotwyes_unsigned(int n)

/* Pack and unpack unknotted wyes with n twists. */

{
  pd_idx_t a,m;

  printf("packing/unpacking all unknotted wyes with %d twists...\n",
	 n);

  for(m=0;m<=n;m++) {

    for(a=0;a<=m;a++) {

      printf("\tbuilding %d-%d-%d unknotted wye (crossings unset)...",a,m-a,n-m);
      pd_code_t *pd;
      pd = pd_build_unknot_wye(a,m-a,n-m);
      printf("done\n");

      int i;
      for(i=0;i<pd->ncross;i++) { pd->cross[i].sign = PD_UNSET_ORIENTATION; }
      
      unsigned int *packform;
      unsigned int packlength;
      pd_code_t *newpd;
      
      packform = pd_pack(pd,&packlength);
      
      printf("\tpacked %d byte pdcode to %d byte bitstring (%g%% compression)...done\n",
	     (int)(pd_size(pd)),(int)(sizeof(unsigned int))*packlength,
	     100.0 - 100.0*(double)(sizeof(unsigned int)*packlength)/(double)(pd_size(pd)));
      
      printf("\tunpacking to new pd code...");
      newpd = pd_unpack(packform);
      printf("done\n");
      
      printf("\tchecking that new and old pd are isomorphic...");
      if (!pd_isomorphic(newpd,pd)) {
	printf("FAIL (not isomorphic)\n");
	return false;
      }
      printf("pass\n");
      
      printf("\tchecking that new and old pd are identical...");
      if (!pd_identical(newpd,pd)) {
	printf("FAIL (not identical)\n");
	return false;
      }
      printf("pass\n");
      
      printf("\thousekeeping...");
      pd_code_free(&pd); pd_code_free(&newpd);
      free(packform);
      printf("done\n\n");
      
    }

  }

  printf("pass.\n"); /* If we didn't crash, we pass. */
  return true;
  
}


bool test_pdpack() {

  printf("--------------------------------\n"
	 "pd_pack/pd_unpack test suite\n"
	 "--------------------------------\n");

  if (!pack_torusknots(30)) { return false; }

  if (!pack_torusknots_unsigned(30)) { return false; }

  int i;
  for(i=0;i<14;i++) { 

    if (!pack_unknotwyes(i)) { return false; }

  }

  for(i=0;i<14;i++) { 

    if (!pack_unknotwyes_unsigned(i)) { return false; }

  }
 
  printf("-----------------------------------\n"
	 "pd_pack/pd_unpack test suite:  PASS\n"
	 "-----------------------------------\n\n");

  return true;

}

bool test_random_bitstring(unsigned int len)
{
  int i;
  
  //printf("testing random bitstring of length %d...",len);
  
  unsigned int *t = calloc(len,sizeof(unsigned int));
  for(i=0;i<len;i++) { t[i] = rand() % 2; }
  unsigned int *packed, packints;

  packed = pd_packbitstring(t,len,&packints);
  unsigned int *unpacked = pd_unpackbitstring(packed,len);
  
  for(i=0;i<len;i++) {
    
    if (t[i] != unpacked[i]) {
      
      printf("FAIL (length %d bitstring different at position %d, "
	     "where test[%d] = %d and unpacked[%d] = %d)\n",
	     len,i,i,t[i],i,unpacked[i]);
      return false;
      
    }
    
  }
  
  //  printf("pass\n");
  free(packed);
  free(unpacked);
  free(t);

  return true;
}

bool test_bitstring() {

  printf("----------------------------------------------\n"
	 "pd_packbitstring/pd_unpackbitstring test suite\n"
	 "----------------------------------------------\n");

  printf("packing bitstring 100110101010...");
  unsigned int t1[12] = {1,0,0,1,1,0,1,0,1,0,1,0};
  unsigned int *packed, packints;

  packed = pd_packbitstring(t1,12,&packints);
  printf("done (size %d bytes)\n",packints*(int)(sizeof(unsigned int)));

  printf("unpacking bitstring...");
  unsigned int *unpacked = pd_unpackbitstring(packed,12);
  printf("done\n");

  printf("comparing strings...");

  int i;
  for(i=0;i<12;i++) {

    if (t1[i] != unpacked[i]) {

      printf("FAIL (different at position %d, where test[%d] = %d and unpacked[%d] = %d)\n",
	     i,i,t1[i],i,unpacked[i]);
      return false;

    }

  }

  printf("pass (strings match)\n");
  free(packed);
  free(unpacked);

  printf("testing random bitstrings of length 1..257 ...");
  
  for(i=1;i<257;i++) {

    if (!test_random_bitstring(i)) { return false; }

  }

  printf("pass (packed and unpacked correctly)\n");
  
  printf("----------------------------------------------\n"
	 "pd_packbitstring/pd_unpackbitstring: PASS     \n"
	 "----------------------------------------------\n\n");

  return true;

}

bool test_random_pdx()
{
  pd_idx_t idx = (pd_idx_t)(rand())%4000;
  
  unsigned int *packed = pd_idx_bitencode(idx,2000);
  pd_idx_t unpacked = pd_idx_bitdecode(packed,2000);

  if (idx != unpacked) {

    printf("FAIL (packed %d, but unpacked %d)\n",(int)(idx),(int)(unpacked));
    return false;

  }

  free(packed);
  return true;
}

bool test_pdx() {

  printf("----------------------------------------------\n"
	 "pd_idx_bitencode/pd_idx_bitdecode test suite  \n"
	 "----------------------------------------------\n");

  printf("packing/unpacking numbers 0..9 with ncross=5 ...");
  int i;

  for(i=0;i<10;i++) {

    unsigned int *packed = pd_idx_bitencode((pd_idx_t)(i),5);
    pd_idx_t unpacked = pd_idx_bitdecode(packed,5);

    if ((pd_idx_t)(i) != unpacked) {

      printf("FAIL (packed %d, unpacked %d)\n",(int)((pd_idx_t)(i)),(int)(unpacked));
      return false;

    }

    free(packed);

  }

  printf("pass\n");
  
  
  printf("packing/unpacking numbers 0..255 with ncross=128...");

  for(i=0;i<256;i++) {

    unsigned int *packed = pd_idx_bitencode((pd_idx_t)(i),128);
    pd_idx_t unpacked = pd_idx_bitdecode(packed,128);

    if ((pd_idx_t)(i) != unpacked) {

      printf("FAIL (packed %d, unpacked %d)\n",(int)((pd_idx_t)(i)),(int)(unpacked));
      return false;

    }

    free(packed);

  }

  printf("pass\n");

  printf("packing/unpacking 5,000 random numbers ...");
  for(i=0;i<5000;i++) {

    if (!test_random_pdx()) { return false; }

  }

  printf("pass\n");
  
  printf("----------------------------------------------\n"
	 "pd_idx_bitencode/pd_idx_bitdecode: PASS       \n"
	 "----------------------------------------------\n");

  return true;

}

bool pd_stor_test(unsigned int n)
/* Load, pack, and unpack everything in a pdstor... */
{
  printf("----------------------------------------------\n"
	 "pd_stor %d.pdstor test\n"
	 "----------------------------------------------\n",n);
  
  char pdname[257];
  sprintf(pdname,"%d.pdstor",n);

  printf("looking for %s...",pdname);
  FILE *pdstor = fopen(pdname,"r");

  if (pdstor == NULL) {

    sprintf(pdname,"../data/pdstors/%d.pdstor",n);
    pdstor = fopen(pdname,"r");

    if (pdstor == NULL ) {
    
      printf("abort (couldn't find test file)\n");
      return true;

    }

  }

  printf("done\n");
  printf("parsing header...");
  
  unsigned int nelts_claimed,nelts_actual,nhashes;

  if (fscanf(pdstor,"pdstor nelts %d/%d (claimed/actual) nhashes %d ",
	     &nelts_claimed, &nelts_actual, &nhashes) != 3) {

    printf("FAIL (file must be in unfamiliar format)\n");
    return false;

  }

  printf("pass (%d pdcodes/%d hashes)\n",nelts_claimed,nhashes);

  printf("reading, packing, and unpacking pd codes...");

  int i;
  
  for(i=0;i<nelts_claimed;i++) {

    pd_code_t *pd = pd_read(pdstor);
    if (pd == NULL) {

      printf("FAIL (couldn't read code %d)\n",i);
      return false;
      
    }
    
    unsigned int *packform;
    unsigned int packlength;
    pd_code_t *newpd;
    
    packform = pd_pack(pd,&packlength);
    newpd = pd_unpack(packform);
    
    if (!pd_isomorphic(newpd,pd)) {
      printf("FAIL (code %d packed/unpacked not isomorphic)\n",i);
      return false;
    }

    if (!pd_identical(newpd,pd)) {
      printf("FAIL (code %d packed/unpacked not identical)\n",i);
      return false;
    }


    pd_code_free(&pd);
    pd_code_free(&newpd); free(packform);

  }

  printf("pass\n");
  fclose(pdstor);

  return true;
  
}

bool test_pdstors(unsigned int maxn)
{
  int i;

  for(i=3;i<maxn;i++) {

    if (!pd_stor_test(i)) { return false; }

  }

  return true;

}

  
int main() {

  printf("test_pdpack (%s)\n",PACKAGE_STRING);
  printf("---------------------------------------\n"
	 "Unit tests for pd_pack.c\n"
	 "=======================================\n");

  if (!test_bitstring() || !test_pdx() || !test_pdpack() || !test_pdstors(7)) {

    printf("=====================================\n");
    printf("test_pdpack:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================\n");
    printf("test_pdpack:  PASS.\n");
    exit(0);

  }

  return 0;

}
