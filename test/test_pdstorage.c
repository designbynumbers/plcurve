/* 

   test_pdstorage.c : Unit tests for the code in pd_storage.c.


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

#include<assert.h>

#include<ordie.h>
#include<plcTopology.h>

#include<pd_multidx.h>
#include<pd_perm.h>
#include<pd_dihedral.h>

#include<pd_isomorphisms.h>

#include<pd_storage.h>

#include"argtable2.h" /* We use a local copy of argtable */



bool insert_torusknots(int nknots,pd_stor_t **pdstor,pd_equivalence_t eq)

/* Insert nknots (distinct) torusknots into the pdstor, using equivalence */
/* type eq (NONE, ISOMORPHISM, or DIAGRAM_ISOTOPY) */
  
{
  int i;

  (*pdstor) = pd_new_pdstor();

  printf("\ninserting %d (distinct) torus knots in pdstor with %d elements\n"
	 "equivalence relation is %s...",
	 nknots,pd_stor_nelts(*pdstor),
	 eq == NONE ? "none" : (eq == ISOMORPHISM ? "isomorphism" : "diagram isotopy"));

  if (PD_VERBOSE) { printf("\n"); }

  for(i=0;i<nknots;i++) {

    pd_code_t *pd;
    pd = pd_build_torus_knot(2,i+3);
    pd_addto_pdstor(*pdstor,pd,eq); /* Do the insert. */
    pd_code_free(&pd);

  }

  printf("pass.\n"); /* If we didn't crash, we pass. */
  printf("checking stats...");

  unsigned int nhashes, nelts;
  pd_stor_stats(*pdstor,&nhashes,&nelts);
  
  if (nhashes == nknots && nelts == nknots) {

    printf("pass (nhashes %u == nelts %u == nknots %d).\n",nhashes,nelts,nknots);

  } else {

    printf("FAIL (nhashes %u, nelts %u, nknots %d).\n",nhashes,nelts,nknots);
    return false;

  }

  if (nknots < 20) { 

    printf("writing pdstor contents...\n");
    pd_display_pdstor(stdout,*pdstor);
    printf("done.\n\n");

  }

  return true;
  
}

bool insert_torusknots_with_dups(int nknots,pd_stor_t **pdstor,pd_equivalence_t eq)

/* Insert nknots (distinct) torusknots into the pdstor, 
   along with (many) duplicate entries. */

{
  int i;

  (*pdstor) = pd_new_pdstor();

  printf("\ninserting %d torus knots w/dups in pdstor with %d elements\n"
	 "equivalence relation is %s...",
	 nknots,pd_stor_nelts(*pdstor),
	 eq == NONE ? "none" : (eq == ISOMORPHISM ? "isomorphism" : "diagram isotopy"));

  if (PD_VERBOSE) { printf("\n"); }

  int added = 0;

  for(i=0;i<nknots;i++) {

    pd_code_t *pd;

    pd = pd_build_torus_knot(2,i+3);
    pd_addto_pdstor(*pdstor,pd,eq); /* Do the insert. */
    pd_code_free(&pd);
    added++;

    pd = pd_build_torus_knot(2,i+3);
    pd_addto_pdstor(*pdstor,pd,eq); /* Do a repeat insert. */
    pd_code_free(&pd);
    added++;

  }

  for(i=nknots-1;i>=0;i--) { /* Try again in backwards order. */

    pd_code_t *pd;

    pd = pd_build_torus_knot(2,i+3);
    pd_addto_pdstor(*pdstor,pd,eq); /* Do the insert. */
    pd_code_free(&pd);
    added++;
    
    pd = pd_build_torus_knot(2,i+3);
    pd_addto_pdstor(*pdstor,pd,eq); /* Do a repeat insert. */
    pd_code_free(&pd);
    added++;

  }

  printf("pass.\n"); /* If we didn't crash, we pass. */
  printf("checking stats...");

  unsigned int nhashes, nelts;
  pd_stor_stats(*pdstor,&nhashes,&nelts);

  if (eq == ISOMORPHISM || eq == DIAGRAM_ISOTOPY) { 
  
    if (nhashes == nknots && nelts == nknots) {
      
      printf("pass (nhashes %u == nelts %u == nknots %d).\n",nhashes,nelts,nknots);
      
    } else {
      
      printf("FAIL (nhashes %u, nelts %u, nknots %d).\n",nhashes,nelts,nknots);
      return false;
      
    }

  } else if (eq == NONE) {

    if (nhashes != nknots) {

      printf("fail (number of hashes %d != number of knot types %d)\n",nhashes,nknots);
      return false;

    }

    if (nelts != added) {

      printf("fail (added %d knots using no equivalence relation, but have %d in pdstor)\n",
	     added,nelts);
      return false;

    }

  }
	    

  if (nknots < 20) { 

    printf("writing pdstor contents...\n");
    pd_display_pdstor(stdout,*pdstor);
    printf("done.\n\n");

  }

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

unsigned int count_unknotwye_diagram_isotopy_types(int n)

/* # diagram_isotopy types == # of (cyclically ordered) triples which sum to ntwists */

/* 
   This is a Mathematica calculation; see diagram_wye_isotopy_types.nb (this directory)

*/
{
  int dwye_isotopy[16] = {1,1,2,4,5,7,10,12,15,19,22,26,31,35,40,46};

  assert(n <= 15);
  return dwye_isotopy[n];  
}

bool insert_unknotwyes(int n,pd_stor_t **pdstor,pd_equivalence_t eq)

/* Insert all of the unknotted wyes with n twists into the pdstor. */
/* This is designed to generate some hash collisions. */

{
  pd_idx_t a,m;
  int added = 0;

  (*pdstor) = pd_new_pdstor();

  printf("\ninserting all unknotted wyes with %d twists in pdstor with %d elements\n"
	 "equivalence relation is %s...",
	 n,pd_stor_nelts(*pdstor),
	 eq == NONE ? "none" : (eq == ISOMORPHISM ? "isomorphism" : "diagram isotopy"));

  if (PD_VERBOSE) { printf("\n"); }
  if (n < 4) { printf("\n"); }

  for(m=0;m<=n;m++) {

    for(a=0;a<=m;a++) {

      pd_code_t *pd;
      pd = pd_build_unknot_wye(a,m-a,n-m);
      pd_addto_pdstor(*pdstor,pd,eq); /* Do the insert. */
      pd_code_free(&pd);
      added++;

      if (n < 4) {

	printf("\t added (%d-%d-%d), nelts = %d\n",a,m-a,n-m,pd_stor_nelts(*pdstor));

      }	

    }

  }

  printf("pass.\n"); /* If we didn't crash, we pass. */
  printf("checking stats...");

  unsigned int nhashes, nelts;
  pd_stor_stats(*pdstor,&nhashes,&nelts);

  unsigned int expected_nelts;

  if (eq == ISOMORPHISM) {

    expected_nelts = count_unknotwye_isomorphism_types(n);

  } else {

    expected_nelts = count_unknotwye_diagram_isotopy_types(n);

  }

  if (eq == ISOMORPHISM || eq == DIAGRAM_ISOTOPY) { 
  
    if (nelts == expected_nelts) {
      
      printf("pass (nhashes %u, nelts %u == expected_nelts %u).\n",nhashes,nelts,expected_nelts);
      
    } else {
      
      printf("FAIL (nhashes %u, nelts %u != expected_nelts %u).\n",nhashes,nelts,expected_nelts);
      return false;

    }

  } else if (eq == NONE) {

    if (nelts != added) {

      printf("FAIL (no eq relation, but nelts = %d != added = %d)\n",nelts,added);
      return false;

    }

    if (nhashes > expected_nelts) {

      printf("FAIL (have more hashes (%d) than isomorphism types (%d))\n",
	     nhashes,expected_nelts);
      return false;

    }

  }

  if (nelts < 20) { 

    printf("writing pdstor contents...\n");
    pd_display_pdstor(stdout,*pdstor);
    printf("done.\n\n");

  }

  return true;
  
}

    
bool free_test(pd_stor_t **pdstor) 
{
  printf("deleting all memory associated with %d element pdstor...",pd_stor_nelts(*pdstor));

  pd_free_pdstor(pdstor);

  if (*pdstor != NULL) { 

    printf("FAIL. (didn't free pointer)\n");
    return false;

  }

  printf("pass.\n");
  printf("(should be run under valgrind to check that no memory has been lost)\n");

  return true;

}

bool test_copyinto() {

  pd_stor_t *pdstor;

  printf("pd_copyinto_pdstor test suite\n"
	 "--------------------------------\n");

  if (!insert_torusknots(10,&pdstor,NONE)) { return false; }
  if (!free_test(&pdstor)) { return false; }

  if (!insert_torusknots_with_dups(10,&pdstor,NONE)) { return false; }
  if (!free_test(&pdstor)) { return false; }

  if (!insert_unknotwyes(1,&pdstor,NONE)) { return false; }
  if (!free_test(&pdstor)) { return false; }

  if (!insert_unknotwyes(2,&pdstor,NONE)) { return false; }
  if (!free_test(&pdstor)) { return false; }
  
  if (!insert_unknotwyes(3,&pdstor,NONE)) { return false; }
  if (!free_test(&pdstor)) { return false; }

  if (!insert_unknotwyes(10,&pdstor,NONE)) { return false; }
  if (!free_test(&pdstor)) { return false; }

  if (!insert_torusknots(10,&pdstor,NONE)) { return false; }
  if (!free_test(&pdstor)) { return false; }

  /** */
  
  if (!insert_torusknots_with_dups(10,&pdstor,ISOMORPHISM)) { return false; }
  if (!free_test(&pdstor)) { return false; }

  if (!insert_unknotwyes(1,&pdstor,ISOMORPHISM)) { return false; }
  if (!free_test(&pdstor)) { return false; }

  if (!insert_unknotwyes(2,&pdstor,ISOMORPHISM)) { return false; }
  if (!free_test(&pdstor)) { return false; }
  
  if (!insert_unknotwyes(3,&pdstor,ISOMORPHISM)) { return false; }
  if (!free_test(&pdstor)) { return false; }

  if (!insert_unknotwyes(10,&pdstor,ISOMORPHISM)) { return false; }
  if (!free_test(&pdstor)) { return false; }

  if (!insert_torusknots(10,&pdstor,ISOMORPHISM)) { return false; }
  if (!free_test(&pdstor)) { return false; }

  /*** */

  if (!insert_torusknots_with_dups(10,&pdstor,DIAGRAM_ISOTOPY)) { return false; }
  if (!free_test(&pdstor)) { return false; }

  if (!insert_unknotwyes(1,&pdstor,DIAGRAM_ISOTOPY)) { return false; }
  if (!free_test(&pdstor)) { return false; }

  if (!insert_unknotwyes(2,&pdstor,DIAGRAM_ISOTOPY)) { return false; }
  if (!free_test(&pdstor)) { return false; }
  
  if (!insert_unknotwyes(3,&pdstor,DIAGRAM_ISOTOPY)) { return false; }
  if (!free_test(&pdstor)) { return false; }

  if (!insert_unknotwyes(10,&pdstor,DIAGRAM_ISOTOPY)) { return false; }
  if (!free_test(&pdstor)) { return false; }
 
  printf("-----------------------------------\n"
	 "pd_copyinto_pdstor test suite        PASS\n");

  return true;

}

bool should_fail_iso_search(pd_code_t *pd,pd_stor_t *pdstor) 

{
  pd_code_t   *retpd;
  pd_iso_t   **isos = (pd_iso_t **)(1234);
  unsigned int nisos = 54;
 
  retpd = pd_search_pdstor_by_isomorphism(pdstor,pd,&isos,&nisos);
  
  if (retpd != NULL) {
    
    printf("FAIL (failed search didn't return NULL).\n");
    return false;
    
  }
  
  if (isos != NULL) {
    
    printf("FAIL (failed search didn't clear isos).\n");
    return false;
    
  }
  
  if (nisos != 0) { 
    
    printf("FAIL (failed search didn't set nisos to zero).\n");
    return false;
    
  }
  
  return true;
  
}

bool should_succeed_iso_search(pd_code_t *pd,pd_stor_t *pdstor) 

{
  pd_code_t   *retpd;
  pd_iso_t   **isos = (pd_iso_t **)(1234);
  unsigned int nisos = 54;

  pd_iso_t   **check_isos;
  unsigned int check_nisos;

  /* First, figure out which automorphisms should be present. */

  check_isos = pd_build_isos(pd,pd,&check_nisos);

  /* Now do the actual search */

  retpd = pd_search_pdstor_by_isomorphism(pdstor,pd,&isos,&nisos);
  
  if (retpd == NULL) {
    
    printf("FAIL (expected to find pd, but didn't).\n");
    return false;
    
  }

  if (!pd_ok(retpd)) { 

    printf("FAIL (returned pd not ok).\n");
    return false;

  }

  if (!pd_isomorphic(retpd,pd)) { 

    printf("FAIL (returned pd not actually isomorphic to original).\n");
    return false;

  }
  
  if (isos == NULL) {
    
    printf("FAIL (didn't allocate any isos).\n");
    return false;
    
  }
  
  if (nisos != check_nisos) { 
    
    printf("FAIL. (should have found %d isos, found %d).\n",check_nisos,nisos);
    return false;
    
  }

  int i;

  for(i=0;i<nisos;i++) {

    if (pd_iso_cmp(&(isos[i]),&(check_isos[i])) != 0) { 

      printf("FAIL. (iso %d doesn't match check_iso %d).\n",i,i);
      return false;

    }

  }

  pd_free_isos(&nisos,&isos);
  pd_free_isos(&check_nisos,&check_isos);
  pd_code_free(&retpd); /* This was a new allocation */
  return true;
  
}

bool torusknot_insert_and_search(int nknots,pd_stor_t **pdstor)

/* Insert torusknots into the pdstor, along with (many)
   duplicate entries, then attempt various isomorphism
   searches. */

{
  int i;
  pd_code_t *pd;

  printf("testing copyinto and search for torusknots with up to %d crossings\n",nknots);

  (*pdstor) = pd_new_pdstor();

  printf("testing searches on empty pdstor ... ");

  for(i=0;i<nknots;i++) {

    pd = pd_build_torus_knot(2,i+3); 
    if (!should_fail_iso_search(pd,*pdstor)) { return false; }
    pd_code_free(&pd);

  }

  printf("pass.\n");
  printf("now inserting collection of torus knots ... ");

  for(i=0;i<nknots;i++) {

    pd_code_t *pd;

    pd = pd_build_torus_knot(2,i+3);
    pd_addto_pdstor(*pdstor,pd,ISOMORPHISM); /* Do the insert. */
    pd_code_free(&pd);

  }

  printf("done.\n"); /* If we didn't crash, we pass. */
  printf("checking stats ... ");

  unsigned int nhashes, nelts;
  pd_stor_stats(*pdstor,&nhashes,&nelts);
  
  if (nhashes == nknots && nelts == nknots) {

    printf("pass (nhashes %u == nelts %u == nknots %d).\n",nhashes,nelts,nknots);

  } else {

    printf("FAIL (nhashes %u, nelts %u, nknots %d).\n",nhashes,nelts,nknots);
    return false;

  }

  printf("now searching for torus knots which should be present ... ");
  fflush(stdout);

  for(i=0;i<nknots;i++) {

    pd = pd_build_torus_knot(2,i+3); 
    if (!should_succeed_iso_search(pd,*pdstor)) { return false; }
    pd_code_free(&pd);
    
  }

  printf("pass.\n");
  
  printf("now searching for various knots/links which shouldn't be present ... ");

  for(i=2;i<14;i++) {

    pd_code_t *pd;
    pd = pd_build_twist_knot(i);
    if (!should_fail_iso_search(pd,*pdstor)) {
      printf("FAIL (found twist_knot %d)\n",i);
      return false; }
    pd_code_free(&pd);

  }

  for(i=3;i<5;i++) {

    pd_code_t *pd;
    pd = pd_build_simple_chain(i);
    if (!should_fail_iso_search(pd,*pdstor)) {
      printf("FAIL (found simple_chain %d)\n",i);
      return false; }
    pd_code_free(&pd);

  }

  for(i=2;i<9;i++) {

    pd_code_t *pd;
    pd = pd_build_unknot(i);
    if (!should_fail_iso_search(pd,*pdstor)) {
      printf("FAIL (found unknot %d)\n",i);
      return false; }
    pd_code_free(&pd);

  }

  printf("pass\n");

  return true;
     
}


bool test_search() {

  pd_stor_t *pdstor;

  printf(
	 "\n"
	 "pd_search_pdstor_by_isomorphism test suite\n"
	 "--------------------------------\n"
	 );

  if (!torusknot_insert_and_search(10,&pdstor)) { return false; }
  if (!free_test(&pdstor)) { return false; }
 
  printf("-----------------------------------\n"
	 "pd_copyinto_pdstor test suite        PASS\n");

  return true;

}

bool read_write_test()

/* Insert a bunch of crap into the pdstor, along with (many)
   duplicate entries, write it to a test file, read it from 
   the test file, and compare. */

{
  int i;
  pd_stor_t *pdstor;
  pd_code_t *pd;

  printf("testing read and write for pdstor\n");

  pdstor = pd_new_pdstor();

  printf("inserting unknotwyes and torusknots ... ");

  if (PD_VERBOSE) { printf("\n"); }

  for(i=0;i<14;i++) {

    pd = pd_build_torus_knot(2,i+3); 

    if (PD_VERBOSE) {

      printf("\t torus knot (2,%d) hash -> %s \n",i+3,pd->hash);

    }

    pd_addto_pdstor(pdstor,pd,ISOMORPHISM);
    pd_code_free(&pd);

  }

  unsigned int j,k;

  for(i=0;i<4;i++) { 

    for(j=0;j<4;j++) { 

      for(k=0;k<4;k++) {

	if (i+j+k <= 14-3) { 

	  pd = pd_build_unknot_wye(i,j,k);
	  
	  if (PD_VERBOSE) {
	    
	    printf("\t unknotwye (%d,%d,%d) hash -> %s \n",i,j,k,pd->hash);
	    
	  }
	  
	  pd_addto_pdstor(pdstor,pd,ISOMORPHISM);
	  pd_code_free(&pd);

	}

      }

    }

  }

  if (PD_VERBOSE) { printf("\n"); }

  printf("done.\n");

  unsigned int nelts,nhashes;
  pd_stor_stats(pdstor,&nhashes,&nelts);

  printf("created mixed pdstor with %u hashes and %u elts\n",nhashes,nelts);

  /* Now we go ahead and open a file for writing */

  FILE *outfile;
  outfile = fopen("./test.pdstor","w");
  assert(outfile != NULL);

  printf("writing pdstor to ./test.pdstor ... ");
  fflush(stdout);

  pd_write_pdstor(outfile,pdstor);
  fclose(outfile);

  printf("done.\n");

  pd_stor_t *new_pdstor;
  unsigned int new_nhashes,new_nelts;

  printf("now reopening test.pdstor and attempting read ... ");

  FILE *infile;
  infile = fopen("./test.pdstor","r");

  if (infile == NULL) { 

    printf("FAIL (couldn't open test.pdstor).\n");
    return false;

  }
  
  new_pdstor = pd_read_pdstor(infile,ISOMORPHISM);

  if (new_pdstor == NULL) {

    printf("FAIL (couldn't parse test.pdstor)\n");
    return false;

  }

  printf("pass.\n");
  printf("comparing # hashes in recovered version ... ");

  pd_stor_stats(new_pdstor,&new_nhashes,&new_nelts);

  if (new_nhashes != nhashes) {

    printf("FAIL (new_nhashes == %u != nhashes == %u)\n",new_nhashes,nhashes);
    return false;

  }
 
  printf("pass (new_nhashes == %u == nhashes == %u)\n",new_nhashes,nhashes);

  printf("compared # elts in recovered version ... ");

  if (new_nelts != nelts) { 

    printf("FAIL (new_nelts == %u != nelts == %u)\n",new_nelts,nelts);
    return false;

  }
 
  printf("pass (new_nelts == %u == nelts == %u)\n",new_nelts,nelts);

  printf("comparing elements of pdstors ... ");
  fflush(stdout);
  
  for(pd = pd_stor_firstelt(pdstor);pd != NULL;pd = pd_stor_nextelt(pdstor)) {

    if (!pd_ok(pd)) { 

      printf("FAIL. (element of original pd not ok).\n");
      return false; 

    }

    pd_code_t    *cmp_pd;
    pd_iso_t     **isos;
    unsigned int nisos;
    
    cmp_pd = pd_search_pdstor_by_isomorphism(new_pdstor,pd,&isos,&nisos);
    pd_free_isos(&nisos,&isos);

    if (cmp_pd == NULL) {

      printf("FAIL. (Couldn't find pd in recovered pdstor).\n");
      return false;

    }

    if (!pd_isomorphic(pd,cmp_pd)) {

      printf("FAIL. (Recovered pd not isomorphic to search).\n");
      return false;

    }

    pd_code_free(&cmp_pd); pd_code_free(&pd); /* These are new-memory allocations from the 
			       search functions. */

  }
   
  printf("pass.\n");

  pd_free_pdstor(&pdstor);
  pd_free_pdstor(&new_pdstor);
  
  return true;
     
}

bool test_io() {

  printf(
	 "\n"
	 "pd_write_pdstor and pd_read_pdstor test suite\n"
	 "--------------------------------\n"
	 );

  if (!read_write_test()) { return false; }
 
  printf("-----------------------------------\n"
	 "pd_copyinto_pdstor test suite        PASS\n");

  return true;

}


int main() {

  printf("test_pdstorage (%s)\n",PACKAGE_STRING);
  printf("---------------------------------------\n"
	 "Unit tests for pd_storage.c\n"
	 "=======================================\n");

  if (!test_copyinto() || !test_search() || !test_io()) {

    printf("=====================================\n");
    printf("test_pdstorage:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================\n");
    printf("test_pdstorage:  PASS.\n");
    exit(0);

  }

  return 0;

}
