/* 

   test_pdcode.c : Unit tests for the code in pdcode.c


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

#include<plcTopology.h>

#include<pd_multidx.h>
#include<pd_perm.h>
#include<pd_dihedral.h>

#include<pd_isomorphisms.h>

int PD_VERBOSE=50;

bool ok_test(pd_code_t *pd) {

  printf("\t pd_ok ... ");
 
  if (pd_ok(pd)) {
    
    printf("pass\n"); return true;
  
  } else {

    printf("FAIL (pd not ok).\n");
    return false;

  }

}

bool ncomps_test(pd_code_t *pd,pd_idx_t ncomps_expected) {

  printf("\t ncomps == %d test ... ",ncomps_expected);

  if (pd->ncomps == ncomps_expected) { 
    
    printf("pass\n");
    return true;

  } else {

    pd_printf("FAIL. (ncomps == %d != %d) in pd %PD",pd,pd->ncomps,ncomps_expected);
    return false;

  }

}

int cmp_nedges(const void *A,const void *B)
{
  pd_idx_t *pA = (pd_idx_t *)(A);
  pd_idx_t *pB = (pd_idx_t *)(B);

  return *pB - *pA;
}

bool face_edgevec_test(pd_code_t *pd,pd_idx_t *face_edges_expected)

/* Checks whether the number of edges in each face matches the expected number. */

{
  printf("\t checking # of edges on each face ... ");
  qsort(face_edges_expected,pd->nfaces,sizeof(pd_idx_t),cmp_nedges);

  pd_idx_t face;

  for(face=0;face<pd->nfaces;face++) {

    if (pd->face[face].nedges != face_edges_expected[face]) { 

      pd_printf("FAIL. (%FACE has %d edges != expected # == %d) in pd %PD",pd,
		face,pd->face[face].nedges,face_edges_expected[face]);
      return false;

    }

  }

  printf("pass.\n");
  return true;

}

bool comp_edgevec_test(pd_code_t *pd,pd_idx_t *comp_edges_expected)

/* Checks whether the number of edges on each component matches the expected number. */

{
  printf("\t checking # of edges on each component ... ");
  qsort(comp_edges_expected,pd->ncomps,sizeof(pd_idx_t),cmp_nedges);

  pd_idx_t comp;

  for(comp=0;comp<pd->ncomps;comp++) {

    if (pd->comp[comp].nedges != comp_edges_expected[comp]) { 

      pd_printf("FAIL. (%COMP has %d edges != expected # == %d) in pd %PD",pd,
		comp,pd->comp[comp].nedges,comp_edges_expected[comp]);
      return false;

    }

  }

  printf("pass.\n");
  return true;

}



bool test_twistn(pd_idx_t n) {

  pd_code_t *pd;

  printf("testing twist knot with %d twists...\n",n);
  pd = pd_build_twist_knot(n);

  if (!ok_test(pd)) { return false; }
  if (!ncomps_test(pd,1)) { return false; }

  /* Now test the number of edges on each comp */
  
  pd_idx_t *cee;
  cee = calloc(pd->MAXCOMPONENTS,sizeof(pd_idx_t));

  cee[0] = pd->nedges;

  if (!comp_edgevec_test(pd,cee)) { return false; }

  /* Now test the number of edges on each face */

  pd_idx_t *fee,face;
  fee = calloc(pd->MAXFACES,sizeof(pd_idx_t));

  fee[0] = 2; 
  fee[1] = 3;
  
  for(face=2;face<n+1;face++) {

    fee[face] = 2;

  } 

  fee[face++] = 3;
  fee[face++] = n+1;
  fee[face++] = n+1;

  if (!face_edgevec_test(pd,fee)) { return false; }

  printf("\t %d-twist knot ... pass.\n\n",n);

  pd_code_free(&pd);
  free(cee);
  free(fee);

  return true;

}

bool test_twist() {

  if (!test_twistn(1)) { return false; }
  if (!test_twistn(2)) { return false; }
  if (!test_twistn(3)) { return false; }
  if (!test_twistn(4)) { return false; }
  if (!test_twistn(5)) { return false; }
  if (!test_twistn(6)) { return false; }
  
  return true;
}

bool test_unknot_wye_abcsum(pd_idx_t n) {

  pd_code_t *pd;
  char      test_hash[2*PD_HASHSIZE];

  printf("testing unknot wyes with a + b + c = %d twists ...\n\n",n);
  
  pd = pd_build_unknot_wye(n,0,0);
  strncpy(test_hash,pd->hash,PD_HASHSIZE);
  printf("everything should have hash %s.\n",test_hash);
  pd_code_free(&pd);

  pd_idx_t A,B;

  for(A=0;A<=n;A++) {

    for (B=A;B<=n;B++) { 

      pd_idx_t a,b,c;

      a = A;
      b = B-A;
      c = n-B;

      printf("building unknot wye with %d, %d, %d twists ... \n",a,b,c);
      
      pd = pd_build_unknot_wye(a,b,c);

      if (!ok_test(pd)) { printf("FAIL (not ok).\n"); return false; }
      if (!ncomps_test(pd,1)) { printf("FAIL (wrong # of components)\n"); return false; }

      /* Now test the number of edges on each comp */
      
      pd_idx_t *cee;
      cee = calloc(pd->MAXCOMPONENTS,sizeof(pd_idx_t));
      
      cee[0] = pd->nedges;
      
      if (!comp_edgevec_test(pd,cee)) { printf("FAIL (wrong # of edges)\n"); return false; }
      
      /* Now test the number of edges on each face */
      
      pd_idx_t *fee,face;
      fee = calloc(pd->MAXFACES,sizeof(pd_idx_t));

      fee[0] = 1; 
      fee[1] = 1;
      fee[2] = 1;
      
      for(face=3;face<n+3;face++) {
	
	fee[face] = 2;
	
      } 
     
      fee[face++] = 3;
      fee[face++] = 2*(n+3);
      
      if (!face_edgevec_test(pd,fee)) { printf("FAIL (wrong # of edges per face)\n"); return false; }

      printf("\t checking hash ... ");
      if (strcmp(pd->hash,test_hash)) { printf("FAIL (hash %s != expected %s)\n",pd->hash,test_hash); return false; }
      printf("pass.\n\n");

      pd_code_free(&pd);
      free(cee);
      free(fee);

    }

  }

  printf("unknot wyes with a + b + c = %d twists ... pass.\n\n",n);
  return true;

}

bool test_unknotwye() {

  if (!test_unknot_wye_abcsum(1)) { return false; }
  if (!test_unknot_wye_abcsum(2)) { return false; }
  if (!test_unknot_wye_abcsum(3)) { return false; }
  if (!test_unknot_wye_abcsum(4)) { return false; }
  if (!test_unknot_wye_abcsum(5)) { return false; }
  if (!test_unknot_wye_abcsum(6)) { return false; }
  
  return true;
}

bool test_torus2q(pd_idx_t q) {

  pd_code_t *pd;

  printf("testing (2,%d) torus knot...\n",q);
  pd = pd_build_torus_knot(2,q);

  if (!ok_test(pd)) { return false; }

  if (!ncomps_test(pd,(q%2)==0 ? 2:1)) { return false; }

  /* Now test the number of edges on each comp */

  pd_idx_t *cee;
  cee = calloc(pd->MAXCOMPONENTS,sizeof(pd_idx_t));
  
  if ((q%2)==0) { /* link case */
    
    cee[0] = (pd_idx_t)(pd->nedges/2);
    cee[1] = (pd_idx_t)(pd->nedges/2);

  } else {
    
    cee[0] = pd->nedges;

  }

  if (!comp_edgevec_test(pd,cee)) { return false; }

  /* Now test the number of edges on each face */

  pd_idx_t *fee,face;
  fee = calloc(pd->MAXFACES,sizeof(pd_idx_t));

  for(face=0;face<q;face++) {

    fee[face] = 2;

  } 

  fee[face++] = q;
  fee[face++] = q;

  if (!face_edgevec_test(pd,fee)) { return false; }

  printf("\t (2,%d)-torus knot ... pass.\n\n",q);
  
  pd_code_free(&pd);
  free(cee);
  free(fee);

  return true;

}

bool test_torus() {

  if (!test_torus2q(2)) { return false; }
  if (!test_torus2q(3)) { return false; }
  if (!test_torus2q(4)) { return false; }
  if (!test_torus2q(5)) { return false; }
  if (!test_torus2q(6)) { return false; }
  
  return true;

}

bool test_nsimplechain(pd_idx_t n) {

  pd_code_t *pd;

  printf("testing %d-link simple chain ...\n",n);
  pd = pd_build_simple_chain(n);

  if (!ok_test(pd)) { return false; }

  if (!ncomps_test(pd,n)) { return false; }

  /* Now test the number of edges on each comp */

  pd_idx_t *cee,comp=0;
  cee = calloc(pd->MAXCOMPONENTS,sizeof(pd_idx_t));

  cee[0] = cee[1] = 2;

  for(comp=2;comp<n;comp++) { cee[comp] = 4; }

  if (!comp_edgevec_test(pd,cee)) { return false; }

  /* Now test the number of edges on each face */

  pd_idx_t *fee,face,i;
  fee = calloc(pd->MAXFACES,sizeof(pd_idx_t));

  fee[0] = fee[1] = 2; /* end faces */
  face = 2;
  
  for(i=0;i<n-1;i++,face++) {  /* n-1 "Middle" bigons */

    fee[face] = 2;

  } 

  for(i=0;i<n-2;i++,face++) {  /* n-2 internal 4-gons */

    fee[face] = 4;

  } 

  fee[face++] = 2*(n-2)+2; 
  /* One gigantic exterior face-- 2 edges from each interior link, 1 from each end link */

  if (!face_edgevec_test(pd,fee)) { return false; }

  printf("\t %d-component simple chain ... pass.\n\n",n);

  pd_code_free(&pd);
  free(fee);
  free(cee);

  return true;

}

bool test_simple_chain() {

  if (!test_nsimplechain(2)) { return false; }
  if (!test_nsimplechain(3)) { return false; }
  if (!test_nsimplechain(4)) { return false; }
  if (!test_nsimplechain(5)) { return false; }
  if (!test_nsimplechain(6)) { return false; }
  
  return true;

}


bool test_nunknot(pd_idx_t n) {

  pd_code_t *pd;

  printf("testing %d-crossing unknot ...\n",n);
  pd = pd_build_unknot(n);

  if (!ok_test(pd)) { return false; }

  if (!ncomps_test(pd,1)) { return false; }

  if (n > 0) { 

    /* Now test the number of edges on each comp */

    pd_idx_t *cee;
    cee = calloc(pd->MAXCOMPONENTS,sizeof(pd_idx_t));

    cee[0] = 2*n;
    if (!comp_edgevec_test(pd,cee)) { return false; }

    /* Now test the number of edges on each face */

    pd_idx_t *fee,face,i;
    fee = calloc(pd->MAXFACES,sizeof(pd_idx_t));

    fee[0] = fee[1] = 1; /* end faces */
    face = 2;
  
    for(i=0;i<n-1;i++,face++) {  /* n-1 "Middle" bigons */

      fee[face] = 2;

    } 

    fee[face++] = 2*(n-1)+2; 
    /* One gigantic exterior face-- 2 edges from each interior bigon, 1 from each end link */

    if (!face_edgevec_test(pd,fee)) { return false; }
    
    free(fee);
    free(cee);

  } 

  printf("\t testing pd_code_free on this unknot ...");
  pd_code_free(&pd);
  printf("pass (didn't crash)\n");
  printf("\t %d-crossing unknot ... pass.\n\n",n);
  
  return true;

}

bool test_unknot() {

  if (!test_nunknot(0)) { return false; }
  if (!test_nunknot(1)) { return false; }
  if (!test_nunknot(2)) { return false; }
  if (!test_nunknot(3)) { return false; }
  if (!test_nunknot(4)) { return false; }
  if (!test_nunknot(5)) { return false; }
  if (!test_nunknot(6)) { return false; }
  
  return true;

}

bool test_component_tag_0crossing() { 

  printf("------------------------------------------------------\n"
	 "tests to recover component tag from 0-crossing unknot \n"
	 "------------------------------------------------------\n");

  /* ---------- */

  printf("building temp_file_name...");
  char template[4096] = "/tmp/pdcodeXXXXXX";
  int outfile_fd = mkstemp(template);

  if (outfile_fd == -1) { 

    printf("fail (couldn't open %s)\n",template);
    return false;

  }

  FILE *outfile = fdopen(outfile_fd,"w");
  
  if (outfile == NULL) { 

    printf("fail (couldn't convert filedescriptor %d to stream)\n",outfile_fd);
    return false;

  }

  printf("done (%s)\n",template);
  
  /* ------------ */

  printf("creating 0-crossing unknot...");
  pd_code_t *pd = pd_build_unknot(0);

  if (pd_ok(pd)) { 

    printf("pass (passes pd_ok)\n");

  } else { 

    printf("FAIL (doesn't pass pd_ok)\n");
    return false;

  }

  printf("assigning tag 'P' to component...");
  pd->comp[0].tag = 'P';
  printf("done\n");

  printf("writing to file...");
  pd_write(outfile,pd);
  printf("pass (didn't crash)\n");

  /* ------------ */

  printf("closing file...");
  fclose(outfile);
  printf("done\n");
  
  printf("reopening file...");
  FILE *infile = fopen(template,"r");
  if (infile == NULL) { 
    printf("fail (couldn't reopen %s)\n",template);
    remove(template);
    return false;
  }
  printf("pass\n");

  /* ------------- */

  printf("reading new_pd from file...");
  pd_code_t *new_pd = pd_read(infile);
  if (new_pd == NULL) { 
    printf("fail (couldn't parse file)\n");
    remove(template);
    return false;
  }

  if (!pd_ok(new_pd)) { 
    pd_printf("fail (input %PD doesn't pass pd_ok)\n",new_pd);
    remove(template);
    return false;
  }

  printf("pass (read %d cross, %d component pd)\n",new_pd->ncross,new_pd->ncomps);

  /* ------------- */

  printf("testing isomorphic to original...");

  if (!pd_isomorphic(pd,new_pd)) { 
    pd_printf("fail (read pd %PD \n is not isomorphic to original",new_pd);
    pd_printf("%PD)",pd);
    remove(template);
    return false;
  }

  printf("pass\n");
  
  /* -------------- */

  printf("testing component 0 tag is 'P'...");
  if (new_pd->comp[0].tag == 'P') { 

    printf("pass\n");

  } else {

    printf("FAIL (new component 0 tag is %c)\n",new_pd->comp[0].tag);
    return false;

  }

  printf("housecleaning...");
  fclose(infile);
  remove(template);
  pd_code_free(&new_pd);
  printf("done\n");

  printf("------------------------------------------------------\n"
	 "recovering component tag from 0-crossing unknot: PASS \n"
	 "------------------------------------------------------\n");
  

  return true;
}
  

bool test_readwrite_pd(pd_code_t *pd,char *name) { 

  printf("----------------------------------------\n"
	 "tests for pd_write and pd_read of %s \n"
	 "----------------------------------------\n",name);

  /* ---------- */

  printf("building temp_file_name...");
  char template[4096] = "/tmp/pdcodeXXXXXX";
  int outfile_fd = mkstemp(template);

  if (outfile_fd == -1) { 

    printf("fail (couldn't open %s)\n",template);
    return false;

  }

  FILE *outfile = fdopen(outfile_fd,"w");
  
  if (outfile == NULL) { 

    printf("fail (couldn't convert filedescriptor %d to stream)\n",outfile_fd);
    return false;

  }

  printf("done (%s)\n",template);
  
  /* ------------ */

  printf("writing to file...");
  pd_write(outfile,pd);
  printf("pass (didn't crash)\n");

  /* ------------ */

  printf("closing file...");
  fclose(outfile);
  printf("done\n");
  
  printf("reopening file...");
  FILE *infile = fopen(template,"r");
  if (infile == NULL) { 
    printf("fail (couldn't reopen %s)\n",template);
    remove(template);
    return false;
  }
  printf("pass\n");

  /* ------------- */

  printf("reading new_pd from file...");
  pd_code_t *new_pd = pd_read(infile);
  if (new_pd == NULL) { 
    printf("fail (couldn't parse file)\n");
    remove(template);
    return false;
  }

  if (!pd_ok(new_pd)) { 
    pd_printf("fail (input %PD doesn't pass pd_ok)\n",new_pd);
    remove(template);
    return false;
  }

  printf("pass (read %d cross, %d component pd)\n",new_pd->ncross,new_pd->ncomps);

  /* ------------- */

  printf("testing isomorphic to original...");

  if (!pd_isomorphic(pd,new_pd)) { 
    pd_printf("fail (read pd %PD \n is not isomorphic to original",new_pd);
    pd_printf("%PD)",pd);
    remove(template);
    return false;
  }

  printf("pass\n");
  
  /* -------------- */

  printf("housecleaning...");
  fclose(infile);
  remove(template);
  pd_code_free(&new_pd);
  printf("done\n");

  printf("-------------------------------------------\n"
	 "tests for pd_write and pd_read of %s: PASS \n"
	 "-------------------------------------------\n",name);
  

  return true;
}
  
bool test_rw() {

  pd_code_t *pd;
  pd_idx_t cr;

  pd = pd_build_torus_knot(2,5);
  if (!test_readwrite_pd(pd,"(2,5) torus knot (crossings set)")) { return false; }  
  for(cr=0;cr<pd->ncross;cr++) { pd->cross[cr].sign = PD_UNSET_ORIENTATION; }
  if (!test_readwrite_pd(pd,"(2,5) torus knot (crossings unset)")) { return false; }
  pd_code_free(&pd);

  pd = pd_build_torus_knot(2,6);
  if (!test_readwrite_pd(pd,"(2,6) torus link (crossings set)")) { return false; }  
  for(cr=0;cr<pd->ncross;cr++) { pd->cross[cr].sign = PD_UNSET_ORIENTATION; }
  if (!test_readwrite_pd(pd,"(2,6) torus link (crossings unset)")) { return false; }
  pd_code_free(&pd);

  pd = pd_build_simple_chain(3);
  if (!test_readwrite_pd(pd,"3 link chain (crossings set)")) { return false; }  
  for(cr=0;cr<pd->ncross;cr++) { pd->cross[cr].sign = PD_UNSET_ORIENTATION; }
  if (!test_readwrite_pd(pd,"3 link chain (crossings unset)")) { return false; }
  pd_code_free(&pd);

  pd = pd_build_unknot(0);
  if (!test_readwrite_pd(pd,"0 crossing unknot diagram")) { return false; }  
  pd_code_free(&pd);
    
  return true;

}

bool test_read_pd(char *pd_code_file, char *testname,bool xpass) { 

  printf("----------------------------------------\n"
	 "tests for pd_read of %s \n"
	 "----------------------------------------\n",testname);

  /* ---------- */

  printf("building temp_file_name...");
  char template[4096] = "/tmp/pdcodeXXXXXX";
  int outfile_fd = mkstemp(template);

  if (outfile_fd == -1) { 

    printf("fail (couldn't open %s)\n",template);
    return false;

  }

  FILE *outfile = fdopen(outfile_fd,"w");
  
  if (outfile == NULL) { 

    printf("fail (couldn't convert filedescriptor %d to stream)\n",outfile_fd);
    return false;

  }

  printf("done (%s)\n",template);
  
  /* ------------ */

  printf("writing to file...");
  fprintf(outfile,"%s",pd_code_file);
  printf("pass (didn't crash)\n");

  /* ------------ */

  printf("closing file...");
  fclose(outfile);
  printf("done\n");
  
  printf("reopening file...");
  FILE *infile = fopen(template,"r");
  if (infile == NULL) { 
    printf("fail (couldn't reopen %s)\n",template);
    remove(template);
    return false;
  }
  printf("pass\n");

  /* ------------- */

  pd_code_t *new_pd;

  printf("reading new_pd from file ");
  if (xpass) { 
    printf("(expect pass)...");
    PD_VERBOSE = 50;
    new_pd = pd_read(infile); /* We WANT errors to stop the test. */
  } else {
    printf("(expect fail)...");
    PD_VERBOSE = 0;
    new_pd = pd_read(infile); /* We don't want errors to stop the test. */
    PD_VERBOSE = 50;
  }

  if (new_pd == NULL && xpass) { 
    printf("test FAILS (expected to parse file, couldn't)\n");
    remove(template);
    return false;
  }

  if (new_pd != NULL && !xpass) {

    printf("test FAILS (expected to NOT parse file, but could)\n");
    remove(template);
    return false;

  }

  printf("pass (did as expected)\n");

  if (xpass) { /* We only check the pd if, y'know, we read one */

    if (!pd_ok(new_pd)) {
 
      pd_printf("fail (input %PD doesn't pass pd_ok)\n",new_pd);
      remove(template);
      return false;

    }

    printf("pass (read %d cross, %d component pd)\n",new_pd->ncross,new_pd->ncomps);

  }
  
  /* -------------- */

  printf("housecleaning...");
  fclose(infile);
  remove(template);
  pd_code_free(&new_pd);  /* Should work even if new_pd is NULL */
  printf("done (didn't crash) \n");

  printf("-------------------------------------------\n"
	 "test for pd_read of %s: PASS \n"
	 "-------------------------------------------\n",testname);
  

  return true;
}

bool test_rw_altforms() { 

  /* This uses 

     bool test_read_pd(char *pd_code_file, char *testname,bool xpass) 

     together with some input pd codes to test various cases of the 
     pd_read code and make sure that bad input doesn't cause segfaults
     and that the various forms of good input are all ok. 
  */ 

  bool xpass = true;
  bool xfail = false;

  if (!test_read_pd(
		    "pd \n"
		    "nv 3\n"
		    "0 0 5 1 -\n"
		    "1 4 2 5 +\n"
		    "2 3 3 4 +\n"
		    "ne 6 \n"
		    "0,0 -> 0,1 \n" 
		    "0,3 -> 1,0 \n"
		    "1,2 -> 2,0 \n"
		    "2,2 -> 2,1 \n"
		    "2,3 -> 1,1 \n"
		    "1,3 -> 0,2 \n"
		    "nc 1 \n"
		    "6 : 0 1 2 3 4 5 tag A \n"
		    "nf 5 \n"
		    "4 : - 1 - 5 + 2 + 4 \n"
		    "3 : - 0 + 1 + 5 \n"
		    "3 : - 2 - 4 + 3 \n"
		    "1 : + 0 \n"
		    "1 : - 3 \n","no hash/uid test",xpass)) { return false; }

  if (!test_read_pd(
		    "pd pdunsethash 0\n"
		    "nv 3\n"
		    "0 0 5 1 -\n"
		    "1 4 2 5 +\n"
		    "2 3 3 4 +\n"
		    "ne 6 \n"
		    "0,0 -> 0,1 \n" 
		    "0,3 -> 1,0 \n"
		    "1,2 -> 2,0 \n"
		    "2,2 -> 2,1 \n"
		    "2,3 -> 1,1 \n"
		    "1,3 -> 0,2 \n"
		    "nc 1 \n"
		    "6 : 0 1 2 3 4 5 tag A \n"
		    "nf 5 \n"
		    "4 : - 1 - 5 + 2 + 4 \n"
		    "3 : - 0 + 1 + 5 \n"
		    "3 : - 2 - 4 + 3 \n"
		    "1 : + 0 \n"
		    "1 : - 3 \n","short hash/uid test",xpass)) { return false; }

  if (!test_read_pd(
		    "pd pdunsethash \n"
		    "nv 3\n"
		    "0 0 5 1 -\n"
		    "1 4 2 5 +\n"
		    "2 3 3 4 +\n"
		    "ne 6 \n"
		    "0,0 -> 0,1 \n" 
		    "0,3 -> 1,0 \n"
		    "1,2 -> 2,0 \n"
		    "2,2 -> 2,1 \n"
		    "2,3 -> 1,1 \n"
		    "1,3 -> 0,2 \n"
		    "nc 1 \n"
		    "6 : 0 1 2 3 4 5 tag A \n"
		    "nf 5 \n"
		    "4 : - 1 - 5 + 2 + 4 \n"
		    "3 : - 0 + 1 + 5 \n"
		    "3 : - 2 - 4 + 3 \n"
		    "1 : + 0 \n"
		    "1 : - 3 \n","hash, no uid test",xfail)) { return false; }

  if (!test_read_pd(
		    "pd \n"
		    "nv 3\n"
		    "0 0 5 1 -\n"
		    "1 4 2 5 +\n"
		    "2 3 3 4 +\n"
		    "ne 6 \n"
		    "0,0 -> 0,1 \n" 
		    "0,3 -> 1,0 \n"
		    "1,2 -> 2,0 \n"
		    "2,2 -> 2,1 \n"
		    "2,3 -> 1,1 \n"
		    "1,3 -> 0,2 \n"
		    "nc 1 \n"
		    "6 : 0 1 2 3 4 5 \n"
		    "nf 5 \n"
		    "4 : - 1 - 5 + 2 + 4 \n"
		    "3 : - 0 + 1 + 5 \n"
		    "3 : - 2 - 4 + 3 \n"
		    "1 : + 0 \n"
		    "1 : - 3 \n","no component tags",xpass)) { return false; }

  if (!test_read_pd(
		    "pd \n"
		    "nv 3\n"
		    "0 0 5 1 -\n"
		    "1 4 2 5 +\n"
		    "2 3 3 4 +\n"
		    "ne 6 \n"
		    "0,0 -> 0,1 \n" 
		    "0,3 -> 1,0 \n"
		    "1,2 -> 2,0 \n"
		    "2,2 -> 2,1 \n"
		    "2,3 -> 1,1 \n"
		    "1,3 -> 0,2 \n"
		    "nc 1 \n"
		    "6 : 0 1 2 3 4 5 tag 12\n"
		    "nf 5 \n"
		    "4 : - 1 - 5 + 2 + 4 \n"
		    "3 : - 0 + 1 + 5 \n"
		    "3 : - 2 - 4 + 3 \n"
		    "1 : + 0 \n"
		    "1 : - 3 \n","bad component tags",xfail)) { return false; }

  if (!test_read_pd(
		    "pd \n"
		    "nv 3\n"
		    "0 0 5 1 \n"
		    "1 4 2 5 \n"
		    "2 3 3 4 \n"
		    "ne 6 \n"
		    "0,0 -> 0,1 \n" 
		    "0,3 -> 1,0 \n"
		    "1,2 -> 2,0 \n"
		    "2,2 -> 2,1 \n"
		    "2,3 -> 1,1 \n"
		    "1,3 -> 0,2 \n"
		    "nc 1 \n"
		    "6 : 0 1 2 3 4 5 tag A\n"
		    "nf 5 \n"
		    "4 : - 1 - 5 + 2 + 4 \n"
		    "3 : - 0 + 1 + 5 \n"
		    "3 : - 2 - 4 + 3 \n"
		    "1 : + 0 \n"
		    "1 : - 3 \n","no crossing signs",xpass)) { return false; }

  if (!test_read_pd(
		    "pd \n"
		    "nv 3\n"
		    "0 0 5 1 +\n"
		    "1 4 2 5 +\n"
		    "2 3 3 4 +\n"
		    "ne 6 \n"
		    "0,0 -> 0,1 \n" 
		    "0,3 -> 1,0 \n"
		    "1,2 -> 2,0 \n"
		    "2,2 -> 2,1 \n"
		    "2,3 -> 1,1 \n"
		    "1,3 -> 0,2 \n"
		    "nc 1 \n"
		    "6 : 0 1 2 3 4 5 tag A\n"
		    "nf 5 \n"
		    "4 : - 1 - 5 + 2 + 4 \n"
		    "3 : - 0 + 1 + 5 \n"
		    /*	"3 : - 2 - 4 + 3 \n" */
		    "1 : + 0 \n"
		    "1 : - 3 \n","missing face",xfail)) { return false; }

  if (!test_read_pd(
		    "pd \n"
		    "nv 3\n"
		    "0 0 5 1 +\n"
		    "1 4 2 5 +\n"
		    "2 3 3 4 +\n"
		    "ne 6 \n"
		    "0,0 -> 0,1 \n" 
		    "0,3 -> 1,0 \n"
		    "1,2 -> 2,0 \n"
		    "2,2 -> 2,1 \n"
		    "2,3 -> 1,1 \n"
		    "nc 1 \n"
		    "6 : 0 1 2 3 4 5 tag A\n"
		    "nf 5 \n"
		    "4 : - 1 - 5 + 2 + 4 \n"
		    "3 : - 0 + 1 + 5 \n"
		    "3 : - 2 - 4 + 3 \n"
		    "1 : + 0 \n"
		    "1 : - 3 \n","missing edge",xfail)) { return false; }

  return true;

}

  

  


int main() {

  printf("test_pdcode (%s)\n",PACKAGE_STRING);
  printf("---------------------------------------\n"
	 "Unit tests for pdcode.c\n"
	 "=======================================\n");

  if (!test_component_tag_0crossing() || !test_rw_altforms() || !test_simple_chain() || !test_unknot() || !test_rw() || !test_twist() || !test_torus() ||  !test_unknotwye()) {

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
