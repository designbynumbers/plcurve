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

  printf("\t %d-crossing unknot ... pass.\n\n",n);

  pd_code_free(&pd);
  free(fee);
  free(cee);

  return true;

}

bool test_unknot() {

  if (!test_nunknot(2)) { return false; }
  if (!test_nunknot(3)) { return false; }
  if (!test_nunknot(4)) { return false; }
  if (!test_nunknot(5)) { return false; }
  if (!test_nunknot(6)) { return false; }
  
  return true;

}


int main() {

  printf("test_pdcode (%s)\n",PACKAGE_STRING);
  printf("---------------------------------------\n"
	 "Unit tests for pdcode.c\n"
	 "=======================================\n");

  if (!test_twist() || !test_torus() || !test_simple_chain() || !test_unknot() || !test_unknotwye()) {

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
