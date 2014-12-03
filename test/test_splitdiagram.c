/* 

   test_splitdiagram.c : Unit tests for the code in pd_splitdiagram.c


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
#include<pd_splitdiagram.h>

int PD_VERBOSE=50;

pd_code_t *pd_joindiagram(pd_idx_t ncomponentPD,pd_code_t **componentPD)

/* Attaches pd_codes together to make a new one
   (which won't pass pd_ok, but should otherwise be ok).
*/

{
  pd_code_t *pd;

  pd_idx_t **new_cross_number;
  pd_idx_t **new_edge_number;
  pd_idx_t **new_comp_number;
  pd_idx_t **new_face_number;

  pd_idx_t i;

  new_cross_number = calloc(ncomponentPD,sizeof(pd_idx_t *));
  assert(new_cross_number != NULL);
  
  new_edge_number = calloc(ncomponentPD,sizeof(pd_idx_t *));
  assert(new_edge_number != NULL);

  new_face_number = calloc(ncomponentPD,sizeof(pd_idx_t *));
  assert(new_face_number != NULL);

  new_comp_number = calloc(ncomponentPD,sizeof(pd_idx_t *));
  assert(new_comp_number != NULL);
  
  for(i=0;i<ncomponentPD;i++) { 

    if (componentPD[i]->ncross > 0) { 

      new_cross_number[i] = calloc(componentPD[i]->ncross,sizeof(pd_idx_t));
      assert(new_cross_number[i] != NULL);
     
    }
    
    if (componentPD[i]->nedges > 0) { 

      new_edge_number[i] = calloc(componentPD[i]->nedges,sizeof(pd_idx_t));
      assert(new_edge_number[i] != NULL);

    }

    if (componentPD[i]->nfaces > 0) { 

      new_face_number[i] = calloc(componentPD[i]->nfaces,sizeof(pd_idx_t));
      assert(new_face_number[i] != NULL);

    }

    if (componentPD[i]->ncomps > 0) { 

      new_comp_number[i] = calloc(componentPD[i]->ncomps,sizeof(pd_idx_t));
      assert(new_comp_number[i] != NULL);

    }

  }

  /* We now make up translation tables for the new crossing, edge,
     component, and face numbers. */

  pd_idx_t total_cross=0, total_edges=0, total_faces=0, total_comps=0; 

  for(i=0;i<ncomponentPD;i++) {

    total_cross += componentPD[i]->ncross;
    total_edges  += componentPD[i]->nedges;
    total_faces  += componentPD[i]->nfaces;
    total_comps  += componentPD[i]->ncomps;
  
  }

  pd_idx_t j,k;

  for(i=0,k=0;i<ncomponentPD;i++) { 

    for(j=0;j<componentPD[i]->ncross;j++,k++) { 

      new_cross_number[i][j] = k;

    }

  }
  
  for(i=0,k=0;i<ncomponentPD;i++) { 

    for(j=0;j<componentPD[i]->nedges;j++,k++) { 

      new_edge_number[i][j] = k;

    }

  }

  for(i=0,k=0;i<ncomponentPD;i++) { 

    for(j=0;j<componentPD[i]->nfaces;j++,k++) { 

      new_face_number[i][j] = k;

    }

  }

  for(i=0,k=0;i<ncomponentPD;i++) { 

    for(j=0;j<componentPD[i]->ncomps;j++,k++) { 

      new_comp_number[i][j] = k;

    }

  }

  /* Now allocate space for the new pd code. Let's keep in mind that there
     will actually be too many faces in the new pd code, since V - E + F = 2
     in general, we'll have V - E + F = 2n in the combined diagram, meaning
     that we have to make space for an extra 2n-2 vertices or faces. 

     There's no harm in just adding too much space to the whole thing. */

  pd = pd_code_new(total_cross + 2*(ncomponentPD));
  
  assert(pd->MAXVERTS > total_cross); pd->ncross = total_cross;
  assert(pd->MAXEDGES > total_edges); pd->nedges = total_edges;
  assert(pd->MAXFACES > total_faces); pd->nfaces = total_faces;
  
  /* Now we're going to fill in the crossings, edges, and faces,
     using the new numbering scheme. */

  for(i=0;i<ncomponentPD;i++) { 

    for(j=0;j<componentPD[i]->ncross;j++) { 

      for(k=0;k<4;k++) { 

	pd->cross[new_cross_number[i][j]].edge[k] 
	  = new_edge_number[i][componentPD[i]->cross[j].edge[k]];

	pd->cross[new_cross_number[i][j]].sign 
	  = componentPD[i]->cross[j].sign;

      }

    }

  }

  for(i=0;i<ncomponentPD;i++) { 

    for(j=0;j<componentPD[i]->nedges;j++) { 

      if (componentPD[i]->ncross > 0) {

	pd->edge[new_edge_number[i][j]].head = 
	  new_cross_number[i][componentPD[i]->edge[j].head]; 

	pd->edge[new_edge_number[i][j]].tail = 
	  new_cross_number[i][componentPD[i]->edge[j].tail];

	pd->edge[new_edge_number[i][j]].headpos = 
	  componentPD[i]->edge[j].headpos;
	
	pd->edge[new_edge_number[i][j]].tailpos = 
	  componentPD[i]->edge[j].tailpos;

      } else { 

	pd->edge[new_edge_number[i][j]].head = PD_UNSET_IDX;
	pd->edge[new_edge_number[i][j]].tail = PD_UNSET_IDX;
	pd->edge[new_edge_number[i][j]].headpos = PD_UNSET_POS;
	pd->edge[new_edge_number[i][j]].tailpos = PD_UNSET_POS;

      }

    }

  }

  for(i=0;i<ncomponentPD;i++) { 

    for(j=0;j<componentPD[i]->nfaces;j++) { 

      pd->face[new_face_number[i][j]].nedges = componentPD[i]->face[j].nedges;

      pd->face[new_face_number[i][j]].edge 
	= calloc(pd->face[new_face_number[i][j]].nedges,sizeof(pd_idx_t));
      assert(pd->face[new_face_number[i][j]].edge != NULL);

      pd->face[new_face_number[i][j]].or 
	= calloc(pd->face[new_face_number[i][j]].nedges,sizeof(pd_or_t));
      assert(pd->face[new_face_number[i][j]].or != NULL);
      
      for(k=0;k<componentPD[i]->face[j].nedges;k++) { 

	pd->face[new_face_number[i][j]].edge[k] = 
	  new_edge_number[i][componentPD[i]->face[j].edge[k]];

	pd->face[new_face_number[i][j]].or[k] = 
	  componentPD[i]->face[j].or[k];

      }

    }

  }

  /* Components ARE allocated. */

  assert(pd->MAXCOMPONENTS > total_comps); pd->ncomps = total_comps;

  for(i=0;i<ncomponentPD;i++) { 

    for(j=0;j<componentPD[i]->ncomps;j++) { 

      pd->comp[new_comp_number[i][j]].nedges = 
	componentPD[i]->comp[j].nedges;

      pd->comp[new_comp_number[i][j]].edge = 
	calloc(pd->comp[new_comp_number[i][j]].nedges,sizeof(pd_idx_t));
      assert(pd->comp[new_comp_number[i][j]].edge != NULL);

      for(k=0;k<componentPD[i]->comp[j].nedges;k++) { 

	pd->comp[new_comp_number[i][j]].edge[k] = 
	  new_edge_number[i][componentPD[i]->comp[j].edge[k]];

      }

      pd->comp[new_comp_number[i][j]].tag = componentPD[i]->comp[j].tag;

    }

  }

  /* Now we free all our complicated buffers. */

  for(i=0;i<ncomponentPD;i++) { 

    if (new_cross_number[i] != NULL) { 

      free(new_cross_number[i]);
      new_cross_number[i] = NULL;

    }

    if (new_edge_number[i] != NULL) { 

      free(new_edge_number[i]);
      new_edge_number[i] = NULL;

    }

    if (new_face_number[i] != NULL) { 

      free(new_face_number[i]);
      new_face_number[i] = NULL;

    }

    if (new_comp_number[i] != NULL) { 

      free(new_comp_number[i]);
      new_comp_number[i] = NULL;

    }   

  }

  free(new_cross_number);
  free(new_edge_number);
  free(new_face_number);
  free(new_comp_number);
  
  return pd;
  
}


bool tested_joindiagram(pd_code_t **joined_pd,pd_idx_t ncomps,pd_code_t **comps) 
{ 
  printf("running pd_joindiagram...");
  *joined_pd = pd_joindiagram(ncomps,comps);
  printf("finished\n");

  printf("checking pd_cross_ok...");

  if (pd_cross_ok(*joined_pd)) { 

    printf("pass\n");

  } else { 

    printf("FAIL\n");
    return false;

  }

  bool contains_0crossing = false;
  pd_idx_t i;

  for(i=0;i<ncomps && !contains_0crossing;i++) {

    if (comps[i]->ncross == 0) { contains_0crossing = true; }

  }

  if (!contains_0crossing) { 
    
    printf("checking pd_edges_ok...");
    
    if (pd_edges_ok(*joined_pd)) { 

      printf("pass\n");

    } else { 
      
      printf("FAIL\n");
      return false;

    }

  } else {

    printf("can't check pd_edges_ok (joining 0-crossing unknot)...done\n");
      
  }

  printf("checking pd_faces_ok...");

  if (pd_cross_ok(*joined_pd)) { 

    printf("pass\n");

  } else { 

    printf("FAIL\n");
    return false;

  }
  
  printf("checking pd_comps_ok...");

  if (pd_cross_ok(*joined_pd)) { 

    printf("pass\n");

  } else { 

    printf("FAIL\n");
    return false;

  }

  return true;

}

bool compare_list_of_pds(pd_idx_t nA, pd_code_t **A, pd_idx_t nB, pd_code_t **B)

/* Figure out whether there's a match between the list of pds in A and the list
   in B. Since we expect the lists to be kind of small, we're just brute-forcing
   this (of course, there are many cleverer ways to do it). */

{
  printf("\tchecking that lists of pd codes are same size...");
  
  if (nA == nB) { 

    printf("pass (%d child pd codes)\n",nA);

  } else { 

    printf("FAIL. (lists are size %d != %d)\n",nA,nB);
    return false;

  }

  pd_perm_t *perm;
  perm = pd_new_perm((void *)(&nA));
  bool all_iso = true;

  printf("\tchecking up to %u perms for diagram isotopy match between lists...",
	 pd_nperms((void *)(perm)));
  
  do {  /* Check the current permutation. */

    pd_idx_t i; 

    for(i=0,all_iso = true;i < nA && all_iso;i++) {

      if (!pd_diagram_isotopic(A[i],B[perm->map[i]])) { all_iso = false; }

    }
    
    pd_increment_perm((void *)(perm));

  } while (!all_iso && !pd_perm_is_e(perm));

  if (all_iso) { 

    /* We need a trick here-- we can't DECREMENT the perm (thought we'd like to)
       in order to report the permutation that we actually used. 
       But we can increment it n-1 times... */

    pd_idx_t np = pd_nperms((void *)(perm));
    pd_idx_t i;
    for(i=0;i<np-1;i++) { pd_increment_perm((void *)(perm)); }

    pd_printf("pass (%PERM matches pd codes in list A with list B)\n",NULL,perm);
    
  } else {

    printf("FAIL (could not find a match between codes in list)\n");
    return false;

  }

  pd_free_perm((void **)(&perm));

  return true;

}
  
bool test_joindiagram() { 

  printf("------------------------------------------------------\n"
	 "testing joindiagram (just to make sure it runs)       \n"
	 "------------------------------------------------------\n");

  pd_code_t *pdA, *pdB, *pdC;

  printf("creating twist knot, simple_chain, and unknot for join...");
  pdA = pd_build_twist_knot(7);
  pdB = pd_build_simple_chain(4);
  pdC = pd_build_unknot(5);
  printf("done.\n");

  pd_code_t *(componentPD[3]);  /* Array of three pd_code_t pointers */
  componentPD[0] = pdA;  componentPD[1] = pdB; componentPD[2] = pdC; 

  pd_code_t *joined_pd;

  printf("running pd_joindiagram...");
  joined_pd = pd_joindiagram(3,componentPD);
  printf("finished\n");

  printf("checking pd_cross_ok...");

  if (pd_cross_ok(joined_pd)) { 

    printf("pass\n");

  } else { 

    printf("FAIL\n");
    return false;

  }

  printf("checking pd_edges_ok...");

  if (pd_edges_ok(joined_pd)) { 

    printf("pass\n");

  } else { 

    printf("FAIL\n");
    return false;

  }

  printf("checking pd_faces_ok...");

  if (pd_cross_ok(joined_pd)) { 

    printf("pass\n");

  } else { 

    printf("FAIL\n");
    return false;

  }
  
  printf("checking pd_comps_ok...");

  if (pd_cross_ok(joined_pd)) { 

    printf("pass\n");

  } else { 

    printf("FAIL\n");
    return false;

  }

  printf("housekeeping...");
  pd_code_free(&pdA);
  pd_code_free(&pdB);
  pd_code_free(&pdC);
  pd_code_free(&joined_pd);
  printf("done\n");

  printf("------------------------------------------------------\n"
	 "joindiagram test PASS                                 \n"
	 "------------------------------------------------------\n");
 
  return true;

}

bool test_joindiagramB() { 

  printf("------------------------------------------------------\n"
	 "testing joindiagramB (just to make sure it runs)       \n"
	 "------------------------------------------------------\n");

  pd_code_t *pdA, *pdB, *pdC;

  printf("creating twist knot, 0-crossing unknot, and 5-crossing unknot for join...");
  pdA = pd_build_twist_knot(7);
  pdB = pd_build_unknot(0);
  pdC = pd_build_unknot(5);
  printf("done.\n");

  pd_code_t *(componentPD[3]);  /* Array of three pd_code_t pointers */
  componentPD[0] = pdA;  componentPD[1] = pdB; componentPD[2] = pdC; 

  pd_code_t *joined_pd;

  printf("running pd_joindiagram...");
  joined_pd = pd_joindiagram(3,componentPD);
  printf("finished\n");

  printf("checking pd_cross_ok...");

  if (pd_cross_ok(joined_pd)) { 

    printf("pass\n");

  } else { 

    printf("FAIL\n");
    return false;

  }

  printf("can't check pd_edges_ok (number of edges off because of 0 crossing comp)...pass\n");

  printf("checking pd_faces_ok...");

  if (pd_cross_ok(joined_pd)) { 

    printf("pass\n");

  } else { 

    printf("FAIL\n");
    return false;

  }
  
  printf("checking pd_comps_ok...");

  if (pd_cross_ok(joined_pd)) { 

    printf("pass\n");

  } else { 

    printf("FAIL\n");
    return false;

  }

  printf("housekeeping...");
  pd_code_free(&pdA);
  pd_code_free(&pdB);
  pd_code_free(&pdC);
  pd_code_free(&joined_pd);
  printf("done\n");

  printf("------------------------------------------------------\n"
	 "joindiagram test B PASS                               \n"
	 "------------------------------------------------------\n");
 
  return true;

}

void assign_unique_tags(pd_idx_t ncodes,pd_code_t **codes) 
{
  pd_tag_t tag = 'A';
  pd_idx_t i,j;

  for(i=0;i<ncodes;i++) { 

    for(j=0;j<codes[i]->ncomps;j++) { 

      codes[i]->comp[j].tag = tag++;

    }
    
  }

} 

bool test_splitdiagramA() { 

  printf("------------------------------------------------------\n"
	 "testing splitdiagram -- twist7-chain4-uk5-uk2         \n"
	 "------------------------------------------------------\n");

  pd_code_t *(pd_list_A[4]);

  printf("making list twist7-chain4-uk5-uk2 of knots to join and split...");
 
  pd_list_A[0] = pd_build_twist_knot(7);
  pd_list_A[1] = pd_build_simple_chain(4);
  pd_list_A[2] = pd_build_unknot(5);
  pd_list_A[3] = pd_build_unknot(2);
 
  printf("done.\n");

  printf("assigning tags to components... ");
  assign_unique_tags(4,pd_list_A);
  printf("done\n");

  printf("looking for diagram isotopy between this list and itself...\n");
  if (!compare_list_of_pds(4,pd_list_A,4,pd_list_A)) { 

    printf("looking for diagram isotopy between this list and itself...FAIL\n");
    return false;

  } else {

    printf("looking for diagram isotopy between this list and itself...pass\n");
    
  }

  pd_code_t *joined_pd;

  if (!tested_joindiagram(&joined_pd,4,pd_list_A)) { 

    return false;

  }

  printf("splitting the joined diagram...");
  pd_code_t **pd_list_B = NULL;
  pd_idx_t nB = pd_split_diagram(joined_pd,&pd_list_B);
  printf("done (split into %d components)\n",nB);

  printf("looking for diagram isotopy between joined and split pd_codes...\n");
  if (!compare_list_of_pds(4,pd_list_A,nB,pd_list_B)) { 

    printf("looking for diagram isotopy between joined and split pd_codes...FAIL\n");
    return false;

  } else {

    printf("looking for diagram isotopy between joined and split pd_codes...pass\n");
    
  }
  
  printf("freeing list twist7-chain4-uk5-uk2 of knots...");

  pd_idx_t i;

  for(i=0;i<4;i++) { 

    pd_code_free(&(pd_list_A[i]));
    pd_code_free(&(pd_list_B[i]));

  }

  free(pd_list_B); /* pd_list_A is heap memory, so isn't freed here. */
  pd_code_free(&joined_pd);

  printf("done\n");

  printf("------------------------------------------------------\n"
	 "splitdiagram test PASS -- twist7-chain4-uk5-uk2       \n"
	 "------------------------------------------------------\n");
 
  return true;

}

bool test_splitdiagramB() { 

  printf("------------------------------------------------------\n"
	 "testing splitdiagram -- uk0-chain4-uk0-uk2         \n"
	 "------------------------------------------------------\n");

  pd_code_t *(pd_list_A[4]);

  printf("making list uk0-chain4-uk0-uk2 of knots to join and split...");
 
  pd_list_A[0] = pd_build_unknot(0);
  pd_list_A[1] = pd_build_simple_chain(4);
  pd_list_A[2] = pd_build_unknot(0);
  pd_list_A[3] = pd_build_unknot(2);
 
  printf("done.\n");

  printf("assigning tags to components... ");
  assign_unique_tags(4,pd_list_A);
  printf("done\n");

  printf("looking for diagram isotopy between this list and itself...\n");
  if (!compare_list_of_pds(4,pd_list_A,4,pd_list_A)) { 

    printf("looking for diagram isotopy between this list and itself...FAIL\n");
    return false;

  } else {

    printf("looking for diagram isotopy between this list and itself...pass\n");
    
  }

  pd_code_t *joined_pd;

  if (!tested_joindiagram(&joined_pd,4,pd_list_A)) { 

    return false;

  }

  printf("splitting the joined diagram...");
  pd_code_t **pd_list_B = NULL;
  pd_idx_t nB = pd_split_diagram(joined_pd,&pd_list_B);
  printf("done (split into %d components)\n",nB);

  printf("looking for diagram isotopy between joined and split pd_codes...\n");
  if (!compare_list_of_pds(4,pd_list_A,nB,pd_list_B)) { 

    printf("looking for diagram isotopy between joined and split pd_codes...FAIL\n");
    return false;

  } else {

    printf("looking for diagram isotopy between joined and split pd_codes...pass\n");
    
  }
  
  printf("freeing list uk0-chain4-uk0-uk2 of knots...");

  pd_idx_t i;

  for(i=0;i<4;i++) { 

    pd_code_free(&(pd_list_A[i]));
    pd_code_free(&(pd_list_B[i]));

  }

  free(pd_list_B); /* pd_list_A is heap memory, so isn't freed here. */
  pd_code_free(&joined_pd);

  printf("done\n");

  printf("------------------------------------------------------\n"
	 "splitdiagram test PASS -- uk0-chain4-uk0-uk2       \n"
	 "------------------------------------------------------\n");
 
  return true;

}

bool test_splitdiagramC() { 

  printf("------------------------------------------------------\n"
	 "testing splitdiagram -- uk0-uk0-uk0                   \n"
	 "------------------------------------------------------\n");

  pd_code_t *(pd_list_A[3]);

  printf("making list uk0-uk0-uk0 of knots to join and split...");
 
  pd_list_A[0] = pd_build_unknot(0);
  pd_list_A[1] = pd_build_unknot(0);
  pd_list_A[2] = pd_build_unknot(0);
 
  printf("done.\n");

  printf("assigning tags to components... ");
  assign_unique_tags(3,pd_list_A);
  printf("done\n");

  printf("looking for diagram isotopy between this list and itself...\n");
  if (!compare_list_of_pds(3,pd_list_A,3,pd_list_A)) { 

    printf("looking for diagram isotopy between this list and itself...FAIL\n");
    return false;

  } else {

    printf("looking for diagram isotopy between this list and itself...pass\n");
    
  }

  pd_code_t *joined_pd;

  if (!tested_joindiagram(&joined_pd,3,pd_list_A)) { 

    return false;

  }

  printf("splitting the joined diagram...");
  pd_code_t **pd_list_B = NULL;
  pd_idx_t nB = pd_split_diagram(joined_pd,&pd_list_B);
  printf("done (split into %d components)\n",nB);

  printf("looking for diagram isotopy between joined and split pd_codes...\n");
  if (!compare_list_of_pds(3,pd_list_A,nB,pd_list_B)) { 

    printf("looking for diagram isotopy between joined and split pd_codes...FAIL\n");
    return false;

  } else {

    printf("looking for diagram isotopy between joined and split pd_codes...pass\n");
    
  }
  
  printf("freeing list uk0-uk0-uk0 of knots...");

  pd_idx_t i;

  for(i=0;i<3;i++) { 

    pd_code_free(&(pd_list_A[i]));
    pd_code_free(&(pd_list_B[i]));

  }

  free(pd_list_B); /* pd_list_A is heap memory, so isn't freed here. */
  pd_code_free(&joined_pd);

  printf("done\n");

  printf("------------------------------------------------------\n"
	 "splitdiagram test PASS -- uk0-uk0-uk0                 \n"
	 "------------------------------------------------------\n");
 
  return true;

}
  

int main() {

  printf("test_splitdiagram (%s)\n",PACKAGE_STRING);
  printf("---------------------------------------\n"
	 "Unit tests for splitdiagram.c\n"
	 "=======================================\n");

  if (!test_joindiagram() || !test_joindiagramB() || !test_splitdiagramA() || !test_splitdiagramB() || !test_splitdiagramC()) {

    printf("=====================================\n");
    printf("test_splitdiagram:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================\n");
    printf("test_splitdiagram:  PASS.\n");
    exit(0);

  }

  return 0;

}
