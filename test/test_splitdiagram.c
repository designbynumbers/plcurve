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

#include<plcTopology.h>

#include<pd_multidx.h>
#include<pd_perm.h>
#include<pd_dihedral.h>

#include<pd_isomorphisms.h>

int PD_VERBOSE=50;

pd_code_t *pd_joindiagram(pd_idx_t ncomponentPD,pd_code_t **componentPD)

/* Attaches pd_codes together to make a new one
   (which won't pass pd_ok, but should otherwise be ok).
*/

{
  pd_code_t *pd;

  pd_idx_t **new_crossing_number;
  pd_idx_t **new_edge_number;
  pd_idx_t **new_comp_number;
  pd_idx_t **new_face_number;

  pd_idx_t i;

  new_crossing_number = calloc(ncomponentPD,sizeof(pd_idx_t *));
  assert(new_crossing_number != NULL);
  
  new_edge_number = calloc(ncomponentPD,sizeof(pd_idx_t *));
  assert(new_edge_number != NULL);

  new_face_number = calloc(ncomponentPD,sizeof(pd_idx_t *));
  assert(new_face_number != NULL);

  new_comp_number = calloc(ncomponentPD,sizeof(pd_idx_t *));
  assert(new_comp_number != NULL);
  
  for(i=0;i<ncomponentPD;i++) { 

    if (componentPD[i]->ncross > 0) { 

      new_crossing_number[i] = calloc(componentPD[i]->ncross,sizeof(pd_idx_t));
      assert(new_crossing_number[i] != NULL);
    
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

  pd_idx_t total_cross=0, total_edge=0, total_face=0, total_comp=0; 

  for(i=0;i<ncomponentPD;i++) {

    total_cross += componentPD[i]->ncross;
    total_edge  += componentPD[i]->nedges;
    total_face  += componentPD[i]->nfaces;
    total_comp  += componentPD[i]->ncomp;
  
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

  /************** Everything below this line is unfinished ***********/
  
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
