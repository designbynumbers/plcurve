/*
  pd_operations.c : Code to perform various topological operations on 
                    a pd-code. Note that this C file basically contains
		    primitive operations-- more sophisticated operations
		    are implemented in Python. 

*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
#endif

#ifdef HAVE_SYS_STAT_H
  #include<sys/stat.h>
#endif

#ifdef HAVE_SYS_TYPES_H
  #include<sys/types.h>
#endif

#ifdef HAVE_ASSERT_H
  #include<assert.h>
#endif

#ifdef HAVE_STDINT_H
  #include<stdint.h>
#endif

#ifdef HAVE_STDIO_H
  #include<stdio.h>
#endif

#ifdef HAVE_STDINT_H
  #include<stdint.h>
#endif

#ifdef HAVE_STDBOOL_H
  #include<stdbool.h>
#endif

#ifdef HAVE_STDLIB_H
  #include<stdlib.h>
#endif

#ifdef HAVE_STRING_H
  #include<string.h>
#endif

#include<plcTopology.h>
#include<pd_multidx.h>
#include<pd_perm.h>
#include<pd_isomorphisms.h>


void pd_reorient_component(pd_code_t *pd, pd_idx_t cmp,pd_or_t or) 
{
  
  pd_check_notnull(SRCLOC,"pd",pd);
  pd_check_cmp(SRCLOC,pd,cmp);

  if (or == PD_NEG_ORIENTATION) {

    /* Several things need to happen in a proper component reversal:

     1) Run through the list of edges in the component, 
        reversing orientation on each one, and flipping the 
	sign on each head crossing (note that if a crossing 
	is between two strands on the same component, this will
	happen twice, and that's correct).

     2) Renumber the edges involved (because they are supposed
        to be consecutive and increasing along the component).
	Update the edge buffer and the crossings with the new
	edge numbering scheme.
	
     3) Regenerate faces (to make sure the +/- data is correct 
        when the edge is referred to on the face).
	
    */

    /* Step 1. Travel around component, flipping crossing signs. */

    pd_idx_t i;
    
    for(i=0;i<pd->comp[cmp].nedges;i++) {
      
      pd_idx_t e = pd->comp[cmp].edge[i];
      pd_idx_t headcross = pd->edge[e].head;
      
      if (pd->cross[headcross].sign != PD_UNSET_ORIENTATION) {
	
	pd->cross[headcross].sign = 
	  (pd->cross[headcross].sign == PD_POS_ORIENTATION) 
	? PD_NEG_ORIENTATION : PD_POS_ORIENTATION; 
      
      }
      
    }
    
    /* Step 2. Build and apply the transformation of edges. */

    /* We'll keep the first edge of the current component in place, then
       renumber all the other edges to match.

       We build this as an "edgemap", and use "pd_apply_edgemap"
       in order to make use of some of the self-checking code 
       for edgemaps. */
    
    pd_edgemap_t *emap = pd_new_edgemap(&pd->nedges);
    pd_idx_t e;
    
    /* First, set the edgemap to the identity. This will be the base for
       edges not affected by this transformation. */
    
    for(e=0;e<pd->nedges;e++) { 

      emap->or[e] = PD_POS_ORIENTATION;
      emap->perm->map[e] = e;

    }

    /* Now renumber and reorient the edges of THIS COMPONENT ONLY. */
    
    pd_idx_t cmp_edges = pd->comp[cmp].nedges;
    pd_idx_t first_edge, last_edge;

    first_edge = pd->comp[cmp].edge[0];
    last_edge  = pd->comp[cmp].edge[cmp_edges-1];
    
    assert(last_edge - first_edge + 1 == cmp_edges); 
    /* This is, we're assuming that the edges are consecutive. */
    
    emap->or[first_edge] = PD_NEG_ORIENTATION;
    emap->perm->map[first_edge] = first_edge;
 
    for(e=1;e<cmp_edges;e++) { 
      
      emap->or[first_edge+e] = PD_NEG_ORIENTATION;
      emap->perm->map[first_edge+e] = (last_edge+1) - e;

    }
    
    /* In case n is small, we need this to pass selftest. */
    pd_regenerate_pcidx(emap->perm);
    assert(pd_edgemap_ok(emap));
    
    pd_apply_edgemap(pd,emap);
    pd_free_edgemap(&emap);
    
    /* Step 3. Regenerate faces. */

    pd_regenerate_faces(pd); 
  
  }

}
