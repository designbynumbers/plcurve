/* 

   Code for splitting diagrams into connected components.
   This isn't part of the exposed API because it operates 
   on diagrams that wouldn't pass pd_ok (on the other hand,
   it does _produce_ diagrams which should pass pd_ok).

   Jason Cantarella, November 2014.

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

#ifdef HAVE_STDLIB_H
   #include<stdlib.h>
#endif

#ifdef HAVE_STDINT_H
   #include<stdint.h>
#endif

#ifdef HAVE_STDBOOL_H
   #include<stdbool.h>
#endif

#ifdef HAVE_STDARG_H
   #include<stdarg.h>
#endif

#ifdef HAVE_ASSERT_H
   #include<assert.h>
#endif

#ifdef HAVE_CTYPE_H
   #include<ctype.h>
#endif

#ifdef HAVE_MATH_H
  #include<math.h>
#endif

#include<plcTopology.h>
//#include<libcassie/cassie.h>
//#include"/usr/local/include/thrift/Thrift.h"
//#include<python2.7/Python.h>
#include<pd_multidx.h>
#include<pd_dihedral.h>
#include<pd_perm.h>
#include<pd_orientation.h>
#include<pd_isomorphisms.h>
#include<pd_sortedbuf.h>
#include<pd_splitdiagram.h>

pd_idx_t pd_split_diagram(pd_code_t *pd,pd_code_t **pd_children)

/* Splits a disconnected diagram with valid (though maybe disordered)
   components, edges, and crossings into a collection of child
   diagrams.  The child diagrams WILL have valid faces (and
   canonicalized ordering for everything), and should pass pd_ok. Some
   of the children may be 0-crossing unknots. Component tags are
   transferred by this operation. The case where there is only one
   child is not special-- we just return 1 and the buffer of pd_children
   contains a single pointer.

   The array pd_children will be allocated in this procedure (and is the 
   caller's responsibility to dispose of), and the number of children 
   is returned. 
*/

{
  /* Basically, this is the usual "make a graph and split it into
     connected components" type of thing. The only clever idea is 
     that we actually partition the edges of the graph (and then 
     identify components from there) as it's faster to determine
     when edges interact than when components do. */

  pd_idx_t *edge_code;
  pd_idx_t ncodes = 1;

  edge_code = calloc(pd->nedges,sizeof(pd_idx_t));
  assert(edge_code != NULL);

  pd_idx_t i;
  for(i=0;i<pd->nedges;i++) { edge_code[i] = PD_UNSET_IDX; }

  pdint_split_worker(pd,&ncodes,edge_code,0);

  /* We should now have assigned a code to each edge, and be ready 
     to partition the diagram. This is going to be kind of interesting.
     For each code we have to do a lot of work:

     Identify the number of edges involved and allocate space.
     Identify the components involved, and copy the components into
     place in the new pd code. 
     Sort the components, and copy into place with tags.
     Copy edges into place according to their new numbering.
     Copy crossings into place, updating their edge records as we go.
     Canonicalize crossings.
     Regenerate faces. 

  */

  *pd_children = calloc(ncodes,sizeof(pd_code_t *));
  assert(pd_children != NULL);
     
  for(i=0;i<ncodes;i++) { 

    pd_idx_t j;
    pd_idx_t nedges = 0;
    pd_component_t ncomps = 0;

    for(j=0;j<pd->nedges;j++) { 
    
      if(edge_code[j] == i) { nedges++; }

    }

    for(j=0;j<pd->ncomps;j++) { 

      if (edge_code[pd->comp[j].edge[0]] == i) { ncomps++; }

    }

    /* We can now allocate space for the child pd code using 
       the fact that the number of vertices is always the number
       of edges divided by 2. */

    pd_children[i] = pd_code_new(nedges/2);

    pd_children[i].ncomps = 0;
    pd_children[i].comp = calloc(ncomps,sizeof(pd_component_t));

    for(j=0;j<pd->ncomps;j++) { 

      if (edge_code[pd->comp[j].edge[0]] == i) { 

	pd_children[i]->comp[pd_children[i]->ncomps].nedges = 
	  pd->comp[j].nedges;
	
	pd_children[i]->comp[pd_children[i]->ncomps].tag = 
	  pd->comp[j].tag;

	if (pd->comp[j].nedges > 0) { 

	  pd_children[i]->comp[pd_children[i]->ncomps].edge = 
	    calloc(pd_children[i]->comp[pd_children[i]->ncomps].nedges,
		   sizeof(pd_idx_t));

	  assert(pd_children[i]->comp[pd_children[i]->ncomps].edge != NULL);

	  memcpy(pd_children[i]->comp[pd_children[i]->ncomps].edge,
		 pd->comp[j].edge,pd->comp[j].nedges * sizeof(pd_idx_t));

	} else { /* There aren't any edges here. */

	  pd_children[i]->comp[pd_children[i]->ncomps].edge = NULL;

	}

	pd_children[i].ncomps++;

      }

    }

    /* We've now copied the new components into place. Our job is now to 
       sort them and prepare to assign new edge numbers. */

    qsort(pd_children[i]->comp,pd_children[i]->ncomps,
	  sizeof(pd_component_t),pd_component_cmp);
    
    /* Copy edge records into place in this new ordering and
       set component edge records in the child to the new 
       ordering. */

    pd_idx_t edgenum = 0;
    pd_idx_t *new_edge_num;
    new_edge_num = calloc(pd_children[i]->MAXEDGES,sizeof(pd_idx_t));
    assert(new_edge_num != NULL);
    
    for(j=0;j<pd_children[i]->ncomps;j++) { 

      for(k=0;k<pd_children[i]->comp[j].nedges;k++) { 

	new_edge_num[pd_children[i]->comp[j].edge[k]] = edgenum;
	pd_children[i]->edge[edgenum] = pd->edge[pd_children[i]->comp[j].edge[k]];
	pd_children[i]->comp[j].edge[k] = edgenum;

	edgenum++;

      }
      
    }

    /* We now need to make a list of crossings to copy. The fastest
       way to get it is to go ahead and scan the list of edges for head
       and tail crossings. Each crossing that we need will appear 4 times
       on the resulting list, but we can then sort it and only look at 
       the indices which are 0 mod 4. */

    pd_idx_t *crossings_to_copy;
    crossings_to_copy = calloc(2*pd_children[i]->nedges,sizeof(pd_idx_t));
    assert(crossings_to_copy != NULL);

    pd_idx_t *new_crossing_num;
    new_crossing_num = calloc(pd->nverts,sizeof(pd_idx_t));
    assert(new_crossing_num != NULL);
    
    for(j=0;j<pd_children[i]->nedges;j++) { 

      crossings_to_copy[2*j] = pd_children[i]->edge[j].head;
      crossings_to_copy[2*j + 1] = pd_children[i]->edge[j].tail;

    }

    qsort(crossings_to_copy,2*pd_children[i]->nedges,
	  sizeof(pd_idx_t),pd_idx_cmp);
    
    for(j=0;j<2*pd_children[i]->nedges;j+=4) { 

      pd_children[i]->cross[j/4].sign = pd->cross[crossings_to_copy[j]].sign;
      new_crossing_num[crossings_to_copy[j]] = j/4;
      
      for(k=0;k<4;k++) { 

	pd_children[i]->cross[j/4].edge[k] =    
	  new_edge_num[pd->cross[crossings_to_copy[j]].edge[k]];
	
      }

    }

    /* Now update the edges in the child with the new crossing numbers. */

    for(j=0;j<pd_children[i]->nedges;j++) { 

      pd_children[i]->edge[j].head = new_crossing_num[pd_children[i]->edge[j].head];
      pd_children[i]->edge[j].tail = new_crossing_num[pd_children[i]->edge[j].tail];

    }

    /* Ok, we've now moved components and edges over to the child. We can 
       free our working memory and canonicalize the crossings in the child. */

    free(new_edge_num);
    free(new_crossing_num);
    free(crossings_to_copy);

    pd_regenerate_crossings(pd_children[i]);
    pd_regenerate_faces(pd_children[i]);

    if (!pd_ok(pd_children[i])) { 

      pd_error(SRCLOC,"pd_split_diagram: child diagram %d given by %PD doesn't pass pd_ok\n",pd_children[i],i);
      exit(1);

    }

  }

  /* We've generated the children, so we can free the codes buffer. */

  free(edge_code);  
      
  /* We've generated the child pd codes in place, so there's not much left
     to do. However, we'll do some reasonable self-checking-- the total number
     of crossings, edges, and components in the children should add up to the 
     original total. */

  pd_idx_t total_crossings=0, total_edges=0, total_components=0;

  for(i=0;i<ncodes;i++) { 

    total_crossings += pd_children[i]->nverts;
    total_edges += pd_children[i]->nedges;
    total_components += pd_children[i]->ncomp;

  }

  if (total_crossings != pd->nverts) { 

    pd_printf("pd_split_diagram: The number of crossings in the %d child diagrams\n"
	      "of the pd %PD don't add up to the number of crossings (%d) in the pd\n"
	      "The child diagrams are:\n\n",
	      pd,ncodes,pd->nverts);

    for(j=0;j<ncodes;j++) { pd_printf("diagram %d:\n %PD ",pd_children[i],i); }

    exit(1);

  }

  if (total_edges != pd->nedges) { 

    pd_printf("pd_split_diagram: The number of edges in the %d child diagrams\n"
	      "of the pd %PD don't add up to the number of edges (%d) in the pd\n"
	      "The child diagrams are:\n\n",
	      pd,ncodes,pd->nedges);

    for(j=0;j<ncodes;j++) { pd_printf("diagram %d:\n %PD ",pd_children[i],i); }

    exit(1);

  }

  if (total_components != pd->ncomps) { 

    pd_printf("pd_split_diagram: The number of components in the %d child diagrams\n"
	      "of the pd %PD don't add up to the number of components (%d) in the pd\n"
	      "The child diagrams are:\n\n",
	      pd,ncodes,pd->ncomps);

    for(j=0;j<ncodes;j++) { pd_printf("diagram %d:\n %PD ",pd_children[i],i); }

    exit(1);

  }

  /* This is as much checking as we're willing to do on each run */
  
  return ncodes;
      	
} 
  


 
