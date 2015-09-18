/* 

   Code for performing various topology-preserving moves on pd-codes.
   These will eventually include all the Reidemeister moves as well as some 
   generalized moves such as the Flype and the "clump slide". 

   Jason Cantarella, June 2014.

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

/* There are some unexposed primitives here which are going to 
   be useful in all of the simplifying operations: */

int  pd_idx_cmp(const void *a, const void *b)
{
  const pd_idx_t *idxa = (const pd_idx_t *) a;
  const pd_idx_t *idxb = (const pd_idx_t *) b;
  
  return (*idxa > *idxb) - (*idxa < *idxb);
}

void pd_compacting_copy(void *source, size_t obj_size, size_t nobj, 
			pd_idx_t ndeletions, pd_idx_t *deletions,
			void **target_var,
			pd_idx_t **target_idx,
			pd_idx_t *ntarget) /* Returns the size of the target buffer. */

/* This primitive copies (in order) all the elements of source_buf EXCEPT those
   in the buffer "deletions" to target buf. We allocate "target_var" ourselves,
   and return in the target_idx a buffer of size nobj so that

   target_idx[index in source] = index of copy, if the object was copied
                               = PD_UNSET_IDX, if the object was deleted

   The basic idea is that we're going to sort the buffer of deletions, and 
   then read through while copying to the target buffer, incrementing the 
   deletions buffer when we find one of the elements to delete. 

   This is linear time, and most importantly, we only need to write (and debug)
   it once. If we delete EVERYTHING (this can happen), then we set target_var to 
   NULL, ntarget to zero, and all the entries of target_idx to PD_UNSET_IDX.

*/
		   
{
  assert(source != NULL);
  assert(obj_size != 0);
  assert(target_var != NULL);

  pd_idx_t i,j;

  if (ndeletions == 0) { 

    /* Don't touch the deletions buffer, which may be null. */
    /* Just do a straight memcpy, and return a buffer telling */
    /* you that the indices haven't changed. */

    *target_var = malloc(nobj*obj_size);
    assert(*target_var != NULL);
    memcpy(*target_var,source,nobj*obj_size);
    
    *target_idx = malloc(nobj*sizeof(pd_idx_t));
    assert(*target_idx != NULL);
    for(i=0;i<nobj;i++) { (*target_idx)[i] = i; }
    *ntarget = nobj;

    return;

  }

  /* Now we know that there is at least one deletion to do. */
  /* We start by sorting the buffer of deletions. */
  
  assert(deletions != NULL);
  qsort(deletions,ndeletions,sizeof(pd_idx_t),pd_idx_cmp);

  /* Now we scan and make sure that there are no duplicate deletions. */
  /* This _could_ happen, depending on how careful the caller is, and */
  /* it's cheaper to check for it here than to worry about some corner */
  /* case later. */

  pd_idx_t *uniquedeletions;
  uniquedeletions = calloc(ndeletions,sizeof(pd_idx_t));
  assert(uniquedeletions != NULL);
  pd_idx_t nuniquedeletions = 0;
  
  uniquedeletions[0] = deletions[0];

  for(i=1,j=1;i<ndeletions;i++) { 

    if (uniquedeletions[j-1] != deletions[i]) { 

      uniquedeletions[j] = deletions[i];
      j++;

    }

  }
  
  /* We now scan to make sure that the deletions are in-range. */

  for(i=0;i<j;i++) { 

    assert(uniquedeletions[i] < nobj);

  }

  nuniquedeletions = j;

  if (nuniquedeletions == nobj) { /* Special case: we're deleting everything. */

    *target_var = NULL;
    *ntarget = 0;
    pd_idx_t *tidx = calloc(nobj,sizeof(pd_idx_t));
    assert(tidx != NULL);
    *target_idx = tidx;
    for(i=0;i<nobj;i++) { 
      tidx[i] = PD_UNSET_IDX;
    }
    return;
  }

  assert(nuniquedeletions < nobj); /* We've already taken care of the = case, 
				      and something is very wrong if we have 
				      MORE than nobj to delete! */

  /* Now we prepare space for the target array and index buffer. */

  void *target = calloc(nobj-nuniquedeletions,obj_size);
  assert(target != NULL);
  *target_var = target;
  
  pd_idx_t *tidx = calloc(nobj,sizeof(pd_idx_t));
  assert(tidx != NULL);
  *target_idx = tidx;

  *ntarget = nobj-nuniquedeletions;

  /* We're ready to run the main loop now. */

  pd_idx_t source_i,target_i,deletions_i;
  bool     all_deletions_used = false;

  for(source_i = target_i = 0, deletions_i = 0;source_i < nobj;source_i++) { 
    
    if (!all_deletions_used && source_i == uniquedeletions[deletions_i]) { 
      /* We're not copying this guy. */
	
	tidx[source_i] = PD_UNSET_IDX;
	deletions_i++;
	all_deletions_used = (deletions_i >= nuniquedeletions);
	
    } else {

      /* We're not at the current deletion, or there are no deletions left. */
      /* So we're definitely copying THIS guy. We need to do the copy with */
      /* memcpy because we're just copying bytes of memory-- we don't know */
      /* the data type involved. */

      memcpy(target + target_i*obj_size,source + source_i*obj_size,obj_size);
      tidx[source_i] = target_i;
      target_i++;

    }
    
  }

  /* Now do the requisite housecleaning. */

  free(uniquedeletions);
   
}

void pd_transform_and_compact_indices(pd_idx_t *source, 
				      pd_idx_t nobj, 
				      pd_idx_t *transform,
				      pd_idx_t **target,
				      pd_idx_t *ntarget)

/* This variant of pd_compacting_copy copies everything in the source
   buffer of indices to a compacted position in the newly allocated target buffer,
   transforming each index by transform, UNLESS the transform is set to 
   PD_UNSET_IDX, in which case we delete the element.

   We use the pd_compacting_copy API. */

{
  pd_idx_t i, ndeletions=0;
  for(i=0;i<nobj;i++) { 
    if (transform[source[i]] == PD_UNSET_IDX) { 
      ndeletions++; 
    } 
  }

  /*  void pd_compacting_copy(void *source, size_t obj_size, size_t nobj, 
                              pd_idx_t ndeletions, pd_idx_t *deletions,
			      void **target_var,
			      pd_idx_t **target_idx,
			      pd_idx_t *ntarget) */

  pd_idx_t *target_idx;
  pd_idx_t *deletions;

  if (ndeletions == 0) { 
    
    pd_compacting_copy((void *)(source),sizeof(pd_idx_t),(size_t)(nobj),ndeletions,NULL,(void **)(target),&target_idx,ntarget);

  } else {

    deletions = calloc(ndeletions,sizeof(pd_idx_t));
    assert(deletions != NULL);

    ndeletions = 0;
    
    for(i=0;i<nobj;i++) { 
      if (transform[source[i]] == PD_UNSET_IDX) { 
	deletions[ndeletions++] = i; 
      } 
    }

    pd_compacting_copy((void *)(source),sizeof(pd_idx_t),(size_t)(nobj),ndeletions,deletions,
		       (void **)(target),&target_idx,ntarget);

    free(deletions);  /* We're done with this */
    
  }

  free(target_idx); /* We don't need it. */

  /* Now do the transform step. */
  
  for(i=0;i<*ntarget;i++) { 

    (*target)[i] = transform[(*target)[i]];

  }

  /* And we're done. */

}

void pd_compact_components(pd_code_t *pd, pd_code_t *outpd,pd_idx_t *new_edge_idx)

/* Performs a compacting copy on each component of pd, eliminating edges
   whose new_edge_idx is PD_UNSET_IDX. Then resorts components by length.
   Afterwards, we apply an edgemap to return the edges to canonical order. */

/* 
   This will use the transform_and_compact_indices api:
 
     void pd_transform_and_compact_indices(pd_idx_t *source, pd_idx_t nobj, 
				      pd_idx_t *transform,
				      pd_idx_t **target)
*/

{
  pd_idx_t i;

  /* First, we make sure that outpd doesn't have anything in the components
     data, and then we resize the buffer to make sure it's big enough to work in. */

  if (outpd->ncomps > 0) { 

    for(i=0;i<outpd->ncomps;i++) { 

      if (outpd->comp[i].edge != NULL) { 

	free(outpd->comp[i].edge);
	outpd->comp[i].edge = NULL;
	outpd->comp[i].nedges = 0;

      }

    }

  }

  if (outpd->comp != NULL) { free(outpd->comp); }
  outpd->ncomps = 0;
  outpd->comp = calloc(pd->ncomps,sizeof(pd_component_t));
  assert(outpd->comp != NULL);
  outpd->MAXCOMPONENTS = pd->ncomps;

  /* Now we're ready to do the transformation. */
  
  for(i=0;i<pd->ncomps;i++) { 

    pd_transform_and_compact_indices(pd->comp[i].edge,pd->comp[i].nedges,new_edge_idx,
				     &(outpd->comp[i].edge),&(outpd->comp[i].nedges));
    outpd->comp[i].tag = pd->comp[i].tag;

  }

  /* Keep in mind that some components may be entirely NOT PRESENT in 
     outpd. These will have survived "transform_and_compact_indices", but
     will have NULL edgebuffers and nedges == 0 */

  outpd->ncomps = pd->ncomps;
  qsort(outpd->comp,outpd->ncomps,sizeof(pd_component_t),pd_component_cmp);

  for(i=0,outpd->ncomps=0;i<pd->ncomps && outpd->comp[i].nedges > 0;i++) { outpd->ncomps++; }
  
  /* 
     We now need to transform the edge indices again to make the edges
     consecutive along the new components. Having done so, we'll scramble
     the edge buffer in the corresponding way, which will cascade down to the
     crossing buffer (because we've changed the edge number references in 
     the crossings, we have to resort the crossings as well). 
  */ 

  pd_idx_t *resort_edge_idx;
  resort_edge_idx = calloc(outpd->nedges,sizeof(pd_idx_t));
  assert(resort_edge_idx != NULL);

  pd_edge_t *new_edge_buf = calloc(outpd->nedges,sizeof(pd_edge_t));
  assert(new_edge_buf != NULL);

  pd_idx_t j,edge=0;

  for(i=0;i<outpd->ncomps;i++) { 

    for(j=0;j<outpd->comp[i].nedges;j++,edge++) { 

      resort_edge_idx[outpd->comp[i].edge[j]] = edge;
      new_edge_buf[edge] = outpd->edge[outpd->comp[i].edge[j]];

    }

  }

  /* Now swap in the new edge buffer. */

  free(outpd->edge);
  outpd->edge = new_edge_buf;

  /* Finally, we apply the transformation to the crossings and components. 
     Faces we'll just regenerate. */

  for(i=0;i<outpd->ncross;i++) { 

    for(j=0;j<4;j++) { 

      outpd->cross[i].edge[j] = resort_edge_idx[outpd->cross[i].edge[j]];

    }

  }

  for(i=0;i<outpd->ncomps;i++) { 

    for(j=0;j<outpd->comp[i].nedges;j++) { 

      outpd->comp[i].edge[j] = resort_edge_idx[outpd->comp[i].edge[j]];

    }

  }

  /* Finally, we'll regenerate crossings and faces. */

  pd_regenerate_crossings(outpd);
  pd_regenerate_faces(outpd);

  free(resort_edge_idx);

}


  

pd_code_t *pd_R1_loopdeletion(pd_code_t *pd,pd_idx_t cr)

/* Performs a loop deletion Reidemeister 1 move at the crossing cr. 
   (A loop addition is a really different move, computationally speaking.) 

            \->e1->-/
             \     /
              \   /
               \ /
            cr  \            ---->    +----e0->---+  
               / \                    |           |
              /   \                   |           |
             e0    e2                 |           |
            /       \
*/

{
  pd_check_cr(SRCLOC,pd,cr);

  /* First, we find and label the edges e0, e1, and e2. */
  /* and figure out their positions ep[0], ep[1], and ep[2] */
  /* in the crossing. We have to keep track of orientation */
  /* here for everything else to make sense. */

  pd_idx_t e[3];  
  pd_pos_t ep[3];
  pd_idx_t i,j;
  bool found_flag = false;

  for(i=0;i<4 && !found_flag;i++) { 
    
    if (pd->cross[cr].edge[i] == pd->cross[cr].edge[(i+1)%4]) { 
      
      e[1] = pd->cross[cr].edge[i];

      ep[1] = pd->edge[e[1]].headpos;
      ep[0] = (pd->edge[e[1]].tailpos+2)%4; /* Across from the tail */
      ep[2] = (pd->edge[e[1]].headpos+2)%4; /* Across from the head */
      
      pd_idx_t j;
      for(j=0;j<3;j++) e[j] = pd->cross[cr].edge[ep[j]];
      
      found_flag = true;

    }
    
  }
  
  if (!found_flag) { 

    pd_error(SRCLOC, 
	     "crossing %CROSS of %PD is not a valid location\n"
	     "for loop deletion since it does not have a \n"
	     "repeated edge.\n",pd,cr);
    exit(1);

  }

  if (e[0] == e[2]) { 

    /* This is a one-crossing figure-8, or something very bad has happened. */

    if (pd->ncomps == 1 && pd->ncross == 1 && pd->nedges == 2 && pd->nfaces == 3) { 

      return pd_build_unknot(0);

    } else {

      pd_error(SRCLOC,
	       "only two edges are incident to %CROSS of %PD, but the pd_code"
	       "does not match the code for a 1-crossing figure-8.",pd,cr);
      exit(1);

    }

  }
  
  /* We now store the component number of this crossing for future reference. */

  pd_idx_t working_comp,working_pos;
  pd_component_and_pos(pd,e[0],&working_comp,&working_pos);
  
  /* We can now allocate space for the new pdcode, which has 
     one fewer crossing and two fewer edges than this one. */

  pd_code_t *outpd = pd_code_new(pd->MAXVERTS-1);

  /* Our first move is to copy the crossing buffer from the 
     old pd code into the new one, without the crossing we're 
     eliminating, of course! */

  /* void pd_compacting_copy(void *source, size_t obj_size, size_t nobj, 
                             pd_idx_t ndeletions, pd_idx_t *deletions,
			     void **target_var,
			     pd_idx_t **target_idx,
			     pd_idx_t *ntarget) */

  /* We're using the compacting_copy api here. */

  void *source = pd->cross;
  size_t obj_size = sizeof(pd_crossing_t);
  size_t nobj = pd->ncross;

  pd_idx_t ndeletions = 1;
  pd_idx_t deletions[1] = {cr};
  
  free(outpd->cross);  /* We need to free the previously allocated buffer, as we're replacing it. */
  outpd->cross = NULL; /* Just to mark it as freed. */

  void **target_var = (void **)(&(outpd->cross));
  pd_idx_t *new_crossing_idx;
  pd_idx_t *ntarget = &(outpd->ncross);

  pd_compacting_copy(source,obj_size,nobj,
		     ndeletions,deletions,
		     target_var,&(new_crossing_idx),ntarget);

  outpd->MAXVERTS = outpd->ncross; /* We trashed the buffer, so this is it! */

  /* We have now created a new buffer of crossings, and kept track 
     of the numberings associated to them with target_idx. */

  /* Our next move is to replace the reference to the old e[2]
     at the head of the old e[2] with a reference to e[0], since 
     the new e[0] will sew the tail of the old e[0] onto the head
     of the old e[2], eliminating the crossing in between.... */

  /* This is the slyest line of code in the whole thing-- basically
     everything else is housekeeping. Some pointers:

     1) The crossings have been renumbered, so the crossing we're
        looking to modify is NOT the head of edge 2 but the new 
	new_crossing_idx numbering of that crossing.

     2) edge e[2], like edge e[1] is gonna go away entirely. It just
        hasn't happened yet. */ 

  outpd->cross[new_crossing_idx[pd->edge[e[2]].head]].edge[pd->edge[e[2]].headpos] = e[0];

  /* We now make a compacting copy of the edge buffer which eliminates 
     edges e[1] and e[2]. */

  /* void pd_compacting_copy(void *source, size_t obj_size, size_t nobj, 
                             pd_idx_t ndeletions, pd_idx_t *deletions,
			     void **target_var,
			     pd_idx_t **target_idx,
			     pd_idx_t *ntarget) */

  source = (void *)(pd->edge);
  obj_size = sizeof(pd_edge_t);
  nobj = pd->nedges;

  ndeletions = 2;
  pd_idx_t newdeletions[2] = {e[1],e[2]};
  
  free(outpd->edge); /* We're replacing it. */
  outpd->edge = NULL;

  target_var = (void **)(&(outpd->edge));
  ntarget = &(outpd->nedges);
  
  pd_idx_t *new_edge_idx;

  pd_compacting_copy(source,obj_size,nobj,
		     ndeletions,newdeletions,
		     target_var,&(new_edge_idx),ntarget);

  outpd->MAXEDGES = outpd->nedges; /* Again, we trashed the buffer, so this is it. */

  /* Now we've eliminated edges e[1] and e[2] from the outpd
     edge buffer, which is good. But the crossing numbers in 
     the remaining edges still refer to crossing numbers in pd
     (not in outpd), so we need to translate. */

  for(i=0;i<outpd->nedges;i++) { 

    outpd->edge[i].head = new_crossing_idx[outpd->edge[i].head];
    outpd->edge[i].tail = new_crossing_idx[outpd->edge[i].tail];

  }

  /* Now we've just set the head of e[0] to PD_UNSET_IDX, since
     it used to point to the crossing, which we've just eliminated. 
     We need to change it to point to the (new numbering of) the 
     head of e[2], and remember that it must show up in the same position.*/

  outpd->edge[new_edge_idx[e[0]]].head = new_crossing_idx[pd->edge[e[2]].head];
  outpd->edge[new_edge_idx[e[0]]].headpos = pd->edge[e[2]].headpos;

  /* We've fixed the fact that the edges in outpd referred to the crossing
     numbers in pd. But we haven't fixed the fact that the crossings in outpd
     (though they've been compacted) still refer to the edge numbers in pd. 
     Now is a good time to change that. */

  for(i=0;i<outpd->ncross;i++) { 

    for(j=0;j<4;j++) { 

      outpd->cross[i].edge[j] = new_edge_idx[outpd->cross[i].edge[j]];

    }
  
  }

  /* We've now translated the new pd to have correct crossings and
     edges, which should be mutually consistent. Unfortunately, the
     crossings are out of order, so we first make sure that they are
     canonically ordered. This actually involves yet another
     translation of the crossing references in the edge data, but
     there's nothing to be done about that right now! */

  pd_regenerate_crossings(outpd);
  pd_compact_components(pd,outpd,new_edge_idx);  
  /* Transform, compact, and resort the components. (This fixes faces, too.) */
  /* The new pd now has crossings, edges, components, and faces which make sense. */

  pd_regenerate_hash(outpd);

  /* Finally, we check that we survived all this. */

  assert(pd_ok(outpd));
  
  /* Housecleaning */
  
  free(new_crossing_idx);
  free(new_edge_idx);

  /* Now the sanity check that we'd like to do is just to check numbers of things:
     We should have 

     1 fewer face
     1 fewer crossing
     2 fewer edges

  */

  assert(outpd->nfaces == pd->nfaces-1);
  assert(outpd->ncross == pd->ncross-1);
  assert(outpd->nedges == pd->nedges-2);

  return outpd;

}
      
  
void pd_R2_bigon_elimination(pd_code_t *pd,pd_idx_t cr[2],
			     pd_idx_t  *noutpd,
			     pd_code_t ***outpd)

/* Performs a bigon-elimination Reidemeister 2 move. 

|                                         |                      |                             |
|                                         |                      |        outA  pd code        |
|                                         |                      v                             ^
|                                         |                      |                             |
|eA[0]                                    |eA[2]                 |                             |
|                  eB[1]                  |                      |            eA[0]            |
|         +-------->------------+         |                      +-------------->--------------+
v         |                     |         ^                                                     
|         |                     |         |                                                     
|         |                     |         |                                                     
|         |       eA[1]         |         |                                   eB[0]                    
+-------------->---------->---------------+                       +------------->--------------+
     cr[0]|                     | cr[1]                           |                            |
          |                     |                                 |                            |
          |                     |                                 |                            |
          ^                     V                                 |       outB pd code         |
          |eB[0]                |eB[2]                            ^      (if different)        v
          |                     |                                 |                            |
          |                     |                                 |                            |

  Input is a pd code and two crossings defining the bigon. Output is a pair of 
  pointers to child pd codes. There are several possibilities. 

  noutpd counts the number of connected components of the output pd. 
  
  outpd[0] is the pd code containing the component "on top" in the bigon. It is NULL if 
  this component is 0-crossing diagram. 

  outpd[1] is always NULL if noutpd == 1 (in this case, both components are part of the same code).
  
  if noutpd == 2, then 

    outpd[1] is NULL if the component is a 0-crossing diagram (split unknot)
    outpd[1] otherwise contains the pd code of the component "on the bottom" 
    in the bigon. 

  There are several possible complications and special cases here.

  #1) You might entirely disconnect the diagram, meaning that you need to return 2 pd codes instead of 1.

      #1a) Either (or both) of these might be 0-crossing diagrams. 

  #2) The edges coming out of the diagram might not be distinct. 

  #3) The crossings must be checked for the correct over/under pattern.

  #4) We can ensure that the crossings cr[0] and cr[1] are in order on 
  component A, but they might occur in either order on component B
  since we can't predict the orientation of component B.

  #5) We're eliminating four edges (eB[1], eB[2], eA[1], and eA[2]) and 
  two crossings (cr[0] and cr[1]) from the diagram, and must do 
  compactions and relabelings accordingly. 

*/
 
{
  /* Step 1: Sanity check the input. */
  
  assert(pd_ok(pd));
  pd_check_cr(SRCLOC,pd,cr[0]);
  pd_check_cr(SRCLOC,pd,cr[1]);
  assert(noutpd != NULL);
  assert(outpd != NULL);

  /* Step 2: Try to verify that we're actually in 
     position for the move. */

  pd_idx_t eA[3],eB[3];

  pd_idx_t in_over,out_over;
  pd_overstrand(pd,cr[0],&in_over,&out_over);

  if (pd->edge[out_over].head != cr[1]) { 
    /* Crossings are not in order along cmpA */
    /* Try reversing order of crossings. */

    pd_idx_t swap; 
    swap = cr[0]; cr[0] = cr[1]; cr[1] = swap;

  }

  pd_overstrand(pd,cr[0],&in_over,&out_over);
  
  if (pd->edge[out_over].head != cr[1]) { 
    /* The crossings just aren't consecutive 
       along an overstrand. We have to quit. */

    pd_error(SRCLOC,
	     "crossings %CROSS and %CROSS aren't consecutive \n"
	     "along an overstrand in %PD\n"
	     "this is not a potential site for an R2 bigon elimination.\n",
	     pd,cr[0],cr[1]);
    exit(1);

  }

  eA[0] = in_over;
  eA[1] = out_over;
  eA[2] = pd_next_edge(pd,eA[1]);

  /* We're now sure that cr[1] follows cr[0] along the overstrand 
     from cr[0]. But is this the overstand at cr[1]? */

  pd_overstrand(pd,cr[1],&in_over,&out_over);
  
  if (eA[1] != in_over || eA[2] != out_over) { 

    pd_error(SRCLOC,
	     "edges eA[1] = %EDGE and eA[2] = %EDGE aren't over at crossing %CROSS in %PD\n"
	     "meaning that we can't do an R2 bigon elimination at %CROSS and %CROSS because\n"
	     "the crossing signs are wrong",pd,
	     eA[1],eA[2],cr[1],cr[0],cr[1]);
    exit(1);

  }

  /* Ok, we're now sure that the "overstrand" sequence eA[0] -> eA[2] is right. */
  /* Let's try to identify the "understrand" sequence. */

  pd_idx_t in_under,out_under;

  pd_understrand(pd,cr[0],&in_under,&out_under);
  
  if (pd->edge[out_under].head == cr[1]) { 
  /*
      (we're in this case)
   
     A                    A                   
     |   +---->------+    |  
     |   |           |    |                       
     +-cr[0]--->---cr[1]--+                      
         |           |         
   F1    B           B      F2   

  */
	     
    eB[0] = in_under;
    eB[1] = out_under;
    eB[2] = pd_next_edge(pd,eB[1]);

    /* We need to check we're going under at cr[1]. */

    pd_idx_t cr1_under_in,cr1_under_out;
    pd_understrand(pd,cr[1],&cr1_under_in,&cr1_under_out);

    if (eB[1] != cr1_under_in || eB[2] != cr1_under_out) { 

      pd_error(SRCLOC,
	     "edges eB[1] = %EDGE and eB[2] = %EDGE aren't under at crossing %CROSS in %PD\n"
	     "meaning that we can't do an R2 bigon elimination at %CROSS and %CROSS because\n"
	     "the crossing signs are wrong",pd,
	     eB[1],eB[2],cr[1],cr[0],cr[1]);
      exit(1);

    }

  } else if (pd->edge[in_under].tail == cr[1]) { 
    
 /*
      (we're in this case)
   
     A                    A                   
     |   +----<------+    |  
     |   |           |    |                       
     +-cr[0]--->---cr[1]--+                      
         |           |         
         B           B      

  */

    eB[2] = out_under;
    eB[1] = in_under;
    eB[0] = pd_previous_edge(pd,eB[1]);
    
    /* We need to check we're going under at cr[1]. */

    pd_idx_t cr1_under_in,cr1_under_out;
    pd_understrand(pd,cr[1],&cr1_under_in,&cr1_under_out);

    if (eB[0] != cr1_under_in || eB[1] != cr1_under_out) { 

      pd_error(SRCLOC,
	     "edges eB[0] = %EDGE and eB[1] = %EDGE aren't under at crossing %CROSS in %PD\n"
	     "meaning that we can't do an R2 bigon elimination at %CROSS and %CROSS because\n"
	     "the crossing signs are wrong",pd,
	     eB[0],eB[1],cr[1],cr[0],cr[1]);
      exit(1);

    }

  } else { /* Crossings just aren't consecutive along strand B, period. */

     pd_error(SRCLOC,
	     "crossings %CROSS and %CROSS aren't consecutive \n"
	     "along an understrand in %PD\n"
	     "this is not a potential site for an R2 bigon elimination.\n",
	     pd,cr[0],cr[1]);
     exit(1);

  }

  /* We now eliminate some special cases: 

+--------------+
|              |
|   +------+   |
|   |      |   |      -> 2 0-crossing diagrams.
+--------------+
    |      |    
    |      |    
    +------+    

or

            +----<---+    
            |        |    
            v        |    
     +--->------->------->--+   -> 1 0-crossing diagram. 
     |      |        |      |
     |      |        |      |
     +---<--+        +---<--+


  */

  if (pd->ncross == 2) { 

    if (pd->ncomps == 2) { 

      *noutpd = 2;
      *outpd = calloc(2,sizeof(pd_code_t *));
      assert(*outpd != NULL);

      (*outpd)[0] = pd_build_unknot(0);
      (*outpd)[1] = pd_build_unknot(0);

      /* Now we need to figure out how to set the tags. */

      pd_idx_t over_in, over_out;
      pd_idx_t under_in, under_out;

      pd_overstrand(pd,cr[0],&over_in,&over_out);
      pd_understrand(pd,cr[0],&under_in,&under_out);

      pd_idx_t over_comp, over_pos;
      pd_idx_t under_comp, under_pos;

      pd_component_and_pos(pd,over_in,&over_comp,&over_pos);
      pd_component_and_pos(pd,under_in,&under_comp,&under_pos);

      /* And finally we can set them. */

      (*outpd)[0]->comp[0].tag = pd->comp[over_comp].tag;
      (*outpd)[1]->comp[0].tag = pd->comp[under_comp].tag;

      return;

    } else if (pd->ncomps == 1) { 

      *noutpd = 1;
      *outpd = calloc(1,sizeof(pd_code_t *));
      assert(*outpd != NULL);
      (*outpd)[0] = pd_build_unknot(0);
      (*outpd)[0]->comp[0].tag = pd->comp[0].tag;

      return;

    }

  } 

  /* Step 3. Begin the actual work of the bigon elimination. */
  /* We basically bifurcate here into the cases where we split */
  /* and we don't split. The key is whether 

     |                    |     |               |
     A                    A     |               |
     v   +-----------+    ^     +------A--------+
     |   |           |    | ->      merged F1 and F2               
     +-cr[0]-->----cr[1]--+                      
     F1  |           |  F2      +------B--------+
         B           B          |               |

  the faces marked F1 and F2 are the same or different. 
  If they are the same, then we've split the diagram. 
  If different, we haven't. The course of action is 
  really different in each case.

  Determining which faces ARE F1 and F2 is not trivial. 

  |                 |           |              |
  v  +----------+   ^           |   +-A->--+   |
  |  |          |   |           |   |      |   |
  +------A->--------+     vs    +---|------|---+
     |          |                   ^      v 
     |          |                   |      |

  In the case at left, we're looking for the negfaces
  of eA[0] and eA[2], while in the case at right, we're 
  looking for the posfaces, so orientation isn't a tell.

  Instead, we should be looking for whether B[0] (or B[2])
  is to the left or right at the crossings.

  */

  pd_idx_t F1, F2;
  bool use_negfaces_of_eA0_and_eA2;

  /* To figure this out, we need to identify the bigon face. */
  
  pd_idx_t pfA,pfpA,nfA,nfpA;
  pd_idx_t pfB,pfpB,nfB,nfpB;

  pd_face_and_pos(pd,eA[1],&pfA,&pfpA,&nfA,&nfpA);
  pd_face_and_pos(pd,eB[1],&pfB,&pfpB,&nfB,&nfpB);

  if (pfA == nfB || pfB == nfA) { 

    /* Orientations AGREE: 

       |                 |           |              |
       v  +---B>-----+   ^           |   +-A->--+   |
       |  |          |   |           |   |      |   |
       +------A->--------+     vs    +---|--B>--|---+
          |          |                   ^      v 
	  |          |                   |      |
	  
    */

    use_negfaces_of_eA0_and_eA2 = 
      (pd->cross[cr[0]].edge[(pd->edge[eA[0]].headpos+1)%4] == eB[0]);

    /* That is, we use negfaces <=> the next (counterclockwise)
       edge from eA[0] at cr[0] is eB[0] */

  } else if (pfA == pfB || nfA == nfB) { 

     /* Orientations DISagree: 

       |                 |           |              |
       v  +---B<-----+   ^           |   +-A->--+   |
       |  |          |   |           |   |      |   |
       +------A->--------+     vs    +---|--B<--|---+
          |          |                   ^      v 
	  |          |                   |      |
	  
    */
    
    use_negfaces_of_eA0_and_eA2 =
      (pd->cross[cr[0]].edge[(pd->edge[eA[0]].headpos+1)%4] == eB[2]);

    /* That is, we use negfaces <=> the next (counterclockwise)
       edge from eA[0] at cr[0] is eB[2] */

  } else { 

    pd_error(SRCLOC,
	     "edges %EDGE and %EDGE don't share a bigon face in %PD"
	     "which means they cannot be the 'middle' edges in an R2\n"
	     "bigon elimination\n\n",pd,eA[1],eB[1]);
    exit(1);

  }
       
  /* We're going to use pd_face_and_pos to identify the faces. */
  
  pd_idx_t pf0,nf0,pfp0,nfp0;
  pd_idx_t pf2,nf2,pfp2,nfp2;

  pd_face_and_pos(pd,eA[0],&pf0,&pfp0,&nf0,&nfp0);
  pd_face_and_pos(pd,eA[2],&pf2,&pfp2,&nf2,&nfp2);
  
  if (use_negfaces_of_eA0_and_eA2) { 

    F1 = nf0; F2 = nf2;

  } else {

    F1 = pf0; F2 = pf2;

  }

  if (F1 != F2) { 
    /* This is the conventional case where we output a single PD code */

    *noutpd = 1;
    *outpd = calloc(1,sizeof(pd_code_t *));
    assert(*outpd != NULL);
    
    pd_code_t *out = pd_code_new(pd->ncross-2);
    (*outpd)[0] = out;

    /* Now we start by copying (compacted) the crossings. */

    free(out->cross);
    out->cross = NULL;
    out->MAXVERTS = out->ncross = 0;

    pd_idx_t crossing_deletions[2] = {cr[0],cr[1]};
    pd_idx_t *new_crossing_idx;
    
    pd_compacting_copy((void *)(pd->cross),sizeof(pd_crossing_t),(size_t)(pd->ncross),
		       2,crossing_deletions,
		       (void **)(&out->cross),
		       &new_crossing_idx,
		       &out->ncross);

    out->MAXVERTS = out->ncross;

    /* For the edges, we'll have a different set of deletions... */

    free(out->edge);
    out->edge = NULL;
    out->MAXEDGES = out->nedges = 0;

    pd_idx_t edge_deletions[4] = {eA[1],eA[2],eB[1],eB[2]};
    pd_idx_t *new_edge_idx;

    pd_compacting_copy((void *)(pd->edge),sizeof(pd_edge_t),(size_t)(pd->nedges),
		       4,edge_deletions,
		       (void **)(&out->edge),
		       &new_edge_idx,
		       &out->nedges);

    out->MAXEDGES = out->nedges;

    /* Now we need to actually sew together the strands to do the 
       move. As usual, these few lines of code are the most subtle
       bits of the whole move. On the left hand side, we're modifying
       positions in the new pd code. So these are addressed by their 
       NEW indices. On the right, we refer to OLD crossing and edge
       indices because we're about to translate absolutely everything
       to new indices and we don't want to translate twice.
    */

   

    if (eB[0] == eA[2]) {

      /*  
                +---eA[1]---+    
                v           ^    
		|           |
                |cr[1]      |cr[0]
      +--eB[0]--|->-eB[1]---|-->-eB[2]->
      |         |           |    
      |         |           |    
      +--eA[2]--+           +----eA[0]-<

      We're going to sew the head of eA[0] directly onto the head of eB[2] */

      out->cross[new_crossing_idx[pd->edge[eB[2]].head]].edge[pd->edge[eB[2]].headpos] = eA[0];      
      out->edge[new_edge_idx[eA[0]]].head = pd->edge[eB[2]].head;
      out->edge[new_edge_idx[eA[0]]].headpos = pd->edge[eB[2]].headpos;

    } else if (eB[2] == eA[0]) {

      /*
                +---eA[1]---+    
                v           ^    
		|           |
                |cr[1]      |cr[0]
      >--eB[0]--|->-eB[1]---|-->-eB[2]-+
                |           |          |
                |           |          |
      <--eA[2]--+           +----eA[0]-+

      We're going to sew the head of eB[0] directly onto the head of eA[2].

      */
	       

      out->cross[new_crossing_idx[pd->edge[eA[2]].head]].edge[pd->edge[eA[2]].headpos] = eB[0];
      out->edge[new_edge_idx[eB[0]]].head = pd->edge[eA[2]].head;
      out->edge[new_edge_idx[eB[0]]].headpos = pd->edge[eA[2]].headpos;

    } else { 

       /*
	 This is the generic case. 

                +---eA[1]---+    
                v           ^    
		|           |
                |cr[1]      |cr[0]
      >--eB[0]--|->-eB[1]---|-->-eB[2]->
                |           |          
                |           |          
      <--eA[2]--+           +----eA[0]-<

      We're going to sew the head of eA[0] directly onto the head of eA[2]
      and the head of eB[0] directly onto the head of eB[2].

      */

      out->cross[new_crossing_idx[pd->edge[eA[2]].head]].edge[pd->edge[eA[2]].headpos] = eA[0]; 
      out->edge[new_edge_idx[eA[0]]].head = pd->edge[eA[2]].head;
      out->edge[new_edge_idx[eA[0]]].headpos = pd->edge[eA[2]].headpos;

      out->cross[new_crossing_idx[pd->edge[eB[2]].head]].edge[pd->edge[eB[2]].headpos] = eB[0]; 
      out->edge[new_edge_idx[eB[0]]].head = pd->edge[eB[2]].head;
      out->edge[new_edge_idx[eB[0]]].headpos = pd->edge[eB[2]].headpos;

    }

    pd_idx_t i,j;

    for(i=0;i<out->ncross;i++) { 
      
      for(j=0;j<4;j++) { 
	
	out->cross[i].edge[j] = new_edge_idx[out->cross[i].edge[j]];

      }

    }

    for(i=0;i<out->nedges;i++) { 

      out->edge[i].head = new_crossing_idx[out->edge[i].head];
      out->edge[i].tail = new_crossing_idx[out->edge[i].tail];

    }

    /* We now need to regenerate crossings. */

    pd_regenerate_crossings(out);
     
    /* Components are interesting. We have the same set of 
       components (and want tags to persist). So we'll use 
       the pd_compact_components function...
    */
  
    pd_compact_components(pd,out,new_edge_idx);
    
    /* The new pd now has crossings, edges, and components which make sense. */
    /* It's ok to just regenerate the faces. */

    pd_regenerate_faces(out);
    pd_regenerate_hash(out);

    /* Finally, we check that we survived all this. */

    assert(pd_ok(out));
  
    /* Housecleaning */
  
    free(new_crossing_idx);
    free(new_edge_idx);

    /* Now the sanity check that we'd like to do is just to check numbers of things:
       We should have 

       2 fewer faces
       2 fewer crossings
       4 fewer edges
     
    */

    assert(out->nfaces == pd->nfaces-2);
    assert(out->ncross == pd->ncross-2);
    assert(out->nedges == pd->nedges-4);

  } else { 

    /* F1 == F2.

       This means that we're going to separate everything into two new
       PD codes.  This isn't going to be particuarly easy, because we
       don't even know which part of the diagram belongs to which
       child pd code until we figure out how the components are
       connected to one another at crossings.

       The straightforward strategy here is to generate a graph of
       connections between components and then find connected
       components of the graph. The problem is that we're not in
       Python, where that would be a reasonably easy thing to do. In
       C, we're going to be forced to import a heavyweight graph
       library, or introduce a weird dependency.

       So the idea is basically that we're going to mark the
       components as "A_unwalked", "B_unwalked", "A_walked",
       "B_walked", or "Unknown". The initial labels are "A_unwalked"
       and "B_unwalked" for the components containing arcs A and B,
       and "Unknown_unwalked" for everything else. 

       We're then going to make passes through the list of components,
       walking each "(A|B)_unwalked" component, marking every unknown
       component we see (except at the forbidden crossings cr[0] and
       cr[1]) as "(A|B)_unwalked" and ignoring every component which
       is "(A|B)_(un)walked". Finally, we mark this component as
       "walked", and keep searching for an unwalked component. 

       When we make a complete pass without finding any unwalked components,
       we're done. We then use the A/B classification to define the child 
       pd_codes.

    */

    typedef enum {A_unwalked,A_walked,B_unwalked,B_walked,Unknown} component_state_t;
    component_state_t *compstate;
    compstate = calloc(pd->ncomps,sizeof(component_state_t));
    assert(compstate != NULL);
    pd_idx_t i,j; 

    for(i=0;i<pd->ncomps;i++) { compstate[i] = Unknown; }
    
    pd_idx_t compA, compB;
    pd_idx_t cposA, cposB;

    pd_component_and_pos(pd,eA[0],&compA,&cposA);
    pd_component_and_pos(pd,eB[0],&compB,&cposB);
    assert(compA != compB);

    compstate[compA] = A_unwalked;
    compstate[compB] = B_unwalked;

    bool all_walked = false;
    pd_idx_t safety = 0;

    for(;!all_walked;safety++) { 

      all_walked = true; /* Hope this is the case */

      for(i=0;i<pd->ncomps;i++) { 

	if (compstate[i] == A_unwalked) { 

	  for(j=0;j<pd->comp[i].nedges;j++) { 

	    pd_idx_t thiscr,thispos,cross_cmp,cross_pos;

	    thiscr = pd->edge[pd->comp[i].edge[j]].head;
	    thispos = pd->edge[pd->comp[i].edge[j]].headpos;

	    if (thiscr != cr[0] && thiscr != cr[1]) { 

	      pd_component_and_pos(pd,pd->cross[thiscr].edge[(thispos+1)%4],&cross_cmp,&cross_pos);	    
	      if (compstate[cross_cmp] == Unknown) { compstate[cross_cmp] = A_unwalked; }

	    }

	  }

	  compstate[i] = A_walked;

	} else if (compstate[i] == B_unwalked) { 

	  for(j=0;j<pd->comp[i].nedges;j++) { 

	    pd_idx_t thiscr,thispos,cross_cmp,cross_pos;

	    thiscr = pd->edge[pd->comp[i].edge[j]].head;
	    thispos = pd->edge[pd->comp[i].edge[j]].headpos;

	    if (thiscr != cr[0] && thiscr != cr[1]) { 

	      pd_component_and_pos(pd,pd->cross[thiscr].edge[(thispos+1)%4],&cross_cmp,&cross_pos);	    
	      if (compstate[cross_cmp] == Unknown) { compstate[cross_cmp] = B_unwalked; }

	    }

	  }

	  compstate[i] = B_walked;
	  
	}

      }

      /* We've been through the list. Now scan for any unwalked/unknown components. */

      for(all_walked=true,i=0;i<pd->ncomps && all_walked;i++) { 
	if (compstate[i] == Unknown || compstate[i] == A_unwalked || compstate[i] == B_unwalked) { all_walked = false; } 
      }

    }

    /* Now we've either walked all the components or something has gone really wrong. */
    assert(all_walked);

    /* We are now ready to partition the diagram. We know the components involved in each 
       child diagram, so we'll have to get the edges and crossings as well now. */

    pd_idx_t *A_crossings, *B_crossings, *A_edges, *B_edges;

    A_crossings = calloc(2*pd->ncross,sizeof(pd_idx_t));
    assert(A_crossings != NULL);
    B_crossings = calloc(2*pd->ncross,sizeof(pd_idx_t));
    assert(B_crossings != NULL);

    A_edges = calloc(pd->nedges,sizeof(pd_idx_t));
    assert(A_edges != NULL);
    B_edges = calloc(pd->nedges,sizeof(pd_idx_t));
    assert(B_edges != NULL);

    pd_idx_t A_nedges = 0, B_nedges = 0, A_ncross = 0, B_ncross = 0;

    for(i=0;i<pd->ncomps;i++) { 

      if (compstate[i] == A_walked) { 

	for(j=0;j<pd->comp[i].nedges;j++) { 

	  pd_idx_t e;
	  e = pd->comp[i].edge[j];
	  
	  A_crossings[A_ncross++] = pd->edge[e].head;
	  A_edges[A_nedges++] = e;

	}

      } else if (compstate[i] == B_walked) { 

	for(j=0;j<pd->comp[i].nedges;j++) { 

	  pd_idx_t e;
	  e = pd->comp[i].edge[j];
	  
	  B_crossings[B_ncross++] = pd->edge[e].head;
	  B_edges[B_nedges++] = e;

	}

      } else {

	printf("compstate[%d] == %d, which is not A_walked (%d) or B_walked (%d)\n"
	       "suspect memory corruption bug, but your guess is as good as mine\n",
	       i,compstate[i],A_walked,B_walked);
	exit(1);

      }

    }

    /* We now expect the crossings in A and B to have been entered
       twice into our record-keeping scheme, since they will have 
       been encounted as heads of TWO edges. 
	 
       The edges, on the other hand, have only been encountered once. 
       
       We are now going to make up new pd codes, and make compacting
       copies, defining the crossings of A to be what's left when we 
       DELETE the list of B crossings (and vice versa) and the same 
       with the edges (we'll use B_edges as DELETIONS to define the
       actual edges of A implicitly).

       Note that we have already included the common crossings cr[0] 
       and cr[1] on both lists of deletions. But we haven't included
       the edges to be deleted. */

    A_edges[A_nedges++] = eB[1]; A_edges[A_nedges++] = eB[2];
    B_edges[B_nedges++] = eA[1]; B_edges[B_nedges++] = eA[2];

    pd_code_t *pdA, *pdB;
    pdA = pd_code_new(pd->ncross); /* Way too big, but that's ok right now */
    pdB = pd_code_new(pd->ncross); /* Again, too big, but we're going to trim the buffers anyway */

    pd_idx_t *new_crossing_idx_A, *new_edge_idx_A, *new_crossing_idx_B, *new_edge_idx_B;

    /* Handle crossings of A, remembering that we get A's crossings by DELETING B's */

    free(pdA->cross); pdA->cross = NULL; pdA->ncross = pdA->MAXVERTS = 0; 
    pd_compacting_copy((void *)(pd->cross),sizeof(pd_crossing_t),(size_t)(pd->ncross),
		       B_ncross,B_crossings,
		       (void *)(&(pdA->cross)),
		       &new_crossing_idx_A,
		       &(pdA->ncross));
    pdA->MAXVERTS = pdA->ncross;

    /* Handle edges of A, remembering that we get A's edges by DELETING B's */

    free(pdA->edge); pdA->edge = NULL; pdA->nedges = pdA->MAXEDGES = 0;
    pd_compacting_copy((void *)(pd->edge),sizeof(pd_edge_t),(size_t)(pd->nedges),
		       B_nedges,B_edges,
		       (void *)(&(pdA->edge)),
		       &new_edge_idx_A,
		       &(pdA->nedges));
    pdA->MAXEDGES = pdA->nedges;

    /* Handle crossings of B, remembering we get B's crossings by DELETING A's */
      
    free(pdB->cross); pdB->cross = NULL; pdB->ncross = pdB->MAXVERTS = 0; 
    pd_compacting_copy((void *)(pd->cross),sizeof(pd_crossing_t),(size_t)(pd->ncross),
		       A_ncross,A_crossings,
		       (void *)(&(pdB->cross)),
		       &new_crossing_idx_B,
		       &(pdB->ncross));
    pdB->MAXVERTS = pdB->ncross;

    /* Handle edges of B, remembering that we get B's edges by DELETING A's */

    free(pdB->edge); pdB->edge = NULL; pdB->nedges = pdB->MAXEDGES = 0;
    pd_compacting_copy((void *)(pd->edge),sizeof(pd_edge_t),(size_t)(pd->nedges),
		       A_nedges,A_edges,
		       (void *)(&(pdB->edge)),
		       &new_edge_idx_B,
		       &(pdB->nedges));
    pdB->MAXEDGES = pdB->nedges;

    /* Now we do the clever sewing bit, keeping in mind that on the left we 
       use NEW indices because we're in a child pd, and on the right we use
       OLD indices because we're referring to the parent pd. 

       We have to deal with a special case here when one or the other of the 
       child pds is a 0-crossing diagram.
    */

    if (pdA->ncross == 0) { 

      pd_code_free(&pdA);
      pdA = pd_build_unknot(0);
      pdA->comp[0].tag = pd->comp[compA].tag;

    } else {
      
      pdA->cross[new_crossing_idx_A[pd->edge[eA[2]].head]].edge[pd->edge[eA[2]].headpos] = eA[0]; 
      /* Since eA[2] was deleted, splice in eA[0]. */
      pdA->edge[new_edge_idx_A[eA[0]]].head = pd->edge[eA[2]].head;    /* Reset the head of edge eA[0] */
      pdA->edge[new_edge_idx_A[eA[0]]].headpos = pd->edge[eA[2]].headpos;

      /* Now we do the usual update of referred edge indices in the new crossings,
	 and referred crossing indices in the new edges. */

      for(i=0;i<pdA->ncross;i++) { 

	for(j=0;j<4;j++) { 

	  pdA->cross[i].edge[j] = new_edge_idx_A[pdA->cross[i].edge[j]];

	}

      }

      for (i=0;i<pdA->nedges;i++) { 

	pdA->edge[i].head = new_crossing_idx_A[pdA->edge[i].head];
	pdA->edge[i].tail = new_crossing_idx_A[pdA->edge[i].tail];

      }

    }

    /* Aaand for child diagram B.... */

    if (pdB->ncross == 0) { 

      pd_code_free(&pdB);
      pdB = pd_build_unknot(0);
      pdB->comp[0].tag = pd->comp[compB].tag;

    } else {

      pdB->cross[new_crossing_idx_B[pd->edge[eB[2]].head]].edge[pd->edge[eB[2]].headpos] = eB[0]; /* Since eB[2] was deleted, splice in eB[0]. */
      pdB->edge[new_edge_idx_B[eB[0]]].head = pd->edge[eB[2]].head;    /* Reset the head of edge eB[0] */
      pdB->edge[new_edge_idx_B[eB[0]]].headpos = pd->edge[eB[2]].headpos;

      for(i=0;i<pdB->ncross;i++) { 

	for(j=0;j<4;j++) { 

	  pdB->cross[i].edge[j] = new_edge_idx_B[pdB->cross[i].edge[j]];

	}

      }

      for (i=0;i<pdB->nedges;i++) { 

	pdB->edge[i].head = new_crossing_idx_B[pdB->edge[i].head];
	pdB->edge[i].tail = new_crossing_idx_B[pdB->edge[i].tail];

      }

    }

    /* Ok, now we need to copy components into the child pd's appropriately. 
       This takes a little thought, because we want the component TAGS to 
       travel appropriately. We need to exempt 0-crossing children from 
       the proceedings here because we've already copied their tags above. */

    if (pdA->ncross != 0) { pd_compact_components(pd,pdA,new_edge_idx_A); }
    if (pdB->ncross != 0) { pd_compact_components(pd,pdB,new_edge_idx_B); }

    /* Ok, at this point the children have crossings, edges, and
      components.  We can regenerate the crossings, and the
      faces. Again, if we're a 0-crossing component, we're exempt from
      this kind of thing-- at best it's unneeded, and at worst
      buggy. */

    if (pdA->ncross != 0) { 

      qsort(pdA->comp,pdA->ncomps,sizeof(pd_component_t),pd_component_cmp);
      pd_regenerate_crossings(pdA); 
      pd_regenerate_faces(pdA);
      pd_regenerate_hash(pdA);

    }

    if (pdB->ncross != 0) { 

      qsort(pdB->comp,pdB->ncomps,sizeof(pd_component_t),pd_component_cmp);

      pd_regenerate_crossings(pdB);
      pd_regenerate_faces(pdB);
      pd_regenerate_hash(pdB);

    }

    /* We can now do some counting to make sure things are sane: */

    assert(pdA->ncomps + pdB->ncomps == pd->ncomps);
    assert(pdA->ncross + pdB->ncross == pd->ncross - 2);

    /* These checks get screwed up for 0 crossing diagrams. */

    if (pdA->ncross != 0 && pdB->ncross != 0) { 
      assert(pdA->nedges + pdB->nedges == pd->nedges - 4);
    }

    if (pdA->ncross != 0) { 
      assert(pdA->ncross - pdA->nedges + pdA->nfaces == 2);
    }

    if (pdB->ncross != 0) { 
      assert(pdB->ncross - pdB->nedges + pdB->nfaces == 2);
    }

    assert(pd_ok(pdA));
    assert(pd_ok(pdB));

    /* Then we do some housecleaning. */

    free(compstate);
    free(new_crossing_idx_A); free(new_crossing_idx_B);
    free(new_edge_idx_A); free(new_edge_idx_B);
    
    /* Finally, we return the results. */

    *noutpd = 2;
    *outpd = calloc(2,sizeof(pd_code_t *));
    (*outpd)[0] = pdA;
    (*outpd)[1] = pdB;
    
  }

}


pd_code_t *pd_clump_slide(pd_code_t *pd,pd_idx_t e[3])

  /*                                                                
     A "clump slide" moves a connect summand K of a pd_code_t across an edge.
     The user is required to specify three edges along a component of pd in order
     to define the desired slide, as below.
     								    
                 f[0]         |                 |                      
	      +--------+      |                 |   +---------+        
              |        |      |                 |   |         |        
      --e[0]--+    K   +-e[1]---e[2]---  >  --------+    K    +------
              |        |      |                 |   |         |        
	      +--------+      |                 |   +---------+        
                 f[1]         |                 |                      
                   
      Notes: We don't know the orientation of the component, so it might be the
      case that e[0], e[1], and e[2] are positively or negatively orientated 
      with respect to the component orientation. 

      Also, the "clump" K may include additional components. However, there
      must be a pair of faces f[0] and f[1] which share e[0] and e[1] to isolate 
      the "clump". 

  */
{
  assert(pd != NULL);
  pd_check_edge(SRCLOC,pd,e[0]);
  pd_check_edge(SRCLOC,pd,e[1]);
  pd_check_edge(SRCLOC,pd,e[2]);

  if (e[0] == e[1] || e[1] == e[2] || e[0] == e[2]) { 

    pd_error(SRCLOC,"edges e[0] = %d, e[1] = %d, and e[2] = %d have to be different in order to clump slide.\n",
	     pd,e[0],e[1],e[2]);
    exit(1);
  }

  /* We need to check that e[1] and e[2] share a vertex */

  pd_or_t order_of_ei;

  if (pd->edge[e[1]].head == pd->edge[e[2]].tail) { 

    order_of_ei = PD_POS_ORIENTATION;

  } else if (pd->edge[e[1]].tail == pd->edge[e[2]].head) {

    order_of_ei = PD_NEG_ORIENTATION;

  } else {

    pd_error(SRCLOC,"clump_slide input error: edges e[1] == %EDGE and e[2] == %EDGE don't meet at a crossing\n",
	     pd,e[1],e[2]);
    exit(1);

  }

  /* Now we make sure that everybody is on the same component. */

  pd_idx_t comp[3], pos[3];
  pd_idx_t i,j;

  for(i=0;i<3;i++) { 
     
    pd_component_and_pos(pd,e[i],&comp[i],&pos[i]);

  }

  if (!(comp[0] == comp[1] && comp[1] == comp[2])) { 

    pd_error(SRCLOC,"clump_slide input error: edges e[0] = %EDGE, e[1] = %EDGE, e[2] = %EDGE are on components %d, %d, %d, which are not the same\n",
	     pd,e[0],e[1],e[2],comp[0],comp[1],comp[2]);
    exit(1);

  }

  /* Now we check for the existence of the faces f[0] and f[1] defining the clump.
      By definition, f[0] is going to be the face on which e[0] and e[1] occur positively
      and f[1] is going to be the face on which they occur negatively. 

                 f[0]
              
              +--------+        
              |        |             
      e[0]->--+        +-->-e[1]
              |        |          
              +--------+     
 
                 f[1]
   
   */
   
  pd_idx_t posface[2], negface[2], pfp[2], nfp[2];
  for(i=0;i<2;i++) { pd_face_and_pos(pd,e[i],&posface[i],&pfp[i],&negface[i],&nfp[i]); }
  pd_idx_t f[2];

  if (posface[0] == posface[1]) { f[0] = posface[0]; } 
  else {
    pd_error(SRCLOC,"edges e[0] = %EDGE and e[1] = %EDGE are positive on faces %FACE and %FACE which are not the same\n",
	     pd,e[0],e[1],posface[0],posface[1]);
    exit(1);
  }

  if (negface[0] == negface[1]) { f[1] = negface[0]; } 
  else {
    pd_error(SRCLOC,"edges e[0] = %EDGE and e[1] = %EDGE are negative on faces %FACE and %FACE which are not the same\n",
	     pd,e[0],e[1],negface[0],negface[1]);
    exit(1);
  }

  if (f[0] == f[1]) {
    pd_error(SRCLOC,"faces f[0] = %FACE and f[1] = %FACE defining a clump are the same face. %PD is corrupted",
	     pd,f[0],f[1]);
  }
    
   
  pd_code_t *out;
  out = pd_code_new(pd->ncross);
  pd_idx_t *new_edge_idx = calloc(pd->nedges,sizeof(pd_idx_t));
  assert(new_edge_idx != NULL);

  /* We've now checked that we are in position for the operation. The operation 
     itself turns out to be surprisingly simple. 
  */

  /* We start by just making a copy of the existing edges and 
     crossings. */
   
  for(i=0;i<pd->nedges;i++) { 
    out->edge[i] = pd->edge[i];
    new_edge_idx[i] = i; 
  } 

  for(i=0;i<pd->ncross;i++) { 
    out->cross[i] = pd->cross[i];
  }

  for(i=0;i<pd->ncomps;i++) { 

    out->comp[i].nedges = pd->comp[i].nedges;
    out->comp[i].tag = pd->comp[i].tag;

    out->comp[i].edge = calloc(out->comp[i].nedges,sizeof(pd_idx_t)); 
    assert(out->comp[i].edge != NULL);
     
    for(j=0;j<out->comp[i].nedges;j++) { 

      out->comp[i].edge[j] = pd->comp[i].edge[j];

    }

  }

  /* Note: Faces are going to get regenerated later, because they really change. */

  /*
    If the orientations AGREE, we're
    doing:

                              |                     |                      
	      +--------+      |                     |   +---------+        
              |        |      |                     |   |         |  e[2] (new indices)      
      --e[0]->+    K   +-e[1]>--e[2]->-  -->  ->e[0]-->-+    K    +-->---------------------
              |        |      |                     |   |         |        
	      +--------+      |                     |   +---------+        
                              |                     |                      
			      
			      which is to say that first, the edges between e[0]+1 and e[1] need to be shifted up 
			      by one index. Then we can splice a new edge in immediately after e[0] at the critical
			      crossing. The crossings are not going to get renumbered (except by resorting).

  */

  if (order_of_ei == PD_POS_ORIENTATION) { 

    /* Note the use of pd_next_edge as an iterator (instead of just
       incrementing i). The point is that we're looping around a
       component, whose edges very well might be consecutive, but
       somewhere in the middle of a large collection of indices. So
       figuring out the index of the "next in line" guy is not just
       some modular arithmetic, but something that's actually a
       pain. */

    for(i=pd_next_edge(pd,e[0]); i != e[2]; i = pd_next_edge(pd,i)) { 
       
      out->edge[pd_next_edge(pd,i)] = pd->edge[i];
      new_edge_idx[i] = pd_next_edge(pd,i);
       
    }

    /* At this point, we've trashed e[2] (in the copy), and duplicated
       the edge after e[0]. 

       We now do the actual sewing, which as usual is the cleverest
       lines in the code. The idea is that we're going to insert the 
       crossing that used to be between e[1] and e[2] in between 
       e[0] and the next edge to e[0], then go into the K stuff. */

    pd_idx_t e0next = pd_next_edge(pd,e[0]);

    out->edge[new_edge_idx[e0next]].head    = pd->edge[e[0]].head;
    out->edge[new_edge_idx[e0next]].headpos = pd->edge[e[0]].headpos;

    out->edge[new_edge_idx[e0next]].tail    = pd->edge[e[2]].tail;
    out->edge[new_edge_idx[e0next]].tailpos = pd->edge[e[2]].tailpos;

    out->edge[new_edge_idx[e[0]]].head    = pd->edge[e[1]].head;
    out->edge[new_edge_idx[e[0]]].headpos = pd->edge[e[1]].headpos;

    out->edge[new_edge_idx[e[1]]].head    = pd->edge[e[2]].head;
    out->edge[new_edge_idx[e[1]]].headpos = pd->edge[e[2]].headpos;

    out->cross[pd->edge[e[1]].head].edge[pd->edge[e[1]].headpos] = e[0]; /* Old edge number. */
    out->cross[pd->edge[e[1]].head].edge[pd->edge[e[2]].tailpos] = e0next; /* Old edge number. */

    /* We're going to translate crossing numbers (and regenerate faces)
       after the branch comes back together. */
     
  } else { 

    /*
      If the orientations DISagree, we're
      doing:

                              |                     |                      
	      +--------+      |                     |            +---------+        
              |        |      |                     |            |         |  e[2] (new indices)      
      --e[0]-<+    K   +-e[1]<--e[2]-<-  -->  -<e[0]--<--e0prev--+    K    +--<---------------------
              |        |      |                     |            |         |        
	      +--------+      |                     |            +---------+        
	      |                     |                      
			      
	      which is to say that first, the edges between e[0]-1 and e[1] need to be shifted DOWN 
	      by one index. Then we can splice a new edge in immediately BEFORE e[0] at the critical
	      crossing. The crossings are not going to get renumbered (except by resorting).

	      This is disgustingly like the previous code, except that we're swapping a lot of
	      heads and tails and prevs and nexts. Still, we'll try to persevere through it.

    */

    for(i=pd_previous_edge(pd,e[0]); i != e[2]; i = pd_previous_edge(pd,i)) { 
       
      out->edge[pd_previous_edge(pd,i)] = pd->edge[i];
      new_edge_idx[i] = pd_previous_edge(pd,i);
       
    }

    /* At this point, we've trashed e[2] (in the copy), and duplicated
       the edge after e[0]. 

       We now do the actual sewing, which as usual is the cleverest
       lines in the code. The idea is that we're going to insert the 
       crossing that used to be between e[1] and e[2] in between 
       e[0] and the previous edge to e[0], then go into the K stuff. */

    pd_idx_t e0prev = pd_previous_edge(pd,e[0]);

    out->edge[new_edge_idx[e0prev]].tail    = pd->edge[e[0]].tail;
    out->edge[new_edge_idx[e0prev]].tailpos = pd->edge[e[0]].tailpos;
     
    out->edge[new_edge_idx[e0prev]].head    = pd->edge[e[2]].head;
    out->edge[new_edge_idx[e0prev]].headpos = pd->edge[e[2]].headpos;

    out->edge[new_edge_idx[e[0]]].tail    = pd->edge[e[1]].tail;
    out->edge[new_edge_idx[e[0]]].tailpos = pd->edge[e[1]].tailpos;

    out->edge[new_edge_idx[e[1]]].tail    = pd->edge[e[2]].tail;
    out->edge[new_edge_idx[e[1]]].tailpos = pd->edge[e[2]].tailpos;

    out->cross[pd->edge[e[1]].tail].edge[pd->edge[e[1]].tailpos] = e[0]; /* Old edge number. */
    out->cross[pd->edge[e[1]].tail].edge[pd->edge[e[2]].headpos] = e0prev; /* Old edge number. */

    /* We're going to translate crossing numbers (and regenerate faces)
       after the branch comes back together. */
  }

  for(i=0;i<out->ncross;i++) { 

    for(j=0;j<4;j++) { 

      out->cross[i].edge[j] = new_edge_idx[out->cross[i].edge[j]];

    }

  }

  /* Now there's no need to do anything other than copy components,
     (which we've already done) because the same edges will still
     make up each component (and they remain in order). For the same
     reason, there's no need to resort the components, as they still
     all have the same number of edges. 

     However, the shift in edge numbers might resort the crossings,
     and the faces and hash are definitely different. 
  */

  pd_regenerate_crossings(out);
  pd_regenerate_faces(out);
  pd_regenerate_hash(out);

  assert(pd_ok(out));

  assert(out->ncross == pd->ncross);
  assert(out->nedges == pd->nedges);
  assert(out->nfaces == pd->nfaces);
  assert(out->ncomps == pd->ncomps);

  /* Now we do some housecleaning */

  free(new_edge_idx);

  return out;
   
}


pd_code_t *pd_connect_sum(pd_code_t *pdA, pd_idx_t edgeA,
			  pd_code_t *pdB, pd_idx_t edgeB)

/* 
  	  +----+	       +-------+	   
       	  |    |       	   +-<-|---+   |       	   
       	  |    |  edgeA	   |   |   |   |	   
       	+-|------<---+ 	   |   |   |   |       	   
      	| |    |     | 	   v   +-------+    	   
      	| +------>---+ 	   +-->----+   	       	   
       	+------+       	   edgeB	  
   				  
            pdA	       	       pdB

                      ||
                      vv

  	  +----+      	       +-------+	   
       	  |    |       	   +-<-|---+   |       	   
       	  |    |       	   |   |   |   |	   
       	+-|------<---------+   |   |   |       	   
      	| |    |       	       +-------+    	   
      	| +------>------------>----+   	       	   
       	+------+       	       		  
       	       	       	       	  
                  (output pd)  	  
		      	     
*/	


{
  assert(pd_ok(pdA));
  assert(pd_ok(pdB));
  pd_check_edge(SRCLOC,pdA,edgeA);
  pd_check_edge(SRCLOC,pdB,edgeB);

  /* The algorithm here is straightforward; allocate space, 
     copy crossings and edges, then reconnect edges. */

  pd_code_t *pd = pd_code_new(pdA->ncross + pdB->ncross);

  /* We need to actually generate an edge set for the connect
     summed guys here, so we can't just get the crossings right.
     Our scheme for crossing numbering is going to go like this:

     <all the crossings of A>, <all the crossings of B>

  */

  pd_idx_t *crossAtt = calloc(pdA->ncross,sizeof(pd_idx_t));
  pd_idx_t *crossBtt = calloc(pdB->ncross,sizeof(pd_idx_t));

  pd_idx_t i;

  for(i=0;i<pdA->ncross;i++) {

    crossAtt[i] = i; /* Crossing i in pdA -> crossing i in pd */

  }

  for(i=0;i<pdB->ncross;i++) {

    crossBtt[i] = i + pdA->ncross; /* Crossing i in pdB -> crossing i + pdA->ncross in pd */

  }

  /* The edge numbering scheme is more complex, because 
     we are connecting (in general) one of many components of A
     to one of many components of B. Let's call the components
     which are joined joinA and joinB.

     <edges of joinA <= edgeA> 
     <edges of joinB> 
     <edges of joinA > edgeA>

     <edges of all other components of A>
     <edges of all other components of B>
  */

  pd_idx_t joinA, joinB;
  pd_idx_t posA, posB;

  pd_component_and_pos(pdA,edgeA,&joinA,&posA);
  pd_component_and_pos(pdB,edgeB,&joinB,&posB);

  pd_idx_t comp,pos,nen=0;
  pd_idx_t *edgeAtt = calloc(pdA->nedges,sizeof(pd_idx_t));

  for(i=0;i<pdA->nedges;i++) { edgeAtt[i] = PD_UNSET_IDX; }

  for(pos=0;pos<=posA;pos++) {

    edgeAtt[pdA->comp[joinA].edge[pos]] = nen++;

  }

  nen += pdB->comp[joinB].nedges; /* Leave room for the compB stuff */

  for(;pos<pdA->comp[joinA].nedges;pos++) {

    edgeAtt[pdA->comp[joinA].edge[pos]] = nen++;

  }
  
  for(comp=0;comp<pdA->ncomps;comp++) {

    if (comp == joinA) { continue; } /* Don't run the loop for the join component */

    for(pos=0;pos<pdA->comp[comp].nedges;pos++) {

      edgeAtt[pdA->comp[comp].edge[pos]] = nen++;  /* New edge number */

    }

  }

  /* We now deal with the edges of B */

  pd_idx_t *edgeBtt = calloc(pdB->nedges,sizeof(pd_idx_t));
  for(i=0;i<pdB->nedges;i++) { edgeBtt[i] = PD_UNSET_IDX; }

  pd_component_t *compB = NULL;
  compB = &(pdB->comp[joinB]);
  nen = 0;
  
  for(pos=(1+posB);pos<(1+posB)+compB->nedges;pos++) {

    edgeBtt[compB->edge[pos%compB->nedges]] = (1+posA) + nen++;

    /* This is a tricky piece of code: the logic is that we're going
       to let edgeA be the edge "posA" in the final pd, and then splice
       in the edges in the join component (compB) of pdB.

       So we walk compB starting at (1 + posB) and wrapping
       around (by the mod operation) until we get back to the start.
       Since we've already used edge numbers 0..posA (inclusive), we
       assign edge numbers posA + 1...posA+1+pdB->comp[joinB].nedges.

    */

  }


  /* Now we skip-- we've already placed all the edges in pdA, and all
     the edges in compB. So we can set */

  nen = pdA->nedges + compB->nedges;

  /* ...and then walk the components of B, skipping the one we've
     already taken care of. */

  for(comp=0;comp<pdB->ncomps;comp++) {

    if (comp == joinB) { continue; }

    for(pos=0;pos<pdB->comp[comp].nedges;pos++) {

      edgeBtt[pdB->comp[comp].edge[pos]] = nen++;

    }

  }

  /* We now check the edge translation tables to make sure that 
     we haven't duplicated any indices in pd->edge and that we've 
     made all assignments */

  for(i=0;i<pdA->nedges;i++) {

    assert(edgeAtt[i] != PD_UNSET_IDX);

  }

  for(i=0;i<pdB->nedges;i++) {

    assert(edgeBtt[i] != PD_UNSET_IDX);

  }

  pd_idx_t *checkbuf = calloc(pdA->nedges + pdB->nedges,sizeof(pd_idx_t));
  pd_idx_t  cb_idx;
  
  for(i=0,cb_idx=0;i<pdA->nedges;i++,cb_idx++) {

    checkbuf[cb_idx] = edgeAtt[i];

  }

  for(i=0;i<pdB->nedges;i++,cb_idx++) {

    checkbuf[cb_idx] = edgeBtt[i];

  }

  qsort(checkbuf,pdA->nedges+pdB->nedges,sizeof(pd_idx_t),pd_idx_cmp);

  for(i=0;i<pdA->nedges+pdB->nedges;i++) {

    assert(checkbuf[i] == i);

  }

  free(checkbuf);  /* If we made it here, we pass */

  /* Now we're ready to copy crossing and edge data from pdA and pdB
     into pd; remember that we have to translate edge references (when
     copying the crossings) and crossing references (when copying the
     edges). */

  for(i=0;i<pdA->ncross;i++) {

    pd_crossing_t *pdcross = NULL;
    pdcross = &(pd->cross[crossAtt[i]]);

    pd_crossing_t *pdAcross = NULL;
    pdAcross = &(pdA->cross[i]);

    pd_idx_t j;

    for(j=0;j<4;j++) {

      pdcross->edge[j] = edgeAtt[pdAcross->edge[j]];

    }

    pdcross->sign = pdAcross->sign;

  }

  for(i=0;i<pdB->ncross;i++) {

    pd_crossing_t *pdcross = NULL;
    pdcross = &(pd->cross[crossBtt[i]]);

    pd_crossing_t *pdBcross = NULL;
    pdBcross = &(pdB->cross[i]);

    pd_idx_t j;

    for(j=0;j<4;j++) {

      pdcross->edge[j] = edgeBtt[pdBcross->edge[j]];

    }

    pdcross->sign = pdBcross->sign;

  }

  pd->ncross = pdA->ncross + pdB->ncross;

  /* We now translate the edge records */

  for(i=0;i<pdA->nedges;i++) {

    pd_edge_t *pdedge = NULL;
    pdedge = &(pd->edge[edgeAtt[i]]);

    pd_edge_t *pdAedge = NULL;
    pdAedge = &(pdA->edge[i]);

    pdedge->head = crossAtt[pdAedge->head];
    pdedge->headpos = pdAedge->headpos;

    pdedge->tail = crossAtt[pdAedge->tail];
    pdedge->tailpos = pdAedge->tailpos;

    
  }

  for(i=0;i<pdB->nedges;i++) {

    pd_edge_t *pdedge = NULL;
    pdedge = &(pd->edge[edgeBtt[i]]);

    pd_edge_t *pdBedge = NULL;
    pdBedge = &(pdB->edge[i]);

    pdedge->head = crossBtt[pdBedge->head];
    pdedge->headpos = pdBedge->headpos;

    pdedge->tail = crossBtt[pdBedge->tail];
    pdedge->tailpos = pdBedge->tailpos;
    
  }

  pd->nedges = pdA->nedges + pdB->nedges;

  /* Finally, we get to the actual connect sum operation,
     which involves joining the two edges edgeA and edgeB. 
     Since we've already translated these two to the new
     numbering scheme, we're going to do this in pd. */

  pd_edge_t *pdedgeA = NULL;
  pdedgeA = &(pd->edge[edgeAtt[edgeA]]);

  pd_edge_t *pdedgeB = NULL;
  pdedgeB = &(pd->edge[edgeBtt[edgeB]]);

  /*  tail	  head 	  tail	     head
       +  	   +       +  	      +	 
       |   edgeA   |       |  	      |	 
       +----->-----+  	   |edgeA     |edgeB	 
       	       	      =>   v   	      ^	 
       +----<------+  	   |  	      |	 
       |   edgeB   |  	   |  	      |	 
       +  	   +  	   +  	      +	 
      head     	  tail	  head	    tail 
  */   	  				 
  					 
  pd_edge_t swapB;
  swapB = *pdedgeB;

  pdedgeB->head = pdedgeA->head;
  pdedgeB->headpos = pdedgeA->headpos;

  pdedgeA->head = swapB.head;
  pdedgeA->headpos = swapB.headpos;

  /* Note that we also have to update the crossing records */

  pd->cross[pdedgeA->head].edge[pdedgeA->headpos] = edgeAtt[edgeA];
  pd->cross[pdedgeB->head].edge[pdedgeB->headpos] = edgeBtt[edgeB];

  /* Finally, we should be able to regenerate pd */

  pd_regenerate(pd);
  assert(pd_ok(pd));
  
  /* Housekeeping */

  free(edgeAtt); free(edgeBtt);
  free(crossAtt); free(crossBtt);

  return pd;
					 
}					 
