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
   in the buffer "deletions" to target buf. We allocate "target_buf" ourselves,
   and return in the target_idx a buffer of size nobj so that

   target_idx[index in source] = index of copy, if the object was copied
                               = PD_UNSET_IDX, if the object was deleted

   The basic idea is that we're going to sort the buffer of deletions, and 
   then read through while copying to the target buffer, incrementing the 
   deletions buffer when we find one of the elements to delete. 

   This is linear time, and most importantly, we only need to write (and debug)
   it once. We should throw an error if EVERYTHING is being deleted.

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

  for(i=0;i<nuniquedeletions;i++) { 

    assert(uniquedeletions[i] < nobj);

  }

  assert(nuniquedeletions != nobj);
    
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
    
    if (!all_deletions_used) { /* There's a chance we might not copy this guy */
    
      if (source_i == uniquedeletions[deletions_i]) { /* We're not copying this guy. */
	
	tidx[source_i] = PD_UNSET_IDX;
	deletions_i++;
	all_deletions_used = (deletions_i >= nuniquedeletions);
	
      }
    
    }

    /* We're not at the current deletion, or there are no deletions left. */
    /* So we're definitely copying THIS guy. We need to do the copy with */
    /* memcpy because we're just copying bytes of memory-- we don't know */
    /* the data type involved. */

    memcpy(target + target_i*obj_size,source + source_i*obj_size,obj_size);
    tidx[source_i] = target_i;
    target_i++;
    
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
    if (transform[i] == PD_UNSET_IDX) { 
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
    for(i=0;i<nobj;i++) { if (transform[i] == PD_UNSET_IDX) { deletions[ndeletions++] = i; } }

    pd_compacting_copy((void *)(source),sizeof(pd_idx_t),(size_t)(nobj),ndeletions,deletions,
		       (void **)(target),&target_idx,ntarget);
    
  }

  free(target_idx); /* We don't need it. */

  /* Now do the transform step. */
  
  for(i=0;i<*ntarget;i++) { 

    (*target)[i] = transform[(*target)[i]];

  }

  /* And we're done. */

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
  /* in the crossing. */

  pd_idx_t e[3];  
  pd_pos_t ep[3];
  pd_idx_t i,j;
  bool found_flag = false;

  for(i=0;i<4 && !found_flag;i++) { 
    
    if (pd->cross[cr].edge[i] == pd->cross[cr].edge[(i+1)%4]) { 

      ep[0] = (i+1)%4; ep[1] = i; ep[2] = (i+2)%4;
      
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

    pd_error(SRCLOC,
	     "crossing %CROSS of %PD is not a valid location\n"
	     "for loop deletion because it belongs to a one-crossing\n"
	     "figure-8 diagram (or something is corrupted!).\n",pd,cr);
    exit(1);

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
  
  free(outpd->cross); /* We need to free the previously allocated buffer, as we're replacing it. */

  void **target_var = (void **)(&(outpd->cross));
  pd_idx_t *target_idx;
  pd_idx_t *ntarget = &(outpd->ncross);

  pd_compacting_copy(source,obj_size,nobj,
		     ndeletions,deletions,
		     target_var,&(target_idx),ntarget);

  /* We have now created a new buffer of crossings, and kept track 
     of the numberings associated to them with target_idx. */

  /* Our next move is to replace the reference to the old e[2]
     at the head of the old e[2] with a reference to e[0], since 
     the new e[0] will sew the tail of the old e[0] onto the head
     of the old e[2], eliminating the crossing in between.... */

  outpd->cross[pd->edge[e[2]].head].edge[pd->edge[e[2]].headpos] = e[0];

  /* We now make a copy of the old edge buffer, in which the crossing
     references refer to the NEW numbering scheme. */

  pd_edge_t *working_edges;
  working_edges = calloc(pd->nedges,sizeof(pd_edge_t));
  
  for(i=0;i<pd->nedges;i++) { 

    working_edges[i].head = target_idx[pd->edge[i].head];
    working_edges[i].headpos = pd->edge[i].headpos;

    working_edges[i].tail = target_idx[pd->edge[i].tail];
    working_edges[i].tailpos = pd->edge[i].tailpos;

  }

  /* Of course, in this blind copy, we're set head of e[0], both tail
     and head of e[1], and the tail of e[2] to PD_UNSET_IDX, since
     these all referred to the deleted crossing. We needn't worry
     about e[1] and e[2] which are simply going away, but the head of
     e[0] should be fixed. */

  working_edges[e[0]].head = target_idx[pd->edge[e[2]].head];
  working_edges[e[0]].headpos = target_idx[pd->edge[e[2]].headpos];

  /* Our goal now is to make a compacting copy from the working_edges
     buffer into the edge buffer of the new pd code. */

  /* void pd_compacting_copy(void *source, size_t obj_size, size_t nobj, 
                             pd_idx_t ndeletions, pd_idx_t *deletions,
			     void **target_var,
			     pd_idx_t **target_idx,
			     pd_idx_t *ntarget) */

  /* We're using the compacting_copy api here. */

  source = working_edges;
  obj_size = sizeof(pd_edge_t);
  nobj = pd->nedges;

  ndeletions = 2;
  pd_idx_t newdeletions[2] = {e[1],e[2]};
  
  free(outpd->edge); /* We need to free the previously allocated buffer, as we're replacing it. */

  target_var = (void **)(&(outpd->edge));
  ntarget = &(outpd->nedges);
  
  free(target_idx); /* Again, we're not going to use the old crossing map again, so we can free it. */

  pd_compacting_copy(source,obj_size,nobj,
		     ndeletions,newdeletions,
		     target_var,&(target_idx),ntarget);

  /* We now have the problem that the new buffer of crossings involves the OLD edge numbers. */

  for(i=0;i<outpd->ncross;i++) { 

    for(j=0;j<4;j++) { 

      outpd->cross[i].edge[j] = target_idx[outpd->cross[i].edge[j]];

    }
  
  }

  /* We've now translated the new pd to have correct crossings and
     edges, which should be mutually consistent. Unfortunately, the
     crossings are out of order, so we first make sure that they are
     canonically ordered. This actually involves yet another
     translation of the crossing references in the edge data, but
     there's nothing to be done about that right now! */

  pd_regenerate_crossings(outpd);
  
  /* We now need to compact and reorder the components, keeping in 
     mind that the edge indices have now changed. This will use the 
     transform_and_compact_indices api:
 
     void pd_transform_and_compact_indices(pd_idx_t *source, pd_idx_t nobj, 
				      pd_idx_t *transform,
				      pd_idx_t **target)
  */
  
  for(i=0;i<pd->ncomps;i++) { 

    pd_transform_and_compact_indices(pd->comp[i].edge,pd->comp[i].nedges,target_idx,
				     &(outpd->comp[i].edge),&(outpd->comp[i].nedges));
    outpd->comp[i].tag = pd->comp[i].tag;

  }

  outpd->ncomps = pd->ncomps;
  qsort(outpd->comp,outpd->ncomps,sizeof(pd_component_t),pd_component_cmp);

  /* The new pd now has crossings, edges, and components which make sense. */
  /* It's ok to just regenerate the faces. */

  pd_regenerate_faces(outpd);

  /* Finally, we check that we survived all this. */

  assert(pd_ok(outpd));
  
  /* Housecleaning */
  
  free(target_idx);
  free(working_edges);

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
      
  
  
  
  
  
