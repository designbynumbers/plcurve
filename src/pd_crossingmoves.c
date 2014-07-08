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

  for(i=0;i<nuniquedeletions;i++) { 

    assert(uniquedeletions[i] < nobj);

  }

  nuniquedeletions = j;
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
     edges e[0] and e[1]. */

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
  
  /* We now need to compact and reorder the components, keeping in 
     mind that the edge indices have now changed. This will use the 
     transform_and_compact_indices api:
 
     void pd_transform_and_compact_indices(pd_idx_t *source, pd_idx_t nobj, 
				      pd_idx_t *transform,
				      pd_idx_t **target)
  */
  
  for(i=0;i<pd->ncomps;i++) { 

    pd_transform_and_compact_indices(pd->comp[i].edge,pd->comp[i].nedges,new_edge_idx,
				     &(outpd->comp[i].edge),&(outpd->comp[i].nedges));
    outpd->comp[i].tag = pd->comp[i].tag;

  }

  outpd->ncomps = pd->ncomps;
  qsort(outpd->comp,outpd->ncomps,sizeof(pd_component_t),pd_component_cmp);

  /* The new pd now has crossings, edges, and components which make sense. */
  /* It's ok to just regenerate the faces. */

  pd_regenerate_faces(outpd);
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
         B           B      

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
  really different in each case. */

  pd_idx_t F1, F2;

  /* We're going to use pd_face_and_pos to identify the faces. */
  /* F1 should be the negface of eA[0] and of eA[2] */

  pd_idx_t pf0,nf0,pfp0,nfp0;
  pd_idx_t pf2,nf2,pfp2,nfp2;

  pd_face_and_pos(pd,eA[0],&pf0,&pfp0,&nf0,&nfp0);
  pd_face_and_pos(pd,eA[2],&pf2,&pfp2,&nf2,&nfp2);
  
  F1 = nf0; F2 = nf2;

  if (F1 != F2) { /* This is the conventional case where we output a single PD code */

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
		       (void *)(out->cross),
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
		       (void *)(out->edge),
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

    out->cross[new_crossing_idx[pd->edge[eA[2]].head]].edge[pd->edge[eA[2]].headpos] = eA[0]; 
    out->edge[new_edge_idx[eA[0]]].head = pd->edge[eA[2]].head;
    out->edge[new_edge_idx[eA[0]]].headpos = pd->edge[eA[2]].headpos;

    out->cross[new_crossing_idx[pd->edge[eB[2]].head]].edge[pd->edge[eB[2]].headpos] = eB[0]; 
    out->edge[new_edge_idx[eB[0]]].head = pd->edge[eB[2]].head;
    out->edge[new_edge_idx[eB[0]]].headpos = pd->edge[eB[2]].headpos;

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
       components (and want tags to persist). So we'll have
       to do the same edge-deletion procedure that we did 
       for the R1 move. 

       This will use the 
       transform_and_compact_indices api:
 
       void pd_transform_and_compact_indices(pd_idx_t *source, pd_idx_t nobj, 
                                             pd_idx_t *transform,
				             pd_idx_t **target)
    */
  
    for(i=0;i<pd->ncomps;i++) { 

      pd_transform_and_compact_indices(pd->comp[i].edge,pd->comp[i].nedges,new_edge_idx,
				       &(out->comp[i].edge),&(out->comp[i].nedges));
      out->comp[i].tag = pd->comp[i].tag;

    }
    
    out->ncomps = pd->ncomps;
    qsort(out->comp,out->ncomps,sizeof(pd_component_t),pd_component_cmp);
    
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
		       (void *)(pdA->cross),
		       &new_crossing_idx_A,
		       &(pdA->ncross));
    pdA->MAXVERTS = pdA->ncross;

    /* Handle edges of A, remembering that we get A's edges by DELETING B's */

    free(pdA->edge); pdA->edge = NULL; pdA->nedges = pdA->MAXEDGES = 0;
    pd_compacting_copy((void *)(pd->edge),sizeof(pd_edge_t),(size_t)(pd->nedges),
		       B_nedges,B_edges,
		       (void *)(pdA->edge),
		       &new_edge_idx_A,
		       &(pdA->nedges));
    pdA->MAXEDGES = pdA->nedges;

    /* Handle crossings of B, remembering we get B's edges by DELETING A's */
      
    free(pdB->cross); pdB->cross = NULL; pdB->ncross = pdB->MAXVERTS = 0; 
    pd_compacting_copy((void *)(pd->cross),sizeof(pd_crossing_t),(size_t)(pd->ncross),
		       A_ncross,A_crossings,
		       (void *)(pdB->cross),
		       &new_crossing_idx_B,
		       &(pdB->ncross));
    pdB->MAXVERTS = pdB->ncross;

    /* Handle edges of B, remembering that we get B's edges by DELETING A's */

    free(pdB->edge); pdB->edge = NULL; pdB->nedges = pdB->MAXEDGES = 0;
    pd_compacting_copy((void *)(pd->edge),sizeof(pd_edge_t),(size_t)(pd->nedges),
		       A_nedges,A_edges,
		       (void *)(pdB->edge),
		       &new_edge_idx_B,
		       &(pdB->nedges));
    pdB->MAXEDGES = pdB->nedges;

    /* Now we do the clever sewing bit, keeping in mind that on the left we 
       use NEW indices because we're in a child pd, and on the right we use
       OLD indices because we're referring to the parent pd. */

    pdA->cross[new_crossing_idx_A[pd->edge[eA[2]].head]].edge[pd->edge[eA[2]].headpos] = eA[0]; /* Since eA[2] was deleted, splice in eA[0]. */
    pdA->edge[new_edge_idx_A[eA[0]]].head = pd->edge[eA[2]].head;    /* Reset the head of edge eA[0] */
    pdA->edge[new_edge_idx_A[eA[0]]].headpos = pd->edge[eA[2]].headpos;

    pdB->cross[new_crossing_idx_B[pd->edge[eB[2]].head]].edge[pd->edge[eB[2]].headpos] = eB[0]; /* Since eB[2] was deleted, splice in eB[0]. */
    pdB->edge[new_edge_idx_B[eB[0]]].head = pd->edge[eB[2]].head;    /* Reset the head of edge eB[0] */
    pdB->edge[new_edge_idx_B[eB[0]]].headpos = pd->edge[eB[2]].headpos;

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

    /* Aaand for child diagram B.... */

    for(i=0;i<pdB->ncross;i++) { 

      for(j=0;j<4;j++) { 

	pdB->cross[i].edge[j] = new_edge_idx_B[pdB->cross[i].edge[j]];

      }

    }

    for (i=0;i<pdB->nedges;i++) { 

      pdB->edge[i].head = new_crossing_idx_B[pdB->edge[i].head];
      pdB->edge[i].tail = new_crossing_idx_B[pdB->edge[i].tail];

    }

    /* Ok, now we need to copy components into the child pd's appropriately. 
       This takes a little thought, because we want the component TAGS to 
       travel appropriately. */

    pdA->ncomps = 0; pdB->ncomps = 0;

    for(i=0;i<pd->ncomps;i++) { 

      if (compstate[i] == A_walked) { 

	pd_transform_and_compact_indices(pd->comp[i].edge,pd->comp[i].nedges,
					 new_edge_idx_A,
					 &(pdA->comp[pdA->ncomps].edge),
					 &(pdA->comp[pdA->ncomps].nedges));
	pdA->comp[pdA->ncomps].tag = pd->comp[i].tag;
	pdA->ncomps++;

      } else if (compstate[i] == B_walked) { 

	pd_transform_and_compact_indices(pd->comp[i].edge,pd->comp[i].nedges,
					 new_edge_idx_B,
					 &(pdB->comp[pdB->ncomps].edge),
					 &(pdB->comp[pdB->ncomps].nedges));
	pdB->comp[pdB->ncomps].tag = pd->comp[i].tag;
	pdB->ncomps++;

      } else {

	printf("compstate[%d] == %d, which is not A_walked (%d) or B_walked (%d)\n"
	       "suspect memory corruption bug, but your guess is as good as mine\n",
	       i,compstate[i],A_walked,B_walked);
	exit(1);

      }

    }

    qsort(pdA->comp,pdA->ncomps,sizeof(pd_component_t),pd_component_cmp);
    qsort(pdB->comp,pdB->ncomps,sizeof(pd_component_t),pd_component_cmp);
    
    /* Ok, at this point the children have crossings, edges, and components. */
    /* We can regenerate the crossings, and the faces. */

    pd_regenerate_crossings(pdA); pd_regenerate_crossings(pdB);
    pd_regenerate_faces(pdA); pd_regenerate_faces(pdB);

    /* We can now do some counting to make sure things are sane: */

    assert(pdA->ncomps + pdB->ncomps == pd->ncomps);
    assert(pdA->ncross + pdB->ncross == pd->ncross - 2);
    assert(pdA->nedges + pdB->nedges == pd->nedges - 4);
    assert(pdA->ncross - pdA->nedges + pdA->nfaces == 2);
    assert(pdB->ncross - pdB->nedges + pdB->nfaces == 2);

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


