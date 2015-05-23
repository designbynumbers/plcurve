/* 

   pd_invariants.c: This is a collection of more-or-less experimental code
   for computing the Arnol'd invariants of pd codes, and some related quantities
   like the number of interlaced crossings. 

*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
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

#ifdef HAVE_STRING_H
  #include<string.h>
#endif

#ifdef HAVE_STDLIB_H
  #include<stdlib.h>
#endif

#ifdef HAVE_STDBOOL_H
  #include<stdbool.h>
#endif

#ifdef HAVE_MATH_H
   #include<math.h>
#endif

#include<plcTopology.h>
#include<pd_container.h>

#include<pd_multidx.h>
#include<pd_dihedral.h>
#include<pd_perm.h>
  
#include<pd_isomorphisms.h>
#include<pd_storage.h>
#include<pd_invariants.h>

int int_cmp (const void *a, const void *b)
{
  const int *ia = (const int *) a;
  const int *ib = (const int *) b;

  return (*ia > *ib) - (*ia < *ib);
}

unsigned int pdint_interlaced_comp(pd_code_t *pd,pd_idx_t comp)
/* 
   Counts the number of interlaced crossings on one component of pd.
   The basic algorithm is simple: first we compute a Gauss code where
   we record the number of each crossing in order around the
   component, for example

   0 2 0 1 2 

   Some crossings (with other components) will occur only
   one. Crossings of this component with itself should occur twice in
   the list. We then scan to compute frequency of each crossing in the
   list, then eliminate crossings that occur only once.

   0 2 0 2

   We then sort the remaining crossings and map to indices in 1..k
   (we're going to use these to hash the pairs of crossings, and don't
   want them to be spread out).

   0 1 0 1

   Now we're going to make a little machine. For each crossing, we
   have two states: ACTIVE and INACTIVE, together with a buffer of
   crossings with smaller numbers, each associated with a tag in the
   form NONE, ONCE, TWICE. Initially, all crossings have state
   INACTIVE and the buffer is filled with NONE.

   We then read from left to right:

   Toggle the state of the current crossing (ACTIVE <-> INACTIVE).
   
   Look at all strictly larger crossings, and for each one that's 
   ACTIVE, find the entry of this crossing, and increment the tag
   NONE -> ONCE -> TWICE. If the tag is already TWICE, abort: 
   something's corrupted.

   At this point the number of interlaced crossings is the total 
   number of ONCE tags in all the crossings. Scan through the state
   for each crossing and count them. 

   Ok, I wrote it out. :) 

   Now let's get a little more clever. We'll start with the same
   Gauss code.

   0 2 0 1 2 

   Some crossings (with other components) will occur only
   one. Crossings of this component with itself should occur twice in
   the list. We then scan to compute frequency of each crossing in the
   list, then eliminate crossings that occur only once. We can rotate 
   the list (also in linear time) so that the lowest numbered crossing
   happens first. For a reason that's about to be clear, we number
   starting with one:

   1 2 1 2

   We now change sign so that each crossing occurs once with -
   (leftmost) and once with + (rightmost). By definition, the 
   first entry in the list is going to be -1. We'll get something
   like

   -1 -2 +1 +2

   We're basically going to bubblesort the list now, keeping track of
   certain switches, but not others. 

   and so forth. The swap rules are obvious: continually scan from left
   to right. Suppose we have p[i] and p[i+1]. If abs(p[i+1]) < abs(p[i]),
   we do the swap. Otherwise we don't. We'll never swap two things with 
   equal absolute value, so in the end the order is:

   -1 +1 -2 +2 -3 +3 ...

   There are couple of options for each pair, which we will call xX
   and yY. Let's assume wlog that X < Y. Then we can enumerate the
   options, and there are really only six. We can arrange these in
   clever way, along with the swaps that get them there
  	   
		 (yYxX)
	       	    |  	   
                 Yx-->xY  
	      	    |	
    YX-->XY    --(yxYX)---   yx-->xy  
              /   	  \ 
             / 	       	   \
          (yxXY)         (xyYX)
   	     \ 	       	   /   	    
      	      \     	  / 	    
     yx-->xy   --(xyXY)---   YX-->XY
       	       	    |
                 yX-->Xy
		    |
                 (xXyY)

   We notice that among these, only yxYX and xyXY are linked; these
   are at distance 1 and 3 from the root; the others are at distance
   0, 2, or 4. 

   So we just need to odd/even the number of swaps. In practice, this
   looks like the following. We'll keep a record for every
   pair. Whenever we make a swap, we'll increment nswaps[i][j]. When
   we're done, we'll sum the remainder (mod 2) over the array to get
   the total number of interlaced pairs.

*/
{
  /* Phase 1. Compute the Gauss code, converting crossings to 1-indexed. */

  int *gausscode = calloc(pd->comp[comp].nedges,sizeof(int));
  assert(gausscode != NULL);
  
  pd_idx_t i;

  for(i=0;i<pd->comp[comp].nedges;i++) {

    gausscode[i] = (int)(pd->edge[pd->comp[comp].edge[i]].head) + 1;

  }

  /* Phase 2. Eliminate crossings that occur only once, 
     assign -+ tags, & renumber */

  int *nseen = calloc(pd->ncross+1,sizeof(int));
  assert(nseen != NULL);

  for(i=0;i<pd->comp[comp].nedges;i++) {

    nseen[gausscode[i]]++;
    
    if (nseen[gausscode[i]] == 1) {

      /* The first time we see it, make it - */

      gausscode[i] *= -1;

    }

  }

  free(nseen); nseen = NULL;

  int *crossingspresent
    = calloc(pd->ncross,sizeof(int)); /* This is a max */
  assert(crossingspresent != NULL);

  pd_idx_t cpsize = 0;

  for(i=0;i<pd->comp[comp].nedges;i++) {

    if (gausscode[i] > 0) { /* This is the SECOND time it occurs. */
    
      crossingspresent[cpsize++] = gausscode[i];

    }

    assert(cpsize <= pd->ncross);

  }

  /* Now it's entirely possible to eliminate ALL the crossings from 
     this list (cf. the (2,q) torus link). We don't want to continue
     at this point! (There are definitely NO interlacements.) */

  if (cpsize == 0) {

    free(crossingspresent);
    free(gausscode);
    
    return 0;

  }
  
  /* We now have a list of only the (+1'd) crossing indices which 
     occur TWICE in the list. We want to use this to translate old
     indices to new ones in 1..gcsize. */

  qsort(crossingspresent,cpsize,sizeof(int),int_cmp);

  int *newcrossinglabel = calloc(pd->ncross+1,sizeof(int));

  for(i=0;i<cpsize;i++) {

    newcrossinglabel[crossingspresent[i]] = i+1;

  }

  free(crossingspresent); crossingspresent = NULL;

  /* The newcrossinglabel is order-preserving. We're now
     going to build the augmented Gauss code, with new labels
     and -+ coding in place. Since every crossing left occurs
     exactly twice, this code will have length 2*cpsize. */

  int *augmentedgausscode = calloc(2*cpsize,sizeof(int));
  int agc_idx = 0;

  for(i=0;i<pd->comp[comp].nedges;i++) {

    /* Recall that "abs" is the INTEGER absolute value function */

    if (newcrossinglabel[abs(gausscode[i])] > 0) { /* This guy occurs twice */
  
      augmentedgausscode[agc_idx++] = ((gausscode[i] > 0) ? 1 : -1) *
	(newcrossinglabel[abs(gausscode[i])]);

    }

  }

  assert(agc_idx == 2*cpsize);
  free(gausscode); gausscode = NULL;
  free(newcrossinglabel); newcrossinglabel = NULL;

  /* Phase 3. Build the matrix to record swaps, and do the bubble sort. */
  
  unsigned int **nswaps = calloc((cpsize+1),sizeof(int *));
  for(i=2;i<=cpsize;i++) {
    nswaps[i] = calloc(i,sizeof(unsigned int));
  }

  bool madeswap;
  int  safety_count = 0;

  do {

    for(i=0,madeswap=false;i<2*cpsize-1;i++) {

      if (abs(augmentedgausscode[i]) > abs(augmentedgausscode[i+1])) {

	int swap;
	swap = augmentedgausscode[i];
	augmentedgausscode[i] = augmentedgausscode[i+1];
	augmentedgausscode[i+1] = swap;

	nswaps[abs(augmentedgausscode[i+1])][abs(augmentedgausscode[i])]++;
	/* nswaps[larger][smaller]++ */

	madeswap = true;
	
      }

    }

    safety_count++;
    assert(safety_count < 2*cpsize+2);

  } while ( madeswap );
    
  /* Phase 4. Count (mod 2) the number of swaps in each bin. */

  unsigned int interlacements = 0;
  int j;

  for(i=2;i<=cpsize;i++) {

    for(j=1;j<i;j++) {

      interlacements += nswaps[i][j] % 2;

    }

  }

  /* Phase 5. Housekeeping. */

  for(i=2;i<=cpsize;i++) {

    free(nswaps[i]);

  }

  free(nswaps);
  free(augmentedgausscode);

  return interlacements;

}


unsigned int *pd_interlaced_crossings(pd_code_t *pd)
/* Returns an array, pd->ncomps long, counting number of crossings
   which occur in the order ABAB along the component. */
{
  unsigned int *interlacements = calloc(pd->ncomps,sizeof(unsigned int));
  pd_idx_t i;

  for(i=0;i<pd->ncomps;i++) {

    interlacements[i] = pdint_interlaced_comp(pd,i);

  }

  return interlacements;

}
