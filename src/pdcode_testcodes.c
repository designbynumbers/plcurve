/* 

   pdcode_testcodes.c : Code to generate standard, valid pd codes for testing purposes. 

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

#include<plcTopology.h>

#include<pd_multidx.h>
#include<pd_dihedral.h>
#include<pd_perm.h>
  
#include<pd_isomorphisms.h>

pd_code_t *pd_build_twist_knot(pd_idx_t n) 

/* Generate a twist knot diagram (scanned documentation in data/pd_twist_knot_doc.pdf) */

{
  pd_code_t *pd = NULL;

  /* Make sure that the number of verts is ok */

  if (n<1) { 

    pd_error(SRCLOC,"Can't generate pd_twist_knot with %d (< 1) twists.\n",pd,n);
    exit(1);

  }
  
  /* Now generate crossings. */

  pd = pd_code_new(n+2); /* Allocate just enough space */
  assert(pd != NULL); 
  pd->ncross = n+2;

  if (pd->ncross == 3) { /* This case is non-generic: hard code it */

    pd->cross[0] = pd_build_cross(0,5,3,1);
    pd->cross[1] = pd_build_cross(0,1,2,4);
    pd->cross[2] = pd_build_cross(2,3,5,4);

  } else { /* The generic case (at least one "internal" twist) */
    
    pd->cross[0] = pd_build_cross(0,1,2,(2*pd->ncross)-2);
    pd->cross[1] = pd_build_cross(0,(2*pd->ncross)-1,3,1);
  
    pd_idx_t i;
    
    for(i=1;i<n+1;i++) {
      
      pd->cross[i+1] = pd_build_cross((2*i),(2*i)+1,(2*i)+3,(2*i)+2);
      
    }

  }

  /* We don't have edges yet, so we can't assign crossing signs. */
    
  pd_regenerate(pd);

  /* Now we have an edge list and can assign crossing signs. It turns
     out that this is trivial as all the crossings have positive sign. */

  pd_idx_t k;
  for(k=0;k<pd->ncross;k++) { 

    pd->cross[k].sign = PD_POS_ORIENTATION;

  }

  assert(pd_ok(pd));

  return pd;
}

pd_code_t *pd_build_torus_knot(pd_idx_t p, pd_idx_t q) 

/*
   q odd: the strand passes through all q crossings before coming back to this one
          for the first time, when it's on the other strand at this crossing, and 
	  passes through all q crossings twice before coming back to this crossing
	  on the original strand.

   +------------------------------------------------------------------------------------------------------+
   |                                                                                                      |
   +-- q-1----\    /---0---\   /--q+1--\    /-(2k+q-1)%2q-\    /-----2k-------\    /----(2k+q+1)%2q--\   /            
               \0 /         \1/         \  /               \2k/                \ (2k+1)               \
	        \            \           ..                 /                    /                     ..
	       / \          / \         /  \               /  \                 / \                   /  \
   +-- 2q-1---/   \----q--/    \---1---/    \---(2k-1)%2q-/    \--(2k+q)%(2q)--/   \----(2k+1)-------/    \      
   |                                                                                                       |
   +-------------------------------------------------------------------------------------------------------+

   q even: the strand passes through all q crossings before coming back, but this time it comes
           back on the SAME strand and closes. The other component starts at edge q and continues to 2q-1.
           Basically, this looks exactly the same, EXCEPT at crossing 0.

   +------------------------------------------------------------------------------------------------------+
   |                                                                                                      |
   +-- 2q-1---\    /---0---\   /--q+1--\    /-(2k+q-1)%2q-\    /-----2k-------\    /----(2k+q+1)%2q--\   /            
               \0 /         \1/         \  /               \2k/                \ (2k+1)               \
	        \            \           ..                 /                    /                     ..
	       / \          / \         /  \               /  \                 / \                   /  \
   +--  q-1---/   \----q--/    \---1---/    \---(2k-1)%2q-/    \--(2k+q)%(2q)--/   \----(2k+1)-------/    \      
   |                                                                                                       |
   +-------------------------------------------------------------------------------------------------------+

*/

{
  pd_code_t *pd = NULL;

  /* Make sure that the number of verts is ok */

  if (p != 2) {

    pd_error(SRCLOC,"pd_torus_knot can only generate (2,q) torus knots (called with (%d,%d)).\n",pd,
	     p,q);
    exit(1);

  }

  if (q < 2) { 
    
    pd_error(SRCLOC,"pd_torus_knot can only generate (2,q) torus knots with q >= 2 (called with (%d,%d))\n",pd,
	     p,q);
    exit(1);

  }
    
  /* Now generate crossings. */

  pd = pd_code_new(q*(p-1));
  assert(pd != NULL); 

  pd->ncross = q; /* Note that this assumes that p = 2 */

  if (q%2 != 0) { /* q odd */

    pd->cross[0] = pd_build_cross(0,q-1,2*q-1,q);

  } else { /* q even */

    pd->cross[0] = pd_build_cross(0,2*q-1,q-1,q);
 
  }
  
  pd_idx_t i,k;

  for(i=1;i<q;i++) {

    if (i%2 == 0) { /* even numbered crossing */

      k = i/2;
      pd->cross[2*k] = pd_build_cross(2*k,(2*k+q-1)%(2*q),(2*k-1)%(2*q),(2*k+q)%(2*q));

    } else { /* odd numbered crossing */

      k = (i-1)/2;
      pd->cross[2*k+1] = pd_build_cross((2*k+q+1)%(2*q),2*k,(2*k+q)%(2*q),(2*k+1)%(2*q));
      
    }
  
  }
  
  /* Regenerating should build edges for the first time.  It doesn't
     make any sense to try to establish crossing signs until this
     happens. */

  pd_regenerate(pd);

  /* All crossings have positive sign. */

  for(i=0;i<pd->ncross;i++) { 

    pd->cross[i].sign = PD_POS_ORIENTATION;

  }

  assert(pd_ok(pd));

  return pd;

}

pd_code_t *pd_build_simple_chain(pd_idx_t n) 
/* Build a simple chain of n links (n > 0) */

{
  pd_code_t *pd = NULL;

  /* First, check the size of the diagram. */

  if (n == 0) { 

    pd_error(SRCLOC,"Can't build simple chain with 0 links.\n",NULL);
    exit(1);

  }

  if (n == 1) { 

    return pd_build_unknot(0);

  } 

  if (n == 2) { 

    pd = pd_code_new(2);

    pd->ncross = 2; 
    pd->cross[0] = pd_build_cross(0,2,1,3);
    pd->cross[1] = pd_build_cross(0,3,1,2);
    pd->cross[0].sign = PD_POS_ORIENTATION;
    pd->cross[1].sign = PD_POS_ORIENTATION;
    
    pd->nedges = 4;
    
    pd->edge[0].head = 1; pd->edge[0].headpos = 0;
    pd->edge[0].tail = 0; pd->edge[0].tailpos = 0;
    
    pd->edge[1].head = 0; pd->edge[1].headpos = 2;
    pd->edge[1].tail = 1; pd->edge[1].tailpos = 2;

    pd->edge[2].head = 0; pd->edge[2].headpos = 1;
    pd->edge[2].tail = 1; pd->edge[2].tailpos = 3;

    pd->edge[3].head = 1; pd->edge[3].headpos = 1;
    pd->edge[3].tail = 0; pd->edge[3].tailpos = 3;
    
    pd_regenerate_comps(pd);
    pd_regenerate_faces(pd);
    pd_regenerate_hash(pd);

    assert(pd_ok(pd));
    return pd;

  }

  pd = pd_code_new(2*(n-1));
  assert(pd != NULL);
  pd->ncross = 2*(n-1);

  if (n == 3) { /* This is a nongeneric case; fill in by hand. */

    pd->cross[0] = pd_build_cross(0,2,1,5);
    pd->cross[1] = pd_build_cross(0,5,1,4);
    pd->cross[2] = pd_build_cross(2,6,3,7);
    pd->cross[3] = pd_build_cross(3,6,4,7);

  } else {

    pd->cross[0]       = pd_build_cross(0,1,5,4);
    pd->cross[1]       = pd_build_cross(0,6,5,1);

    pd_idx_t i;  /* The "interior" components are numbered i = 1, ... , n-2 */

    for(i=1;i<(n-1);i++) {

      pd->cross[(2*i)]   = pd_build_cross((4*i),(4*i)+3,(4*i)+5,(4*i)+4);
      pd->cross[(2*i)+1] = pd_build_cross((4*i)+2,(4*i)+6,(4*i)+5,(4*i)+3);

    }

    pd->cross[(2*n)-4] = pd_build_cross(2,3,(4*n)-8,(4*n)-5);
    pd->cross[(2*n)-3] = pd_build_cross(2,(4*n)-5,(4*n)-6,3);

  }

  /* We now call regenerate to build edges and faces. */

  pd_regenerate(pd);

  /* Now that we have edges, it makes sense to assign
     crossing orientations. */

  pd_idx_t k;
  for(k=0;k<pd->ncross;k++) { 

    pd->cross[k].sign = PD_POS_ORIENTATION;

  }


  assert(pd_ok(pd));

  return pd;

}

pd_code_t *pd_build_unknot(pd_idx_t n) 

/* An n-crossing unknot diagram in the form  


       +--<--+       +--<--+            +-<--+      
       |     |       |     |                 |      
       |     |       |     |                 |      
+------|->-----------|->-------+  ...    +----->---+
|      |     |       |     |                 |     |
|      |     |       |     |                 ^     |
+--<---+     +---<---+     +---+             +-----+

with all crossings positive. This generates the 
case with n == 0 and n == 1 correctly. 

*/

{
  pd_code_t *pd = NULL;
  
  if (n == 0) { 

    pd = pd_code_new(2); /* Make room for 2, even though we won't use them. */
    
    pd->ncross = 0;
    pd->nedges = 1;
    pd->nfaces = 2;
    pd->ncomps = 1;

    pd->edge[0].head = PD_UNSET_IDX; pd->edge[0].headpos = PD_UNSET_POS;
    pd->edge[0].tail = PD_UNSET_IDX; pd->edge[0].tailpos = PD_UNSET_POS;
    
    pd->face[0].nedges = 1;
    
    pd->face[0].edge = calloc(1,sizeof(pd_idx_t));
    assert(pd->face[0].edge != NULL);
    pd->face[0].edge[0] = 0;
    
    pd->face[0].or = calloc(1,sizeof(pd_or_t));
    assert(pd->face[0].or != NULL);
    pd->face[0].or[0] = PD_POS_ORIENTATION;

    pd->face[1].nedges = 1;
    
    pd->face[1].edge = calloc(1,sizeof(pd_idx_t));
    assert(pd->face[1].edge != NULL);
    pd->face[1].edge[0] = 0;
    
    pd->face[1].or = calloc(1,sizeof(pd_or_t));
    assert(pd->face[1].or != NULL);
    pd->face[1].or[0] = PD_NEG_ORIENTATION;

    pd->comp[0].nedges = 1;
    pd->comp[0].edge = calloc(1,sizeof(pd_idx_t));
    assert(pd->comp[0].edge != NULL);
    pd->comp[0].edge[0] = 0;
    pd->comp[0].tag = 'A';

    assert(pd_ok(pd));
    pd_regenerate_hash(pd);
    return pd;

  }

  if (n == 1) { 

    /* 

      Face 0  +--<--+
              | F1  |
           cr0|     0
       +--->--------+
       1  F2  |      
       |      1      
       +--<---+      

    */

    pd = pd_code_new(2); /* Make room for 2, even though we won't use them. */
    
    pd->ncross = 1;
    pd->nedges = 2;
    pd->nfaces = 3;
    pd->ncomps = 1;

    pd->cross[0].edge[0] = 0; pd->cross[0].edge[1] = 0; 
    pd->cross[0].edge[2] = 1; pd->cross[0].edge[3] = 1;
    pd->cross[0].sign = PD_POS_ORIENTATION;

    pd->edge[0].head = 0; pd->edge[0].headpos = 1;
    pd->edge[0].tail = 0; pd->edge[0].tailpos = 0;

    pd->edge[1].head = 0; pd->edge[1].headpos = 2;
    pd->edge[1].tail = 0; pd->edge[1].tailpos = 3;

    pd->face[0].nedges = 2;
    
    pd->face[0].edge = calloc(2,sizeof(pd_idx_t));
    assert(pd->face[0].edge != NULL);
    pd->face[0].edge[0] = 0;
    pd->face[0].edge[1] = 1;
    
    pd->face[0].or = calloc(2,sizeof(pd_or_t));
    assert(pd->face[0].or != NULL);
    pd->face[0].or[0] = PD_NEG_ORIENTATION;
    pd->face[0].or[1] = PD_POS_ORIENTATION;

    pd->face[1].nedges = 1;
    
    pd->face[1].edge = calloc(1,sizeof(pd_idx_t));
    assert(pd->face[1].edge != NULL);
    pd->face[1].edge[0] = 0;
    
    pd->face[1].or = calloc(1,sizeof(pd_or_t));
    assert(pd->face[1].or != NULL);
    pd->face[1].or[0] = PD_POS_ORIENTATION;

    pd->face[2].nedges = 1;
    
    pd->face[2].edge = calloc(1,sizeof(pd_idx_t));
    assert(pd->face[2].edge != NULL);
    pd->face[2].edge[0] = 1;
    
    pd->face[2].or = calloc(1,sizeof(pd_or_t));
    assert(pd->face[2].or != NULL);
    pd->face[2].or[0] = PD_NEG_ORIENTATION;

    pd->comp[0].nedges = 2;
    pd->comp[0].edge = calloc(2,sizeof(pd_idx_t));
    assert(pd->comp[0].edge != NULL);
    pd->comp[0].edge[0] = 0;
    pd->comp[0].edge[1] = 1;
    pd->comp[0].tag = 'A';

    assert(pd_ok(pd));
    pd_regenerate_hash(pd);
    return pd;

  }
   
  /* Now we can generate the general case. */

  pd = pd_code_new(n);
  assert(pd != NULL);
  pd->ncross = n;
  
  pd->cross[0]     = pd_build_cross(0,0,3,2);
  pd->cross[(n-1)] = pd_build_cross(1,1,(2*n)-2,(2*n)-1);

  pd_idx_t i;

  for(i=1;i<n-1;i++) {

    pd->cross[i] = pd_build_cross((2*i),(2*i)+1,(2*i)+3,(2*i)+2);
    
  }

  pd_regenerate(pd);  /* This will unset any crossing information, 
			 because we don't have edges yet. */

  for(i=0;i<pd->ncross;i++) { 

    pd->cross[i].sign = PD_POS_ORIENTATION;

  }

  assert(pd_ok(pd));
  return pd;

}

void pdint_build_tendril(pd_code_t *pd,pd_idx_t cross_start,pd_idx_t edge_start, pd_idx_t ntwists)

/* An "tendril" is a succession of twists (like a plectoneme) which starts with a 
   given edge number (going into the first crossing), and crossing number (first crossing),
   and has ntwists + 1 crossings total-- each twist represents a face with two edges
   along the tendril */

{
  assert(cross_start + ntwists < pd->ncross); /* Assert that there's room for the crossings. */

  /* Classify each (leftmost) crossing by twist number (starting at "twist number 1") */

  pd_idx_t twist_num;
  
  for(twist_num=1;twist_num<=ntwists;twist_num++) { 

    pd_idx_t low_edge, high_edge;

    low_edge  = edge_start + twist_num;
    high_edge = edge_start + 2 * ntwists + 2 - twist_num;

    if (twist_num % 2 == 1) { /* On an ODD-numbered twist */

      pd->cross[cross_start+twist_num-1] = pd_build_cross(low_edge,high_edge,low_edge-1,(high_edge+1)%pd->nedges);

    } else { /* On an EVEN-numbered twist */

      pd->cross[cross_start+twist_num-1] = pd_build_cross(high_edge,low_edge,(high_edge+1)%pd->nedges,low_edge-1);

    }

  }
  
  /* Now we need to finish with the last crossing */

  if (ntwists % 2 == 1) { 

    pd->cross[cross_start+ntwists] = pd_build_cross(edge_start+ntwists+1,edge_start+ntwists+1,
						    (edge_start+ntwists+2)%pd->nedges,edge_start+ntwists);
  } else {

    pd->cross[cross_start+ntwists] = pd_build_cross(edge_start+ntwists+1,edge_start+ntwists+1,
						    edge_start+ntwists,(edge_start+ntwists+2)%pd->nedges);
  }
  
}
  
pd_code_t *pd_build_unknot_wye(pd_idx_t a, pd_idx_t b, pd_idx_t c) 

/* This generates a "wye" with three tendrils, containing "a", "b", and "c" twists, respectively. */
/* The number of twists can be zero. All n-crossing unknot-wyes should hash to the same value. */

{
  pd_code_t *pd = NULL;
  pd_idx_t   n = a + b + c + 3; /* The total wye should contain a+b+c+3 crossings. */

  if (n < 2) { 

    pd_error(SRCLOC,"Can't generate %d-crossing unknot-wye diagram.",pd,n);
    exit(1);

  }

  /* Now we can generate */

  pd = pd_code_new(n);
  assert(pd != NULL);
  pd->ncross = n;
  pd->nedges = 2*n;
  
  pdint_build_tendril(pd,0,0,a);
  pdint_build_tendril(pd,a+1,2*a+2,b);
  pdint_build_tendril(pd,a+b+2,2*a+2*b+4,c);

  /* Now build the edges and faces with regenerate. */
  
  pd_regenerate(pd);

  /* Once we have these, it makes sense to assign crossing signs. */

  pd_idx_t k;

  for(k=0;k<pd->ncross;k++) { 

    pd->cross[k].sign = PD_POS_ORIENTATION;

  }


  assert(pd_ok(pd));

  return pd;

}
  
