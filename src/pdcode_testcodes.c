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

#include<pdcode.h>

#include<pd_multidx.h>
#include<pd_dihedral.h>
#include<pd_perm.h>
  
#include<pd_isomorphisms.h>

pd_code_t *pd_build_twist_knot(pd_idx_t n) 

/* Generate a twist knot diagram (scanned documentation in data/pd_twist_knot_doc.pdf) */

{
  pd_code_t *pd;

  /* Make sure that the number of verts is ok */

  if (n+2 > PD_MAXVERTS) {

    pd_error(SRCLOC,"Can't generate pd_twist_knot diagram with %d (> PD_MAXVERTS == %d) crossings.\n",pd,
	     n+2,PD_MAXVERTS);
    exit(1);

  }

  if (n<1) { 

    pd_error(SRCLOC,"Can't generate pd_twist_knot with %d (< 1) twists.\n",pd,n);
    exit(1);

  }
  
  /* Now generate crossings. */

  pd = calloc(1,sizeof(pd_code_t));
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
    
  pd_regenerate(pd);
  assert(pd_ok(pd));

  return pd;
}

pd_code_t *pd_build_torus_knot(pd_idx_t p, pd_idx_t q) 

/* Generate a (2,q) torus knot diagram (scanned documentation in data/pd_torus_knot_doc.pdf) */

{
  pd_code_t *pd;

  /* Make sure that the number of verts is ok */

  if (p != 2) {

    pd_error(SRCLOC,"pd_torus_knot can only generate (2,q) torus knots (called with (%d,%d)).\n",pd,
	     p,q);
    exit(1);

  }
    
  if (q*(p-1) > PD_MAXVERTS) {

    pd_error(SRCLOC,"Can't generate pd_torus_knot diagram with %d (> PD_MAXVERTS == %d) crossings.\n",pd,
	     q*(p-1),PD_MAXVERTS);
    exit(1);

  }

  /* Now generate crossings. */

  pd = calloc(1,sizeof(pd_code_t));
  assert(pd != NULL); 

  pd->ncross = q;
  pd->cross[0] = pd_build_cross(0,(2*q)-2,(2*q)-1,1);

  pd_idx_t i;

  for(i=1;i<q;i++) {

    pd->cross[i] = pd_build_cross((2*i)-2,(2*i)-1,(2*i)+1,(2*i));

  }

  pd_regenerate(pd);
  assert(pd_ok(pd));

  return pd;

}

pd_code_t *pd_build_simple_chain(pd_idx_t n) 
/* Build a simple chain of n links (n > 2) */

{
  pd_code_t *pd;

  /* First, check the size of the diagram. */

  if (n < 3) { 

    pd_error(SRCLOC,"Can't build a simple chain with %d (< 3) links",pd,n);
    exit(1);

  } 

  if (2*(n-1) > PD_MAXVERTS) {

    pd_error(SRCLOC,"Simple chain with %d links has %d crossings (> PD_MAXVERTS == %d)",pd,
	     n,2*(n-1));
    exit(1);

  }

  pd = calloc(1,sizeof(pd_code_t));
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

  pd_regenerate(pd);
  assert(pd_ok(pd));

  return pd;

}

pd_code_t *pd_build_unknot(pd_idx_t n) 

/* An n-crossing unknot diagram. */

{
  pd_code_t *pd;

  if (n > PD_MAXVERTS) { 
    
    pd_error(SRCLOC,"Can't generate unknot diagram with %d ( > PD_MAXVERTS = %d ) crossings.",pd,n,PD_MAXVERTS);
    exit(1);
    
  }
  
  if (n < 2) { 

    pd_error(SRCLOC,"Can't generate %d-crossing unknot diagram.",pd,n);
    exit(1);

  }

  /* Now we can generate */

  pd = calloc(1,sizeof(pd_code_t));
  assert(pd != NULL);

  pd->ncross = n;
  
  pd->cross[0]     = pd_build_cross(0,0,3,2);
  pd->cross[(n-1)] = pd_build_cross(1,1,(2*n)-2,(2*n)-1);

  pd_idx_t i;

  for(i=1;i<n-1;i++) {

    pd->cross[i] = pd_build_cross((2*i),(2*i)+1,(2*i)+3,(2*i)+2);

  }

  pd_regenerate(pd);
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
  pd_code_t *pd;
  pd_idx_t   n = a + b + c + 3; /* The total wye should contain a+b+c+3 crossings. */

  if (n > PD_MAXVERTS) { 
    
    pd_error(SRCLOC,"Can't generate unknot-wye diagram with %d ( > PD_MAXVERTS = %d ) crossings.",pd,n,PD_MAXVERTS);
    exit(1);
    
  }
  
  if (n < 2) { 

    pd_error(SRCLOC,"Can't generate %d-crossing unknot-wye diagram.",pd,n);
    exit(1);

  }

  /* Now we can generate */

  pd = calloc(1,sizeof(pd_code_t));
  assert(pd != NULL);
  pd->ncross = n;
  pd->nedges = 2*n;
  
  pdint_build_tendril(pd,0,0,a);
  pdint_build_tendril(pd,a+1,2*a+2,b);
  pdint_build_tendril(pd,a+b+2,2*a+2*b+4,c);

  pd_regenerate(pd);

  assert(pd_ok(pd));

  return pd;

}
  
