/*

  kascentpermutation.c : This file provides a fast generator for 
                         random k ascent permutations, using the 
			 nonrecursive method in k_ascent_permutation.nb.

			 We try to be efficient, so we're using various 
			 tricks to save bits of time as we go.

*/

#include <config.h>
#include <plCurve.h>

#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#ifdef HAVE_MATH_H
#include <math.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef HAVE_STDARG_H
#include <stdarg.h>
#endif
#ifdef HAVE_CTYPE_H
#include <ctype.h>
#endif
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#ifdef HAVE_ASSERT_H
#include <assert.h>
#endif
#ifdef HAVE_FLOAT_H
#include <float.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#ifdef HAVE_GSL_GSL_RNG_H
#include <gsl/gsl_rng.h>
#endif

#ifdef HAVE_GSL_GSL_RANDIST_H
#include <gsl/gsl_randist.h>
#endif

#include<naptable.h>

double new_ascent_probability(int n,int k)
/*
   Computes the probability 

   (n-k) Eulerian(n-1,k-1)/Eulerian(n,k)

   using the lookup table naptable. Since naptable stores
   data from 1 <= n <= 1000 and 1 <= k <= n-1, there's
   a somewhat complicated offset function required to unpack
   the data on the fly. We work out the proof in 

   k_ascent_permutation_generation.nb.

*/

{
  int ofs = (2 - 3*n + n*n)/2; // Int division ok.
  return naptable[ofs + k-1];
}
  
typedef struct choicestruct {
  int  n;
  int  k;
} perm_choice_t;

void build_choice_vector(int N, int K,
			 int *nchoices,perm_choice_t **cv,
			 gsl_rng *r)
/*
  To generate a permutation of n letters with k ascents, we work essentially
  recursively: such a permutation comes from a permutation of (n-1) letters
  with (k-1) ascents (with probability new_ascent_probability(n,k) = X) or 
  from a permutation of (n-1) letters with k ascents (with probability 1-X). 

  We record the list of choices of n and k from the initial n,k to the point
  where k = 1 or n-1, returning a buffer recording these pairs. 

  Since this is a probabilistic function, we need a gsl random number generator
  to work. 
*/
  
{
  *cv = malloc(N*sizeof(perm_choice_t));
  assert(*cv != NULL);

  int n=N,k=K;
  int i=0;

  do {

    (*cv)[i].n = n;
    (*cv)[i].k = k;
    i++;
      
    if(gsl_rng_uniform(r) < new_ascent_probability(n,k)) { k--; }
    /* This generates a random number in [0,1). */
    n--;
    
  } while (k != n-1 && k != 0);

  (*cv)[i].n = n;
  (*cv)[i].k = k;

  *nchoices = i+1;

}

/*
  We're going to store the permutation in a linked list while we 
  work. In order to prevent cache misses (remember, these arrays 
  are tiny!), we'll fake the linked list by allocating all the 
  storage in advance and keep track of the storage used.
  
*/

typedef struct permnode_struct {

  int s;
  int next;
  
} perm_node_t;


typedef struct workperm_struct {

  perm_node_t *node;
  int start;
  int max;
  int used;

} working_perm_t;

working_perm_t wperm_new(int n)
{
  working_perm_t wp;
  wp.max = n;
  wp.used = 0;
  wp.start = 0;
  wp.node = calloc(n,sizeof(perm_node_t));
  assert(wp.node != NULL);
  
  return wp;
}

void make_descending_wperm(working_perm_t *wp,int n)
{
  int i;
  assert(n < wp->max);
  
  for(i=0;i<n;i++) {

    wp->node[i].s = n-i;
    wp->node[i].next = i+1;

  }

  wp->node[i-1].next = -1;
  wp->start = 0;
  wp->used = n;
}

void make_ascending_wperm(working_perm_t *wp,int n)
{
  int i;
  assert(n < wp->max);
  
  for(i=0;i<n;i++) {

    wp->node[i].s = i+1;
    wp->node[i].next = i+1;

  }

  wp->node[i-1].next = -1;
  wp->start = 0;
  wp->used = n;
}

void convert_wperm_to_buffer(working_perm_t *wp,int *buf)
/* 
   Traverses the linked list, writing entries to the 
   buffer. We assume that "buf" is allocated, and is 
   large enough to store the entries of wp.
*/

{
  int i = 0;
  int loc = wp->start;

  for(;loc != -1;i++,loc = wp->node[loc].next) {

    buf[i] = wp->node[loc].s;

  }
}

void wperm_free(working_perm_t *wp)
/* 
   Frees memory associated with working perm 
*/
{
  if (wp->node != NULL) {

    free(wp->node);
    wp->node = NULL;

  }

  wp->start = 0;
  wp->max = 0;
  wp->used = 0;
}

void wperm_insert_at_ascent(working_perm_t *wp,
			    int toinsert,
			    int ascent)

/* Walks the linked list wp until we find ascent number "ascent", then
   inserts "toinsert" in the middle of this ascent. */

{
  int ascentcount = -1; /* We haven't seen any ascents yet. */
  int loc;

  for(loc = wp->start;
      wp->node[loc].next != -1;
      loc = wp->node[loc].next) {

    if (wp->node[loc].s < wp->node[wp->node[loc].next].s) {

      ascentcount++;

      if (ascentcount == ascent) {
 
	/* We're at the ascent in question! */

	wp->node[wp->used].s = toinsert;
	wp->node[wp->used].next = wp->node[loc].next;
	wp->node[loc].next = wp->used;
	wp->used++;

	/* Check to make sure we don't use too much memory */

	assert(wp->used <= wp->max);

	return;

      }

    }

  }

  /* We've traversed the entire linked list without finding 
     the correct ascent. At this point, we should just bail. */

  pd_error(SRCLOC,
	   "attempted to insert at ascent %d, but the perm\n"
	   "does not have this many ascents.",NULL,
	   ascent);
  exit(1);

} 

void wperm_insert_at_descent(working_perm_t *wp,
			     int toinsert,
			     int descent)
{
/* Walks the linked list wp until we find descent number "descent", then
   inserts "toinsert" in the middle of this descent. */

  int descentcount = -1; /* We haven't seen any descents yet. */
  int loc;

  for(loc = wp->start;
      wp->node[loc].next != -1;
      loc = wp->node[loc].next) {

    if (wp->node[loc].s > wp->node[wp->node[loc].next].s) {

      descentcount++;

      if (descentcount == descent) {
 
	/* We're at the descent in question! */

	wp->node[wp->used].s = toinsert;
	wp->node[wp->used].next = wp->node[loc].next;
	wp->node[loc].next = wp->used;
	wp->used++;

	/* Check to make sure we don't use too much memory */

	assert(wp->used <= wp->max);

	return;

      }

    }

  }

  /* We've traversed the entire linked list without finding 
     the correct descent. At this point, we should just bail. */

  pd_error(SRCLOC,
	   "attempted to insert at descent %d, but the perm\n"
	   "does not have this many descents.",NULL,
	   descent);
  exit(1);

} 

void wperm_insert_at_start(working_perm_t *wp,
			   int toinsert)
{
  wp->node[wp->used].s = toinsert;
  wp->node[wp->used].next = wp->start;
  wp->start = wp->used;
  wp->used++;
  assert(wp->used <= wp->max);
}

void wperm_insert_at_end(working_perm_t *wp,
			 int toinsert)
{
  /* This is really different, because we 
     don't know where the last node is. We
     can search for it by searching for the
     node with wp->node[i].next = -1. We 
     know this will be faster than traversing
     the list in order. */

  int i;
  for(i=0;i<wp->used;i++) {

    if (wp->node[i].next == -1) { /* This is the last node! */

      wp->node[i].next = wp->used;
      wp->node[wp->used].s = toinsert;
      wp->node[wp->used].next = -1;
      wp->used++;

      assert(wp->used <= wp->max);
      return;
    }

  }

  pd_error(SRCLOC,"searched linked list for last element, but "
	   "didn't find one. List must be corrupt.\n",NULL);
  exit(1);
    
}

int  ascent_count(int *perm,int n)
/* 
  Count the number of ascents in the permutation perm. 
*/
{
  int i, ascent_count = 0;

  for(i=0;i<n-1;i++) {

    if (perm[i+1] > perm[i]) { ascent_count++; }

  }

  return ascent_count;
}



void random_k_ascent_permutation(int n,int k,
				 int *perm,
				 gsl_rng *rng)
  
/* 
   This version generates random k ascent permutations with 
   1 <= n <= 1000, 1 <= k <= n-1. We use a "recursive" algorithm.
   We assume that "perm" is already allocated, and contains space
   for n ints. When we're done, it will be overwritten with a 
   permutation of 1..n containing k ascents.
*/
{
  working_perm_t wp = wperm_new(n);

  /*  
    We can use the subprocedure build_choice_vector to generate
    a sequence in the form 

    {n,k}, {n-1,X}, ...., {nend,kend}

    where kend = 0 or kend = nend - 1. In either case, we know the 
    initial permutation: if there are no ascents, it must be the 
    permutation nend, nend-1, ..., 1 and if there are nend-1 ascents,
    it must be the permutation 1, 2, ..., nend.

  */

  perm_choice_t *cv;
  int nchoices;

  build_choice_vector(n,k,&nchoices,&cv,rng);

  if (cv[nchoices-1].k == 0) {

    make_descending_wperm(&wp,cv[nchoices-1].n);

  } else {

    make_ascending_wperm(&wp,cv[nchoices-1].n);

  }

  /* Now we're going to traverse the choice vector, doing inserts
     into the working permutation as we go. */

  int i;

  for(i=nchoices-1;i>0;i--) {

    int toinsert = cv[i-1].n;

    if (cv[i-1].k != cv[i].k) {
      
      /* We should insert "toinsert" at a point which 
         creates a new ascent. Since we're inserting something
	 larger than all current elements of the permutation,
	 this means we're inserting at the end, or in the middle
	 of a current descent. 

	 Now there are cv[i].k ascents in the current permutation 
	 of cv[i].n letters, so there are (cv[i].n-1-cv[i].k) descents.

	 This means that there are cv[i].n-cv[i].k possible locations
	 in which to do this insert. We need to pick one at random
	 and then insert at that position. */

      int descent = gsl_rng_uniform_int(rng,cv[i].n - cv[i].k);

      if (descent == cv[i].n-cv[i].k-1) {

	wperm_insert_at_end(&wp,toinsert);

      } else {

	wperm_insert_at_descent(&wp,toinsert,descent);

      }

    } else {

      /* In this case, we're not creating a new ascent. Which 
	 means that we're inserting in the middle of an existing
	 ascent or at the start. There are 

	 cv[i].k

	 existing ascents, so there are cv[i].k + 1 possibilities.
      */

      int ascent = gsl_rng_uniform_int(rng,cv[i].k+1);

      if (ascent == cv[i].k) {

	wperm_insert_at_start(&wp,toinsert);

      } else {

	wperm_insert_at_ascent(&wp,toinsert,ascent);

      }

    }

  }

  /* We should have done all the inserts, and the resulting working
     perm should be a complete k ascent permutation. We read it out
     into the buffer, which should already have been allocated, and
     then free the wperm memory. */

  free(cv);
  convert_wperm_to_buffer(&wp,perm);
  wperm_free(&wp);
}
  


      

  
  



