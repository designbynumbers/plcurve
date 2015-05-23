/*

   pd_invariants.h : This is a collection of more-or-less experimental code
   for computing the Arnol'd invariants of pd codes, and some related quantities
   like the number of interlaced crossings. 

*/

#ifndef __PD_INVARIANTS_H 
#define __PD_INVARIANTS_H 1

unsigned int *pd_interlaced_crossings(pd_code_t *pd);
/* Returns an array, pd->ncomps long, counting number of crossings
   which occur in the order ABAB along the component. */

#endif
