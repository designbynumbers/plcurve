/*

   pd_invariants.h : This is a collection of more-or-less experimental code
   for computing the Arnol'd invariants of pd codes, and some related quantities
   like the number of interlaced crossings. 

*/

#ifndef __PD_INVARIANTS_H 
#define __PD_INVARIANTS_H 1

int *pd_interlaced_crossings(pd_code_t *pd);
/* Returns an array, pd->ncomps long, counting signed sum of crossings
   which occur in the order ABAB along the component according to Polyak's
   sign convention. This is (-1/2) *(2 St + J^+) by Polyak's theorem. */

unsigned int *pd_interlaced_crossings_unsigned(pd_code_t *pd);
/* Returns an array, pd->ncomps long, counting the raw number of crossings
   which occur in the order ABAB along the component. */

#endif
