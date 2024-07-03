/* 

   pd_simplify.c: Simplify pd codes by performing Reidemeister moves.

*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
#endif

#include<assert.h>

#ifdef HAVE_STDINT_H
  #include<stdint.h>
#endif

#include<stdio.h>
#include<string.h>
#include<stdlib.h>

#ifdef HAVE_STDBOOL_H
  #include<stdbool.h>
#endif

#include<math.h>
#include<assert.h>

#include"plcTopology.h"

pd_code_t *pd_simplify(pd_code_t *pd)
/* 
   This is a very simple temporary implementation of pd_simplify which just deloops the pdcode.
*/
{
  bool found_monogon;
  pd_code_t *workingpd = pd_copy(pd);

  if (pd->ncross == 0) { /* This is already a 0-crossing unknot */

    return workingpd;

  }
  
  do {

    pd_idx_t f;
    found_monogon = false;
    
    for(f=0;f<workingpd->nfaces && !found_monogon;f++) {

      if (workingpd->face[f].nedges == 1) {

	pd_code_t *newpd = pd_R1_loopdeletion(workingpd,workingpd->edge[workingpd->face[f].edge[0]].head);
	pd_code_free(&workingpd);
	workingpd = newpd;

	found_monogon = true;

      }
      
    }
    
  } while (found_monogon && workingpd->ncross > 0);

  return workingpd;

}

	
  
  

