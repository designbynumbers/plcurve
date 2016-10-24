/* 

   pd_simplify.c: Simplify pd codes by performing Reidemeister moves.

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

#ifdef HAVE_ASSERT_H
   #include<assert.h>
#endif

#include"plcTopology.h"

pd_code_t *pd_simplify(pd_code_t *pd)
/* 
   This is a very simple temporary implementation of pd_simplify which just deloops the pdcode.
*/
{
  bool found_monogon;
  pd_code_t *workingpd = pd_copy(pd);
  
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

	
  
  

