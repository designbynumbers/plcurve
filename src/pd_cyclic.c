/* 

   pd_cyclic.c : The basic data type is pd_cyclic_t, 
   defined in pd_cyclic.h. This is a basic implementation 
   of the cyclic group.

*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
#endif

#ifdef HAVE_SYS_STAT_H
  #include<sys/stat.h>
#endif

#ifdef HAVE_SYS_TYPES_H
  #include<sys/types.h>
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

#ifdef HAVE_STDINT_H
  #include<stdint.h>
#endif

#ifdef HAVE_STDBOOL_H
  #include<stdbool.h>
#endif

#ifdef HAVE_STDLIB_H
  #include<stdlib.h>
#endif

#ifdef HAVE_STRING_H
  #include<string.h>
#endif

#include"plcTopology.h"
#include"pd_multidx.h"
#include"pd_cyclic.h"

pd_iterops_t cyclic_ops = {pd_new_cyclic,pd_free_cyclic,
			   pd_print_cyclic,pd_copy_cyclic,
			   pd_increment_cyclic,pd_ncyclics,
			   pd_cyclic_ok,pd_cyclic_cmp};

void *pd_new_cyclic(void *np)
{
  pd_cyclic_t *c;
  pd_idx_t n = *(pd_idx_t *)(np);

  c = calloc(1,sizeof(pd_cyclic_t));
  assert(c != NULL);

  c->n = n;
  c->map = calloc(n,sizeof(pd_idx_t));
  assert(c->map != NULL);

  pd_idx_t i;
  for(i=0;i<n;i++) { c->map[i] = i; }
  
  return c;
}

void pd_free_cyclic(void **cyclicP)
{
  pd_cyclic_t *c = *(pd_cyclic_t **)(cyclicP);

  if (c == NULL) { return; }
  if (c->map != NULL) { free(c->map); c->map = NULL; }
  c->n = 0;
  
  free(c);
  *cyclicP = NULL;
}

char  *pd_print_cyclic(void *cyclicP) 

/* Allocate and print a cyclic */

{
  pd_cyclic_t *c = (pd_cyclic_t *)(cyclicP);
  bool can_id = true;
  char *buf = NULL;
  pd_idx_t i;
  
  /* This is a rotation: make sure everything goes up */

  for(i=0;i<(c->n)-1;i++) {
	
    if (c->map[i] >= c->map[i+1] && !(c->map[i] == (c->n)-1 && c->map[i+1] == 0)) { can_id = false; }
      
  }

  if (can_id) {

    buf = calloc(56,sizeof(char)); assert(buf != NULL);
    snprintf(buf,56,"rot %d",c->map[0]);
    
  }
    
  if (!can_id) { 

    /* Just print the map itself (this code borrowed from pd_print_perm). */

    size_t bufsize = 256;
    buf = calloc(bufsize,sizeof(char)); assert(buf != NULL); 
    
    pd_idx_t i;
    size_t chars_ptd = 0;

    chars_ptd = snprintf(buf,bufsize-chars_ptd,"[BAD CYCLIC] (");

    for(i=0;i<c->n && chars_ptd < bufsize;i++) {

      chars_ptd += snprintf(buf + chars_ptd,bufsize-chars_ptd,"%d",c->map[i]); 
      if (i != c->n-1) { chars_ptd += snprintf(buf + chars_ptd,bufsize-chars_ptd," "); }

    }

    if (chars_ptd > bufsize) { 

      snprintf(buf,bufsize,"[CYCLIC > %d CHARS]",(int)(bufsize));

    }
    
  }

  return buf;

}

void  *pd_copy_cyclic(void *cyclicP) 

/* Make a new-memory copy of cyclic. */

{
  pd_cyclic_t *c = (pd_cyclic_t *)(cyclicP);
  pd_cyclic_t *new_c;

  new_c = pd_new_cyclic(&(c->n));
  memcpy(new_c->map,c->map,c->n*sizeof(pd_idx_t));
  return (void *)(new_c);
} 

void pd_increment_cyclic(void *cyclicP) 

/* Generates (in place) the next value of d. To do this, we use the algebra of the cyclic group. */
/* We have rotation elements 

   R_0, R_1, ..., R_{n-1}

*/
{
  pd_cyclic_t *c = (pd_cyclic_t *)(cyclicP);

  /* The 1-element cyclic group is weird, and needs to be handled separately */

  if (c->n == 1) { 
    
    return; /* There's nothing TO do. */

  } else { /* Post-compose with rotation. */
    
    pd_idx_t swap,i;
    for(swap = c->map[0],i=0;i<(c->n)-1;i++) { c->map[i] = c->map[i+1]; }
    c->map[(c->n)-1] = swap;
    
  }    

}


unsigned int pd_ncyclics(void *cyclicP) 

/* Count the number of unique values that d can take. */

{
  pd_cyclic_t *c = (pd_cyclic_t *)(cyclicP);
  return c->n;
}


bool pd_cyclic_ok(void *cyclicP) 

/* Verifies that this is a meaningful element of the cyclic group. */

{
  pd_cyclic_t *c = (pd_cyclic_t *)(cyclicP);
  pd_idx_t i,j;

  for(i=0;i<c->n;i++) {

    if (c->map[i] >= c->n) { 

      return pd_error(SRCLOC,"%CYCLIC contains illegal index %d at position %d",NULL,
		      c,c->map[i],i);

    }

  }

  for(i=0;i<c->n;i++) {

    for(j=i+1;j<c->n;j++) {

      if (c->map[i] == c->map[j]) {

	return pd_error(SRCLOC,"%CYCLIC contains repeated index (%d) at positions %d and %d\n",NULL,
			c,c->map[i],i,j);

      }

    }

  }

  /* Now we know that d contains a valid permutation of 0..n-1 */

  for(i=0;i<(c->n)-1;i++) {

    if (c->map[i] >= c->map[i+1] && !(c->map[i] == (c->n)-1 && c->map[i+1] == 0)) { 

      return pd_error(SRCLOC,"%CYCLIC has isn't a rotation.\n",NULL,c);

    }
    
  }
    
  return true;

}

int          pd_cyclic_cmp(const void *cyclicAp,const void *cyclicBp)
{
  pd_cyclic_t *cyclicA = *(pd_cyclic_t **)(cyclicAp);
  pd_cyclic_t *cyclicB = *(pd_cyclic_t **)(cyclicBp);

  pd_idx_t i,n;

  assert(cyclicA->n == cyclicB->n);
  n = cyclicA->n;

  for(i=0;i<n;i++) {

    if (cyclicA->map[i] != cyclicB->map[i]) {

      return cyclicA->map[i] - cyclicB->map[i];

    }

  }

  return 0;

}
 

bool pd_cyclics_unique(unsigned int ncyclics,pd_cyclic_t **cyclic_buf)

{
  assert(cyclic_buf != NULL);

  /* Check to see if the cyclics are unique */

  qsort(cyclic_buf,ncyclics,sizeof(pd_cyclic_t *),pd_cyclic_cmp);

  unsigned int i;

  for(i=0;i<ncyclics-1;i++) {

    if (pd_cyclic_cmp(&(cyclic_buf[i]),&(cyclic_buf[i+1])) == 0) {

      return pd_error(SRCLOC,"cyclic_buf contains %CYCLIC == %CYCLIC at positions %d, %d \n",NULL,
		      cyclic_buf[i],cyclic_buf[i+1],i,i+1);

    }

  }

  return true;

}
