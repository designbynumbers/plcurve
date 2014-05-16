/* 

   pd_orientation.c : The basic data type is pd_orientation_t, 
   defined in pd_orientation.h. This is a basic implementation 
   of the group (Z/2Z)^n, which collects possible orientations 
   of n objects.

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

#include<plcTopology.h>
#include<pd_multidx.h>
#include<pd_orientation.h>

pd_iterops_t orientation_ops = {pd_new_orientation,pd_free_orientation,
				pd_print_orientation,pd_copy_orientation,
				pd_increment_orientation,pd_norientations,
				pd_orientation_ok,pd_orientation_cmp};

void *pd_new_orientation(void *np)
{
  pd_orientation_t *d;
  pd_idx_t n = *(pd_idx_t *)(np);

  d = calloc(1,sizeof(pd_orientation_t));
  assert(d != NULL);

  d->n = n;
  d->or = calloc(n,sizeof(pd_or_t));
  assert(d->or != NULL);

  pd_idx_t i;
  for(i=0;i<n;i++) { d->or[i] = PD_POS_ORIENTATION; }
  
  return d;
}

void pd_free_orientation(void **orientationP)
{
  pd_orientation_t *d = *(pd_orientation_t **)(orientationP);

  if (d == NULL) { return; }
  if (d->or != NULL) { free(d->or); d->or = NULL; }
  d->n = 0;
  
  free(d);
  *orientationP = NULL;
}

char  *pd_print_orientation(void *orientationP) 

/* Allocate space for and print a orientation */

{
  pd_orientation_t *d = (pd_orientation_t *)(orientationP);
  int bufsize = 2*d->n + 5, total_used=0, this_used=0;

  char *buf = calloc(bufsize,sizeof(char));
  assert(buf != NULL);
  char *bufptr = buf;
  pd_idx_t i;

  this_used = snprintf(bufptr,bufsize-total_used,"or ");
  total_used += this_used;
  bufptr += this_used;

  for(i=0;i<d->n;i++) { 
    
    this_used = snprintf(bufptr,bufsize-total_used,"%c",pd_print_or(d->or[i]));
    total_used += this_used;
    bufptr += this_used;

  }
  
  return buf;

}

void  *pd_copy_orientation(void *orientationP) 

/* Make a new-memory copy of orientation. */

{
  pd_orientation_t *d = (pd_orientation_t *)(orientationP);
  pd_orientation_t *new_d;

  new_d = pd_new_orientation(&(d->n));
  memcpy(new_d->or,d->or,d->n*sizeof(pd_or_t));
  return (void *)(new_d);
} 

void pd_increment_orientation(void *orientationP) 

/* Generates (in place) the next value of d. To do this, we 
   implement (basically) a binary add: given a string of 
   binary digits, read from left to right converting 1s to 0s
   until you reach the first 0. Convert it to a 1.*/

{
  pd_orientation_t *d = (pd_orientation_t *)(orientationP);
  pd_idx_t i;

  for(i=0;i<d->n && d->or[i] == PD_POS_ORIENTATION;i++) { 
    
    d->or[i] = PD_NEG_ORIENTATION;

  } 

  if (i < d->n && d->or[i] == PD_NEG_ORIENTATION) { 

    d->or[i] = PD_POS_ORIENTATION;

  }
    
}


unsigned int pd_norientations(void *orientationP) 

/* Count the number of unique values that d can take.

   And yes, I realize that you might be able to do this faster with
   bit shifting. But that's often buggy at the compiler level, and
   anyway this will never be called enough for you to notice 
   any difference.

*/
{
  pd_orientation_t *d = (pd_orientation_t *)(orientationP);
  unsigned int result = 1;
  pd_idx_t i;

  for(i=0;i<d->n;i++) { result *= 2; }
  return result ;
}


bool pd_orientation_ok(void *orientationP) 

/* Verifies that this is a meaningful element of the orientation group. */

{
  pd_orientation_t *d = (pd_orientation_t *)(orientationP);
  pd_idx_t i;

  for(i=0;i<d->n;i++) {

    if (!pd_or_ok(d->or[i])) { 

      return pd_error(SRCLOC,
		      "%ORIENTATION contains illegal orientation %d at position %d",
		      NULL,d,d->or[i],i);

    }

  }

  return true;

}

int pd_orientation_cmp(const void *orientationAp,const void *orientationBp)
{
  pd_orientation_t *orientationA = *(pd_orientation_t **)(orientationAp);
  pd_orientation_t *orientationB = *(pd_orientation_t **)(orientationBp);

  pd_idx_t i,n;

  assert(orientationA->n == orientationB->n);
  n = orientationA->n;

  for(i=0;i<n;i++) {

    if (orientationA->or[i] != orientationB->or[i]) {

      return (orientationA->or[i] == PD_POS_ORIENTATION) ? 1 : -1;

    }

  }

  return 0;

}
 

bool pd_orientations_unique(unsigned int norientations,pd_orientation_t **orientation_buf)

{
  assert(orientation_buf != NULL);

  /* Check to see if the orientations are unique */

  qsort(orientation_buf,norientations,sizeof(pd_orientation_t *),pd_orientation_cmp);

  unsigned int i;

  for(i=0;i<norientations-1;i++) {

    if (pd_orientation_cmp(&(orientation_buf[i]),&(orientation_buf[i+1])) == 0) {

      return pd_error(SRCLOC,"orientation_buf contains %ORIENTATION == %ORIENTATION at positions %d, %d \n",NULL,
		      orientation_buf[i],orientation_buf[i+1],i,i+1);

    }

  }

  return true;

}
