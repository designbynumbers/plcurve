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

#include<assert.h>

#ifdef HAVE_STDINT_H
  #include<stdint.h>
#endif

#include<stdio.h>

#ifdef HAVE_STDBOOL_H
  #include<stdbool.h>
#endif

#include<stdlib.h>
#include<string.h>

#include"plcTopology.h"
#include"pd_multidx.h"
#include"pd_orientation.h"

pd_iterops_t orientation_ops = {pd_new_orientation,pd_free_orientation,
				pd_print_orientation,pd_copy_orientation,
				pd_increment_orientation,pd_norientations,
				pd_orientation_ok,pd_orientation_cmp};

void *pd_new_orientation(void *nodata)
{
  pd_orientation_t *d;
  d = calloc(1,sizeof(pd_orientation_t));
  assert(d != NULL);
  d->orient = PD_POS_ORIENTATION;
 
  return d;
}

void pd_free_orientation(void **orientationP)
{
  pd_orientation_t *d = *(pd_orientation_t **)(orientationP);

  if (d == NULL) { return; }
  free(d);

  *orientationP = NULL;
}

char  *pd_print_orientation(void *orientationP) 

/* Allocate space for and print a orientation */

{
  pd_orientation_t *d = (pd_orientation_t *)(orientationP);
  size_t bufsize = 16;

  char *buf = calloc(bufsize,sizeof(char));
  assert(buf != NULL);

  snprintf(buf,bufsize,"or %c",pd_print_or(d->orient));
  return buf;

}

void  *pd_copy_orientation(void *orientationP) 

/* Make a new-memory copy of orientation. */

{
  pd_orientation_t *d = (pd_orientation_t *)(orientationP);
  pd_orientation_t *new_d;

  new_d = pd_new_orientation(NULL);
  new_d->orient = d->orient;
  return (void *)(new_d);
} 

void pd_increment_orientation(void *orientationP) 

/* Generates (in place) the next value of d. */

{
  pd_orientation_t *d = (pd_orientation_t *)(orientationP);

  if (d->orient == PD_POS_ORIENTATION) { 

    d->orient = PD_NEG_ORIENTATION;

  } else if (d->orient == PD_NEG_ORIENTATION) { 

    d->orient = PD_POS_ORIENTATION;

  } else {

    pd_error(SRCLOC,"can't increment invalid orientation %d",NULL,d->orient);

  }

}


unsigned int pd_norientations(void *orientationP) 

/* Count the number of unique values that d can take.
*/
{
  return 2;
}


bool pd_orientation_ok(void *orientationP) 

/* Verifies that this is a meaningful element of the orientation group. */

{
  pd_orientation_t *d = (pd_orientation_t *)(orientationP);

  if (!pd_or_ok(d->orient)) { 

    return pd_error(SRCLOC,
		    "%ORIENTATION contains illegal orientation %d",
		    NULL,d,d->orient);

  }
  
  return true;

}

int pd_orientation_cmp(const void *orientationAp,const void *orientationBp)
{
  pd_orientation_t *orientationA = *(pd_orientation_t **)(orientationAp);
  pd_orientation_t *orientationB = *(pd_orientation_t **)(orientationBp);

  if (orientationA->orient == orientationB->orient) {

    return 0;

  } 

  if (orientationA->orient == PD_POS_ORIENTATION) {

    return -1;

  } else {

    return 1;

  }

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
