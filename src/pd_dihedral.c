/* 

   pd_dihedral.c : The basic data type is pd_dihedral_t, 
   defined in pd_dihedral.h. This is a basic implementation 
   of the dihedral group.

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
#include"pd_dihedral.h"

pd_iterops_t dihedral_ops = {pd_new_dihedral,pd_free_dihedral,
			     pd_print_dihedral,pd_copy_dihedral,
			     pd_increment_dihedral,pd_ndihedrals,
			     pd_dihedral_ok,pd_dihedral_cmp};

void *pd_new_dihedral(void *np)
{
  pd_dihedral_t *d;
  pd_idx_t n = *(pd_idx_t *)(np);

  d = calloc(1,sizeof(pd_dihedral_t));
  assert(d != NULL);

  d->n = n;
  d->map = calloc(n,sizeof(pd_idx_t));
  assert(d->map != NULL);

  pd_idx_t i;
  for(i=0;i<n;i++) { d->map[i] = i; }
  
  d->orient = PD_POS_ORIENTATION;
  
  return d;
}

void pd_free_dihedral(void **dihedralP)
{
  pd_dihedral_t *d = *(pd_dihedral_t **)(dihedralP);

  if (d == NULL) { return; }
  if (d->map != NULL) { free(d->map); d->map = NULL; }
  d->n = 0;
  
  free(d);
  *dihedralP = NULL;
}

char  *pd_print_dihedral(void *dihedralP) 

/* Allocate and print a dihedral */

{
  pd_dihedral_t *d = (pd_dihedral_t *)(dihedralP);
  bool can_id = true;
  char *buf = NULL;
  pd_idx_t i;
  
  /* First, we try to identify the dihedral as a rotation or reflection. */

  if (d->orient == PD_POS_ORIENTATION) { 

    /* This should match to a rotation: make sure everything goes up */

    for(i=0;i<(d->n)-1;i++) {
	
      if (d->map[i] >= d->map[i+1] && !(d->map[i] == (d->n)-1 && d->map[i+1] == 0)) { can_id = false; }
      
    }

    if (can_id) {

      buf = calloc(32,sizeof(char)); assert(buf != NULL);
      snprintf(buf,32,"rot %d",d->map[0]);
      
    }
    
  } else if (d->orient == PD_NEG_ORIENTATION) { 
    
    for(i=0;i<(d->n)-1;i++) { 
      
      if (d->map[i] <= d->map[i+1] && !(d->map[i] == 0 && d->map[i+1] == (d->n) - 1)) { can_id = false; }
      
    }
    
    if (can_id) {
      
      buf = calloc(32,sizeof(char)); assert(buf != NULL);
      snprintf(buf,32,"ref %d",d->map[0]);
      
    }
    
  }

  if (!can_id) { 

    /* Just print the map itself (this code borrowed from pd_print_perm). */

    size_t bufsize = 256;
    buf = calloc(bufsize,sizeof(char)); assert(buf != NULL); 
    
    pd_idx_t i;
    size_t chars_ptd = 0;

    chars_ptd = snprintf(buf,bufsize-chars_ptd,"[BAD] (");

    for(i=0;i<d->n && chars_ptd < bufsize;i++) {

      chars_ptd += snprintf(buf + chars_ptd,bufsize-chars_ptd,"%d",d->map[i]); 
      if (i != d->n-1) { chars_ptd += snprintf(buf + chars_ptd,bufsize-chars_ptd," "); }

    }

    if (chars_ptd > bufsize) { 

      snprintf(buf,bufsize,"[DIHEDRAL > %d CHARS]",(int)(bufsize));

    }
    
  }

  return buf;

}

void  *pd_copy_dihedral(void *dihedralP) 

/* Make a new-memory copy of dihedral. */

{
  pd_dihedral_t *d = (pd_dihedral_t *)(dihedralP);
  pd_dihedral_t *new_d;

  new_d = pd_new_dihedral(&(d->n));
  memcpy(new_d->map,d->map,d->n*sizeof(pd_idx_t));
  new_d->orient = d->orient;

  return (void *)(new_d);
} 

void pd_increment_dihedral(void *dihedralP) 

/* Generates (in place) the next value of d. To do this, we use the algebra of the dihedral group. */
/* We have rotation elements 

   R_0, R_1, ..., R_{n-1}

   and reflection elements

   S_0, S_1, ..., S_{n-1}

   These obey the relations

   R_i R_j = R_{i+j}

   and 

   S_i R_j = S_{i-j}

   This means that post-composing with R_1 will increment a rotation,
   or decrement a reflection, unless we're supposed to switch from 
   R_{n-1} to S_0 or from S_1 to R_0. */

{
  pd_dihedral_t *d = (pd_dihedral_t *)(dihedralP);

  /* The 1-element dihedral group is weird, and needs to be handled separately */

  if (d->n == 1) { 
    
    if (d->orient == PD_POS_ORIENTATION) { d->orient = PD_NEG_ORIENTATION; }
    else if (d->orient == PD_NEG_ORIENTATION) { d->orient = PD_POS_ORIENTATION; }
    else { assert(d->orient == PD_POS_ORIENTATION || d->orient == PD_NEG_ORIENTATION); }

  } else if (d->orient == PD_POS_ORIENTATION && d->map[0] == d->n - 1) { /* R_{n-1} -> S_0 */

    d->orient = PD_NEG_ORIENTATION;
    pd_idx_t i;
    d->map[0] = 0;
    for(i=1;i<d->n;i++) { d->map[i] = d->n - i; }

  } else if (d->orient == PD_NEG_ORIENTATION && d->map[0] == 1) { /* S_1 -> R_0 */

    d->orient = PD_POS_ORIENTATION;
    pd_idx_t i;
    for(i=0;i<d->n;i++) { d->map[i] = i; }

  } else { /* Neither edge case, so post-compose with rotation. */
    
    pd_idx_t swap,i;
    for(swap = d->map[0],i=0;i<(d->n)-1;i++) { d->map[i] = d->map[i+1]; }
    d->map[(d->n)-1] = swap;
    
  }    

}


unsigned int pd_ndihedrals(void *dihedralP) 

/* Count the number of unique values that d can take. */

{
  pd_dihedral_t *d = (pd_dihedral_t *)(dihedralP);
  return 2*d->n;
}


bool pd_dihedral_ok(void *dihedralP) 

/* Verifies that this is a meaningful element of the dihedral group. */

{
  pd_dihedral_t *d = (pd_dihedral_t *)(dihedralP);
  pd_idx_t i,j;

  for(i=0;i<d->n;i++) {

    if (d->map[i] >= d->n) { 

      return pd_error(SRCLOC,"%DIHEDRAL contains illegal index %d at position %d",NULL,
		      d,d->map[i],i);

    }

  }

  for(i=0;i<d->n;i++) {

    for(j=i+1;j<d->n;j++) {

      if (d->map[i] == d->map[j]) {

	return pd_error(SRCLOC,"%DIHEDRAL contains repeated index (%d) at positions %d and %d\n",NULL,
			d,d->map[i],i,j);

      }

    }

  }

  if (!(d->orient == PD_POS_ORIENTATION || d->orient == PD_NEG_ORIENTATION)) {

    return pd_error(SRCLOC,"%DIHEDRAL contains illegal orientation %d.\n",NULL,
		    d,d->orient);

  }

  /* Now we know that d contains a valid permutation of 0..n-1, and a valid orientation. */

  if (d->orient == PD_NEG_ORIENTATION) { /* This is a reflection and should go down as we read. */

    for(i=0;i<(d->n)-1;i++) { 

      if (d->map[i] <= d->map[i+1] && !(d->map[i] == 0 && d->map[i+1] == (d->n) - 1)) {

	return pd_error(SRCLOC,"%DIHEDRAL has negative orientation but isn't a reflection.\n",NULL,d);

      }

    }
 
  } else { /* This is a rotation and should go up as we read. */

    for(i=0;i<(d->n)-1;i++) {

      if (d->map[i] >= d->map[i+1] && !(d->map[i] == (d->n)-1 && d->map[i+1] == 0)) { 

	return pd_error(SRCLOC,"%DIHEDRAL has positive orientation but isn't a rotation.\n",NULL,d);

      }

    }
    
  }

  return true;

}

int          pd_dihedral_cmp(const void *dihedralAp,const void *dihedralBp)
{
  pd_dihedral_t *dihedralA = *(pd_dihedral_t **)(dihedralAp);
  pd_dihedral_t *dihedralB = *(pd_dihedral_t **)(dihedralBp);

  pd_idx_t i,n;

  assert(dihedralA->n == dihedralB->n);
  n = dihedralA->n;

  if (dihedralA->orient != dihedralB->orient) {

    if (dihedralA->orient == PD_POS_ORIENTATION) { return -1; }
    else { return +1; }

  }

  for(i=0;i<n;i++) {

    if (dihedralA->map[i] != dihedralB->map[i]) {

      return dihedralA->map[i] - dihedralB->map[i];

    }

  }

  return 0;

}
 

bool pd_dihedrals_unique(unsigned int ndihedrals,pd_dihedral_t **dihedral_buf)

{
  assert(dihedral_buf != NULL);

  /* Check to see if the dihedrals are unique */

  qsort(dihedral_buf,ndihedrals,sizeof(pd_dihedral_t *),pd_dihedral_cmp);

  unsigned int i;

  for(i=0;i<ndihedrals-1;i++) {

    if (pd_dihedral_cmp(&(dihedral_buf[i]),&(dihedral_buf[i+1])) == 0) {

      return pd_error(SRCLOC,"dihedral_buf contains %DIHEDRAL == %DIHEDRAL at positions %d, %d \n",NULL,
		      dihedral_buf[i],dihedral_buf[i+1],i,i+1);

    }

  }

  return true;

}
