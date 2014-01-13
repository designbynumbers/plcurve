/* 

   pd_multidx.c : The basic data type is pd_multidx_t, 
   defined in pd_multidx.h.

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

#ifdef HAVE_STDLIB_H
  #include<stdlib.h>
#endif

#ifdef HAVE_STDINT_H
  #include<stdint.h>
#endif

#ifdef HAVE_STDBOOL_H
  #include<stdbool.h>
#endif

#ifdef HAVE_STRING_H
  #include<string.h>
#endif

#include<plcTopology.h>
#include<pd_multidx.h>

bool pd_iterops_eq(pd_iterops_t opsA,pd_iterops_t opsB) {

  if (opsA.new == opsB.new &&
      opsA.free == opsB.free &&
      opsA.print == opsB.print &&
      opsA.copy == opsB.copy &&
      opsA.incr == opsB.incr &&
      opsA.nvals == opsB.nvals &&
      opsA.ok == opsB.ok &&
      opsA.cmp == opsB.cmp ) { return true; }
  else {
    return false;
  }
}

pd_multidx_t *pd_new_multidx(pd_idx_t nobj,void **initialdata,pd_iterops_t ops)

{
  pd_multidx_t *idx;
  pd_idx_t i;

  idx = calloc(1,sizeof(pd_multidx_t)); assert(idx != NULL); 

  idx->nobj = nobj;

  idx->obj  = calloc(nobj,sizeof(void *)); assert(idx->obj != NULL);
  for(i=0;i<nobj;i++) { idx->obj[i] = ops.new(initialdata[i]); }   
  
  idx->i = calloc(nobj,sizeof(unsigned int)); assert(idx->i != NULL);
  for(i=0;i<nobj;i++) { idx->i[i] = 0; }
  
  idx->nvals = calloc(nobj,sizeof(unsigned int)); assert(idx->nvals != NULL);
  for(i=0;i<nobj;i++) {  idx->nvals[i] = ops.nvals(idx->obj[i]); }
  
  idx->ops = ops;
  return idx;
}

void pd_free_multidx(pd_multidx_t **idxp)
{
  pd_multidx_t *idx;
  if (*idxp == NULL) { return; }
  
  idx = *idxp;

  /* First, free the array of "objects" */

  if (idx->obj != NULL) {

    pd_idx_t i;
    for(i=0;i<idx->nobj;i++) { idx->ops.free(&(idx->obj[i])); }
    free(idx->obj);
    
  }

  /* Now free the state indices and maxima (if present) */

  if (idx->i != NULL) { free(idx->i); }
  if (idx->nvals != NULL) { free(idx->nvals); }

  /* Now free "idx" itself */

  free(idx);
  *idxp = NULL;
    
}

void pd_increment_multidx(pd_multidx_t *idx) 

  /* Generates (in place) the next value of idx. */
{

  assert(idx != NULL);

  /* 1. Count # of "ready to rollover" indices in a contiguous block 
        starting at the right of the idx->i vector. */

  pd_idx_t rr;
  
  for(rr=0;rr<idx->nobj;rr++) {

    if (idx->i[(idx->nobj-1) - rr] != idx->nvals[(idx->nobj-1) - rr] - 1) { break; }

  }

  if (rr<idx->nobj) { rr++; } 
  /* The last (non-ready to rollover) index gets incremented, unless EVERY index */
  /* is ready to rollover, in which case incrementing rr would leave us trying to */
  /* increment the "-1"th index. */

  /* 2. Increment these indices and objects (if there is an object incrementer) */

  pd_idx_t i;
  
  for(i=(idx->nobj)-rr;i<idx->nobj;i++) {

    idx->i[i] = (idx->i[i] + 1) % idx->nvals[i];
    idx->ops.incr(idx->obj[i]); /* Again, idx->obj[i] is the _pointer_ to the ith object */

  }
  
}

unsigned int pd_multidx_nvals(pd_multidx_t *idx) 

/* Count the number of unique values that idx can take. */

{
  unsigned int nvals=1;
  pd_idx_t j;

  for(j=0;j<idx->nobj;j++) {
    
    nvals *= idx->nvals[j];

  }

  return nvals;
}


bool  pd_multidx_ok(pd_multidx_t *idx)
{
  /* 0. Make sure there's an object to check. */

  if (idx == NULL) { return true; }

  if (idx->nobj == 0) { /* A trivial case; all buffers should be null. */

    if (idx->obj != NULL || idx->i != NULL || idx->nvals != NULL) { return false; }
    return true;

  }

  if (idx->i == NULL || idx->nvals == NULL) { return false; }

  /* 1. If there are objects, check them using obj_ok. */

  pd_idx_t i;

  for(i=0;i<idx->nobj;i++) {

      if (!idx->ops.ok(idx->obj[i])) { return false; }

  }

  /* 2. Now check that the state indices are within range. */

  for(i=0;i<idx->nobj;i++) {

    if( idx->i[i] >= idx->nvals[i]) { return false; }

  }

  return true;

}

int pd_multidx_cmp(const void *idxAp, const void *idxBp)

/* Assuming that there are objects, sorts in idx->obj_cmp
   order on the object list. Otherwise, lex order on
   indices. */

{
  pd_multidx_t *idxA = *(pd_multidx_t **)(idxAp); 
  pd_multidx_t *idxB = *(pd_multidx_t **)(idxBp); 

  assert(idxA->nobj == idxB->nobj);
  pd_idx_t nobj = idxA->nobj;

  assert(pd_iterops_eq(idxA->ops,idxB->ops));
  pd_iterops_t ops = idxA->ops;

  pd_idx_t i;
  for(i=0;i<nobj;i++) {
      
    int cmp;
    cmp = ops.cmp(&(idxA->obj[i]),&(idxB->obj[i])); /* Compare pointers to pointers to objects */
    if (cmp != 0) { return cmp; }
      
  }
    
  return 0;

}
  
bool pd_multidxs_unique(unsigned int nmultidxs,pd_multidx_t **multidx_buf)

{
  assert(multidx_buf != NULL);

  /* First, we copy the old buffer of pointers. */

  pd_multidx_t **mb_copy;
  mb_copy = calloc(nmultidxs,sizeof(pd_multidx_t *)); assert(mb_copy != NULL);
  memcpy(mb_copy,multidx_buf,nmultidxs*sizeof(pd_multidx_t *));

  /* Now we sort the copy */

  qsort(mb_copy,(size_t)(nmultidxs),sizeof(pd_multidx_t *),pd_multidx_cmp);

  unsigned int i;

  for(i=0;i<nmultidxs-1;i++) {

    if (pd_multidx_cmp(&(mb_copy[i]),&(mb_copy[i+1])) == 0) {

      /* Now search for these pointers in the original multidx_buffer. */

      unsigned int j,k;
      bool ifound = false, iplusfound = false;

      for(j=0;j<nmultidxs;j++) { 

	if (mb_copy[i] == multidx_buf[j]) { ifound = true; break; }

      }

      for (k=0;k<nmultidxs;k++) { 

	if (mb_copy[i+1] == multidx_buf[k]) { iplusfound = true; break; }

      }

      assert(ifound && iplusfound);

      return pd_error(SRCLOC,
		      "multidx_buf contains\n"
		      "%MULTIDX == \n"
		      "%MULTIDX at positions %d, %d \n"
		      "(sorted) mb_copy buf positions are %d, %d\n"
		      ,NULL,
		      mb_copy[i],mb_copy[i+1],j,k,i,i+1);

    }

  }

  /* Housekeeping */

  free(mb_copy);

  return true;
  
}

pd_multidx_t *pd_copy_multidx(pd_multidx_t *idx)

/* Make a new-memory copy of multidx. We can't use the usual constructor */
/* because we don't know that the initial data array means anything anymore. */
/* Instead, we do a custom construction, duplicating each iterator using */
/* ops.copy. */

{
  pd_multidx_t *new_idx;
  pd_idx_t i;

  new_idx = calloc(1,sizeof(pd_multidx_t)); assert(new_idx != NULL);
  new_idx->nobj = idx->nobj;

  new_idx->obj = calloc(new_idx->nobj,sizeof(void *)); assert(new_idx->obj != NULL);
  new_idx->i = calloc(new_idx->nobj,sizeof(unsigned int)); assert(new_idx->i != NULL);
  new_idx->nvals = calloc(new_idx->nobj,sizeof(unsigned int)); assert(new_idx->nvals != NULL);

  new_idx->ops = idx->ops;

  /* Now fill in buffers. */

  for(i=0;i<new_idx->nobj;i++) {

    new_idx->obj[i] = new_idx->ops.copy(idx->obj[i]);
    new_idx->i[i] = idx->i[i];
    new_idx->nvals[i] = idx->nvals[i];

  }

  return new_idx;

}
    
    
    
