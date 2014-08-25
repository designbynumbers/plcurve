/* 

   pd_sortedbuf.c : Contains some utility functions for maintaining buffers
   of pd_idx_t in sorted order. 

*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
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

#ifdef HAVE_STDBOOL_H
  #include<stdbool.h>
#endif

#ifdef HAVE_STRING_H
  #include<string.h>
#endif

#include<plcTopology.h>

/* int pd_idx_cmp(const void *A, const void *B) */
/* { */
/*   pd_idx_t *pdiA = (pd_idx_t *)(A); */
/*   pd_idx_t *pdiB = (pd_idx_t *)(B); */

/*   return (int)(*pdiA) - (int)(*pdiB); */
/* } */

  
bool      pd_sortedbuf_contains(pd_idx_t *buf,pd_idx_t nelems,pd_idx_t i)
/* Returns true if the buffer contains the number "i". */
{
  // This is basically a wrapper for the C library bsearch function.
  pd_idx_t *pos = bsearch(&i,buf,(size_t)(nelems),sizeof(pd_idx_t),pd_idx_cmp);
  if (pos == NULL) { 
    return false;
  } else {
    return true;
  }
}

pd_idx_t *pd_sortedbuf_insert(pd_idx_t *buf,pd_idx_t *nelems,pd_idx_t insert)
/* Inserts "insert" into the buffer, maintaining sort order, 
   and making the buffer bigger if needed. */
{
  pd_idx_t *pos = bsearch(&insert,buf,(size_t)(*nelems),sizeof(pd_idx_t),pd_idx_cmp);

  if (pos != NULL) { /* This is already IN the array. */

    return buf;

  }

  /* We have to move the data anyway, so realloc won't save us much time. */

  pd_idx_t *newbuf = malloc((size_t)((*nelems+1)*sizeof(pd_idx_t)));
  pd_idx_t i,j;

  if (insert < buf[0]) { 

    newbuf[0] = insert; 
    memcpy(&(newbuf[1]),buf,(size_t)(*nelems * sizeof(pd_idx_t)));
    return newbuf;

  } 

  newbuf[0] = buf[0];
  
  for(i=1,j=1;i<*nelems;i++,j++) { 

    if (newbuf[i-1] < insert && insert < newbuf[i]) { 
      
      newbuf[j] = insert;
      j++; 

    }
    
    newbuf[j] = buf[i];

  }

  free(buf);
  (*nelems) = (*nelems) + 1;
  return newbuf;

}
  

bool      pd_buf_contains(pd_idx_t *buf,pd_idx_t size,pd_idx_t i, pd_idx_t *where)
/* Returns true if the buffer contains the number "i" */
{

  pd_idx_t j; 

  for(j=0;j<size;j++) { 

    if (buf[j] == i) { *where = j; return true; }

  } 

  *where = PD_UNSET_IDX;
  return false;

}

    
