/* 

   pd_container.c : A general container type. 

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

#ifdef HAVE_STDLIB_H
  #include<stdlib.h>
#endif

#ifdef HAVE_STDBOOL_H
  #include<stdbool.h>
#endif

#include<ordie.h>
#include<pd_container.h>

pd_container_t *pd_new_container(pd_contidx_t initsize)
/* Allocates space for a new pd_container holding initsize pointers. */

{
  pd_container_t *cont;
  
  cont = calloc(1,sizeof(pd_container_t));
  assert(cont != NULL);

  if (initsize < 100) { 

    cont->initsize = 100;

  } else {

    cont->initsize = initsize;

  }

  cont->bufsize  = cont->initsize;
  cont->nelts    = 0;
  cont->elt      = (pd_contdata_t *) calloc_or_die(cont->bufsize,sizeof(pd_contdata_t));

  return cont;

}

void  pd_free_container(pd_container_t **cont,void eltfree(void **ptr))
/* Frees the memory in a pd_container, calling eltfree on each stored pointer if eltfree != NULL,
   or free on each (non-NULL) stored pointer otherwise. */

{
  pd_contidx_t elt;

  if (*cont == NULL) {  /* Ok to call free as many times as you like on a container. */

    return; 

  }

  if ((*cont)->elt != NULL) { /* There's a buffer to free here. */

    for(elt=0;elt<(*cont)->nelts;elt++) {
      
      if ((*cont)->elt[elt] != NULL) { 
	
	if (eltfree != NULL) {
	  
	  eltfree(&((*cont)->elt[elt]));
	  
	} else {
	  
	  free((*cont)->elt[elt]);
	  (*cont)->elt[elt] = NULL;
	  
	}

      }
      
    }

  }

  free((*cont)->elt); 
  free(*cont);
  *cont = NULL;

  return;

}

void  pd_free_fake(void **ptr)
/* A "fake" free function that shouldn't ever be called. */
{
  printf("pd_free_fake called with pointer %p\n",*ptr);
  exit(1);
}

pd_contidx_t pd_container_nelts(pd_container_t *cont)
/* Returns number of elts in container (used to insulate
   the details of the implementation from the API). */
{
  return cont->nelts;
}

void  pd_addto_container(pd_container_t *cont,void *new_elt)
/* Adds an element to a container. */

{
  pd_contidx_t elt;

  assert(cont != NULL);  /* You should never try to add to a freed container */

  if (cont->nelts == cont->bufsize-1) { /* The container is almost full! */

    pd_contidx_t newsize = (cont->bufsize + cont->initsize)*sizeof(pd_contdata_t);

    cont->elt = (pd_contdata_t *)realloc(cont->elt,newsize);
    assert(cont->elt != NULL);
      
    cont->bufsize += cont->initsize;

    for(elt=cont->nelts;elt<cont->bufsize;elt++) {

      cont->elt[elt] = NULL;

    }

  }

  cont->elt[cont->nelts] = new_elt;
  cont->nelts++;

}

void *pd_pop_container(pd_container_t *cont) 
/* "Pops" the last element added to container, 
   removing the pointer from cont. */
{
  void *lastelt;
  lastelt = cont->elt[cont->nelts-1];
  cont->elt[cont->nelts-1] = NULL;
  cont->nelts--;

  return lastelt;
}

void *pd_container_elt(pd_container_t *cont,pd_contidx_t elt)
/* Reads an element from a container. */

{
  assert(cont != NULL);
  assert(elt < cont->nelts);
  
  return cont->elt[elt];
}

bool pd_container_ok(pd_container_t *cont)
/* Internal sanity checks */
{

  if (PD_VERBOSE > 10) {

    if (cont == NULL) {

      printf("pd_container_ok called on NULL container.\n");
      exit(1);

    }

    if (cont->nelts > cont->bufsize) {

      printf("cont->nelts (%d) > cont->bufsize(%d)\n",cont->nelts,cont->bufsize);
      exit(1);

    }

    if (cont->elt == NULL) {

      printf("cont->elt == NULL\n");
      exit(1);

    }

  } else {

    if (cont == NULL) { return false; }
    if (cont->nelts >= cont->bufsize) { return false; }
    if (cont->elt == NULL) { return false; }

  }

  return true;

}
