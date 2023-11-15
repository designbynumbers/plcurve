/* 

   pd_container.h : 

   This is a general container type (which we can
   elaborate on as needed), used to store items where we
   don't know in advance how many will be generated. They
   are not guaranteed to be particularly quickly
   searchable, as this is intended to be used for
   storage. (If we need fast access, we should use a hash
   or tree.)

*/

#ifndef __PD_CONTAINER_H__
#define __PD_CONTAINER_H__ 1

#ifndef PD_VERBOSE
#define PD_VERBOSE 0
#endif

typedef unsigned int pd_contidx_t;  /* Index to elements in a pd_container_t. */

typedef void * pd_contdata_t; 

typedef struct pd_container_struct {

  pd_contidx_t  nelts;

  pd_contidx_t  bufsize;
  pd_contidx_t  initsize;

  pd_contdata_t *elt;  /* Array of void ptrs */
  
} pd_container_t;

pd_container_t *pd_new_container(pd_contidx_t initsize);
/* Allocates space for a new pd_container holding initsize pointers. */
/* This should be a reasonable guess of the number of elements expected */
/* to go into the container */

void  pd_free_container(pd_container_t **cont,void eltfree(void **ptr));
/* Frees the memory in a pd_container, calling eltfree on each stored pointer if eltfree != NULL,
   or free on each (non-NULL) stored pointer otherwise. */

void  pd_free_fake(void **ptr); 
/* A fake function used to call pd_free_container when the container should be empty: 
   terminates immediately with an error message if called. */

pd_contidx_t pd_container_nelts(pd_container_t *cont); 
/* Returns the number of elements in the container. */

void  pd_addto_container(pd_container_t *cont,void *elt);
/* Adds an element to a container. */

void *pd_pop_container(pd_container_t *cont); 
/* "Pops" the last element added to container, 
   removing the pointer from cont. */

void *pd_container_elt(pd_container_t *cont,pd_contidx_t i);
/* Reads an element from a container. */

bool pd_container_ok(pd_container_t *cont);
/* Internal sanity checks. If we pass, cont is ready to add or read elts. */

#endif
