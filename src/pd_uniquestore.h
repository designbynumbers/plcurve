/*

   pd_storage.h : 

   When we're generating diagrams, we are going to have
   LOTS of diagrams coming in, which may (rarely) repeat
   hash codes and (even more rarely) actually be
   isomorphic. 

   We will need to search these both by hash/uid and by 
   isomorphism. In each case, we'd like the access process
   to be quite fast. 

   We go with a two-layer strategy:

   1) On top, we have a JudySL indexed by hash strings. 
   2) Inside each entry we have a JudyL indexed by uids.

   We emphasize that it's the hash/uid _pair_ which determines
   a pd code uniquely among pd codes with a given number of 
   crossings.

   This is built on top of libJudy.

*/

#ifndef __PD_STORAGE_H 1
#define __PD_STORAGE_H

#include<judy.h>

/* The top data type is a JudySL pointer. */

typedef struct pdc_struct {

  PVoid_t PJSLArray;
  unsigned int npds;

} pdc_t;

/* A pd_ll_t is a linked list of pd codes (with the same hash). */

typedef struct pd_linked_list pd_ll_t;

struct pd_linked_list {

  pd_code_t *pd;
  pd_ll_t   *next;

} pd_ll_t;

/* A pdc_list_state is a state variable used to index the elements of a pdc
   in order to enumerate the entries (primarily for I/O purposes). */

struct pdc_list_state_struct {

  char      hash[32]; /* Current hash */
  pd_ll_t   *hashmem; /* Member of current linked list */

} pdc_liststate_t;

/* The unique container only contains a few primitives. */

pdc_t          *pd_new_pdc();
void            pd_free_pdc(pdc_t **pdc);

unsigned int   *pdc_nelts(pdc_t *pdc);

void            pd_insert_pdc(pdc_t *pdc,pd_code_t *pd);
/* Given a pd with a hash, inserts <=> it is not isomorphic to an existing pd */

/* Access Functions. 

   There are three ways to access elements in a pdc. First, we can simply loop over 
   everything. This is mostly useful for writing to disk or freeing the entries, but
   it's needed for those purposes. 

   Next, we can search by isomorphism type. This returns both a pointer to the entry
   in the pdc (if found) and a list of explicit isomorphisms between the given pd_code and 
   the entry in the pdc. 

   Last, we can search by hash and uid. This method simply locates a particular uid in the
   structure (if we need to refer back to it later). 
*/   

pd_ll_t        *pdc_firstelt(pdc_t *pdc,pdc_liststate_t *pdls); 
/* Finds the first element of pdc, and initializes a liststate to continue enumeration. */
pd_ll_t        *pdc_nextelt(pdc_t *pdc,pdc_liststate_t *pdls);
/* Can be repeatedly called until it returns NULL to iterate over a pdc. */

pd_ll_t        *pd_search_pdc_by_isomorphism(pdc_t *pdc,pd_code_t *pd,pd_iso_t **isos,unsigned int *nisos); 
pd_ll_t        *pd_search_pdc_by_hash_uid(pdc_t *pdc,char *hash, uint_fast64_t uid);  

#endif
