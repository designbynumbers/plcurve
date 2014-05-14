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

   This is built on top of libJudy, which is now rolled
   into the distro. (This may eventually cause trouble, but
   for now it seems to be working relatively well.)

*/

#ifndef __PD_STORAGE_H
#define __PD_STORAGE_H

/* The data type is private b/c it depends on Judy.h */

typedef struct pdstorage_struct pd_stor_t;

/* The pdstor type only contains a few primitives. */

pd_stor_t      *pd_new_pdstor();
void            pd_free_pdstor(pd_stor_t **pdstor);

unsigned int    pd_stor_nelts(pd_stor_t *pdstor);

void            pd_copyinto_pdstor(pd_stor_t *pdstor,pd_code_t *pd);
/*
 Given a pd with a hash, inserts <=> it is not isomorphic
 to an existing pd If the operation succeeds, it makes a
 new-memory copy of pd, so it is always safe to free the
 original pd after this operation. The uid of pd is always
 correct after this operation (regardless of whether we
 found an isomorphic pd or inserted this one).
*/

void            pd_copyinto_cass(pd_code_t *pd);

/* Access Functions.

   There are three ways to access elements in a pdstor. First, we can simply loop over
   everything. This is mostly useful for writing to disk or freeing the entries, but
   it's needed for those purposes.

   Next, we can search by isomorphism type. This returns both a pointer to the entry
   in the pdstor (if found) and a list of explicit isomorphisms between the given pd_code and
   the entry in the pdstor.

   Last, we can search by hash and uid. This method simply locates a particular uid in the
   structure (if we need to refer back to it later).
*/

pd_code_t      *pd_stor_firstelt(pd_stor_t *pdstor);

/* Finds the first element of pdstor, and initializes the internal state of pdstor */
/* Note that this internal state will go stale if an insert or delete operation is performed,
   and so it's reset if we do an insert. */

pd_code_t      *pd_stor_nextelt(pd_stor_t *pdstor);
/* Can be repeatedly called until it returns NULL to iterate over a pdstor */
/* As with pdstor_firstelt, this depends on an internal state which goes stale and
   is reset whenever an insert or delete operation is performed. */

pd_code_t      *pd_search_pdstor_by_isomorphism(pd_stor_t *pdstor,pd_code_t *pd,
						pd_iso_t ***isos,unsigned int *nisos);
/* Return a pointer to the actual stored copy of the unique pd in pdstor which is
   isomorphic to the given pd, along with a buffer of all the isomorphisms from pd
   to the stored copy. Returns NULL if we can't find such a pd in storage. */

pd_code_t      *pd_search_pdstor_by_hash_uid(pd_stor_t *pdstor,char *hash, pd_uid_t uid);

/* Not implemented yet. */

/* File I/O */

void pd_write_pdstor(FILE *stream,pd_stor_t *pdstor);
pd_stor_t *pd_read_pdstor(FILE *stream);

/* Debugging code */

void pd_stor_stats(pd_stor_t *pdstor,unsigned int *nhashes,unsigned int *nelts);
void pd_display_pdstor(FILE *stream,pd_stor_t *pdstor); /* Prints a representation to stream. */




#endif
