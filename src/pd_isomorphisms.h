/* 

   pd_isomorphisms.h : Code to find isomorphisms between pd codes. 


*/

#ifndef __PD_ISOMORPHISMS_H__
#define __PD_ISOMORPHISMS_H__ 1

/* The basic idea here is that everything starts with a permutation 
   of the components. However, this is not an arbitrary permutation:
   it must only permute components with the same number of edges. 

   To that end, we start by grouping components by edge number. */

typedef struct pdint_compgrp_struct {

    pd_idx_t ncomps;
    pd_idx_t comp[PD_MAXCOMPONENTS];
  
} pd_compgrp_t;

bool pd_compgrp_ok(pd_code_t *pd,pd_compgrp_t *compgrp);
pd_compgrp_t *pd_build_compgrps(pd_code_t *pdA,pd_code_t *pdB,pd_idx_t *ngrps); 

/* Once we have grouped the components, we can generate
   the list of possible component permutations. Since
   there could be a large number of such permutations, we
   use an unsigned int to store the number of perms
   generated. */

pd_perm_t **pd_build_compperms(pd_idx_t ngrps, pd_compgrp_t *compgrp,
			       unsigned int *ncomp_perms);

void        pd_free_compperms(unsigned int ncomp_perms,pd_perm_t ***comp_perms);
bool        pd_compperms_ok(unsigned int ncomp_perms,pd_perm_t **comp_perms);

/* For a given component permutation, we can generate all
   of the compatible edgemaps by iterating over a multi-index
   of dihedral groups linking components to components. 

   It's going to prove to be convenient later to store the
   orientation of each edge now (even though we could
   theoretically recover it later from the compmap by
   looking up the component each edge was on). */ 

typedef struct pd_edgemap_struct {
  
  pd_perm_t *perm;
  pd_or_t    or[PD_MAXEDGES];

} pd_edgemap_t;

pd_edgemap_t  *pd_new_edgemap(pd_idx_t *nedges);
void           pd_free_edgemap(pd_edgemap_t **edgemap);

char          *pd_print_edgemap(pd_edgemap_t *edgemap);
void          *pd_copy_edgemap(pd_edgemap_t *edgemap);

bool           pd_edgemap_ok(pd_edgemap_t *edgemap);

void           pd_free_edgemaps(unsigned int nedgemaps,pd_edgemap_t ***edgemaps);
int            pd_edgemap_cmp(const void *A,const void *B); /* A,B are **pd_edgemap_t */
bool           pd_edgemaps_ok(unsigned int nedgemaps,pd_edgemap_t **edgemaps);

pd_edgemap_t **pd_build_edgemaps(pd_code_t *pdA,pd_code_t *pdB,
				 pd_perm_t *comp_perm,unsigned int *nedgemaps);

bool           pd_edgemap_consistent(pd_code_t *pdA,pd_code_t *pdB,pd_edgemap_t *edgemap);
/* A slower, stronger check that edgemap takes components
   to components in an orientation-consistent way. */

pd_edgemap_t  *pd_compose_edgemaps(pd_edgemap_t *edgemapA,pd_edgemap_t *edgemapB);
/* When the edgemap maps a pd to itself, we can iterate. 
   Creates a new memory (A * B)(pd) = A(B(pd)) */
void           pd_stareq_edgemap(pd_edgemap_t *edgemapA,pd_edgemap_t *edgemapB);
/* Compose A with B in-place. */

/* A permutation of the edges MAY induce a (unique)
   permutation of the crossings, either preserving cyclic
   orientation at all crossings or reversing cyclic
   orientation at all of the crossings. Such a crossing
   permutation/global orientation combination is called a
   crossmap. */

typedef struct pd_crossmap_struct {

  pd_perm_t *perm; /* Permutation of the crossings. */
  pd_or_t    or;   /* Effect on (global) orientation of plane */

} pd_crossmap_t;

pd_crossmap_t  *pd_new_crossmap(pd_idx_t *ncross);
void            pd_free_crossmap(pd_crossmap_t **crossmap);

char           *pd_print_crossmap(pd_crossmap_t *crossmap);
void           *pd_copy_crossmap(pd_crossmap_t *crossmap);

bool            pd_crossmap_ok(pd_crossmap_t *crossmap);
int             pd_crossmap_cmp(const void *A,const void *B); 
/* Compares **pd_crossmap_t */

pd_crossmap_t **pd_build_crossmaps(pd_code_t *pdA,pd_code_t *pdB,
				   pd_edgemap_t *edgemap,
				   unsigned int *ncrmaps);

void            pd_free_crossmaps(unsigned int ncrmaps,pd_crossmap_t ***crossmaps);

pd_crossmap_t  *pd_compose_crossmaps(pd_crossmap_t *crossmapA,pd_crossmap_t *crossmapB);
/* When a crossmap maps pd->pd, we can iterate. Makes a new-memory (A*B)(pd) = A(B(pd)). */
void            pd_stareq_crossmap(pd_crossmap_t *crossmapA,pd_crossmap_t *crossmapB);
/* Compose A with B in-place. */

/* An edgemap may induce a map from faces to faces.  The
   idea here is essentially identical to the crossing
   case, but a little easier. Since an edge can only occur
   once on a given face, there is a unique canonical
   ordering for the edges on a face. 

   This means that although the face->face maps are
   technically polygon->polygon maps and hence elements of
   a dihedral group, in reality the group element is
   easily recomputed by putting the image in canonical
   order. 

   Thus, just like for a crossmap, we need to store only a
   (global) orientation and a permutation of the face
   indices. 

   SOFTWARE ENGINEERING NOTE: At this point, you're
   wondering why I didn't just replace both the "crossmap"
   and the "facemap" types with some kind of "oriented
   permutation" type. The reason is that precisely because
   they are functionally identical, there's a real danger
   of mixing up the two types, leading to disaster later
   on. We make them different types in order to get
   compiler warnings if we carelessly mix up the two types
   down the road. The price we pay for this is having to 
   rewrite some (simple) primitives twice. 

*/

typedef struct pd_facemap_struct {

  pd_perm_t    *perm;
  pd_or_t       or;               /* Effect on (global) orientation of plane. */

} pd_facemap_t;

pd_facemap_t  *pd_new_facemap(pd_idx_t *ncross);
void           pd_free_facemap(pd_facemap_t **facemap);

char          *pd_print_facemap(pd_facemap_t *facemap);
void          *pd_copy_facemap(pd_facemap_t *facemap);

bool           pd_facemap_ok(pd_facemap_t *facemap);
int            pd_facemap_cmp(const void *A,const void *B);
/* Compare **pd_facemap_t */ 


pd_facemap_t **pd_build_facemaps(pd_code_t *pdA,pd_code_t *pdB,
				 pd_edgemap_t *edgemap,
				 unsigned int *nfacemaps);

void           pd_free_facemaps(unsigned int nfacemaps,pd_facemap_t ***facemaps);

pd_facemap_t  *pd_compose_facemaps(pd_facemap_t *facemapA,pd_facemap_t *facemapB);
/* When a facemap maps pd->pd, we can iterate. Makes a new-memory (A*B)(pd) = A(B(pd)). */
void           pd_stareq_facemap(pd_facemap_t *facemapA,pd_facemap_t *facemapB);
/* Compose A with B in-place. */

/* An isomorphism of pd codes consists of a collection of
 compatible data. It encodes a combinatorial isomorphism
 between the (oriented) polyhedral determined by the
 pdcode which 

 1. is a bijection from the components of one edge graph
    to the components of the other.

 2. maps each component of pdA to a component of pdB by an
    element of the dihedral group of appropriate size.

 3. is a bijection from the crossings of pdA to the
    crossings of pdB, which either preserves all cyclic
    orderings or reverses all cyclic orderings.

 4. is a bijection from the faces of pdA to the faces of
    pdB, mapping each face->face pair by a dihedral group
    element of the appropriate size, and preserving all
    cyclic orientations or reversing all cyclic
    orientations.

*/

typedef struct pd_iso_struct {

  pd_perm_t     *compperm;
  pd_edgemap_t  *edgemap;
  pd_crossmap_t *crossmap;
  pd_facemap_t  *facemap;

} pd_iso_t;

pd_iso_t *pd_new_iso(pd_code_t *pd); 
void      pd_free_iso(pd_iso_t **iso);

pd_iso_t *pd_copy_iso(pd_iso_t *iso); /* Make a new-memory copy of iso. */
char     *pd_print_iso(pd_iso_t *iso);

pd_iso_t **pd_build_isos(pd_code_t *pdA,pd_code_t *pdB,unsigned int *nisos);
/* Returns a buffer of pointers to pd_iso_t. */

bool      pd_iso_ok(pd_iso_t *iso);
bool      pd_iso_consistent(pd_code_t *pdA,pd_code_t *pdB,pd_iso_t *iso);

int       pd_iso_cmp(const void *A,const void *B);
/* Compare **pd_iso_t. */

bool      pd_isos_unique(unsigned int nisos,pd_iso_t **isobuf);
void      pd_free_isos(unsigned int *nisos,pd_iso_t ***isobuf);

bool      pd_iso_is_e(pd_iso_t *iso); /* Check whether iso is identity (for testing) */

/* When isomorphisms are _auto_morphisms (a particularly
   good case for testing), they have a group structure. At
   the moment, this is largely good for testing purposes
   but identifying these groups will later be important
   for the symmetric diagrams project. */

pd_iso_t    *pd_compose_isos(pd_iso_t *A,pd_iso_t *B); /* product automorphism (A * B)(x) = A(B(x)). */
void         pd_stareq_iso(pd_iso_t *A,pd_iso_t *B);  /* A *= B (updates A in-place) */
unsigned int pd_iso_period(pd_iso_t *A); /* computes period of A in isoutation group. */

/* If we're just interested in finding out whether two pd_codes are isomorphic,
   we don't have to store the isomorphisms. (This also appears in pdcode.h.) */

bool      pd_isomorphic(pd_code_t *pdA,pd_code_t *pdB);
 
#endif

