/*

   pd_pack.h : 

   For large runs, we may need to keep track of more information than
   fits in memory. For this reason, we provide a bit-packed format for
   pd codes which reduces the storage requirements to the absolute
   minimum. This will also help when we're keeping large disk files.

   The basic idea is that we'll need space for a string of n
   crossings. This is going to work in the following way. We're going
   to store an array of "unsigned int"s.

   <first unsigned int> : should be converted to signed int

   +the number of crossings if we are NOT encoding crossing signs
   -the number of crossings if we ARE encoding crossing signs

   Remaining unsigned ints hold a bitstring which packs a string of
   bits encoding crossings in the format:

   (edge0 head/tail bit) (edge0 index) ... (edge3 head/tail bit) (edge3 index)

   if we are NOT encoding crossing signs and

   (sign bit) (edge0 head/tail bit) (edge0 index) ... (edge3 head/tail bit) (edge3 index)

   if we ARE encoding crossing signs.

   and the crossings are going to follow, one after the other, until
   we reach n of them. The number of bits per index is determined by the
   number of edges, and is the smallest number needed to represent the 
   edges we can have.

   At this point we'll round up to the nearest multiple of (8 bits per
   byte) * (sizeof(unsigned int) in bytes) and store the whole
   enchilada as an array of unsigned ints.
*/

#ifndef __PD_PACK_H 
#define __PD_PACK_H 1

typedef enum {SIGNED_CROSSINGS, UNSIGNED_CROSSINGS} pd_signcross_t;

unsigned int pd_size(pd_code_t *pd);
/* Total memory (bytes) used by this pd code */

unsigned int pd_indexbits(unsigned int ncrossings);
/* # of bits to store an edge index in an ncrossings pd code */

unsigned int pd_packbits(unsigned int ncrossm,pd_signcross_t crtype);
/* Returns the length (in bits) of the packed format for n crossing
   codes, not counting the size of the initial unsigned int giving the
   number of crossings. */

unsigned int *pd_packbitstring(unsigned int *intform,
			       unsigned int nints,unsigned int *packints);
/* Convert a string of unsigned ints, each 1 or 0, into a
   corresponding string of bits occupying a buffer of packints
   unsigned ints. */

unsigned int *pd_unpackbitstring(unsigned int *bitform,unsigned int nbits);
/* Convert a bitstring with nbits bits into a buffer of 
   nbits unsigned ints, each set to 1 or 0. */

unsigned int *pd_idx_bitencode(pd_idx_t i,unsigned int ncrossings);
/* Expands the index i (assuming that it's in 0..2*ncrossings-1) 
   into an array of pd_indexbits unsigned ints, each containing 1 or 0 */

pd_idx_t pd_idx_bitdecode(unsigned int *bitcode,unsigned int ncrossings);
/* Recover a pd_idx_t from an array of pd_indexbits(ncrossings) unsigned ints,
   each set to 1 or 0 */

unsigned int *pd_crossing_bitencode(pd_idx_t cr,pd_code_t *pd,pd_signcross_t crtype);
/* Encodes a crossing in the format  

   (head/tail 0 bit) (edge 0 index) ... (head/tail 3 bit) (edge 3 index)

   as a buffer of unsigned ints if the crtype is UNSIGNED_CROSSINGS, and as

   (sign bit) (head/tail 0 bit) (edge 0 index) ... (head/tail 3 bit) (edge 3 index)

   if the crtype is SIGNED_CROSSINGS.

   Note that the head/tail bits are always 1:head, 0:tail.
*/

void pd_crossing_bitdecode(unsigned int *bitcode,pd_idx_t cr,pd_code_t *pd,pd_signcross_t crtype);
/* Decodes a crossing in the format

   (head/tail 0 bit) (edge 0 idx) ... (head/tail 3 bit) (edge 3 index)

   as a buffer of unsigned ints, each holding 0 or 1 if the crtype is UNSIGNED_CROSSINGS,
   and in the format 

   (sign bit) (head/tail 0 bit) (edge 0 index) ... (head/tail 3 bit) (edge 3 index)

   if the crtype is SIGNED_CROSSINGS. Builds the corresponding edge
   and crossing records in pd. */

unsigned int *pd_pack(pd_code_t *pd,unsigned int *packedlength);
/* Returns an array of "packedlength" unsigned ints encoding the pdcode */

pd_code_t    *pd_unpack(unsigned int *bitform);
/* Decodes a pd_code_t from the packed representation in bitform */

#endif
