/* 

   pd_pack.c: Convert pd codes into a very compressed storage format.

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

#ifdef HAVE_STRING_H
  #include<string.h>
#endif

#ifdef HAVE_STDLIB_H
  #include<stdlib.h>
#endif

#ifdef HAVE_STDBOOL_H
  #include<stdbool.h>
#endif

#ifdef HAVE_MATH_H
   #include<math.h>
#endif

#include<plcTopology.h>
#include<pd_container.h>

#include<pd_multidx.h>
#include<pd_dihedral.h>
#include<pd_perm.h>
  
#include<pd_isomorphisms.h>
#include<pd_storage.h>
#include<pd_pack.h>

/* We now define a compressed version of pd codes for storage. The
   basic idea is that we'll need space for a string of n
   crossings. This is going to work in the following way. We're going
   to store an array of "unsigned int"s.

   <first unsigned int> : the number of crossings.

   Remaining unsigned ints hold a bitstring which packs a string of
   bits encoding crossings in the format:

   (edge0 head/tail bit) (edge0 index) ... (edge3 head/tail bit) (edge3 index)

   and the crossings are going to follow, one after the other, until
   we reach n of them.  

   At this point we'll round up to the nearest multiple of (8 bits per
   byte)*(sizeof(unsigned int) bytes) and store the whole enchilada as
   an array of unsigned ints. The number of bits per index is
   determined by the number of edges */

unsigned int pd_indexbits(unsigned int ncrossings)
{
  /* We'll regularly need to know the number of bits in the index variables. */
  return (unsigned int)(ceil(log2((double)(2*ncrossings))));
}

unsigned int pd_packbits(unsigned int ncross,pd_signcross_t crtype)
/* Returns the length (in bits) of the packed format for n crossing codes, not counting
   the size of the initial unsigned int giving the number of crossings. */
{
  return ncross*(4/*sign bits*/ + 4*pd_indexbits(ncross) /* edge indices */ +
		 (crtype == SIGNED_CROSSINGS ? 1 : 0) /* sign bit (if present) */);
}

/* We now give code for unpacking a string of unsigned ints into bits and vice versa */

#define SET_BIT(val, bitIndex)     val |= (1 << bitIndex)
#define CLEAR_BIT(val, bitIndex)   val &= ~(1 << bitIndex)
#define TOGGLE_BIT(val, bitIndex)  val ^= (1 << bitIndex)
#define BIT_IS_SET(val, bitIndex) (val & (1 << bitIndex))

unsigned int *pd_packbitstring(unsigned int *intform,unsigned int nints,unsigned int *packints)
/* Convert a string of unsigned ints, each 1 or 0, into a
   corresponding string of bits occupying a buffer of size
   packints unsigned ints. */
{
  double uintbits = 8.0*(sizeof(unsigned int));
  *packints = (unsigned int)(ceil(nints/uintbits));
  unsigned int *packed = calloc(*packints,sizeof(unsigned int));

  unsigned int thisbit, thisint, thisofs;

  for(thisbit=0,thisint=0;thisint < *packints;thisint++) {

    for(thisofs=0;thisofs < 8*sizeof(unsigned int);thisofs++,thisbit++) {

      if (thisbit == nints) { return packed; }

      assert(intform[thisbit] == 1 || intform[thisbit] == 0);

      if (intform[thisbit] == 1) {

	SET_BIT(packed[thisint],thisofs);

      }

    }

    if (thisbit == nints) { return packed; }

  }

  /* We should never get here */

  assert(1 == 0);
  exit(1);

}

unsigned int *pd_unpackbitstring(unsigned int *bitform,unsigned int nbits)
/* Convert a bitstring with nbits bits into a buffer of 
   nbits unsigned ints, each set to 1 or 0. */
{
  double uintbits = 8.0*(sizeof(unsigned int));
  unsigned int bitformsize = (unsigned int)(ceil(nbits/uintbits));
  unsigned int *unpacked = calloc(nbits,sizeof(unsigned int));

  unsigned int thisbit, thisint, thisofs;

  for(thisbit=0,thisint=0;thisint < bitformsize;thisint++) {

    for(thisofs=0;thisofs < 8*sizeof(unsigned int);thisofs++,thisbit++) {

      if (thisbit == nbits) { return unpacked; }

      if (BIT_IS_SET(bitform[thisint],thisofs)) {

	unpacked[thisbit] = 1;

      } else {

	unpacked[thisbit] = 0;

      }
      
    }

    if (thisbit == nbits) { return unpacked; }

  }

  /* We should never get here */

  assert(1 == 0);
  exit(1);

}
      
unsigned int *pd_idx_bitencode(pd_idx_t i,unsigned int ncrossings)
/* Expands the index i (assuming that it's in 0..2*ncrossings-1) 
   into an array of pd_indexbits unsigned ints, each containing 1 or 0 */
{
  assert(i < 2*ncrossings);
  unsigned int idxbits = pd_indexbits(ncrossings);
  unsigned int k;
  unsigned int *bitcode = calloc(idxbits,sizeof(unsigned int));

  for(k=0;k<idxbits;k++) {

    bitcode[k] = i%2; /* Check odd or even */
    i /= 2;           /* Divide by 2 */         

  }

  /* N.b.: We could have done this by shifting, but we're trying not
     to make assumptions about whether the machine is big-endian or 
     little-endian */

  return bitcode;

}

pd_idx_t pd_idx_bitdecode(unsigned int *bitcode,unsigned int ncrossings)
/* Recover a pd_idx_t from an array of pd_indexbits(ncrossings) unsigned ints,
   each set to 1 or 0 */
{
  pd_idx_t i = 0;
  unsigned int idxbits = pd_indexbits(ncrossings);
  unsigned int k;
  unsigned int placeval = 1;

  for(k=0;k<idxbits;k++) {

    i += bitcode[k]*placeval;
    placeval *= 2;

  }

  return i;
}

unsigned int *pd_crossing_bitencode(pd_idx_t cr,pd_code_t *pd,pd_signcross_t crtype)
/* Encodes a crossing in the format 
   
   (head/tail 0 bit) (edge 0 index) ... (head/tail 3 bit) (edge 3 index)

   as a buffer of unsigned ints if crtype == UNSIGNED_CROSSINGS and in the form

   (sign bit)  (head/tail 0 bit) (edge 0 index) ... (head/tail 3 bit) (edge 3 index)

   if crtype == SIGNED_CROSSINGS.
*/
{
  unsigned int idxbits = pd_indexbits(pd->ncross);
  unsigned int *bitcode = calloc(4 + 4*idxbits,sizeof(unsigned int) + (crtype == SIGNED_CROSSINGS ? 1 : 0));
  unsigned int k,ofs,j;

  if (crtype == SIGNED_CROSSINGS) {

    assert(pd->cross[cr].sign != PD_UNSET_ORIENTATION);
    bitcode[0] = (pd->cross[cr].sign == PD_POS_ORIENTATION ? 1 : 0);
    ofs = 1;
    
  } else {

    ofs = 0;

  }

  for(k=0;k<4;k++,ofs+=1+idxbits) {

    bitcode[ofs] =
      (pd->edge[pd->cross[cr].edge[k]].head == cr && pd->edge[pd->cross[cr].edge[k]].headpos == k) ? 1 : 0;

    unsigned int *idx = pd_idx_bitencode(pd->cross[cr].edge[k],pd->ncross);
    for(j=0;j<idxbits;j++) { bitcode[ofs+1+j] = idx[j]; }
    free(idx);

  }

  return bitcode;

}

void pd_crossing_bitdecode(unsigned int *bitcode,pd_idx_t cr, pd_code_t *pd,pd_signcross_t crtype)
/* Decodes a crossing in the format

   (head/tail 0 bit) (edge 0 idx) ... (head/tail 3 bit) (edge 3 index)

   as a buffer of unsigned ints, each holding 0 or 1 if crtype == UNSIGNED_CROSSINGS
   and in the format 

   (sign bit) (head/tail 0 bit) (edge 0 idx) ... (head/tail 3 bit) (edge 3 index)

   if crtype == SIGNED_CROSSINGS.

   Builds the corresponding edge and crossing records in pd. */
{
  unsigned int idxbits = pd_indexbits(pd->ncross);
  unsigned int k,ofs;

  if (crtype == SIGNED_CROSSINGS) {

    assert(bitcode[0] == 0 || bitcode[0] == 1);
    pd->cross[cr].sign = (bitcode[0] == 1 ? PD_POS_ORIENTATION : PD_NEG_ORIENTATION);
    ofs = 1;
    
  } else {

    ofs = 0;

  }

  for(k=0;k<4;k++,ofs+=idxbits+1) {

    pd_idx_t edge = pd_idx_bitdecode(&(bitcode[ofs+1]),pd->ncross);

    pd->cross[cr].edge[k] = edge;
    
    assert(bitcode[ofs] == 0 || bitcode[ofs] == 1);

    if (bitcode[ofs] == 1) { /* This is the head */

      pd->edge[edge].head = cr;
      pd->edge[edge].headpos = k;

    } else { /* This is the tail */

      pd->edge[edge].tail = cr;
      pd->edge[edge].tailpos = k;

    }

  }
  
}

unsigned int *pd_pack(pd_code_t *pd,unsigned int *packedlength)
/* Returns an array of "packedlength" unsigned ints encoding the pdcode */
{
  /* The first thing we have to check is that all the crossings are 
     signed or all are unset (we don't have a way to bitpack partially 
     set crossings) */

  int k;
  bool crossings_consistent = true;

  if (pd->cross[0].sign == PD_UNSET_ORIENTATION) {

    for(k=1;k<pd->ncross;k++) {

      if (pd->cross[k].sign != PD_UNSET_ORIENTATION) { crossings_consistent = false; }

    }

  } else {

    for(k=1;k<pd->ncross;k++) {

      if (pd->cross[k].sign == PD_UNSET_ORIENTATION) { crossings_consistent = false; }

    }

  }

  assert(crossings_consistent);

  pd_signcross_t crtype;

  if (pd->cross[0].sign == PD_UNSET_ORIENTATION) { crtype = UNSIGNED_CROSSINGS; }
  else { crtype = SIGNED_CROSSINGS; }

  unsigned int  idxbits = pd_indexbits(pd->ncross);
  unsigned int  crossbits = 4*(1+idxbits) + (crtype == SIGNED_CROSSINGS ? 1 : 0); /* Bits per crossing */
  
  unsigned int *intstring = calloc(pd_packbits(pd->ncross,crtype),sizeof(unsigned int));
  int ofs;

  for(ofs=0,k=0;k<pd->ncross;k++,ofs+=crossbits) {

    unsigned int *cr = pd_crossing_bitencode(k,pd,crtype);
    unsigned int j;
    for(j=0;j<crossbits;j++) { intstring[ofs+j] = cr[j]; }
    free(cr);

  }

  unsigned int *packed = pd_packbitstring(intstring,pd->ncross*crossbits,packedlength);
  (*packedlength)++;
  unsigned int *finalcode = calloc(*packedlength,sizeof(unsigned int));

  if (crtype == UNSIGNED_CROSSINGS) {

    finalcode[0] = (unsigned int)(pd->ncross);

  } else { /* crtype == SIGNED_CROSSINGS */

    finalcode[0] = (unsigned int)(-1 * pd->ncross);

  }

  for(k=1;k<*packedlength;k++) { finalcode[k] = packed[k-1]; }

  free(packed); free(intstring);
  return finalcode;
}

int pdint_lowest_edge_unused(bool *edge_used,pd_idx_t nedges)
{
  int k;
  for(k=0;k<nedges;k++) {
    if (!edge_used[k]) { return k; }
  }
  return -1;
}

pd_idx_t pdint_next_edge(pd_code_t *pd,pd_idx_t edge) {

  return pd->cross[pd->edge[edge].head].edge[(pd->edge[edge].headpos + 2) % 4];

}

pd_code_t *pd_unpack(unsigned int *bitform)
/* Decodes a pd_code_t from the packed representation in bitform */
{
  /* If the crossings are unsigned, the (int) casting of bitform[0]
     is positive. If the crossings are signed, the (int) casting 
     of bitform[0] is negative. */
  
  pd_idx_t ncross;
  pd_signcross_t crtype;
  
  if ((int)(bitform[0]) > 0) {

    ncross = (pd_idx_t)(bitform[0]);
    crtype = UNSIGNED_CROSSINGS;

  } else {

    ncross = (pd_idx_t)(-1 * (int)(bitform[0]));
    crtype = SIGNED_CROSSINGS;

  }
  
  pd_code_t *pd = calloc(1,sizeof(pd_code_t));
  pd->ncross = ncross;

  unsigned int nbits = pd_packbits(ncross,crtype);
  unsigned int idxbits = pd_indexbits(ncross);
  unsigned int ofs=0,k;
  unsigned int *intstring = pd_unpackbitstring(&bitform[1],nbits);

  for(ofs=0,k=0;k<pd->ncross;k++,ofs+=4+4*idxbits) {

    pd_crossing_bitdecode(&(intstring[ofs]),k,pd,crtype);

  }

  free(intstring);
  pd->nedges = 2*pd->ncross;
  
  /* We should have (possibly signed) crossings and edges now,
     and it should be ok to call pd_regenerate. */

  pd_regenerate(pd);
  assert(pd_ok(pd));
  
  return pd;
}

    
  
