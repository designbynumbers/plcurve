/* 

   Code for handling pd_codes (link shadows). Part of the
   general COLD (Census of Link Diagrams) project.  

   Jason Cantarella, September 2012.

*/

int PD_VERBOSE;

#ifdef HAVE_CONFIG_H
  #include"config.h"
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

#ifdef HAVE_STDINT_H
   #include<stdint.h>
#endif

#ifdef HAVE_STDBOOL_H
   #include<stdbool.h>
#endif

#ifdef HAVE_STDARG_H
   #include<stdarg.h>
#endif

#ifdef HAVE_ASSERT_H
   #include<assert.h>
#endif

#ifdef HAVE_CTYPE_H
   #include<ctype.h>
#endif

#include<pdcode.h>
//#include<libcassie/cassie.h>
//#include"/usr/local/include/thrift/Thrift.h"
//#include<python2.7/Python.h>
#include<pd_multidx.h>
#include<pd_dihedral.h>
#include<pd_perm.h>

#include<pd_isomorphisms.h>

pd_or_t pd_compose_or(pd_or_t a,pd_or_t b)
/* Returns the composition of the two orientation changes: ++ = -- = +, +- = -+ = - */
{
  if ((a == PD_POS_ORIENTATION && b == PD_POS_ORIENTATION) ||
      (a == PD_NEG_ORIENTATION && b == PD_NEG_ORIENTATION)) { 

    return PD_POS_ORIENTATION;

  } else if ((a == PD_POS_ORIENTATION && b == PD_NEG_ORIENTATION) ||
	     (a == PD_NEG_ORIENTATION && b == PD_POS_ORIENTATION)) {
    
    return PD_NEG_ORIENTATION;

  }

  pd_error(SRCLOC,"pd_compose_or called with orientations %c %c \n",NULL,pd_print_or(a),pd_print_or(b));
  exit(1);
}

bool pd_or_ok(pd_or_t or) /* Check whether or has a legal value. */
{
  return (or == PD_POS_ORIENTATION || or == PD_NEG_ORIENTATION);
}

int pd_or_cmp(const void *A, const void *B)
/* Compare *pd_or_t. */
{
  pd_or_t orA = *(pd_or_t *)(A);
  pd_or_t orB = *(pd_or_t *)(B);

  assert(pd_or_ok(orA) && pd_or_ok(orB));

  if (orA == orB) { return 0; }

  else if (orA == PD_POS_ORIENTATION) { return -1; }

  else if (orB == PD_POS_ORIENTATION) { return 1; }

  /* We should never get here. */

  pd_error(SRCLOC,"pd_or_cmp passed illegal orientation pair orA == %OR, orB == %OR",NULL,
	   orA,orB);
  exit(1);
}
    
/* These functions appear in pdcode.c, but they are not exposed in the header. */

pd_pos_t pdint_cross_rdist(pd_crossing_t *cr);
void pdint_rotate_left(pd_crossing_t *cr,pd_pos_t k);


pd_crossing_t pd_build_cross(pd_idx_t e0,pd_idx_t e1,pd_idx_t e2,pd_idx_t e3) 
  
/* Builds a pd_crossing_t from edge numbers */

{
  pd_crossing_t cr;

  cr.edge[0] = e0; cr.edge[1] = e1; cr.edge[2] = e2; cr.edge[3] = e3;
  
  return cr;
}

void pd_canonorder_cross(pd_crossing_t *cross, pd_or_t or)
/* Reverses (if or == PD_NEG_ORIENTATION) and then rotates
   a crossing into canonical order: cr->edge[0] is the
   lowest edge # */
{
  if (or == PD_NEG_ORIENTATION) { /* Reverse the cyclic order of the edges */

    pd_idx_t rev[4];
    pd_pos_t pos;

    for(pos=0;pos<4;pos++) { rev[pos] = cross->edge[3-pos]; }
    for(pos=0;pos<4;pos++) { cross->edge[pos] = rev[pos]; }

  }
  
  pd_pos_t cross_rotate;

  cross_rotate = pdint_cross_rdist(cross);
  pdint_rotate_left(cross,cross_rotate);

}

void pd_canonorder_face(pd_face_t *face, pd_or_t or)
/* Reverses (if or == PD_NEG_ORIENTATION) and then rotates
   a face into canonical order: face->edge[0] is the
   lowest edge # */
{
  pd_face_t reface;
  pd_idx_t  edge;
  pd_idx_t  nedges;

  /* First we copy (or reverse-copy) the face data into reface. */

  nedges = reface.nedges = face->nedges;

  if (or == PD_POS_ORIENTATION) { 

    memcpy(reface.edge,face->edge,nedges*sizeof(pd_idx_t));
    memcpy(reface.or,face->or,nedges*sizeof(pd_or_t));
    
  } else if (or == PD_NEG_ORIENTATION) { /* There's no library function for reverse-copy */

    for(edge=0;edge < nedges;edge++) { 
      
      reface.edge[edge] = face->edge[(nedges-1) - edge];
      reface.or[edge] = face->or[(nedges-1) - edge];

    }

  }
  
  /* First, search for the position of the lowest edge # */

  pd_idx_t lowE = PD_MAXEDGES+1,lowPos = 0;
  
  for(edge=0;edge<nedges;edge++) { 
    
    if (reface.edge[edge] <= lowE) {
      
      lowE = reface.edge[edge]; lowPos = edge;
      
    }
    
  }
    
  /* Now copy this back into face->edge with the lowPos offset. */

  for(edge=0;edge<nedges;edge++) { 

    face->edge[edge] = reface.edge[(lowPos + edge) % nedges];
    face->or[edge] = reface.or[(lowPos + edge) % nedges];

  }

}

int  pd_cross_cmp(const void *A,const void *B) 

/* Compares crossings in dictionary order. */

{
  pd_crossing_t *crA = (pd_crossing_t *)(A);
  pd_crossing_t *crB = (pd_crossing_t *)(B);
  pd_pos_t pos;

  for(pos=0;pos<4;pos++) {
    
    if (crA->edge[pos] != crB->edge[pos]) {
      
      return crA->edge[pos] - crB->edge[pos];
      
    }
    
  }
  
  return 0;

}

int  pd_face_cmp(const void *A, const void *B) 

/* Compares faces by number of components (largest number first) then dictionary order */

{
  pd_face_t *fA = (pd_face_t *)(A);
  pd_face_t *fB = (pd_face_t *)(B);

  if (fA->nedges != fB->nedges) {

    return fB->nedges - fA->nedges;

  }

  pd_idx_t nedges = fA->nedges,edge;

  for(edge=0;edge<nedges;edge++) {

    if (fA->edge[edge] != fB->edge[edge]) {

      return fA->edge[edge] - fB->edge[edge];

    }

  }

  return 0;
 
}

void pd_component_and_pos(pd_code_t *pd, pd_idx_t edge,
			  pd_idx_t *compP,pd_idx_t *comp_posP)

/* Returns component number and position on component of a
   given edge number (or die) */

{
  pd_idx_t comp, comp_pos;
  
  assert(pd != NULL);
  assert(edge < pd->nedges);
  
  for(comp=0;comp < pd->ncomps;comp++) {

    for(comp_pos=0;comp_pos < pd->comp[comp].nedges;comp_pos++) {

      if (pd->comp[comp].edge[comp_pos] == edge) { 

	if (compP != NULL) { *compP = comp; }
	if (comp_posP != NULL) { *comp_posP = comp_pos; }
	return; 

      }

    }

  }

  fprintf(stderr,"pd_component_and_pos: Couldn't find edge %d anywhere in pd code.\n\n",edge);
  pd_write(stderr,pd);
  exit(1);

}

pd_edge_t pd_oriented_edge(pd_edge_t e,pd_or_t or)
/* Returns original edge if or = PD_POS_ORIENTATION, 
   reversed edge if or = PD_NEG_ORIENTATION */

{ pd_edge_t ret;

  assert(or == PD_POS_ORIENTATION || or == PD_NEG_ORIENTATION);

  if (or == PD_POS_ORIENTATION) {

    ret = e;

  } else if (or == PD_NEG_ORIENTATION) {

    ret.tail = e.head;
    ret.tailpos = e.headpos;

    ret.head = e.tail;
    ret.headpos = e.tailpos;

  } 

  return ret;

}

void pd_reorient_edge(pd_code_t *pd,pd_idx_t edge,pd_or_t or)
/* Flips the edge in pd->edges[] if or == PD_NEG_ORIENTATION */
{
  assert(pd != NULL);
  assert(edge < pd->nedges);
  assert(or == PD_POS_ORIENTATION || or == PD_NEG_ORIENTATION);
  
  if (or == PD_NEG_ORIENTATION) {

    pd_idx_t cross_swap;
    
    cross_swap = pd->edge[edge].head; 
    pd->edge[edge].head = pd->edge[edge].tail;
    pd->edge[edge].tail = cross_swap;

    pd_pos_t pos_swap;

    pos_swap = pd->edge[edge].headpos;
    pd->edge[edge].headpos = pd->edge[edge].tailpos;
    pd->edge[edge].tailpos = pos_swap;

  }

}

/* pd "regeneration" functions. These help a pd code "grow back" 
   auxiliary data after it's "damaged" by some operation such as 
   a Reidemeister move or other simplication. */

pd_code_t *pdint_cross_index_cmp_glob;

int  pdint_cross_index_cmp(const void *A, const void *B)

/* Dictionary order on the edges in the cross records. */

{
  pd_idx_t *idxA = (pd_idx_t *)(A);
  pd_idx_t *idxB = (pd_idx_t *)(B);

   pd_code_t *pd = pdint_cross_index_cmp_glob;

   return pd_cross_cmp(&(pd->cross[*idxA]),&(pd->cross[*idxB]));
}

void pdint_rotate_left(pd_crossing_t *cr,pd_pos_t k) 

/* Rotates the buffer left by k spots */

{
  pd_crossing_t ncr;
  pd_pos_t      pos;

  for(pos=0;pos<4;pos++) {

    ncr.edge[pos] = cr->edge[(pos+k)%4];

  }

  *cr = ncr;
}

pd_pos_t pdint_cross_rdist(pd_crossing_t *cr)
 
/* Detects how many spots to rotate left in order to put 
   lowest edge number first. */

{
  pd_pos_t lowpos, pos;
  pd_idx_t lowedge;

  for(lowpos=0,lowedge=cr->edge[0],pos=1;pos<4;pos++) {
      
    if (cr->edge[pos] < lowedge) {

      lowpos = pos; lowedge = cr->edge[pos];

    }
    
  }
  
  /* Make sure we have the "leftmost" copy. */
  
  if (cr->edge[(lowpos+3)%4] == lowedge) { lowpos = (lowpos + 3) % 4; }

  return lowpos;

}

  
void pd_regenerate_crossings(pd_code_t *pd)

/* Cyclically reorders pd codes so lowest index edge 
   comes first, fixing edgecodes as we go, sorts crossings. */

{
  assert(pd != NULL);

  /* First, decide how far left every crossing needs to be rotated
     and update the crossing records accordingly. */

  pd_idx_t cross;
  pd_pos_t cross_rotate[PD_MAXVERTS];

  for(cross=0;cross<pd->ncross;cross++) {
    
    cross_rotate[cross] = pdint_cross_rdist(&(pd->cross[cross]));
    pdint_rotate_left(&(pd->cross[cross]),cross_rotate[cross]);

  }

  /* Now update the "headpos" and "tailpos" records in the edge
     data structures using the stored cross_rotate values. */

  pd_idx_t edge;

  for(edge=0;edge<pd->nedges;edge++) {

    pd_edge_t *e = &(pd->edge[edge]);
    
    e->headpos = (e->headpos - cross_rotate[e->head] + 4) % 4;
    e->tailpos = (e->tailpos - cross_rotate[e->tail] + 4) % 4;

  }

  /* Now sort the crossings. We'll need to do this indirectly, storing */
  /* a permutation of the crossings and then applying it, because the  */
  /* edge records refer to crossing numbers in their "head" and "tail" */
  /* fields and we'll need to update these as well. */

  pdint_cross_index_cmp_glob = pd;
  pd_idx_t cross_perm[PD_MAXVERTS];
  for(cross=0;cross<pd->ncross;cross++) { cross_perm[cross] = cross; }
  qsort(cross_perm,pd->ncross,sizeof(pd_idx_t),pdint_cross_index_cmp);

  pd_idx_t new_cross_num[PD_MAXVERTS];
  for(cross=0;cross<pd->ncross;cross++) { new_cross_num[cross_perm[cross]] = cross; }


  /* First, we change the references to crossings in the edge records. */

  for(edge=0;edge<pd->nedges;edge++) {
    
    pd->edge[edge].head = new_cross_num[pd->edge[edge].head];
    pd->edge[edge].tail = new_cross_num[pd->edge[edge].tail];

  }

  /* Next, we rearrange the crossings themselves. */

  pd_crossing_t newcross[PD_MAXVERTS];

  for(cross=0;cross<pd->ncross;cross++) { newcross[cross] = pd->cross[cross_perm[cross]]; }
  for(cross=0;cross<pd->ncross;cross++) { pd->cross[cross] = newcross[cross]; }

}

void pdint_next_comp_edge(pd_code_t *pd,pd_idx_t edge,pd_idx_t *nxt,pd_or_t *nxt_or)

/* Assuming that this edge is oriented positively, find the next edge
   (the one across from the head of this edge) and its orientation. */

{
  assert(pd != NULL && nxt != NULL && nxt_or != NULL);
  assert(edge < pd->nedges);

  pd_idx_t cross = pd->edge[edge].head;
  pd_pos_t pos   = pd->edge[edge].headpos;
  pd_pos_t nxtpos = (pos+2)%4;

  *nxt = pd->cross[cross].edge[nxtpos];
  
  if (cross == pd->edge[*nxt].tail && nxtpos == pd->edge[*nxt].tailpos) { 
    *nxt_or = PD_POS_ORIENTATION; 
  } else if (cross == pd->edge[*nxt].head && nxtpos == pd->edge[*nxt].headpos) {
    *nxt_or = PD_NEG_ORIENTATION; 
  } else {
    pd_error(SRCLOC,"Mismatch between %CROSS, which has edge %d at position %d"
	     "and record %EDGE which does not match either head/headpos or tail/tailpos",
	     pd,cross,*nxt,nxtpos,*nxt);
    exit(1);
  }
}

pd_idx_t pdint_find_unused(pd_code_t *pd,bool *edge_used, bool *all_used)

{
  assert(pd != NULL && edge_used != NULL && all_used != NULL);

  pd_idx_t edge;
  
  for(edge=0;edge<pd->nedges;edge++) {

    if (!edge_used[edge]) { 

      *all_used = false; 
      edge_used[edge] = true;
      return edge; 

    }

  }

  *all_used = true;
  return 0;

}

int  pd_component_cmp(const void *A,const void *B) {

  /* Sorts components first by number of edges (longest FIRST),
     then by dictionary order on the edge numbers. */

  pd_component_t *pdcA = (pd_component_t *)(A);
  pd_component_t *pdcB = (pd_component_t *)(B);

  if (pdcA->nedges != pdcB->nedges) {

    return pdcB->nedges - pdcA->nedges;

  } 

  pd_idx_t edge,nedges = pdcA->nedges;

  for(edge=0;edge<nedges;edge++) {

    if (pdcA->edge[edge] != pdcB->edge[edge]) {

      return pdcA->edge[edge] - pdcB->edge[edge];

    }

  }

  return 0;

}
  

void pd_regenerate_comps(pd_code_t *pd) 

/* This one is somewhat complicated, as it involves
   an edge renumbering, and so goes all the way back 
   to changing vertex data. */

{
  assert(pd != NULL);
  assert(pd->ncross <= PD_MAXVERTS);
  
  /* Step 0: Generate SOME list of edges to get off the
     ground with using only crossing information. We don't
     know that ANYTHING about the pd is good except for cross
     and ncross yet, so we reset nedges first. */
  
  pd_idx_t cross,edge;
  bool     hit[PD_MAXEDGES];

  pd->nedges = pd->ncross*2;

  for(edge=0;edge<pd->nedges;edge++) { hit[edge] = false; };
  
  for(cross=0;cross<pd->ncross;cross++) {

    pd_pos_t pos;

    for(pos=0;pos<4;pos++) {

      edge = pd->cross[cross].edge[pos];

      if (!hit[edge]) { 

	pd->edge[edge].tail = cross;
	pd->edge[edge].tailpos = pos;
	hit[edge] = true;

      } else {

	pd->edge[edge].head = cross;
	pd->edge[edge].headpos = pos;

      }

    }

  }

  bool edge_used[PD_MAXEDGES];
  bool all_used = false;

  /* Step 1. Run around the components, assembling the (old) edge
     numbers of the edges, and reorienting edges as we go. */

  for(edge=0;edge<pd->nedges;edge++) { edge_used[edge] = false; }

  edge = 0;                /* Now we start with edge 0 */
  pd->ncomps = 0;
  edge_used[edge] = true;

  for(;!all_used && pd->ncomps < PD_MAXCOMPONENTS+1;
      pd->ncomps++,edge=pdint_find_unused(pd,edge_used,&all_used)) {

    pd_idx_t comp = pd->ncomps;
    pd_or_t  next_or;
    pd_idx_t next_edge;

    /* We are now starting a new component with an unused edge. */

    pd->comp[comp].nedges = 1;
    pd->comp[comp].edge[0] = edge;
    pdint_next_comp_edge(pd,edge,&next_edge,&next_or);

    for(;next_edge != pd->comp[comp].edge[0] && pd->comp[comp].nedges < PD_MAXEDGES + 1;
	pdint_next_comp_edge(pd,edge,&next_edge,&next_or),pd->comp[comp].nedges++) {

      pd_reorient_edge(pd,next_edge,next_or);
      pd->comp[comp].edge[pd->comp[comp].nedges] = next_edge;
      edge_used[next_edge] = true;
      edge = next_edge;

    }

    assert(pd->comp[comp].nedges <= PD_MAXEDGES);

  }

  assert(pd->ncomps <= PD_MAXCOMPONENTS);

  /* Step 2. Sort the components longest-first. */

  qsort(pd->comp,pd->ncomps,sizeof(pd_component_t),pd_component_cmp);

  /* Step 3. Build a translation table for edge reordering and update the
     comp records to standard edge numbers. */

  pd_idx_t old_edge_num[PD_MAXEDGES];
  pd_idx_t new_edge_num[PD_MAXEDGES];

  pd_idx_t comp;
  pd_idx_t i;
  
  for(i=0,comp=0;comp<pd->ncomps;comp++) {
    
    for(edge=0;edge<pd->comp[comp].nedges;edge++,i++) {
      
      old_edge_num[i] = pd->comp[comp].edge[edge];
      new_edge_num[pd->comp[comp].edge[edge]] = i;

      pd->comp[comp].edge[edge] = i; /* This SHOULD be the edge here. */
    
    }
    
  }
  
  /* Step 4. Rearrange the edges buffer accordingly. */
  
  pd_edge_t new_edges[PD_MAXEDGES];
  
  for(edge=0;edge<pd->nedges;edge++) {
    
    new_edges[edge] = pd->edge[old_edge_num[edge]];
    
  }
  
  for(edge=0;edge<pd->nedges;edge++) {

    pd->edge[edge] = new_edges[edge];

  }

  /* Step 5. Rewrite the crossing data with the
     new edge numbers. */

  for(cross=0;cross<pd->ncross;cross++) { 

    pd_pos_t pos;

    for(pos=0;pos<4;pos++) {

      pd->cross[cross].edge[pos] = new_edge_num[pd->cross[cross].edge[pos]];

    }

  }

  /* Step 6. Resort the crossings. */

  pd_regenerate_crossings(pd);
  
}

/*
** Translation Table as described in RFC1113
*/
static const char cb64[]="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

/*
** encodeblock
**
** encode 3 8-bit binary bytes as 4 '6-bit' characters
*/
static void encodeblock( unsigned char *in, char *out, int len )
{
    out[0] = (char) cb64[ (int)(in[0] >> 2) ];
    out[1] = (char) cb64[ (int)(((in[0] & 0x03) << 4) | ((in[1] & 0xf0) >> 4)) ];
    out[2] = (char) (len > 1 ? cb64[ (int)(((in[1] & 0x0f) << 2) | ((in[2] & 0xc0) >> 6)) ] : '=');
    out[3] = (char) (len > 2 ? cb64[ (int)(in[2] & 0x3f) ] : '=');
}

void base64encode(unsigned char data[24],char *encoded) {

  int data_idx,encoded_idx;

  for(data_idx=0,encoded_idx=0;data_idx<24;data_idx+=3,encoded_idx+=4) {
    
    encodeblock(&(data[data_idx]),&(encoded[encoded_idx]),3);
    
  }

} 

void pd_regenerate_hash(pd_code_t *pd)

/* Fills in a (printable) hash string from component, face, crossing information. */
/* The idea is that the hash value should not depend on the numbering of faces, crossings, etc */
/* and that it should be computable from the pd_code (that is, it cannot depend on crossing sign info */

/* So we start with a string of consecutive unsigned chars: 
		     
   nverts
   nedges
   nfaces

   edges around face 1 (face with the largest number of edges)
   ...
   edges around face (2 + nverts) (the last face, with the smallest number of edges)

   ncomponents
   edges in comp 1 (longest component)
   ...
   edges in comp n (shortest component)

   and then base64encode to make it printable. 
   
*/

{
  assert(pd != NULL);

  unsigned char data[24];
  int i,face,comp;
  
  data[0] = (unsigned char)(pd->ncross);
  data[1] = (unsigned char)(pd->nedges);
  data[2] = (unsigned char)(pd->nfaces);

  for(i=3,face=0;face<pd->nfaces && i < 24;i++,face++) {

    data[i] = (unsigned char)(pd->face[face].nedges);

  }

  if(i < 24) { data[i++] = (unsigned char)(pd->ncomps); }
  
  for(comp=0;comp<pd->ncomps && i < 24;i++,comp++) {
    
    data[i] = (unsigned char)(pd->comp[comp].nedges);
    
  }
  
  for(;i < 24;i++) { data[i] = 0; } /* Pad any remaining space with zeros. */
  
  base64encode(data,pd->hash);

  pd->hash[31] = 0; /* Make sure that we end the string with a zero. */
  
}


void pdint_next_edge_on_face(pd_code_t *pd,pd_idx_t ed,pd_or_t or,
			     pd_idx_t *next_edge,pd_or_t *next_or) 
  
/* Given an edge and an orientation (or), if or = +1, travel to the head of the edge,
   turn left, and record the next edge arrived at (with orientation), if or = -1, travel
   to the tail of the edge, turn left, and record next edge (with orientation). */
  
{
  assert(pd != NULL && next_edge != NULL && next_or != NULL);
  assert(ed < pd->nedges);
  assert(or == PD_POS_ORIENTATION || or == PD_NEG_ORIENTATION);

  pd_idx_t cross, this_pos, next_pos;
  
  if (or == PD_POS_ORIENTATION) {
    
    cross    = pd->edge[ed].head;
    this_pos = pd->edge[ed].headpos;
    
  } else { /* or == PD_NEG_ORIENTATION */
    
    cross    = pd->edge[ed].tail;
    this_pos = pd->edge[ed].tailpos;
    
  } 
  
  next_pos = (this_pos + 3) % 4; /* Turn LEFT (clockwise) */
  (*next_edge) = pd->cross[cross].edge[next_pos];

  /* Now we check data integrity */

  pd_edge_t *nxt = &(pd->edge[*next_edge]);

  assert(nxt->head == cross || nxt->tail == cross);

  if (nxt->head == cross && nxt->tail != cross) {

    assert(nxt->headpos == next_pos);

  } else if (nxt->tail == cross && nxt->head != cross) {

    assert(nxt->tailpos == next_pos);

  } else { /* Both must be equal */

    assert(nxt->headpos == next_pos || nxt->tailpos == next_pos);

  }	   
  
  /* We need to determine the orientation of this (next) edge on the face. */
    
  if (nxt->head == cross && nxt->headpos == next_pos) {
    
    *next_or = PD_NEG_ORIENTATION; /* This is the _head_ of the next edge */
    
  } else { 
    
    *next_or = PD_POS_ORIENTATION; /* This was the _tail_ of the next edge */
    
  }
  
}

void pd_regenerate_faces(pd_code_t *pd) 

/* Generates faces of the polyhedron corresponding to the
   crossing and edge data in pd.  Add that data to the
   pd_code_t structure pd. 

   The basic idea here is that we're computing the orbits
   of pdint_next_edge_on_face. To do that, we need to keep track
   of the edges and orientations we've used already. */

{
  assert(pd != NULL);

  bool pos_or_used[PD_MAXEDGES];
  bool neg_or_used[PD_MAXEDGES];

  pd_idx_t edge;

  for(edge=0;edge < pd->nedges;edge++) { 

    pos_or_used[edge] = false;
    neg_or_used[edge] = false;
  
  }
  
  pd->nfaces = 0;

  /* Start the main (face-counting) loop. */

  int used;
 
  for(used=0;used < 2*pd->nedges;pd->nfaces++) { /* Each edge occurs twice in a face list. */

    pd_idx_t this_edge;
    pd_or_t  this_or; 
    
    /* Scan for an unused edge/orientation pair. */
    
    for(this_or=PD_POS_ORIENTATION,this_edge = 0;
	this_edge < pd->nedges;this_edge++) { if (!pos_or_used[this_edge]) { break; }    }
    
    if (this_edge == pd->nedges) { /* Positives are used; scan for a neg orientation unused. */
      
      for(this_or=PD_NEG_ORIENTATION,this_edge = 0;
	  this_edge < pd->nedges;
	  this_edge++) { if (!neg_or_used[this_edge]) { break; } }
      
    }
    
    assert(this_edge < pd->nedges && 
	   (this_or == PD_POS_ORIENTATION || this_or == PD_NEG_ORIENTATION));
    
    /* Now we have an unused edge/orientation pair to start 
       the iteration on. */
    
    int start_edge, start_or;
    
    /* Add initial edge to face and mark as used in edge_checklist */
    
    pd_face_t *face = &(pd->face[pd->nfaces]); /* A convenience variable */
    
    start_edge = this_edge; start_or = this_or;
    
    face->nedges = 1;
    face->edge[0] = start_edge;
    face->or[0] = start_or;
    
    if (start_or == PD_POS_ORIENTATION) { pos_or_used[start_edge] = true; }
    else { neg_or_used[start_edge] = true; }

    used++; /* We have just used an edge/or pair */
    
    /* Now increment to next edge and start working our way around the face. */
    
    int itcount = 0; /* Add a safety to make sure we don't just hang here if something is very wrong. */
    
    for(pdint_next_edge_on_face(pd,start_edge,start_or,&this_edge,&this_or);
	!(this_edge == start_edge && this_or == start_or) && itcount < 4*pd->nedges;
	pdint_next_edge_on_face(pd,this_edge,this_or,&this_edge,&this_or),itcount++) {

      if (this_or == PD_POS_ORIENTATION) { pos_or_used[this_edge] = true; }
      else { neg_or_used[this_edge] = true; }
      used++; /* We have used another edge/or pair. */

      face->or[face->nedges] = this_or;
      face->edge[face->nedges] = this_edge; /* Add to this face. */
      face->nedges += 1;
      
    }
    
    assert(itcount < 4*pd->nedges); /* Put in an assert to make sure we didn't bail out. */

    /* Make sure the face is in canonical order. */

    pd_canonorder_face(face,PD_POS_ORIENTATION);
    
  }
  
  /* Now sort the faces using pd_face_cmp. */
  
  qsort(pd->face,pd->nfaces,sizeof(pd_face_t),pd_face_cmp);

  if (PD_VERBOSE > 10) {

    if (pd->nfaces != pd->ncross + 2) {

      fprintf(stderr,"pd->nfaces (%d) does not equal pd->ncross + 2 (%d) for pd\n\n",
	      pd->nfaces,pd->ncross);

      pd_write(stderr,pd);

      exit(1);

    }

  } else {

    assert(pd->nfaces == pd->ncross + 2);  /* Euler characteristic check */

  }
  
}

void pd_regenerate(pd_code_t *pd)

/* Starting with a valid list of crossings, regenerate 
   everything else. */

{
  assert(pd != NULL);

  pd_regenerate_comps(pd);
  pd_regenerate_faces(pd);
  pd_regenerate_hash(pd);

}


/* pd sanity checking */

bool pd_cross_ok(pd_code_t *pd)
/* Checks to see that the edge numbers referenced in the
   crossings are consecutive and between 0 and
   pd->nedges-1, and that the crossing data is sorted. */
{
  assert(pd != NULL);
  int edge_seen[PD_MAXEDGES];

  pd_idx_t edge, cross;
  pd_pos_t pos;

  for(edge=0;edge<PD_MAXEDGES;edge++) { edge_seen[edge] = 0; }

  for(cross=0;cross<pd->ncross;cross++) {

    for(pos=0;pos<4;pos++) {

      edge = pd->cross[cross].edge[pos];

      if (!(edge < pd->nedges)) {

	return pd_error(SRCLOC,"%CROSS contains illegal edge number %d in pd %PD",pd,
			cross,edge);

      } 

      edge_seen[edge]++;

    }

  }

  /* Now we check that each edge was seen twice. */

  for(edge=0;edge<pd->nedges;edge++) {

    if (edge_seen[edge] != 2) {

      return pd_error(SRCLOC,"%EDGE not seen twice in %PD",pd,edge);

    }

  }

  /* Now check sort order. */

  for(cross=1;cross<pd->ncross;cross++) {

    if (pd_cross_cmp(&(pd->cross[cross]),&(pd->cross[cross-1])) < 0) {

      return pd_error(SRCLOC,"%CROSS and %CROSS are out of order in %PD",
		      pd,cross-1,cross);

    }

  }

  return true;

}

bool pd_edges_ok(pd_code_t *pd) 

{
  assert(pd != NULL);

  pd_idx_t edge;
  
  for(edge=0;edge < pd->nedges;edge++) {

    pd_edge_t *e = &(pd->edge[edge]);

    if (e->head >= pd->ncross || 
	e->tail >= pd->ncross) {

      return pd_error(SRCLOC,"%EDGE contains reference to illegal crossing in pd %PD",pd,edge);

    } 
    
    if (e->headpos >= 4 ||
	e->tailpos >= 4) {

      return pd_error(SRCLOC,"%EDGE contains reference to illegal position in pd %PD",pd,edge);

    }

    pd_crossing_t *hc = &(pd->cross[e->head]), *tc = &(pd->cross[e->tail]);

    if (hc->edge[e->headpos] != edge) {

      return pd_error(SRCLOC,"%EDGE fails check at head %CROSS in pd %PD",pd,edge,e->head);

    }

    if (tc->edge[e->tailpos] != edge) {

      return pd_error(SRCLOC,"%EDGE fails check at tail %CROSS in pd %PD",pd,edge,e->tail);

    }

  }

  return true;

}

bool pd_comps_ok(pd_code_t *pd) 

/* Check that the edges are numbered correctly around the components 
   and that the edges are oriented head-to-tail, also checks that every
   edge is present in a component. */

{
  pd_idx_t comp, edge, next_edge;
  pd_idx_t correct_edgenum;

  assert(pd != NULL);
  
  if (pd->ncomps > PD_MAXCOMPONENTS) {

    return pd_error(SRCLOC,"Number of components %d > PD_MAXCOMPONENTS %d in pd %PD",pd,pd->ncomps,PD_MAXCOMPONENTS);

  }

  for(correct_edgenum=0,comp=0;comp < pd->ncomps;comp++) {

    for(edge=0;edge < pd->comp[comp].nedges;edge++,correct_edgenum++) {

      if (correct_edgenum != pd->comp[comp].edge[edge]) {

	return pd_error(SRCLOC,"Edge %d in %COMP has incorrect number (should be %d) in pd %PD",
			pd,edge,comp,correct_edgenum);

      } 

      next_edge = (edge + 1) % pd->comp[comp].nedges;

      pd_idx_t e_idx, n_idx;

      e_idx = pd->comp[comp].edge[edge];
      n_idx = pd->comp[comp].edge[next_edge];

      if (pd->edge[e_idx].head != pd->edge[n_idx].tail ||
	  (4 + pd->edge[e_idx].headpos - pd->edge[n_idx].tailpos) % 4 != 2) {

	return pd_error(SRCLOC,
			"Edges %EDGE and %EDGE are consecutive on %COMP at\n"
			"positions %d and %d, but don't match head-to-tail\n"
			"in pd %PD",pd,e_idx,n_idx,comp,edge,next_edge);

      }

    }

  }

  if (correct_edgenum != pd->nedges) {

    return pd_error(SRCLOC,"Count of edges %d doesn't match pd->nedges %d in pd %PD",pd,
		    correct_edgenum,pd->nedges);

  }

  return true;

}
 

bool pd_faces_ok(pd_code_t *pd) 

/* Checks face data */

{
  assert(pd != NULL);
  
  pd_idx_t face, edge, nxt_edge;
  pd_edge_t this,next;

  for(face=0;face<pd->nfaces;face++) {

    pd_face_t *f = &(pd->face[face]);

    for(edge=0,nxt_edge= 1 % f->nedges;
	edge<f->nedges;
	edge++,nxt_edge = (edge+1) % f->nedges) { 

      /* Check that f->edge entries are legal. */

      if (f->edge[edge] >= pd->nedges) {

	return pd_error(SRCLOC,"%FACE contains illegal edge %FEDGE in pd %PD",pd,
			face,face,edge);

      } 

      if (f->edge[nxt_edge] >= pd->nedges) {

	return pd_error(SRCLOC,"%FACE contains illegal edge %FEDGE in pd %PD",pd,
			face,face,nxt_edge);

      } 
			
      /* Check that consecutive edges meet at vertex in left turn. */

      this = pd_oriented_edge(pd->edge[f->edge[edge]],f->or[edge]);
      next = pd_oriented_edge(pd->edge[f->edge[nxt_edge]],f->or[nxt_edge]);

      if (this.head != next.tail || (next.tailpos + 1) % 4 != this.headpos) {

	return pd_error(SRCLOC,
			"%FEDGE -> %FEDGE transition not ok\n"
			"Oriented %EDGE and %EDGE don't meet\n"
			"correctly in pd %PD\n",
			pd,face,edge,face,nxt_edge,
			pd->face[face].edge[edge],
			pd->face[face].edge[nxt_edge]);
      }

    }

  }

  /* Check sort order according to pd_face_cmp */

  for(face=1;face<pd->nfaces;face++) {

    if (pd_face_cmp(&(pd->face[face-1]),&(pd->face[face])) > 0) {

      return pd_error(SRCLOC,"%FACE and %FACE out of order in pd %PD",
		      pd,face-1,face);

    }

  }

  return true;

}

bool pd_ok(pd_code_t *pd) {

  return pd_cross_ok(pd) && pd_edges_ok(pd) && pd_faces_ok(pd) && pd_comps_ok(pd);

}

/* pd operations */

pd_code_t *pd_copy(pd_code_t *pd)
/* Make a new-memory copy of pd. */
{
  pd_code_t *pdA;
  pdA = malloc(sizeof(pd_code_t));
  assert(pdA != NULL);
  (*pdA) = *pd;
  return pdA;
}

void pd_write(FILE *of,pd_code_t *pd) 

/* Writes a pd code (including all precomputed information) in human readable ASCII format. */
/* It's true that this seems inefficient, but it's better to compress the resulting */
/* file for storage than to muck around with a binary file format. */

/* Format:

   pd   <hash> <uid>
   nv   <nverts>
   <nv lines of crossing information in the format edge edge edge edge>
   ne   <nedges>
   <ne lines of crossing information in the format tail, tailpos -> head, headpos>
   nc   <ncomps>
   <nc lines, each containing nedges edge edge .... edge> 
   nf   <nfaces>
   <nf lines, each in the format nedges edge edge ... edge giving face information counterclockwise>

*/
   
{
  /* First, we do a bit of sanity checking */

  if (pd == NULL) { return; }
  if (pd->ncross > PD_MAXVERTS || pd->nedges > PD_MAXEDGES || pd->ncomps > PD_MAXCOMPONENTS || pd->nfaces > PD_MAXFACES) {

    fprintf(of,"INVALID PD CODE\n");
    printf("%s (%d): Invalid PD code passed for writing.\n",SRCLOC);
    return;

  }
  
  pd_idx_t cross,edge,face,comp;
  pd_pos_t  pos;

  fprintf(of,"pd %s %ju\n",pd->hash,(uintmax_t)(pd->uid));
  
  /* Crossing data */
  
  fprintf(of,"nv %u\n",(unsigned int)(pd->ncross));
  
  for(cross=0;cross<pd->ncross;cross++) {

    for(pos=0;pos<4;pos++) {

      fprintf(of,"%u ",(unsigned int)(pd->cross[cross].edge[pos]));

    }

    fprintf(of,"\n");

  }

  /* Edge data */

  fprintf(of,"ne %u\n",(unsigned int)(pd->nedges));
  
  for(edge=0;edge<pd->nedges;edge++) {

    fprintf(of,"%u,%u -> %u,%u \n",
	    (unsigned int)(pd->edge[edge].tail),(unsigned int)(pd->edge[edge].tailpos),
	    (unsigned int)(pd->edge[edge].head),(unsigned int)(pd->edge[edge].headpos));

  }

  /* Component data */

  fprintf(of,"nc %u\n",(unsigned int)(pd->ncomps));

  for(comp=0;comp<pd->ncomps;comp++) {

    fprintf(of,"%u : ",(unsigned int)(pd->comp[comp].nedges));

    for(edge=0;edge<pd->comp[comp].nedges;edge++) {

      fprintf(of,"%u ",(unsigned int)(pd->comp[comp].edge[edge]));

    }

    fprintf(of,"\n");
    
  }

  /* Face data */

  fprintf(of,"nf %u\n",(unsigned int)(pd->nfaces));

  for(face=0;face<pd->nfaces;face++) {

    fprintf(of,"%u : ",(unsigned int)(pd->face[face].nedges));

    for(edge=0;edge<pd->face[face].nedges;edge++) {

      if (pd->face[face].or[edge] == PD_POS_ORIENTATION) {
	
	fprintf(of,"+ ");

      } else {

	fprintf(of,"- ");

      }

      fprintf(of,"%u ",(unsigned int)(pd->face[face].edge[edge]));

    }

    fprintf(of,"\n");

  }

}  

bool pd_read(FILE *infile,pd_code_t *pd) 

/* Reads an (ASCII) pd code written by pd_write. Return true if succeed, false if fail. */

/* This is carefully written to check that the pd codes are compatible with this set of */
/* compile-time constants (PD_MAXVERTS, PD_MAXEDGES, PD_MAXFACES, etc) since we set these as small */
/* as possible to minimize run-time memory use, so it's completely plausible that we're  */
/* trying to load a pd code which is valid, but may overrun our buffers. */

{
  /* pd   <hash> <uid> */ /* Remember that the hash is base64 encoded using a-zA-Z;: character set */

  uintmax_t input_temp,input_temp2,input_temp3,input_temp4;

  if (fscanf(infile," pd %32[a-zA-Z0-9;:]s ",pd->hash) != 1) { return false; }

  if (fscanf(infile," %ju ", &input_temp) != 1) { return false; }
  pd->uid = (pd_uid_t)(input_temp);

  /* nv   <nverts> */

  int cross,edge,comp,pos,face;
  
  if (fscanf(infile,"nv %ju ",&input_temp) != 1) { return false; }
  pd->ncross = (pd_idx_t)(input_temp);

  if (pd->ncross > PD_MAXVERTS) {

    fprintf(stderr,
	    "%s (%d): Reading pd code which appears to be valid but has %d crossings. "
	    "         This version of pdcode compiled for PD_MAXVERTS = %d.\n",__FILE__,__LINE__,
	    pd->ncross,PD_MAXVERTS);
    exit(1);

  }

  /* <nv lines of crossing information in the format edge edge edge edge> */

  for(cross=0;cross<pd->ncross;cross++) {

    for(pos=0;pos<4;pos++) {

      if(fscanf(infile," %ju ",&input_temp) != 1) { return false; }
      pd->cross[cross].edge[pos] = (pd_idx_t)(input_temp);

    }

  }

  /* ne   <nedges> */

  if (fscanf(infile,"ne %ju ",&input_temp) != 1) { return false; }
  pd->nedges = (pd_idx_t)(input_temp);

  if (pd->nedges > PD_MAXEDGES) {

    fprintf(stderr,
	    "%s (%d): Reading pd code which appears to be valid but has %d edges. "
	    "         This version of pdcode compiled for PD_MAXEDGES = %d.\n",
	    __FILE__,__LINE__,pd->nedges,PD_MAXEDGES);
    exit(1);

  }

  /* <ne lines of crossing information in the format tail, tailpos -> head, headpos> */

  for(edge=0;edge<pd->nedges;edge++) {

    if(fscanf(infile," %ju,%ju -> %ju,%ju ",
	      &input_temp,&input_temp2,&input_temp3,&input_temp4) != 4) { return false; }

    pd->edge[edge].tail = (pd_idx_t)(input_temp);
    pd->edge[edge].tailpos = (pd_pos_t)(input_temp2);

    pd->edge[edge].head = (pd_idx_t)(input_temp3);
    pd->edge[edge].headpos = (pd_pos_t)(input_temp4);

  }

  /* nc   <ncomps> */

  if (fscanf(infile,"nc %ju ",&input_temp) != 1) { return false; }
  pd->ncomps = (pd_idx_t)(input_temp);

  if (pd->ncomps > PD_MAXCOMPONENTS) {

    fprintf(stderr,
	    "%s (%d): Reading pd code which appears to be valid but has %d components. "
	    "         This version of pdcode compiled for PD_MAXCOMPONENTS = %d.\n",
	    __FILE__,__LINE__,pd->ncomps,PD_MAXCOMPONENTS);
    exit(1);

  }

  /* <nc lines, each containing nedges : +/- edge +/- edge .... +/- edge> */

  for(comp=0;comp<pd->ncomps;comp++) {
     
    if (fscanf(infile," %ju : ",&input_temp) != 1) { return false; }
    pd->comp[comp].nedges = (pd_idx_t)(input_temp);

    if (pd->comp[comp].nedges > PD_MAXEDGES) {

      fprintf(stderr,
	    "%s (%d): Reading component which appears to be valid but has %d edges. "
	    "         This version of pdcode compiled for PD_MAXEDGES = %d.\n",
	      __FILE__,__LINE__,pd->nedges,PD_MAXEDGES);
      exit(1);

    }

    for(edge=0;edge<pd->comp[comp].nedges;edge++) {

      if(fscanf(infile," %ju ",&input_temp) != 1) { return false; }
      pd->comp[comp].edge[edge] = (pd_idx_t)(input_temp);
      

    }

  }

  /*nf   <nfaces> */
 
  if (fscanf(infile,"nf %ju ",&input_temp) != 1) { return false; }
  pd->nfaces = (pd_idx_t)(input_temp);

  if (pd->nfaces > PD_MAXFACES) {

    fprintf(stderr,
	    "%s (%d): Reading pd code which appears to be valid but has %d faces. "
	    "         This version of pdcode compiled for PD_MAXFACES = %d.\n",
	    __FILE__,__LINE__,pd->nfaces,PD_MAXFACES);
    exit(1);

  }

  /* <nf lines, each in the format nedges edge edge
     ... edge giving face info counterclockwise> */

  for(face=0;face<pd->nfaces;face++) {
    
    if (fscanf(infile," %ju : ",&input_temp) != 1) { return false; }
    pd->face[face].nedges = (pd_idx_t)(input_temp);

    if (pd->face[face].nedges > PD_MAXEDGES) {

      fprintf(stderr,
	    "%s (%d): Reading face which appears to be valid but has %d edges. "
	    "         This version of pdcode compiled for PD_MAXEDGES = %d.\n",
	      __FILE__,__LINE__,pd->nedges,PD_MAXEDGES);
      exit(1);

    }

    for(edge=0;edge<pd->face[face].nedges;edge++) {

      char orientation[2];
      if(fscanf(infile," %1[+-]s ",orientation) != 1) {

	return false;

      }

      if (fscanf(infile," %ju ",&input_temp) != 1) { return false; }

      pd->face[face].edge[edge] = (pd_idx_t)(input_temp);
      pd->face[face].or[edge] = (orientation[0] == '+') ? PD_POS_ORIENTATION : PD_NEG_ORIENTATION;

    }

  }

  return true;
  
}


/* Pd human output */

char pd_print_or(pd_or_t or) 
/* Returns a single-character version of or: +, -, U (for unset), or ? (anything else) */
{
  if (or == PD_POS_ORIENTATION) { return '+'; }
  else if (or == PD_NEG_ORIENTATION) { return '-'; }
  else if (or == PD_UNSET_ORIENTATION) { return 'U'; }
  else return '?';
}

void pd_vfprintf(FILE *stream, char *infmt, pd_code_t *pd, va_list ap )

/* This (internal) function is designed to make it easy to
   report information referencing parts of a pd_code. It
   is a specialized printf which takes several new format
   conversions which convert pd_idx_t types in a special
   way

   Tag        pd_idx_t           output

   %FACE      face number        fnum (e1 (or1) -> e2 (or2) -> ... -> e1 (or1))        
   %EDGE      edge number        enum (tail (tailpos) -> head (headpos) )
   %CROSS     cross number       cnum (e0 e1 e2 e3)
   %COMP      comp number        compnum (e1 -> e2 -> e3 -> ..... -> e1 (or1))
   %PD        (no argument)      (\n\n output of pd_write \n\n)

   We also have conversions for pointers. 

   %OR        *pd_or_t           +, -, U (unset), or ? (anything else)
   %MULTIDX   *pd_multidx_t      multidx (i[0] i[1] ... i[n-1])
   %COMPGRP   *pd_compgrp_t      compgrp (comp[0] .. comp[n-1])
 
   %EDGEMAP   *pd_edgemap_t      edgemap (+/- map[0] +/- map[1] ... +/- map[n-1])
   %CROSSMAP  *pd_crossmap_t     crossmap +/- ( permutation )
   %FACEMAP   *pd_facemap_t      facemap +/- ( permutation )
   %CROSSPTR  *pd_crossing_t     cross (e0 e1 e2 e3)

   It also converts %d specifications in the usual way.

   We ignore any other format conversions present in fmt,
   passing them unchanged into the output stream. The intention
   is that these would have already been processed in fmt (using 
   sprintf) before calling pd_vfprintf.

*/

{
  char *fmtptr, *nxtconv;
  char *fmt = strdup(infmt);

  assert(stream != NULL);

  for(fmtptr=fmt,nxtconv = strchr(fmtptr,'%');
      fmtptr != NULL;
      fmtptr = nxtconv, nxtconv = strchr(fmtptr,'%')) {

    /* Start by copying the characters between fmtptr and 
       nxtconv into outbuf. */
    
    if (nxtconv == NULL) { /* No more conversions remain */

      fprintf(stream,"%s",fmtptr);
      break;

    } 

    *nxtconv = 0;
    fprintf(stream,"%s",fmtptr);
    *nxtconv = '%'; /* Restore the previous % character */

    /* Now expand the current format conversion. */

    if (!strncmp(nxtconv,"%FACE",5)) { /* %FACE conversion */

      pd_idx_t face = (pd_idx_t) va_arg(ap,int);
      pd_idx_t edge;

      fprintf(stream,"face %d (",face);
      for(edge=0;edge<pd->face[face].nedges-1 && edge<PD_MAXEDGES;edge++) {

	fprintf(stream," (%c) %d ->",
		pd->face[face].or[edge] == PD_POS_ORIENTATION ? '+':'-',
		pd->face[face].edge[edge]);

      }

      fprintf(stream,"(%c) %d) ",
	      pd->face[face].or[edge] == PD_POS_ORIENTATION ? '+':'-',
	      pd->face[face].edge[edge]);

      nxtconv += 5;

    } else if (!strncmp(nxtconv,"%EDGE ",6)) { /* %EDGE conversion */

      pd_idx_t edge = (pd_idx_t) va_arg(ap,int);

      fprintf(stream,"edge %d (%d,%d -> %d,%d)",edge,
	      pd->edge[edge].tail,pd->edge[edge].tailpos,
	      pd->edge[edge].head,pd->edge[edge].headpos);

      nxtconv += 5;

    } else if (!strncmp(nxtconv,"%COMP ",6)) { /* %COMP conversion */

      pd_idx_t comp = (pd_idx_t) va_arg(ap,int);
      pd_idx_t edge;
      
      fprintf(stream,"comp %d (",comp);
      for(edge=0;edge<pd->comp[comp].nedges-1 && edge<PD_MAXEDGES;edge++) {

	fprintf(stream," %d ->",
		pd->comp[comp].edge[edge]);

      }

      fprintf(stream," %d ) ",
	      pd->comp[comp].edge[edge]);

      nxtconv += 5;

    } else if (!strncmp(nxtconv,"%CROSS ",7)) { /* %CROSS conversion */

      pd_idx_t cross = (pd_idx_t) va_arg(ap,int);

      fprintf(stream,"cross %d (%d %d %d %d)",cross,
	      pd->cross[cross].edge[0],	      
	      pd->cross[cross].edge[1],
	      pd->cross[cross].edge[2],
	      pd->cross[cross].edge[3]);

      nxtconv += 6;

    } else if (!strncmp(nxtconv,"%CROSSPTR ",9)) { /* %CROSSPTR conversion */

      pd_crossing_t *cross = (pd_crossing_t *) va_arg(ap,void *);

      fprintf(stream,"cross (%d %d %d %d)",
	      cross->edge[0],	      
	      cross->edge[1],
	      cross->edge[2],
	      cross->edge[3]);

      nxtconv += 8;

    } else if (!strncmp(nxtconv,"%PD",3)) { /* %PD conversion */

      fprintf(stream,"\n");
      pd_write(stream,pd);
      fprintf(stream,"\n");

      nxtconv += 3;

    } else if (!strncmp(nxtconv,"%FEDGE",6)) { /* FACE/EDGE conversion */

      pd_idx_t face = (pd_idx_t) va_arg(ap,int);
      pd_idx_t edge = (pd_idx_t) va_arg(ap,int);
       
      fprintf(stream,"%d (%c)",pd->face[face].edge[edge],
	      pd->face[face].or[edge] == PD_POS_ORIENTATION ? '+':'-');

      nxtconv += 6;

    } else if (!strncmp(nxtconv,"%MULTIDX",8)) { /* multidx conversion */

      pd_idx_t i;
      pd_multidx_t *idx = (pd_multidx_t *) va_arg(ap,void *);
      
      fprintf(stream,"multidx ( ");
      for(i=0;i<idx->nobj;i++) { 

	char *objp;
	objp = idx->ops.print(idx->obj[i]);
	fprintf(stream,"%d/%d (%s) ",idx->i[i],idx->nvals[i],objp);
	free(objp);

      }
      
      fprintf(stream,")");

      nxtconv += 8;

    } else if (!strncmp(nxtconv,"%COMPGRP",8)) { /* Component group */

      pd_idx_t comp;
      pd_compgrp_t *grp = (pd_compgrp_t *) va_arg(ap,void *);

      fprintf(stream,"compgrp ( ");
      for(comp=0;comp<grp->ncomps && comp<PD_MAXCOMPONENTS;comp++) { fprintf(stream,"%d ",grp->comp[comp]); }
      fprintf(stream,")");

      nxtconv += 8;

    } else if (!strncmp(nxtconv,"%EDGEMAP",8)) { /* Edge map */

      pd_edgemap_t *edgemap = (pd_edgemap_t *) va_arg(ap,void *);
      char *printed_form;

      printed_form = pd_print_edgemap(edgemap);
      fprintf(stream,"%s",printed_form);
      free(printed_form);

      nxtconv += 8;

    } else if (!strncmp(nxtconv,"%OR ",3)) { /* Edge map */

      pd_or_t *or = (pd_or_t *) va_arg(ap,void *);
      fprintf(stream,"or (%c)",pd_print_or(*or));

      nxtconv += 3;

    }  else if (!strncmp(nxtconv,"%CROSSMAP ",10)) { /* Crossing map */

      pd_crossmap_t *crmap = (pd_crossmap_t *) va_arg(ap,void *);
      char *printed;

      printed = pd_print_crossmap(crmap);
      fprintf(stream,"%s",printed);
      free(printed);

      nxtconv += 9;

    } else if (!strncmp(nxtconv,"%ISO ",5)) { /* pd_code -> pd_code isomorphism */

      pd_iso_t *iso = (pd_iso_t *) va_arg(ap,void *);
      char *printed;

      printed = pd_print_iso(iso);
      fprintf(stream,"%s",printed);
      free(printed);

      nxtconv += 4;

    } else if (!strncmp(nxtconv,"%FACEMAP ",9)) { /* Face map */

      pd_facemap_t *facemap = (pd_facemap_t *) va_arg(ap,void *);
      char *printed;

      printed = pd_print_facemap(facemap);
      fprintf(stream,"%s",printed);
      free(printed);

      nxtconv += 8;

    } else if (!strncmp(nxtconv,"%DIHEDRAL ",10)) { /* Dihedral group element */

      pd_dihedral_t *d = (pd_dihedral_t *) va_arg(ap,void *);

      fprintf(stream,"dihedral ");

      char *dprint = pd_print_dihedral(d);
      fprintf(stream,"%s",dprint);
      free(dprint);

      nxtconv += 9;

    } else if (!strncmp(nxtconv,"%PERM ",6)) { /* Permutation group elt */

      pd_idx_t i;
      pd_perm_t *perm = (pd_perm_t *) va_arg(ap,void *);

      fprintf(stream,"perm ( ");

      for(i=0;i<perm->n;i++) { fprintf(stream,"%d ",perm->map[i]); }
      fprintf(stream,") idx %d",perm->pc_idx);
	
      char *pform;
      pform = pd_print_perm(perm);

      if (pform != NULL) { 

	fprintf(stream," %s",pform);
	free(pform);

      }

      nxtconv += 5;

    } else if (!strncmp(nxtconv,"%d",2)) { /* Standard decimal conversion */

      int decimal = (int) va_arg(ap,int);
      fprintf(stream,"%d",decimal);
      nxtconv += 2;

    } else { /* A standard or unrecognized format conversion */

      char *adv;
      for(adv=nxtconv;!isspace(*adv) && *adv != 0;adv++); 
      /* Advance the pointer to the next whitespace or end of string. */
      *adv = 0; /* Terminate the string there */

      fprintf(stream,"%s",nxtconv); /* Copy the format conversion into stream */
      *adv = ' ';
      nxtconv = adv;

    }

  }

  free(fmt);
  
}
      
bool pd_error(char *file, int line, char *fmt, pd_code_t *pd, ...)

{
  va_list ap;
  
  if (PD_VERBOSE < 10) { return false; }
  
  fprintf(stderr,"%s (%d): pd_error \n",file,line);
  
  va_start(ap,pd);
  pd_vfprintf(stderr,fmt,pd,ap);
  va_end(ap);
  
  exit(1);
}

void pd_printf(char *fmt,pd_code_t *pd, ... )

{
  va_list ap;

  va_start(ap,pd);
  pd_vfprintf(stdout,fmt,pd,ap);
  va_end(ap);
}

bool pd_isomorphic_strings(char *pdcodeA, int nA, char*pdcodeB, int nB)

{
  // First check that two shadows have the same number of crossings.
  if(nA != nB) {
    return false;
  }
  
  // Allocate the memory for the two pd code structs.
  pd_code_t *pdA = calloc(1,sizeof(pd_code_t));
  assert(pdA != NULL);

  pdA->ncross = nA;

  pd_code_t *pdB =calloc(1,sizeof(pd_code_t));
  assert(pdA !=NULL);

  pdB->ncross = nB;

  // Loop through the passed in character arrays and build the
  //  crossings. Note that the passed in arrays should contain
  //  characters whose ASCII values are the appropriate integers
  //  plus 65 (so that we start with A).
  int iA;
  int iB;
  int j = 0;
    
  for(iA=0;iA<nA;iA++) {
    pdA->cross[iA] = pd_build_cross((int)pdcodeA[j]-65,(int)pdcodeA[j+1]-65,
				    (int)pdcodeA[j+2]-65,(int)pdcodeA[j+3]-65);
    //printf("%d,%d,%d,%d \n",(int)pdcodeA[j],(int)pdcodeA[j+1],(int)pdcodeA[j+2],(int)pdcodeA[j+3]);
    j=j+4;    
  }

  j=0;

  for(iB=0;iB<nB;iB++) {
    pdB->cross[iB] = pd_build_cross((int)pdcodeB[j]-65,(int)pdcodeB[j+1]-65,
				    (int)pdcodeB[j+2]-65,(int)pdcodeB[j+3]-65);
    //printf("%d,%d,%d,%d \n",(int)pdcodeB[j],(int)pdcodeB[j+1],(int)pdcodeB[j+2],(int)pdcodeB[j+3]);
    j=j+4;
  }

  // Regenerate the pd codes
  pd_regenerate(pdA);
  assert(pd_ok(pdA));

  pd_regenerate(pdB);
  assert(pd_ok(pdB));

  // Return the call to pd_isomorphic
  return pd_isomorphic(pdA,pdB);
}
