/*

   Code for handling pd_codes (link shadows). Part of the
   general COLD (Census of Link Diagrams) project. This fork
   incorporated in plCurve is significantly different, because
   it changes the way that the code handles very large diagrams.

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

#ifdef HAVE_MATH_H
  #include<math.h>
#endif

#include<plcTopology.h>
//#include<libcassie/cassie.h>
//#include"/usr/local/include/thrift/Thrift.h"
//#include<python2.7/Python.h>
#include<pd_multidx.h>
#include<pd_dihedral.h>
#include<pd_cyclic.h>
#include<pd_perm.h>
#include<pd_orientation.h>

#include<pd_isomorphisms.h>

pd_code_t *pd_code_new(pd_idx_t mv) {

  assert(mv != 0);
  pd_idx_t maxverts = (mv == 1) ? 2 : mv;

  pd_code_t *pd = calloc(1,sizeof(pd_code_t));
  assert(pd != NULL);
  pd->uid = 0;

  pd->MAXVERTS = maxverts;
  pd->MAXEDGES = 2*maxverts+1;
  pd->MAXCOMPONENTS = (pd_idx_t)(floor(maxverts/2.0)) + 2;
  pd->MAXFACES = maxverts+2;

  pd->ncross = 0;
  pd->nedges = 0;
  pd->ncomps = 0;
  pd->nfaces = 0;

  sprintf(pd->hash," ");
  pd->edge = calloc(pd->MAXEDGES,sizeof(pd_edge_t));
  assert(pd->edge != NULL);

  pd->comp = calloc(pd->MAXCOMPONENTS,sizeof(pd_component_t));
  assert(pd->comp != NULL);

  pd->cross = calloc(pd->MAXVERTS,sizeof(pd_crossing_t));
  assert(pd->cross != NULL);

  pd->face = calloc(pd->MAXFACES,sizeof(pd_face_t));
  assert(pd->face != NULL);

  return pd;

}

void pd_code_eltfree(void **PD) {

  pd_code_t **pd = (pd_code_t **)(PD);
  pd_code_free(pd);

}

void pd_code_free(pd_code_t **PD) {

  int cmp;

  if (*PD == NULL) { return; }

  pd_code_t *pd = *PD;

  if (pd->MAXCOMPONENTS != 0 && pd->comp != NULL) {

    for(cmp=0;cmp<pd->MAXCOMPONENTS;cmp++) {

      if (pd->comp[cmp].edge != NULL) {

	free(pd->comp[cmp].edge);
	pd->comp[cmp].edge = NULL;
	pd->comp[cmp].nedges = 0;

      }

    }

    free(pd->comp);
    pd->comp = NULL;
    pd->MAXCOMPONENTS = 0;
    pd->ncomps = 0;

  }

  pd_idx_t face;

  if (pd->MAXFACES != 0 && pd->face != NULL) {

    for(face=0;face<pd->MAXFACES;face++) {

      if (pd->face[face].edge != NULL) {

	free(pd->face[face].edge);
	pd->face[face].edge = NULL;

      }

      if (pd->face[face].or != NULL) {

	free(pd->face[face].or);
	pd->face[face].or = NULL;

      }

      pd->face[face].nedges = 0;

    }

    free(pd->face);
    pd->face = NULL;
    pd->MAXFACES = 0;
    pd->nfaces = 0;

  }

  if (pd->edge != NULL) {

    free(pd->edge);
    pd->edge = NULL;
    pd->nedges = 0;
    pd->MAXEDGES = 0;

  }

  if (pd->cross != NULL) {

    free(pd->cross);
    pd->cross = NULL;
    pd->ncross = 0;
    pd->MAXVERTS = 0;

  }

  sprintf(pd->hash," ");
  pd->uid = 0;
  pd->MAXVERTS = 0;
  pd->MAXEDGES = 0;
  pd->ncross = 0;
  pd->nedges = 0;

  free(pd);
  *PD = NULL;

}

/* These are the basic allocate/deallocate functions. */

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

  cr.edge[0] = e0; cr.edge[1] = e1; cr.edge[2] = e2; cr.edge[3] = e3; cr.sign = 0; /* We don't know, so initialize to 0 */

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
  reface.edge = calloc(reface.nedges,sizeof(pd_idx_t));
  reface.or = calloc(reface.nedges,sizeof(pd_or_t));
  assert(reface.edge != NULL && reface.or != NULL);

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

  pd_idx_t lowE = reface.edge[0],lowPos = 0;

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

  /* Finally, free the scratch buffers reface.edge and reface.or. */

  free(reface.edge);
  free(reface.or);

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

void pd_face_and_pos(pd_code_t *pd, pd_idx_t edge,
		     pd_idx_t *posface, pd_idx_t *posface_pos,
		     pd_idx_t *negface, pd_idx_t *negface_pos)

  /* Finds the two faces which edge occurs on, which should include
     one face where edge appears in positive orientation (posface)
     and one where edge appears with negative orientation (negface).
     If we don't find two with this description, die.

     Return the position (index) of the edge on each face. */

{
  bool pos_found = false, neg_found = false;
  pd_idx_t face;
  pd_idx_t fedge;

  for(face=0;face < pd->nfaces;face++) {

    for(fedge=0;fedge < pd->face[face].nedges; fedge++)  {

      if (pd->face[face].edge[fedge] == edge) {

	if (pd->face[face].or[fedge] == PD_POS_ORIENTATION) {

	  if (pos_found == true) { /* We have a problem */

	    pd_error(SRCLOC,
		     "Found edge %EDGE twice with positive orientation \n"
		     "%FACE, position %d\n"
		     "%FACE, position %d\n",
		     pd,edge,*posface,*posface_pos,face,fedge);

	  }

	  pos_found = true;
	  *posface = face; *posface_pos = fedge;

	} else if (pd->face[face].or[fedge] == PD_NEG_ORIENTATION) {

	  if (neg_found == true) { /* We have a problem */

	    pd_error(SRCLOC,
		     "Found edge %EDGE twice with negative orientation \n"
		     "%FACE, position %d\n"
		     "%FACE, position %d\n",
		     pd,edge,*negface,*negface_pos,face,fedge);

	  }

	  neg_found = true;
	  *negface = face; *negface_pos = fedge;

	} else {

	  pd_error(SRCLOC,"pd %PD contains unknown orientation %d for edge %d of face %FACE",pd,
		   pd->face[face].or[fedge],fedge,face);

	}

      }

    }

  }

  if (!pos_found || !neg_found) {

    pd_error(SRCLOC,"couldn't find edge %EDGE on two faces of %PD with positive and negative orientation.\n",
	     pd,edge);

  }

}

bool pd_edge_on_face(pd_code_t *pd, pd_idx_t edge, pd_idx_t face)
  /* Returns true if the edge is on the face (with either sign). */
{

  pd_idx_t i;
  bool found = false;

  pd_check_edge(SRCLOC,pd,edge);
  pd_check_face(SRCLOC,pd,face);

  for(i=0;i<pd->face[face].nedges && !found;i++) {

    if (pd->face[face].edge[i] == edge) {

      found = true;

    }

  }

  return found;

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

  ncr.sign = cr->sign;
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
  pd_pos_t cross_rotate[pd->MAXVERTS];

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
  pd_idx_t cross_perm[pd->MAXVERTS];
  for(cross=0;cross<pd->ncross;cross++) { cross_perm[cross] = cross; }
  qsort(cross_perm,pd->ncross,sizeof(pd_idx_t),pdint_cross_index_cmp);

  pd_idx_t new_cross_num[pd->MAXVERTS];
  for(cross=0;cross<pd->ncross;cross++) { new_cross_num[cross_perm[cross]] = cross; }


  /* First, we change the references to crossings in the edge records. */

  for(edge=0;edge<pd->nedges;edge++) {

    pd->edge[edge].head = new_cross_num[pd->edge[edge].head];
    pd->edge[edge].tail = new_cross_num[pd->edge[edge].tail];

  }

  /* Next, we rearrange the crossings themselves. */

  pd_crossing_t newcross[pd->MAXVERTS];

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

#define PD_UNSET_COMPONENT -1
#define PD_UNSET_POSITION -1

typedef struct edge_assignment_struct {

  /* When building the components from a pd code, we need
     to keep track of the component a given edge has been
     assigned to, and the position of the edge within that
     component.

     This structure stores that information. If the component
     hasn't been assigned yet, the component is set to PD_UNSET_COMPONENT
     and the pos is set to PD_UNSET_POSITION. */

    int comp;
    int pos;

} pdint_edge_assignment;

pd_idx_t pdint_find_unassigned_edge(pd_code_t *pd,pdint_edge_assignment *edge_assigned, bool *all_assigned)

{
  assert(pd != NULL && edge_assigned != NULL && all_assigned != NULL);

  pd_idx_t edge;

  for(edge=0;edge<pd->nedges;edge++) {

    if (edge_assigned[edge].comp == PD_UNSET_COMPONENT) {

      *all_assigned = false;
      return edge;

    }

  }

  *all_assigned = true;
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

typedef struct crossing_marker_struct {

  bool     pos_used[4];
  pd_idx_t num_used;

} pd_crossmark_t;

pd_crossmark_t *pdint_crossmark_new(pd_idx_t ncross) {

  pd_crossmark_t *cm;
  cm = calloc(ncross,sizeof(pd_crossmark_t));
  assert(cm != NULL);
  pd_idx_t i,j;

  for(i=0;i<ncross;i++) {

    for(j=0;j<4;j++) {

      cm[i].pos_used[j] = false;

    }

    cm[i].num_used = 0;

  }

  return cm;

}

bool pdint_find_unused_cross(pd_crossmark_t *cm,pd_idx_t ncross,pd_idx_t *cross,pd_idx_t *pos)
{
  pd_idx_t i,j;

  for(i=0;i<ncross;i++) {

    if (cm[i].num_used != 4) {

      for(j=0;j<4;j++) {

	if (!cm[i].pos_used[j]) {

	  *cross = i; *pos = j;
	  return true;

	}

      }

    }

  }

  return false;

}

void pdint_find_edge_in_crossings(pd_crossing_t *cross,pd_idx_t ncross,
				  pd_idx_t edge,
				  pd_idx_t ends[2], pd_pos_t pos[2])

/* Search for this edge in crossing buffer. */

{
  pd_idx_t nfound = 0;
  pd_idx_t i,j;

  for(i=0;i<ncross;i++) {

    for(j=0;j<4;j++) {

      if (cross[i].edge[j] == edge) {

	ends[nfound] = i; pos[nfound] = j;
	nfound++;

	if (nfound == 2) { return; }

      }

    }

  }

  /* We shouldn't get here. */

  printf("pdint_find_edge_in_crossings:"
	 "Couldn't find edge %d twice"
	 "in list of crossings\n",edge);
  for(i=0;i<ncross;i++) {

    printf("\t %d %d %d %d \n",
	   cross[i].edge[0],
	   cross[i].edge[1],
	   cross[i].edge[2],
	   cross[i].edge[3]);
  }

  exit(1);

}

void pd_regenerate_edges(pd_code_t *pd)

/*
  We do this the slow-but-clean way, generating edges in
  order around components by direct search. The idea is
  to mark used positions on each crossing, then enter
  edge records as we go. It's slow the first time, but
  cached for future accesses.
*/

{

  /* First, mark all the crossings and positions as unused */

  pd_crossmark_t *cm = pdint_crossmark_new(pd->ncross);

  /* Second, flush the current edgebuffer. */

  assert(pd->edge != NULL);
  pd_idx_t i,j;

  for(i=0;i<pd->MAXEDGES;i++) {

    pd->edge[i].head    = PD_UNSET_IDX;
    pd->edge[i].headpos = PD_UNSET_POS;
    pd->edge[i].tail    = PD_UNSET_IDX;
    pd->edge[i].tailpos = PD_UNSET_POS;

  }

  pd->nedges = 2*pd->ncross;

  /* Now we recreate the old code, which simply loops through and generates
     edge records with no thought of consistent orientation. */

  for(i=0;i<pd->nedges;i++) {

    pd_idx_t foundcr[2], foundpos[2];
    pdint_find_edge_in_crossings(pd->cross,pd->ncross,i,foundcr,foundpos);

    pd->edge[i].head = foundcr[0];  pd->edge[i].headpos = foundpos[0];
    pd->edge[i].tail = foundcr[1];  pd->edge[i].tailpos = foundpos[1];

  }

  /* Ok, we don't have any sense of consistent orientation here. So
     this is where we'll try to generate one. */

  pd_idx_t start_edge_idx;
  pd_edge_t start_edge;
  pd_idx_t failsafe = 0;

  while (pdint_find_unused_cross(cm,pd->ncross,&i,&j)) {

    start_edge_idx = pd->cross[i].edge[j];
    start_edge = pd->edge[start_edge_idx];
    /* This edge is by assumption oriented positively. */

    cm[start_edge.head].pos_used[start_edge.headpos] = true;
    cm[start_edge.tail].pos_used[start_edge.tailpos] = true;
    cm[start_edge.head].num_used++;
    cm[start_edge.tail].num_used++;

    pd_idx_t next_edge_idx;
    pd_or_t  next_or;
    pd_idx_t compsafe = 0;

    for (pdint_next_comp_edge(pd,start_edge_idx,
			      &next_edge_idx,&next_or);
	 next_edge_idx != start_edge_idx;
	 pdint_next_comp_edge(pd,next_edge_idx,
			      &next_edge_idx,&next_or),
	   compsafe++) {

      pd->edge[next_edge_idx]
	= pd_oriented_edge(pd->edge[next_edge_idx],
			   next_or);

      pd_edge_t next_edge = pd->edge[next_edge_idx];

      cm[next_edge.head].pos_used[next_edge.headpos] = true;
      cm[next_edge.tail].pos_used[next_edge.tailpos] = true;
      cm[next_edge.head].num_used++;
      cm[next_edge.tail].num_used++;

      if (compsafe >= pd->nedges) {

	pd_error(SRCLOC,"component appears to have %d edges in %d edge pd %PD",pd,
		 compsafe+1,pd->nedges);
	exit(1);

      }

    }

    failsafe++;
    if (failsafe > pd->ncross) {

      pd_error(SRCLOC,
	       "Ran through list of edges %d times looking for components in %PD",pd,failsafe);
      exit(1);

    }

  }

  /* At this point, we should pass pd_edges_ok. */

  assert(pd_edges_ok(pd));

  /* And we can safely do some housecleaning. */

  free(cm);

}


void pd_regenerate_comps(pd_code_t *pd)

/*
   This procedure assumes that edges exist, and possibly that components still
   exist (because we've used this pdcode for something before this point). We
   delete any existing component records and reassemble the edges into components,
   flipping orientations as required. Note that pd_edges_ok is simply a check that
   the vertex indices are valid-- it doesn't attempt to check orientations (although
   the pd_regenerate_edges code should actually produce consistently oriented edges).
*/

{
  assert(pd != NULL);
  assert(pd->ncross <= pd->MAXVERTS);
  assert(pd_cross_ok(pd));
  assert(pd_edges_ok(pd));

  if (pd->ncomps != 0 && pd->ncomps != PD_UNSET_IDX) { /* We have components allocated. In this case, destroy them: */

    pd_idx_t comp;

    for(comp=0;comp<pd->ncomps;comp++) {

      if (pd->comp[comp].edge != NULL) {

	free(pd->comp[comp].edge); pd->comp[comp].edge = NULL;

      }

      pd->comp[comp].nedges = PD_UNSET_IDX;
      pd->comp[comp].tag = PD_UNSET_TAG;

    }

    pd->ncomps = 0; /* There are 0 components NOW! */

  }

  /* Step 0.

     We're going to generate a new collection of components in any
     case.  If we don't already have components, then we'll assign
     tags now.

     We DON'T need to free the pd->comp array itself, which remains a
     fixed size, and was allocated by pd_code_new, regardless of the
     number of components that are actually used.

     But we don't need to reallocate it, either.

  */

  assert(pd->comp != NULL);

  /* Step 1. Run around the components,
     assembling the (old) edge
     numbers of the edges, and
     reorienting edges as we go. */

  pdint_edge_assignment *edge_assigned;
  edge_assigned = calloc(pd->nedges,sizeof(pdint_edge_assignment));
  assert(edge_assigned != NULL);

  bool all_assigned = false;
  pd_idx_t edge;

  for(edge=0;edge<pd->nedges;edge++) {

    edge_assigned[edge].comp = PD_UNSET_COMPONENT;
    edge_assigned[edge].pos = PD_UNSET_POSITION;

  }

  /* We're now going to work our way through the edges, assigning
     them to components and positions. */

  for(edge=0,pd->ncomps=0;            /* This loop is over components */
      !all_assigned;
      pd->ncomps++,edge=pdint_find_unassigned_edge(pd,edge_assigned,&all_assigned)) {

    pd_idx_t comp = pd->ncomps;
    pd_or_t  next_or;
    pd_idx_t next_edge, start_edge;
    pd_idx_t nedges = 1;          /* Keep a counter here to prevent infinite loops */

    /* We are now starting a new component with an unassigned edge. The first step
       is to assign it to this component in position 0. This is the first edge in
       the component, so we keep track of it in order to know when we've completed
       the loop. */

    edge_assigned[edge].comp = comp;
    edge_assigned[edge].pos = 0;
    start_edge = edge;

    pdint_next_comp_edge(pd,edge,&next_edge,&next_or);

    for(;next_edge != start_edge && nedges <= pd->MAXEDGES+1;    /* This loop is over edges in a component */
	pdint_next_comp_edge(pd,edge,&next_edge,&next_or),nedges++) {

      /* We've now discovered that next_edge is in this component.
	 Go ahead and assign it to this component, in the position given
	 by nedges, and change the orientation (if needed). */

      pd_reorient_edge(pd,next_edge,next_or);
      edge_assigned[next_edge].comp = comp;
      edge_assigned[next_edge].pos = nedges;

      /* Now set ourselves up to go to the next step by resetting edge. */

      edge = next_edge;

      /* Finally, we include a safety check to prevent infinite loops. */

      if (nedges >= pd->nedges) {

	fprintf(stderr,
		"pd_regenerate_comp: Found %d edges in component %d of %d edge pdcode.\n"
		"                    Suspect that crossing information is poorly formed.\n",
		nedges+1,comp,pd->nedges);
	exit(1);

      }

    } /* End of loop over this component. */

    if (pd->ncomps > pd->MAXCOMPONENTS) {

      fprintf(stderr,
	      "pd_regenerate_comps: Found %d components in %d edge pdcode without\n"
	      "                     detecting \"all_assigned\". Suspect bug in \n"
	      "                     assignment code?\n",
	      pd->ncomps,pd->nedges);
      exit(1);

    }

    /* We now know that component "comp" has "nedges" total edges.
       Use this to allocate space in the correct record in the pd->comp
       array. */

    pd->comp[comp].nedges = nedges;
    pd->comp[comp].edge   = calloc(nedges,sizeof(pd_idx_t));
    assert(pd->comp[comp].edge != NULL);

  } /* End of loop over all components. */

  /* We now need to actually build the component records from the
     data in the edge_assignment array. We've made space, so this
     is just a matter of cruising throught the assignment array
     again, putting everybody in the right place. */

  for(edge=0;edge<pd->nedges;edge++) {

    assert(edge_assigned[edge].comp >= 0 && edge_assigned[edge].comp < pd->ncomps);
    assert(edge_assigned[edge].pos >= 0 && edge_assigned[edge].pos < pd->comp[edge_assigned[edge].comp].nedges);

    pd->comp[edge_assigned[edge].comp].edge[edge_assigned[edge].pos] = edge;

  }

  /* We can now free the edge_assignment data, since we've assigned
     all edges to components. */

  free(edge_assigned);
  edge_assigned = NULL;

  /* Step 2. Sort the components longest-first. */

  qsort(pd->comp,pd->ncomps,sizeof(pd_component_t),pd_component_cmp);

  /* Step 2a. This is the canonical order for components (remember, this
     only ran if we didn't have components in the first place), so we assign
     tags now. */

  pd_tag_t tag = 'A';
  pd_idx_t i;

  for(i=0;i<pd->ncomps;i++,tag++) {
    pd->comp[i].tag = tag; }

  /* Step 3. Build a translation table for edge reordering and update the
     comp records to standard edge numbers. */

  pd_idx_t old_edge_num[pd->MAXEDGES];
  pd_idx_t new_edge_num[pd->MAXEDGES];

  pd_idx_t comp;

  for(i=0,comp=0;comp<pd->ncomps;comp++) {

    for(edge=0;edge<pd->comp[comp].nedges;edge++,i++) {

      old_edge_num[i] = pd->comp[comp].edge[edge];
      new_edge_num[pd->comp[comp].edge[edge]] = i;

      pd->comp[comp].edge[edge] = i; /* This SHOULD be the edge here. */

    }

  }

  /* Step 4. Rearrange the edges buffer accordingly. */

  pd_edge_t new_edges[pd->MAXEDGES];

  for(edge=0;edge<pd->nedges;edge++) {

    new_edges[edge] = pd->edge[old_edge_num[edge]];

  }

  for(edge=0;edge<pd->nedges;edge++) {

    pd->edge[edge] = new_edges[edge];

  }

  /* Step 5. Rewrite the crossing data with the
     new edge numbers. */

  pd_pos_t j;

  for(i=0;i<pd->ncross;i++) {

    for(j=0;j<4;j++) {

      pd->cross[i].edge[j] =
	new_edge_num[pd->cross[i].edge[j]];

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

  next_pos = (this_pos + 3) % 4; /* Turn LEFT (counterclockwise) */
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

#define PD_UNSET_FACE -1

/* Each edge should occur on two faces, with a
   positive orientation on the "pos_face" and a
   negative orientation on the "neg_face". We also
   keep track of _where_ this edge is going to show
   up on the face. */

typedef struct face_assignment_struct {

  int pos_face;
  pd_idx_t pos_face_position;

  int neg_face;
  pd_idx_t neg_face_position;

} pdint_face_assignment;

void pdint_next_face_unassigned_edge_or(pd_code_t *pd,
					pdint_face_assignment *face_assigned,
					pd_idx_t *start_edge,pd_or_t *start_or,
					bool *all_assigned)

/* Search the array of face assignments for unset pos or neg faces, in the
   order edge 0-pos_face, edge 1-pos_face, ..., edge n-pos_face,
         edge 0-neg_face, edge 1-neg_face, ..., edge n-neg_face.
   When we find one, reset the start_edge/start_or values appropriately.
*/

{

  /* First, the usual defensive asserts designed to make sure the
     input isn't garbage. */

  assert(face_assigned != NULL);
  assert(pd != NULL);
  assert(start_edge != NULL);
  assert(start_or != NULL);

  assert(*start_edge < pd->nedges);
  assert(*start_or == PD_POS_ORIENTATION || *start_or == PD_NEG_ORIENTATION);
  assert(all_assigned != NULL);

  /* Now we can actually start work. */

  pd_idx_t edge;

  if (*start_or == PD_POS_ORIENTATION) {

    for(edge=*start_edge+1;edge < pd->nedges;edge++) {

      if (face_assigned[edge].pos_face == PD_UNSET_FACE) {

	*start_edge = edge; *start_or = PD_POS_ORIENTATION;
	*all_assigned = false;
	return;

      }

    }

    /* Now we're cruising the negative orientations from the
       start. */

    for(edge=0;edge < pd->nedges;edge++) {

      if (face_assigned[edge].neg_face == PD_UNSET_FACE) {

	*start_edge = edge; *start_or = PD_NEG_ORIENTATION;
	*all_assigned = false;
	return;

      }

    }

    /* Now we've actually checked everything. */

    *start_edge = 0; *start_or = PD_UNSET_ORIENTATION;
    *all_assigned = true;
    return;

  } else if (*start_or == PD_NEG_ORIENTATION) { /* We started with a negative orientation */

    for(edge = *start_edge+1;edge < pd->nedges;edge++) {

      if (face_assigned[edge].neg_face == PD_UNSET_FACE) {

	*start_edge = edge; *start_or = PD_NEG_ORIENTATION;
	*all_assigned = false;
	return;

      }

    }

    /* We've actually reached the end of the list. */

    *start_edge = 0; *start_or = PD_UNSET_ORIENTATION;
    *all_assigned = true;
    return;

  }

  /* It should actually be impossible to get here, but it
     means that the input orientation was bogus. */

  fprintf(stderr,"pdint_next_face_unassigned_edge_or: Shouldn't get here.\n");
  exit(1);

}

void pd_regenerate_faces(pd_code_t *pd)

/* Generates faces of the polyhedron corresponding to the
   crossing and edge data in pd.  Add that data to the
   pd_code_t structure pd.

   The basic idea here is that we're computing the orbits
   of pdint_next_edge_on_face. To do that, we need to keep track
   of the edges and orientations we've used already.

   Again, this has to be done kind of carefully, because we
   may or may not have allocated memory for the individual face
   records, and in any case don't know how much memory to allocate
   for each face until we're through the orbit. */

{
  assert(pd != NULL);
  assert(pd->ncross > 0);
  assert(pd->nedges > 1);

  /* Step 0. Look at the face information we've got, and if
     the buffers are assigned, free them. We don't need to
     worry about the pd->face buffer itself because it's
     allocated by pd_code_new. */

  assert(pd->face != NULL);
  pd_idx_t face;

  for(face=0;face<pd->MAXFACES;face++) {

    if (pd->face[face].edge != NULL) { free(pd->face[face].edge); pd->face[face].edge = NULL; }
    if (pd->face[face].or != NULL) { free(pd->face[face].or); pd->face[face].or = NULL; }
    pd->face[face].nedges = 0;

  }

  /* Step 1. Clean the buffer of face assignments, reflecting
     that no edge has yet been assigned to any face, in either
     the positive or negative orientation. */

  pdint_face_assignment *face_assigned;
  face_assigned = calloc(pd->nedges,sizeof(pdint_face_assignment));
  assert(face_assigned != NULL);

  pd_idx_t edge;

  for(edge=0;edge < pd->nedges;edge++) {

    face_assigned[edge].pos_face = PD_UNSET_FACE;
    face_assigned[edge].pos_face_position = pd->nedges+1; /* Set this to something impossible */

    face_assigned[edge].neg_face = PD_UNSET_FACE;
    face_assigned[edge].neg_face_position = pd->nedges+1;

  }

  /* Step 3. Start the main (face-counting) loop. */

  pd_idx_t start_edge;
  pd_or_t  start_or;
  bool     all_assigned = false;

  for(pd->nfaces=0, start_edge = 0, start_or = PD_POS_ORIENTATION;
      /* Loop over faces. */
      !all_assigned;
      pd->nfaces++,pdint_next_face_unassigned_edge_or(pd,face_assigned,
						      &start_edge,&start_or,
						      &all_assigned)) {

    assert(start_edge < pd->nedges &&
	   (start_or == PD_POS_ORIENTATION || start_or == PD_NEG_ORIENTATION));

    /* Now we have an unused edge/orientation pair to start
       the next loop (around the face) on. */

    pd_idx_t this_edge = start_edge;
    pd_or_t  this_or = start_or;
    pd_idx_t nedges = 0; /* Stores the number of edges in this face. */

    do {

      /* Assign this edge/orientation pair to the current face. */

      assert(this_or == PD_POS_ORIENTATION || this_or == PD_NEG_ORIENTATION);

      if (this_or == PD_POS_ORIENTATION) {

	face_assigned[this_edge].pos_face = pd->nfaces;
	face_assigned[this_edge].pos_face_position = nedges;

      } else if (this_or == PD_NEG_ORIENTATION) {

	face_assigned[this_edge].neg_face = pd->nfaces;
	face_assigned[this_edge].neg_face_position = nedges;

      }

      pd_idx_t nedges_in_this_face = nedges;
      assert(nedges_in_this_face < pd->nedges);

      /* Increment the edge to the next edge on the face. */

      pdint_next_edge_on_face(pd,this_edge,this_or,&this_edge,&this_or);
      nedges++;

    } while (!(this_edge == start_edge && this_or == start_or));
    /* End of loop around this face. */

    /* Make sure the number of faces is sane. */

    if (pd->nfaces >= pd->MAXFACES) {

      fprintf(stderr,
	      "pd_regenerate_faces: Found %d faces on %d crossing diagram (maxfaces == %d).\n"
	      "                     Suspect error in face assignment code?\n",
	      pd->nfaces,pd->ncross,pd->MAXFACES);
      exit(1);

    }

    /* We now know the number of edges in this face. Use this to allocate space. */

    pd->face[pd->nfaces].nedges = nedges;
    pd->face[pd->nfaces].edge = calloc(nedges,sizeof(pd_idx_t));
    pd->face[pd->nfaces].or = calloc(nedges,sizeof(pd_or_t));
    assert(pd->face[pd->nfaces].edge != NULL && pd->face[pd->nfaces].or != NULL);

  } /* End of loop over faces. */

  /* Now that we've allocated space in the individual face records, we can take
     another tour through the face_assigned array and put edge numbers and orientations
     in the face records as appropriate. */

  for(edge=0;edge<pd->nedges;edge++) {

    assert(face_assigned[edge].pos_face < pd->nfaces);
    assert(face_assigned[edge].pos_face_position < pd->face[face_assigned[edge].pos_face].nedges);

    pd->face[face_assigned[edge].pos_face].edge[face_assigned[edge].pos_face_position] = edge;
    pd->face[face_assigned[edge].pos_face].or[face_assigned[edge].pos_face_position] = PD_POS_ORIENTATION;

    assert(face_assigned[edge].neg_face < pd->nfaces);
    assert(face_assigned[edge].neg_face_position < pd->face[face_assigned[edge].neg_face].nedges);

    pd->face[face_assigned[edge].neg_face].edge[face_assigned[edge].neg_face_position] = edge;
    pd->face[face_assigned[edge].neg_face].or[face_assigned[edge].neg_face_position] = PD_NEG_ORIENTATION;

  }

  /* We have finished building the faces, so we can release the face_assigned buffer. */

  free(face_assigned);
  face_assigned = NULL;

  /* Now for each face, make sure it is in canonical order. */

  for(face=0;face<pd->nfaces;face++) {

    pd_canonorder_face(&(pd->face[face]),PD_POS_ORIENTATION);

  }

  /* Now sort the faces using pd_face_cmp. */

  qsort(pd->face,pd->nfaces,sizeof(pd_face_t),pd_face_cmp);

  /* Finally, check to make sure the face count is correct. */

  if (pd->nfaces != pd->ncross + 2) {

    pd_error(SRCLOC,"pd->nfaces (%d) does not equal pd->ncross + 2 (%d) for pd %PD\n",
	     pd,pd->nfaces,pd->ncross);
    exit(1);

  }

}

void pd_regenerate(pd_code_t *pd)

/* Starting with a valid list of crossings, regenerate
   everything else. */

{
  assert(pd != NULL);

  pd_regenerate_crossings(pd); /* This can always be run */

  /* Now regenerating the edges will destroy crossing sign information.
     Therefore, we want to run it

     a) definitely, if all the crossings are set to
        PD_UNSET_ORIENTATION.

     b) otherwise, only if the edges are actually corrupted or there
        are no edges yet (in which case, the crossing information
        doesn't mean anything either, so we should lose it).

  */

  bool crossing_signs_all_unset = true;
  pd_idx_t i;

  for(i=0;i<pd->ncross && crossing_signs_all_unset;i++) {

    if (pd->cross[i].sign != PD_UNSET_ORIENTATION) {

      crossing_signs_all_unset = false;

    }

  }

  if (crossing_signs_all_unset) {

    pd_regenerate_edges(pd);

  } else {

    pd_idx_t PD_VERBOSITY_STORE = PD_VERBOSE;
    PD_VERBOSE = 0;

    if (!pd_edges_ok(pd)) {

      for(i=0;i<pd->ncross;i++) {

	pd->cross[i].sign = PD_UNSET_ORIENTATION;

      }

      pd_regenerate_edges(pd);

    }

    PD_VERBOSE = PD_VERBOSITY_STORE;

  }

  pd_regenerate_comps(pd);
  pd_regenerate_faces(pd);
  pd_regenerate_hash(pd);

}


/* pd sanity checking */

bool pd_cross_ok(pd_code_t *pd)
/* Checks to see that the edge numbers referenced in the
   crossings are consecutive and between 0 and
   pd->nedges-1, and that the crossing data is sorted,
   and that the crossing sign has a legal value for an orientation. */
{
  assert(pd != NULL);
  if (pd->ncross == 0) { return true; }
  /* The 0-crossing unknot diagram has nothing to check */

  int *edge_seen;
  edge_seen = calloc(pd->MAXEDGES,sizeof(int));
  assert(edge_seen != NULL);

  pd_idx_t edge, cross;
  pd_pos_t pos;

  for(edge=0;edge<pd->MAXEDGES;edge++) { edge_seen[edge] = 0; }

  for(cross=0;cross<pd->ncross;cross++) {

    if (!(pd->cross[cross].sign == PD_POS_ORIENTATION || pd->cross[cross].sign == PD_NEG_ORIENTATION || pd->cross[cross].sign == PD_UNSET_ORIENTATION)) {

      return pd_error(SRCLOC,"%CROSS contains illegal sign %d in pd %PD",pd,
		      cross,pd->cross[cross].sign);

    }

    for(pos=0;pos<4;pos++) {

      edge = pd->cross[cross].edge[pos];

      if (!(edge < pd->nedges)) {

	return pd_error(SRCLOC,"%CROSS contains illegal edge number %d (not < number of edges %d) in pd %PD",
			pd,
			cross,edge,pd->nedges);

      }

      edge_seen[edge]++;

    }

  }

  /* Now we check that each edge was seen twice. */

  for(edge=0;edge<pd->nedges;edge++) {

    if (edge_seen[edge] != 2) {

      if ((pd->edge[edge].head == PD_UNSET_IDX &&
	   pd->edge[edge].tail == PD_UNSET_IDX &&
	   pd->edge[edge].headpos == PD_UNSET_POS &&
	   pd->edge[edge].tailpos == PD_UNSET_POS) &&
	  edge_seen[edge] == 0) {

	/* It's ok. This isn't an error. If we have split, unknotted
	   components (this occurs only in testing), we can get stray
	   edges which don't interact with any crossing. */

      } else {

	return pd_error(SRCLOC,"%EDGE not seen twice in %PD",pd,edge);

      }
    }

  }

  /* Now check sort order. */

  for(cross=1;cross<pd->ncross;cross++) {

    if (pd_cross_cmp(&(pd->cross[cross]),&(pd->cross[cross-1])) < 0) {

      return pd_error(SRCLOC,"%CROSS and %CROSS are out of order in %PD",
		      pd,cross-1,cross);

    }

  }

  free(edge_seen);
  return true;

}

bool pd_edges_ok(pd_code_t *pd)

{
  assert(pd != NULL);

  pd_idx_t edge;

  if (pd->ncross != 0 && pd->nedges != pd->ncross*2) {
    /* We first check for the correct number of edges! */

    return pd_error(SRCLOC,"number of edges (%d) is not 2 * number of crossings (%d)\n",pd,
		    pd->nedges,pd->ncross);

  }

  for(edge=0;edge < pd->nedges;edge++) {

    pd_edge_t *e = &(pd->edge[edge]);

    if (e->head >= pd->ncross ||
	e->tail >= pd->ncross) {

      if (pd->ncross != 0) { /* In which case this is XFAIL, because
				head and tail will be PD_UNSET_IDX */

	return pd_error(SRCLOC,"%EDGE contains reference to illegal crossing in pd %PD\n",pd,edge);

      } else if (e->head != PD_UNSET_IDX || e->tail != PD_UNSET_IDX) {

	return pd_error(SRCLOC,"%EDGE in 0-crossing diagram has head or tail not PD_UNSET_IDX",pd,edge);

      }

    }

    if (e->headpos >= 4 ||
	e->tailpos >= 4) {

      if (pd->ncross != 0) {

	return pd_error(SRCLOC,"%EDGE contains reference to illegal position in pd %PD\n",pd,edge);

      } else if (e->headpos != PD_UNSET_POS || e->tailpos != PD_UNSET_POS) {

	return pd_error(SRCLOC,"%EDGE in 0-crossing diagram has headpos or tailpos not PD_UNSET_POS",pd,edge);

      }

    }

    if (pd->ncross > 0) { /* The next checks really only make sense
			     if the diagram HAS crossings. */

      pd_crossing_t *hc = &(pd->cross[e->head]), *tc = &(pd->cross[e->tail]);

      if (hc->edge[e->headpos] != edge) {

	return pd_error(SRCLOC,"%EDGE is not in correct position at head %CROSS in pd %PD\n",pd,edge,e->head);

      }

      if (tc->edge[e->tailpos] != edge) {

	return pd_error(SRCLOC,"%EDGE is not in correct position at tail %CROSS in pd %PD\n",pd,edge,e->tail);

      }

      /* We now check for consistent orientations along the components. */

      pd_pos_t nextpos = (e->headpos+2)%4;
      pd_pos_t prevpos = (e->tailpos+2)%4;

      if (pd->edge[hc->edge[nextpos]].tail != e->head ||
	  pd->edge[hc->edge[nextpos]].tailpos != nextpos) {

	return pd_error(SRCLOC,"%EDGE is supposed to be followed by %EDGE at %CROSS, but orientations disagree\n",pd,edge,hc->edge[nextpos]);

      }

      if (pd->edge[tc->edge[prevpos]].head != e->tail ||
	  pd->edge[tc->edge[prevpos]].headpos != prevpos) {

	return pd_error(SRCLOC,"%EDGE is supposed to be followed by %EDGE at %CROSS, but orientations disagree\n",pd,edge,tc->edge[prevpos]);

      }

    }

  }

  return true;

}

int pd_tag_cmp(const void *A,const void *B) {

  pd_tag_t *tagA = (pd_tag_t *)(A);
  pd_tag_t *tagB = (pd_tag_t *)(B);

  return (int)(*tagA - *tagB);

}

bool pd_comps_ok(pd_code_t *pd)

/* Check that the edges are numbered correctly around the components
   and that the edges are oriented head-to-tail, also checks that every
   edge is present in a component. Last, checks that tags are assigned
   to every component in A..Z (and forward) and that the tag set is consecutive
   and contains unique elements only. */

{
  pd_idx_t comp, edge, next_edge;
  pd_idx_t correct_edgenum;

  assert(pd != NULL);

  if (pd->ncomps > pd->MAXCOMPONENTS) {

    return pd_error(SRCLOC,"Number of components %d > pd->MAXCOMPONENTS %d in pd %PD",pd,pd->ncomps,pd->MAXCOMPONENTS);

  }

  if (pd->ncross == 0) {

    if (pd->ncomps != 1) {

      return pd_error(SRCLOC,"Wrong number of components %d != 1 in 0-crossing diagram %PD",pd,pd->ncomps);

    }

    if (pd->comp[0].nedges != 1) {

       return pd_error(SRCLOC,"Wrong number of edges %d != 1 in component 0 of 0-crossing diagram %PD",pd,pd->comp[0].nedges);

    }

    if (pd->comp[0].edge[0] != 0) {

      return pd_error(SRCLOC,"Edge 0 in component 0 of 0-crossing pd %PD is %d != 0",pd,pd->comp[0].edge[0]);

    }

    return true;

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

  pd_tag_t *tagset;
  tagset = calloc(pd->ncomps,sizeof(pd_tag_t));
  assert(tagset != NULL);
  pd_idx_t i;

  for(i=0;i<pd->ncomps;i++) {

    tagset[i] = pd->comp[i].tag;

  }

  qsort(tagset,pd->ncomps,sizeof(pd_tag_t),pd_tag_cmp);

  for(i=1;i<pd->ncomps;i++) {

    if (tagset[i] == tagset[i-1]) {

      pd_idx_t compA, compB;
      bool Afound = false, Bfound = false;

      for(compA=0;compA<pd->ncomps && !Afound;compA++) {

	if (pd->comp[compA].tag == tagset[i]) {

	  Afound = true;

	}

      }

      for(compB=compA+1;compB<pd->ncomps && !Bfound;compB++) {

	if (pd->comp[compB].tag == tagset[i]) {

	  Bfound = true;

	}

      }

      return pd_error(SRCLOC,"components %COMP and %COMP of %PD have duplicate tags",pd,compA,compB);

    }

    if (tagset[i] < 'A') {

      pd_idx_t compA;
      bool Afound = false;

      for(compA=0;compA<pd->ncomps && !Afound;compA++) {

	if (pd->comp[compA].tag == tagset[i]) {

	  Afound = true;

	}

      }

      return pd_error(SRCLOC,"component %COMP has illegal tag (< 'A')",pd,compA);

    }

  }

  free(tagset);

  return true;

}


bool pd_faces_ok(pd_code_t *pd)

/* Checks face data */

{
  assert(pd != NULL);

  if (pd->nfaces > pd->MAXFACES) {

    return pd_error(SRCLOC,"Number of faces %d > pd->MAXFACES %d in pd %PD",pd,pd->nfaces,pd->MAXFACES);

  }

  if (pd->ncross == 0) {

    if (pd->nfaces != 2) {

      return pd_error(SRCLOC,"Wrong number of faces %d != 2 in 0-crossing diagram %PD",pd,pd->nfaces);

    }

    if (pd->face[0].nedges != 1) {

      return pd_error(SRCLOC,"Wrong number of edges %d != 1 in face 0 of 0-crossing diagram %PD",pd,pd->face[0].nedges);

    }

    if (pd->face[0].edge[0] != 0) {

      return pd_error(SRCLOC,"Illegal edge reference %d in face 0 of 0-crossing diagram %PD",pd,pd->face[0].edge[0]);

    }

    if (pd->face[1].nedges != 1) {

      return pd_error(SRCLOC,"Wrong number of edges %d != 1 in face 1 of 0-crossing diagram %PD",pd,pd->face[1].nedges);

    }

    if (pd->face[1].edge[0] != 0) {

      return pd_error(SRCLOC,"Illegal edge reference %d in face 0 of 0-crossing diagram %PD",pd,pd->face[1].edge[0]);

    }

    return true;

  }

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
			"in face %FACE of pd\n"
			"Oriented %EDGE and %EDGE don't meet\n"
			"correctly in pd %PD\n",
			pd,face,edge,face,nxt_edge,
			face,
			pd->face[face].edge[edge],
			pd->face[face].edge[nxt_edge]);
      }

    }

  }

  /* Check that each face is rotated into canonical order */

  for(face=0;face<pd->nfaces;face++) {

    pd_idx_t lowE = pd->face[face].edge[0];
    pd_idx_t i;

    for(i=1;i<pd->face[face].nedges;i++) {

      if (pd->face[face].edge[i] < lowE) {

	lowE = pd->face[face].edge[i];

      }

    }

    if (lowE != pd->face[face].edge[0]) {

      return pd_error(SRCLOC,"face %FACE is not in canonical order in %PD",pd,face);

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

bool pd_is_alternating(pd_code_t *pd)
/* Checks to see if pd is an alternating link. Returns an error if some crossings are not set. */
{
  pd_idx_t cr;

  for(cr=0;cr<pd->ncross;cr++) {

    if (pd->cross[cr].sign == PD_UNSET_ORIENTATION) {

      pd_error(SRCLOC,"can't decide whether pd %PD is alternating since it has unset crossing signs",
	       pd);
      exit(1);

    }

  }

  pd_idx_t comp;
  pd_idx_t compedge;

  for(comp=0;comp<pd->ncomps;comp++) {

    for(compedge=0;compedge<pd->comp[comp].nedges;compedge++) {

      pd_idx_t edge = pd->comp[comp].edge[compedge];

      /* We need to check that this edge goes UNDER->OVER or OVER->UNDER. */

      pd_idx_t incoming_under_pos, outgoing_under_pos;
      pd_idx_t incoming_over_pos, outgoing_over_pos;

      pd_understrand_pos(pd,pd->edge[edge].tail,&incoming_under_pos,&outgoing_under_pos);
      pd_overstrand_pos(pd,pd->edge[edge].tail,&incoming_over_pos,&outgoing_over_pos);

      if (outgoing_under_pos == pd->edge[edge].tailpos) { /* We are coming from an UNDER */

	pd_understrand_pos(pd,pd->edge[edge].head,&incoming_under_pos,&outgoing_under_pos);
	pd_overstrand_pos(pd,pd->edge[edge].head,&incoming_over_pos,&outgoing_over_pos);

	if (incoming_under_pos == pd->edge[edge].headpos) { /* We go TO an UNDER */

	  return false;

	} else if (incoming_over_pos != pd->edge[edge].headpos) { /* Some weird error */

	  pd_error(SRCLOC,"edge %EDGE enters crossing %CROSS but has tailpos\n"
		   "not equal to either the incoming UNDER pos (%d) or \n"
		   "incoming OVER pos (%d) in %PD",
		   pd,edge,pd->edge[edge].head,incoming_under_pos,incoming_over_pos);

	  exit(1);

	}

      } else if (outgoing_over_pos == pd->edge[edge].tailpos) { /* We are coming from an OVER */

	pd_understrand_pos(pd,pd->edge[edge].head,&incoming_under_pos,&outgoing_under_pos);
	pd_overstrand_pos(pd,pd->edge[edge].head,&incoming_over_pos,&outgoing_over_pos);

	if (incoming_over_pos == pd->edge[edge].headpos) { /* We go TO an OVER */

	  return false;

	} else if (incoming_under_pos != pd->edge[edge].headpos) { /* Some weird error */

	  pd_error(SRCLOC,"edge %EDGE enters crossing %CROSS but has tailpos\n"
		   "not equal to either the incoming UNDER pos (%d) or \n"
		   "incoming OVER pos (%d) in %PD",
		   pd,edge,pd->edge[edge].head,incoming_under_pos,incoming_over_pos);

	  exit(1);

	}

      } else { /* There's an error. */

	pd_error(SRCLOC,"edge %EDGE leaves crossing %CROSS but has tailpos\n"
		 "not equal to either the outgoing UNDER pos (%d) or \n"
		 "outgoing OVER pos (%d) in %PD",
		 pd,edge,pd->edge[edge].tail,outgoing_under_pos,outgoing_over_pos);

	exit(1);

      }

    }

  }

  return true;

}

pd_code_t *pd_copy(pd_code_t *pd)
/* Make a new-memory copy of pd. This can require some care, for instance if the face and comp arrays aren't allocated yet. */
{
  if (pd->ncross == 0) {

    return pd_build_unknot(0);

  }

  pd_code_t *pdA;
  pdA = pd_code_new(pd->ncross);

  assert(pdA->MAXVERTS      >= pd->ncross);
  assert(pdA->MAXEDGES      >= pd->nedges);
  assert(pdA->MAXCOMPONENTS >= pd->ncomps);
  assert(pdA->MAXFACES      >= pd->nfaces);

  pd_idx_t i, edge, cmp, face, cr;

  pdA->uid = pd->uid;

  pdA->ncross = pd->ncross;
  pdA->nedges = pd->nedges;
  pdA->ncomps = pd->ncomps;
  pdA->nfaces = pd->nfaces;

  for(i=0;i<PD_HASHSIZE;i++) { pdA->hash[i] = pd->hash[i]; }
  for(edge=0;edge<pd->nedges;edge++) { pdA->edge[edge] = pd->edge[edge]; };

  for(cmp=0;cmp<pd->ncomps;cmp++) {

    assert(   (pd->comp[cmp].nedges == 0 && pd->comp[cmp].edge == NULL)
	   || (pd->comp[cmp].nedges != 0 && pd->comp[cmp].edge != NULL));

    if (pd->comp[cmp].nedges != 0) {

      pdA->comp[cmp].nedges = pd->comp[cmp].nedges;
      pdA->comp[cmp].edge = calloc(pdA->comp[cmp].nedges,sizeof(pd_idx_t));
      assert(pdA->comp[cmp].edge != NULL);

      for(edge=0;edge<pdA->comp[cmp].nedges;edge++) { pdA->comp[cmp].edge[edge] = pd->comp[cmp].edge[edge]; }

      pdA->comp[cmp].tag = pd->comp[cmp].tag;

    }

  }

  for(cr=0;cr<pdA->ncross;cr++) {

    pdA->cross[cr] = pd->cross[cr];

  }

  for(face=0;face<pdA->nfaces;face++) {

    assert((pd->face[face].nedges == 0 && (pd->face[face].edge == NULL && pd->face[face].or == NULL)) ||
	   (pd->face[face].nedges != 0 && (pd->face[face].edge != NULL && pd->face[face].or != NULL)));

    if (pd->face[face].nedges != 0) {

      pdA->face[face].nedges = pd->face[face].nedges;
      pdA->face[face].edge   = calloc(pdA->face[face].nedges,sizeof(pd_idx_t));
      pdA->face[face].or     = calloc(pdA->face[face].nedges,sizeof(pd_or_t));

      assert(pdA->face[face].edge != NULL && pdA->face[face].or != NULL);

      for(edge=0;edge<pdA->face[face].nedges;edge++) {

	pdA->face[face].edge[edge] = pd->face[face].edge[edge];
	pdA->face[face].or[edge] = pd->face[face].or[edge];

      }

    }

  }

  return pdA;
}

void pd_write_KnotTheory(FILE *of, pd_code_t *pd)
/* Writes a pd code to file in the style of KnotTheory */
{
  pd_idx_t cross,under_in,under_out;
  pd_pos_t pos;


  fprintf(of,"PD[");

  for(cross=0;cross<pd->ncross;cross++) {
    pd_understrand_pos(pd, cross, &under_in, &under_out);
    fprintf(of,"X[");
    for(pos=0;pos<4;pos++) {
      fprintf(of,"%u",(pd->cross[cross].edge[(under_in+pos)%4])+1);
      if(pos<3) {
	fprintf(of,",");
      }
    }
    fprintf(of,"]");
    if(cross<(pd->ncross)-1) {
      fprintf(of,",");
    }
  }

  fprintf(of,"]\n");
}



void pd_write(FILE *of,pd_code_t *pd)

/* Writes a pd code (including all precomputed information) in human readable ASCII format. */
/* It's true that this seems inefficient, but it's better to compress the resulting */
/* file for storage than to muck around with a binary file format. */

/* Format:

   pd   <hash> <uid>
   nv   <nverts>
   <nv lines of crossing information in the format edge edge edge edge> <sign +/- optional>
   ne   <nedges>
   <ne lines of crossing information in the format tail, tailpos -> head, headpos>
   nc   <ncomps>
   <nc lines, each containing nedges edge edge .... edge and a final "tag (character)">
   nf   <nfaces>
   <nf lines, each in the format nedges edge edge ... edge giving face information counterclockwise>

   (nts: Add ordering information!)

*/

{
  /* First, we do a bit of sanity checking */

  if (pd == NULL) { return; }
  if (pd->ncross > pd->MAXVERTS || pd->nedges > pd->MAXEDGES || pd->ncomps > pd->MAXCOMPONENTS || pd->nfaces > pd->MAXFACES) {

    fprintf(of,"INVALID PD CODE\n");
    printf("%s (%d): Invalid PD code passed for writing.\n",SRCLOC);
    return;

  }

  pd_idx_t cross,edge,face,comp;
  pd_pos_t  pos;

  fprintf(of,"pd %s %lu\n",pd->hash,(unsigned long int)(pd->uid));

  /* Crossing data */

  fprintf(of,"nv %u\n",(unsigned int)(pd->ncross));

  for(cross=0;cross<pd->ncross;cross++) {

    for(pos=0;pos<4;pos++) {

      fprintf(of,"%u ",(unsigned int)(pd->cross[cross].edge[pos]));

    }

    if (pd->cross[cross].sign != PD_UNSET_ORIENTATION) {

      fprintf(of,"%c",pd_print_or(pd->cross[cross].sign));

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

    fprintf(of,"tag %c",(char)(pd->comp[comp].tag));
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

void pd_write_c(FILE *outfile, pd_code_t *pd, char *pdname)
/* Writes a c procedure which recreates the pd code pd.
   The procedure will be called pd_create_(pdname) */
{
  fprintf(outfile,"pd_code_t *pd_create_%s() { \n\n",pdname);
  fprintf(outfile,
	  "/* This procedure is machine generated by pd_write_c */\n"
	  "/* and probably shouldn't be hand-edited. */\n\n");

  fprintf(outfile,
	  "pd_code_t *pd;\n"
	  "pd = pd_code_new(%d);\n"
	  "assert(pd != NULL);\n",pd->MAXVERTS);

  fprintf(outfile,
	  "pd->ncross = %d;\n"
	  "pd->nedges = %d;\n"
	  "pd->ncomps = %d;\n"
	  "pd->nfaces = %d;\n"
	  "sprintf(pd->hash,\"%%s\",\"%s\");\n",
	  pd->ncross,pd->nedges,pd->ncomps,pd->nfaces,pd->hash);

  pd_idx_t i,j;

  fprintf(outfile,"\n/* Crossing data. */\n\n");

  /* Now rebuild the crossing buffer. */

  for(i=0;i<pd->ncross;i++) {

    for(j=0;j<4;j++) {

      fprintf(outfile,"pd->cross[%d].edge[%d] = %d;\n",i,j,pd->cross[i].edge[j]);

    }

    fprintf(outfile,"pd->cross[%d].sign = %d;\n\n",i,pd->cross[i].sign);

  }

  /* The edge buffer... */

  fprintf(outfile,"\n/* Edge data */\n\n");

  for(i=0;i<pd->nedges;i++) {

    fprintf(outfile,
	    "pd->edge[%d].head = %d;\n"
	    "pd->edge[%d].headpos = %d;\n"
	    "pd->edge[%d].tail = %d;\n"
	    "pd->edge[%d].tailpos = %d;\n\n",
	    i, pd->edge[i].head, i, pd->edge[i].headpos,
	    i, pd->edge[i].tail, i, pd->edge[i].tailpos);

  }

  /* Component data */

  fprintf(outfile,"\n/* Component Data */\n\n");

  for(i=0;i<pd->ncomps;i++) {

    fprintf(outfile,
	    "pd->comp[%d].nedges = %d;\n"
	    "pd->comp[%d].tag = '%c';\n\n",
	    i,pd->comp[i].nedges,i,pd->comp[i].tag);

    fprintf(outfile,
	    "pd->comp[%d].edge = calloc(pd->comp[%d].nedges,sizeof(pd_idx_t));\n"
	    "assert(pd->comp[%d].edge != NULL);\n\n",
	    i,i,i);

    for(j=0;j<pd->comp[i].nedges;j++) {

      fprintf(outfile,
	      "pd->comp[%d].edge[%d] = %d;\n",i,j,pd->comp[i].edge[j]);

    }

    fprintf(outfile,"\n");

  }

  /* Face data. */

  fprintf(outfile,"\n/* Face data */\n\n");

  for(i=0;i<pd->nfaces;i++) {

    fprintf(outfile,
	    "pd->face[%d].nedges = %d;\n"
	    "pd->face[%d].edge = calloc(pd->face[%d].nedges,sizeof(pd_idx_t));\n"
	    "pd->face[%d].or = calloc(pd->face[%d].nedges,sizeof(pd_or_t));\n"
	    "assert(pd->face[%d].edge != NULL);\n"
	    "assert(pd->face[%d].or != NULL);\n\n",
	    i,pd->face[i].nedges,i,i,i,i,i,i);

    for(j=0;j<pd->face[i].nedges;j++) {

      fprintf(outfile,
	      "pd->face[%d].edge[%d] = %d;\n"
	      "pd->face[%d].or[%d] = %d;\n\n",
	      i,j,pd->face[i].edge[j],
	      i,j,pd->face[i].or[j]);

    }

  }

  fprintf(outfile,
	  "\n/* End of data. */\n\n"
	  "assert(pd_ok(pd));\n"
	  "return pd;\n\n"
	  "}\n\n");

}

pd_code_t *pd_read_err(FILE *infile, int *err)

/* Reads an (ASCII) pd code written by pd_write. Return a pointer if we succeed, NULL if we fail. */

{
  /* pd   <hash> <uid> */ /* Remember that the hash is base64 encoded using a-zA-Z;: character set */

  unsigned long int input_temp,input_temp2,input_temp3,input_temp4;
  char hash[PD_HASHSIZE];
  pd_uid_t uid;

  /* The first thing we do is spin through whitespace until we encounter a "p" */
  int read_char;
  for(read_char = fgetc(infile);isspace(read_char);read_char = fgetc(infile));
  ungetc(read_char,infile);

  /* We now read the rest of the line, to see if it contains a hash
     and uid. We're going to bet on hashes (and hence, lines) no
     longer than 4096 characters here. */

  char pd_line[4096];

  if (fgets(pd_line,4096,infile) == NULL) {

    pd_error(SRCLOC,"infile is already at EOF-- can't pd_read from it\n",NULL);
    pd_err_set(err, PD_EOF);
    return NULL;

  }

  /* Now we try the alternatives. If we're going to read a hash, we have
     to do a little dodge to construct the pattern for sscanf, because we
     want the sscanf pattern to reflect PD_HASHSIZE (which might change) */

  char hash_template[1024];
  sprintf(hash_template," pd %%%d[a-zA-Z0-9;:] %%lu ",PD_HASHSIZE);

  int sscanf_result = sscanf(pd_line,hash_template,hash,&input_temp);

  if (sscanf_result != 2) {

    /* Check for an incorrect hash. */

    char bogus_hash[4096];
    if (sscanf(pd_line," pd %4096s %lu ",bogus_hash,&input_temp) == 2) {

      pd_error(SRCLOC,
	       "first line of pdcode file\n"
	       "%s"
	       "is in format\n"
	       "pd <hash> <uid>\n"
	       "but hash == %s, which violates a hash rule\n"
	       "the hash may be at most %d characters long\n"
	       "and should only contain characters [a-zA-Z0-9;:]\n"
	       "\n"
	       "If you don't know the hash, it should be omitted\n"
	       "and this line should read\n"
	       "\n"
	       "pd\n",NULL,pd_line,bogus_hash,PD_HASHSIZE);
      pd_err_set(err, PD_BAD_FORMAT);
      return NULL;

    } else if (sscanf(pd_line," pd %4096s %lu ",bogus_hash,&input_temp) == 1) {

      pd_error(SRCLOC,
	       "first line of pdcode file appears to be in format\n"
	       "pd <hash> <uid>\n"
	       "and includes hash %s, but does not include uid (an unsigned integer)\n",NULL,bogus_hash);
      pd_err_set(err, PD_BAD_FORMAT);
      return NULL;

    }

    /* Well, it looks like we didn't even TRY to provide a hash and uid. */
    /* In this case, we should match "(whitespace)pd(whitespace)" */

    char pd_string[32];

    if (sscanf(pd_line," %32s ",pd_string) != 1) {

      pd_error(SRCLOC,
	       "first line of pdcode file is neither in format\n"
	       "pd <hash> <uid>\n"
	       "nor\n"
	       "pd\n",NULL);
      pd_err_set(err, PD_BAD_FORMAT);
      return NULL;

    }

    if (strcmp(pd_string,"pd") != 0) {

      pd_error(SRCLOC,
	       "first line of pdcode file contains\n"
	       "%s\n"
	       "which is neither in format\n"
	       "pd <hash> <uid>\n"
	       "nor\n"
	       "pd\n",NULL,pd_string);
      pd_err_set(err, PD_BAD_FORMAT);
      return NULL;

    }

    /* We DO actually match this. This means that we should set the
       hash and the uid ourselves. */

    uid = PD_UNSET_UID;
    sprintf(hash,"unset");

  } else { /* We matched the pd <hash> <uid> format */

    uid = (pd_uid_t)(input_temp);

  }

  /* nv   <nverts> */

  int cross,edge,comp,pos,face;
  if (fscanf(infile," nv %lu ",&input_temp) != 1) {

    pd_error(SRCLOC,
	     "second line of pdcode file should be in format\n"
	     "nv <nverts>\n"
	     "where <nverts> == number of crossings in pd_code\n",NULL);
    pd_err_set(err, PD_BAD_FORMAT);
    return NULL;

  }

  /* We now know the number of crossings to allocate space for,
     and we can call pd_new. We make sure that we have space for
     at least two more crossings than we see now (to prevent zero
     crossing diagrams from crashing everything). */

  if ((pd_idx_t)(input_temp) == 0) {

    /* If we're a zero-crossing unknot, there is one piece of
       information that can actually matter, which is the tag
       of the (single) component. We're going to discard the
       rest of the file, but first we'll run through and extract
       the "tag" (if present). */

    char rolling_buffer[4];
    pd_tag_t tag = 'A';

    rolling_buffer[0] = (char)(getc(infile));
    rolling_buffer[1] = (char)(getc(infile));
    rolling_buffer[2] = (char)(getc(infile));
    rolling_buffer[3] = 0;

    for(;!feof(infile);) {

      if (!strcmp(rolling_buffer,"tag")) {

	fscanf(infile," %c ",&tag);

      }

      rolling_buffer[0] = rolling_buffer[1]; rolling_buffer[1] = rolling_buffer[2];
      rolling_buffer[2] = (char)(getc(infile));

    }

    pd_code_t *outcode;
    outcode = pd_build_unknot(0);
    outcode->comp[0].tag = tag;

    return outcode;

  }

  pd_code_t *pd = pd_code_new((pd_idx_t)(input_temp));

  /* We now have to copy the hash into the new pd code. */

  strncpy(pd->hash,hash,PD_HASHSIZE);
  pd->uid = uid;
  pd->ncross = (pd_idx_t)(input_temp);

  if (pd->ncross > pd->MAXVERTS) {

    pd_error(SRCLOC,
	    "Reading pd code which appears to be valid but has %d crossings. "
             "We allocated space for pd->MAXVERTS = %d.\n",pd,
             pd->ncross,pd->MAXVERTS);
    pd_err_set(err, PD_BAD_FORMAT);
    pd_code_free(&pd);
    return NULL;

  }

  /* <nv lines of crossing information in the format edge edge edge edge> */
  /* There's an optional crossing sign at the end of the line (+ or -) */

  for(cross=0;cross<pd->ncross;cross++) {

    for(pos=0;pos<4;pos++) {

      if(fscanf(infile," %lu ",&input_temp) != 1) {

	pd_error(SRCLOC,"error on crossing %d (of %d), in pd (so far) %PD",pd,cross,pd->ncross);
        pd_err_set(err, PD_BAD_FORMAT);
	pd_code_free(&pd);
	return NULL;

      }

      pd->cross[cross].edge[pos] = (pd_idx_t)(input_temp);

    }

    int peek;
    peek = fgetc(infile);
    if (peek == '+') {

      pd->cross[cross].sign = PD_POS_ORIENTATION;

    } else if (peek == '-') {

      pd->cross[cross].sign = PD_NEG_ORIENTATION;

    } else {

      ungetc(peek,infile);
      pd->cross[cross].sign = PD_UNSET_ORIENTATION;

    }

  }

  /* ne   <nedges> */

  if (fscanf(infile," ne %lu ",&input_temp) != 1) {

    pd_error(SRCLOC,
	     "after crossing lines, should have\n"
	     "ne <nedges>\n"
	     "but this file does not\n",pd);
    pd_err_set(err, PD_BAD_FORMAT);
    pd_code_free(&pd);
    return NULL;

  }
  pd->nedges = (pd_idx_t)(input_temp);

  if (pd->nedges > pd->MAXEDGES) {

    pd_error(SRCLOC,
	    "Reading pd code which appears to be valid but has %d edges. "
	    "We allocated space for pd->MAXEDGES = %d.\n",
	    pd,pd->nedges,pd->MAXEDGES);
    pd_err_set(err, PD_BAD_FORMAT);
    pd_code_free(&pd);
    return NULL;

  }

  /* <ne lines of crossing information in the format tail, tailpos -> head, headpos> */

  for(edge=0;edge<pd->nedges;edge++) {

    if(fscanf(infile," %lu,%lu -> %lu,%lu ",
	      &input_temp,&input_temp2,&input_temp3,&input_temp4) != 4) {

      pd_error(SRCLOC,
	       "edge %d is not in the format\n"
	       "<crossing>,<pos> -> <crossing>,<pos>\n"
	       "where <crossing> and <pos> are positive integers\n",pd,edge);
      pd_err_set(err, PD_BAD_FORMAT);
      pd_code_free(&pd);
      return NULL;

    }

    pd->edge[edge].tail = (pd_idx_t)(input_temp);
    pd->edge[edge].tailpos = (pd_pos_t)(input_temp2);

    pd->edge[edge].head = (pd_idx_t)(input_temp3);
    pd->edge[edge].headpos = (pd_pos_t)(input_temp4);

  }

  /* nc   <ncomps> */

  if (fscanf(infile,"nc %lu ",&input_temp) != 1) {

    pd_error(SRCLOC,
	     "after edge lines, pdcode file should have a line\n"
	     "nc <ncomps>\n"
	     "where <ncomps> is a positive integer\n"
	     "but this one doesn't",pd);
    pd_code_free(&pd); return NULL;

  }
  pd->ncomps = (pd_idx_t)(input_temp);

  if (pd->ncomps > pd->MAXCOMPONENTS) {

    pd_error(SRCLOC,
	    "Reading pd code which appears to be valid but has %d components. "
	    "We allocated space for pd->MAXCOMPONENTS = %d in %PD\n",
	    pd,pd->ncomps,pd->MAXCOMPONENTS);
    pd_err_set(err, PD_BAD_FORMAT);
    pd_code_free(&pd);
    return NULL;

  }

  /* <nc lines, each containing nedges : +/- edge +/- edge .... +/- edge> */

  pd_tag_t next_tag = 'A';

  for(comp=0;comp<pd->ncomps;comp++) {

    if (fscanf(infile," %lu : ",&input_temp) != 1) {

      pd_error(SRCLOC,
	       "component %d should start with\n"
	       "<nedges> : \n"
	       "where <nedges> is a positive integer\n"
	       "but this file doesn't\n",pd,comp);
      pd_err_set(err, PD_BAD_FORMAT);
      pd_code_free(&pd);
      return NULL;

    }

    pd->comp[comp].nedges = (pd_idx_t)(input_temp);

    /* We haven't allocated space in the component. */

    pd->comp[comp].edge = calloc(pd->comp[comp].nedges,sizeof(pd_idx_t));
    assert(pd->comp[comp].edge != NULL);

    if (pd->comp[comp].nedges > pd->MAXEDGES) {

      pd_error(SRCLOC,
               "Reading component which appears to be valid but has %d edges. "
               "We expect only a maximum of pd->MAXEDGES = %d.\n",
               pd,pd->nedges,pd->MAXEDGES);
      pd_err_set(err, PD_BAD_FORMAT);
      pd_code_free(&pd);
      return NULL;

    }

    for(edge=0;edge<pd->comp[comp].nedges;edge++) {

      if(fscanf(infile," %lu ",&input_temp) != 1) {

	pd_error(SRCLOC,"edge %d of component %d should be a positive integer\n"
		 "but isn't in this file\n",pd,edge,comp);
        pd_err_set(err, PD_BAD_FORMAT);
	pd_code_free(&pd);
        return NULL;

      }

      pd->comp[comp].edge[edge] = (pd_idx_t)(input_temp);


    }

    /* Now we're at the end, and the next character is either "tag X" or a newline. */
    /* Skip whitespace until we find a character... */

    char testchar;
    fscanf(infile," %c",&testchar);

    if (testchar == 't') { /* for tag */

      ungetc((int)(testchar),infile);
      if (fscanf(infile," tag %c ",&(pd->comp[comp].tag)) != 1) {

	pd_error(SRCLOC,
		 "component %d of pd_code has extra characters after %d edges\n"
		 "which start with t (so we're assuming a tag) but they don't match\n"
		 "tag <tag>\n"
		 "where <tag> is a single ASCII character\n",
		 pd,comp,edge);
        pd_err_set(err, PD_BAD_FORMAT);
	pd_code_free(&pd);
	return NULL;
      }

    } else {

      ungetc((int)(testchar),infile);
      pd->comp[comp].tag = next_tag++;

    }

  }

  /*nf   <nfaces> */

  if (fscanf(infile," nf %lu ",&input_temp) != 1) {

    pd_error(SRCLOC,
	     "after component data, pdfile should contain\n"
	     "nf <nfaces>\n"
	     "where <nfaces> is a positive integer\n"
	     "but this one doesn't\n",pd);
    pd_err_set(err, PD_BAD_FORMAT);
    pd_code_free(&pd);
    return NULL;

  }
  pd->nfaces = (pd_idx_t)(input_temp);

  if (pd->nfaces > pd->MAXFACES) {

    pd_error(SRCLOC,
	    "Reading pd code which appears to be valid but has %d faces. "
	     "We expected a maximum of pd->MAXFACES = %d in %PD\n",pd,
	     pd->nfaces,pd->MAXFACES);
    pd_err_set(err, PD_BAD_FORMAT);
    pd_code_free(&pd);
    return NULL;

  }

  /* <nf lines, each in the format nedges edge edge
     ... edge giving face info counterclockwise> */

  for(face=0;face<pd->nfaces;face++) {

    if (fscanf(infile," %lu : ",&input_temp) != 1) {

      pd_error(SRCLOC,
	       "face record %d (of %d) is expected to begin\n"
	       "<nedges> : \n"
	       "where <nedges> is a positive integer\n"
	       "but this one doesn't in %PD\n",pd,
	       face,pd->nfaces);
      pd_err_set(err, PD_BAD_FORMAT);
      pd_code_free(&pd);
      return NULL;

    }

    pd->face[face].nedges = (pd_idx_t)(input_temp);

    if (pd->face[face].nedges > pd->MAXEDGES) {

      pd_error(SRCLOC,
	       "Reading face which appears to be valid but has %d edges. "
	       "We expected a maximum of pd->MAXEDGES = %d.\n",
	       pd,pd->face[face].nedges,pd->MAXEDGES);
      pd_err_set(err, PD_BAD_FORMAT);
      pd_code_free(&pd);
      return NULL;

    }

    pd->face[face].edge = calloc(pd->face[face].nedges,sizeof(pd_idx_t));
    pd->face[face].or = calloc(pd->face[face].nedges,sizeof(pd_or_t));
    assert(pd->face[face].edge != NULL && pd->face[face].or != NULL);

    for(edge=0;edge<pd->face[face].nedges;edge++) {

      char orientation[2];
      if(fscanf(infile," %1[+-]s ",orientation) != 1) {

	pd_error(SRCLOC,
		 "edge %d of face %d in pdcode file\n"
		 "doesn't have a sign which matches [+-]\n",
		 pd,edge,face);
        pd_err_set(err, PD_BAD_FORMAT);
	pd_code_free(&pd);
	return NULL;

      }

      if (fscanf(infile," %lu ",&input_temp) != 1) {

	pd_error(SRCLOC,
		 "edge %d of face %d in pdcode file\n"
		 "has sign (+/-) but doesn't match\n"
		 "<+-> <edge>\n"
		 "where <edge> is a positive integer\n",
		 pd,edge,face);
        pd_err_set(err, PD_BAD_FORMAT);
	pd_code_free(&pd); return NULL;
      }

      pd->face[face].edge[edge] = (pd_idx_t)(input_temp);
      pd->face[face].or[edge] = (orientation[0] == '+') ? PD_POS_ORIENTATION : PD_NEG_ORIENTATION;

    }

  }

  pd_regenerate_crossings(pd); /* If the crossings were not in canonical order, fix it. */

  pd_err_check(pd_ok(pd), err, PD_NOT_OK, pd);
  if (strstr(pd->hash,"unset") != NULL) { pd_regenerate_hash(pd); }

  return pd;

}
pd_code_t *pd_read(FILE *infile) {
    return pd_read_err(infile, NULL);
}

pd_code_t *pd_read_KnotTheory(FILE *infile)

/*

Reads an (ASCII) pd code written by Mathematica. Return a pointer if we succeed, NULL if we fail.
These files will be single-line files which look like:

PD[X[1, 6, 2, 7], X[3, 8, 4, 9], X[5, 10, 6, 1], X[7, 2, 8, 3], X[9, 4, 10, 5]]

Here each "X" denotes a crossing, with X[i,j,k,l] denoting a crossing in the form

            k
            ^
            |
     l <----------> j
            |
            |
            i

the direction of the l,j strand is determined by the ordering of the edge numbers
l and j. We know that each component is ordered consecutively by edges, but there
is no particular ordering of the components in one of these codes.

We are going to use a large buffer to store the incoming crossings, then discard the
ones we don't need.

*/

{
  pd_crossing_t *crbuf;
  int crbuf_size=10000, crbuf_used = 0;
  crbuf = calloc(crbuf_size,sizeof(pd_crossing_t));
  assert(crbuf != NULL);

  if (fscanf(infile," PD [ X[ %d, %d, %d, %d ]",
	     &(crbuf[0].edge[0]),
	     &(crbuf[0].edge[1]),
	     &(crbuf[0].edge[2]),
	     &(crbuf[0].edge[3])) != 4) {

    printf("pd_read_KnotTheory: Can't read first crossing from KnotTheory format PD code file");
    free(crbuf);
    return NULL;

  }

  crbuf_used++;

  for(;fgetc(infile) == ',';crbuf_used++) {

    if (fscanf(infile," X[ %d, %d, %d, %d ]",
	       &(crbuf[crbuf_used].edge[0]),
	       &(crbuf[crbuf_used].edge[1]),
	       &(crbuf[crbuf_used].edge[2]),
	       &(crbuf[crbuf_used].edge[3])) != 4) {

      printf("pd_read_KnotTheory: Couldn't read crossing %d\n",crbuf_used);

      free(crbuf);
      return NULL;

    }

    if (crbuf_used == crbuf_size-10) {

      printf("pd_read_KnotTheory: Can't read a single PD code with more than %d crossings.\n",crbuf_size-10);
      free(crbuf);
      return NULL;

    }

  }

  /* We have now read all the crossings and won't bother with the file anymore. Our job is
     to extract the information needed to allocate and fill a pd_code_t from this buffer,
     then free it. The first task is to allocate space for the code. */

  pd_code_t *pd = pd_code_new(crbuf_used+2);
  assert(pd != NULL);

  pd_idx_t cr;

  /* Now we're going to scan through and set crossings AND EDGES. */

  for(cr=0;cr<crbuf_used;cr++) {

    pd_idx_t j;
    for(j=0;j<4;j++) {

      pd->cross[cr].edge[j] = crbuf[cr].edge[j]-1;

      if (j==0) { /* We're definitely going IN */

	pd->edge[pd->cross[cr].edge[j]].head = cr;
	pd->edge[pd->cross[cr].edge[j]].headpos = 0;

      } else if (j==2) { /* We're definitely going OUT */

	pd->edge[pd->cross[cr].edge[j]].tail = cr;
	pd->edge[pd->cross[cr].edge[j]].tailpos = 2;

      }

      /* In positions 1 and 3 we don't know yet, and we'll have to
	 use more information to figure out edge orientation. */

    }

    /* Determining sign is a lot trickier. Recall that the
       convention used to determine sign is this:

            ^
            |
       ----------->
            |
            |
	    0

     positive crossing
     (upper tangent vector) x (lower tangent vector) points OUT of screen.
     This is true when the edges in positions 1 and 3 are in the order 3 -> 1.

            ^
            |
      0-----|----->
            |
            |

     negative crossing
     (upper tangent vector) x (lower tangent vector) points INTO screen.
     This is true when the edges in positions 1 and 3 are in the order 1 -> 3.

     The problem (and it's a large one) is that the edge numbers wrap
     around at the end of the component.

     So if pd->cross[cr].edge[3] and pd->cross[cr].edge[1] differ by
     1, we know which way the strand is going -- it's going in the
     "up" direction.

     If they differ by more than one, the strand indices are wrapping
     around, and the strand goes in the "down" direction.

     For a 2 edge component, we're going to need a better algorithm:

        If the two-edge component lies OVER the remainder of the link,
	its' orientation is genuinely unspecified by the KnotTheory
	data. In this case, we can assign either orientation to the
	component, but must do so consistently.

	If the two-edge component is Hopf-linked to the rest of the link,
	it's orientation IS specified, but only by the crossing where it
	goes UNDER.

    */

    if (pd->cross[cr].edge[3] - pd->cross[cr].edge[1] == 1) { // 1 -> 3

      pd->cross[cr].sign = PD_NEG_ORIENTATION;

      pd->edge[pd->cross[cr].edge[1]].head = cr;
      pd->edge[pd->cross[cr].edge[1]].headpos = 1;

      pd->edge[pd->cross[cr].edge[3]].tail = cr;
      pd->edge[pd->cross[cr].edge[3]].tailpos = 3;

    } else if (pd->cross[cr].edge[1] - pd->cross[cr].edge[3] == 1) { // 3 -> 1

      pd->cross[cr].sign = PD_POS_ORIENTATION;

      pd->edge[pd->cross[cr].edge[3]].head = cr;
      pd->edge[pd->cross[cr].edge[3]].headpos = 3;

      pd->edge[pd->cross[cr].edge[1]].tail = cr;
      pd->edge[pd->cross[cr].edge[1]].tailpos = 1;

    } else if (pd->cross[cr].edge[1] > pd->cross[cr].edge[3]) { // 1 -> 3, wrapping

      pd->cross[cr].sign = PD_NEG_ORIENTATION;

      pd->edge[pd->cross[cr].edge[1]].head = cr;
      pd->edge[pd->cross[cr].edge[1]].headpos = 1;

      pd->edge[pd->cross[cr].edge[3]].tail = cr;
      pd->edge[pd->cross[cr].edge[3]].tailpos = 3;

    } else { // 3 -> 1, wrapping

      pd->cross[cr].sign = PD_POS_ORIENTATION;

      pd->edge[pd->cross[cr].edge[3]].head = cr;
      pd->edge[pd->cross[cr].edge[3]].headpos = 3;

      pd->edge[pd->cross[cr].edge[1]].tail = cr;
      pd->edge[pd->cross[cr].edge[1]].tailpos = 1;


    }

  }

  pd->ncross = crbuf_used;
  pd->nedges = 2*crbuf_used;
  pd->nfaces = crbuf_used + 2;

  /* We now search for two-edge components in the list. Unfortunately,
     none of the existing "regenerate" codes are going to do the trick here,
     because they all resort the crossings, losing the information we'll
     need to correct our past errors (if any).
  */

  pd_idx_t crA, crB;

  for(crA=0;crA<pd->ncross-1;crA++) {

    for(crB=crA+1;crB<pd->ncross;crB++) {


      /*       +------<--e2-----+
	       |                |
	       |0               |2
	   1   |    3        3  |   1
	+-------------->---------------+
	 crA   |                |   crB
	       |2               |0
	       |                |
	       +------->--e1----+

      */

      if (pd->cross[crB].edge[2] == pd->cross[crA].edge[0] &&
	  pd->cross[crB].edge[0] == pd->cross[crA].edge[2] &&
	  pd->cross[crB].edge[3] == pd->cross[crA].edge[3]) {

	pd->cross[crB].sign = PD_POS_ORIENTATION;
	pd->cross[crA].sign = PD_NEG_ORIENTATION;

	pd_idx_t e1, e2;

	e1 = pd->cross[crB].edge[0];
	e2 = pd->cross[crB].edge[2];

	if (pd->edge[e1].head != crB) { pd_reorient_edge(pd,e1,PD_NEG_ORIENTATION); }
	if (pd->edge[e2].head != crA) { pd_reorient_edge(pd,e2,PD_NEG_ORIENTATION); }

      }


      /*       +------->-e2-----+
	       |                |
	       |2               |0
	   3   |    1        1  |   3
	+-------------->---------------+
	 crA   |                |   crB
	       |0               |2
	       |                |
	       +-------<-e1-----+

      */

      if (pd->cross[crB].edge[2] == pd->cross[crA].edge[0] &&
	  pd->cross[crB].edge[0] == pd->cross[crA].edge[2] &&
	  pd->cross[crB].edge[1] == pd->cross[crA].edge[1]) {

	pd->cross[crB].sign = PD_NEG_ORIENTATION;
	pd->cross[crA].sign = PD_POS_ORIENTATION;

	pd_idx_t e1, e2;

	e1 = pd->cross[crB].edge[2];
	e2 = pd->cross[crB].edge[0];

	if (pd->edge[e1].head != crA) { pd_reorient_edge(pd,e1,PD_NEG_ORIENTATION); }
	if (pd->edge[e2].head != crB) { pd_reorient_edge(pd,e2,PD_NEG_ORIENTATION); }

      }

      /*       +------->-e2-----+
	       |                |
	       |3               |3
	   0   |    2        0  |   2
	+------|------->--------|------+
	 crA   |                |   crB
	       |1               |1
	       |                |
	       +-------<--e1----+

	   taken as default, instead of

	       +----------------+
	       |                |
	       |3               |3
	   0   |    2        0  |   2
	+------|------->--------|------+
	 crA   |                |   crB
	       |1               |1
	       |                |
	       +------->--------+

	   which has the same pd code.
      */

      if (pd->cross[crB].edge[3] == pd->cross[crA].edge[1] &&
	  pd->cross[crB].edge[1] == pd->cross[crA].edge[3] &&
	  pd->cross[crB].edge[0] == pd->cross[crA].edge[2]) {

	pd->cross[crB].sign = PD_POS_ORIENTATION;
	pd->cross[crA].sign = PD_NEG_ORIENTATION;

	pd_idx_t e1, e2;

	e1 = pd->cross[crB].edge[1];
	e2 = pd->cross[crB].edge[3];

	if (pd->edge[e1].head != crA) { pd_reorient_edge(pd,e1,PD_NEG_ORIENTATION); }
	if (pd->edge[e2].head != crB) { pd_reorient_edge(pd,e2,PD_NEG_ORIENTATION); }

      }

      /*       +-------->-e2----+
	       |                |
	       |1               |1
	   2   |    0        2  |   0
	+------|-------<--------|------+
	 crA   |                |   crB
	       |3               |3
	       |                |
	       +-------<-e1-----+

	   taken as default, instead of

	       +----------------+
	       |                |
	       |1               |1
	   2   |    0        2  |   0
	+------|-------<--------|------+
	 crA   |                |   crB
	       |3               |3
	       |                |
	       +------->--------+

	   which has the same pd code.
      */

      if (pd->cross[crB].edge[3] == pd->cross[crA].edge[1] &&
	  pd->cross[crB].edge[1] == pd->cross[crA].edge[3] &&
	  pd->cross[crB].edge[2] == pd->cross[crA].edge[0]) {

	pd->cross[crB].sign = PD_NEG_ORIENTATION;
	pd->cross[crA].sign = PD_POS_ORIENTATION;

	pd_idx_t e1, e2;

	e1 = pd->cross[crB].edge[3];
	e2 = pd->cross[crB].edge[1];

	if (pd->edge[e1].head != crA) { pd_reorient_edge(pd,e1,PD_NEG_ORIENTATION); }
	if (pd->edge[e2].head != crB) { pd_reorient_edge(pd,e2,PD_NEG_ORIENTATION); }

      }

      /*       +-------<---e2---+
	       |                |
	       |3               |2
	   0   |    2        3  |   1
	+------|------->---------------+
	 crA   |                |   crB
	       |1               |0
	       |                |
	       +------->---e1---+

      */

      if (pd->cross[crB].edge[2] == pd->cross[crA].edge[3] &&
	  pd->cross[crB].edge[0] == pd->cross[crA].edge[1] &&
	  pd->cross[crB].edge[3] == pd->cross[crA].edge[2]) {

	pd->cross[crB].sign = PD_POS_ORIENTATION;
	pd->cross[crA].sign = PD_POS_ORIENTATION;

	pd_idx_t e1, e2;

	e1 = pd->cross[crB].edge[0];
	e2 = pd->cross[crB].edge[2];

	if (pd->edge[e1].head != crB) { pd_reorient_edge(pd,e1,PD_NEG_ORIENTATION); }
	if (pd->edge[e2].head != crA) { pd_reorient_edge(pd,e2,PD_NEG_ORIENTATION); }

      }

      /*       +------<---e2----+
	       |                |
	       |1               |2
	   2   |    0        3  |   1
	+------|-------<---------------+
	 crA   |                |   crB
	       |3               |0
	       |                |
	       +------->--e1----+

      */

      if (pd->cross[crB].edge[2] == pd->cross[crA].edge[1] &&
	  pd->cross[crB].edge[0] == pd->cross[crA].edge[3] &&
	  pd->cross[crB].edge[3] == pd->cross[crA].edge[0]) {

	pd->cross[crB].sign = PD_NEG_ORIENTATION;
	pd->cross[crA].sign = PD_NEG_ORIENTATION;

	pd_idx_t e1, e2;

	e1 = pd->cross[crB].edge[0];
	e2 = pd->cross[crB].edge[2];

	if (pd->edge[e1].head != crB) { pd_reorient_edge(pd,e1,PD_NEG_ORIENTATION); }
	if (pd->edge[e2].head != crA) { pd_reorient_edge(pd,e2,PD_NEG_ORIENTATION); }

      }

      /*       +------->---e2---+
	       |                |
	       |1               |0
	   2   |    0        1  |   3
	+------|-------<---------------+
	 crA   |                |   crB
	       |3               |2
	       |                |
	       +-------<---e1---+

      */

      if (pd->cross[crB].edge[0] == pd->cross[crA].edge[1] &&
	  pd->cross[crB].edge[3] == pd->cross[crA].edge[2] &&
	  pd->cross[crB].edge[1] == pd->cross[crA].edge[0]) {

	pd->cross[crB].sign = PD_POS_ORIENTATION;
	pd->cross[crA].sign = PD_POS_ORIENTATION;

	pd_idx_t e1, e2;

	e1 = pd->cross[crB].edge[2];
	e2 = pd->cross[crB].edge[0];

	if (pd->edge[e1].head != crA) { pd_reorient_edge(pd,e1,PD_NEG_ORIENTATION); }
	if (pd->edge[e2].head != crB) { pd_reorient_edge(pd,e2,PD_NEG_ORIENTATION); }

      }

      /*       +------>---e2----+
	       |                |
	       |3               |0
	   0   |    2        1  |   3
	+------|------->---------------+
	 crA   |                |   crB
	       |1               |2
	       |                |
	       +-------<--e1----+

      */

      if (pd->cross[crB].edge[0] == pd->cross[crA].edge[3] &&
	  pd->cross[crB].edge[2] == pd->cross[crA].edge[1] &&
	  pd->cross[crB].edge[1] == pd->cross[crA].edge[2]) {

	pd->cross[crB].sign = PD_NEG_ORIENTATION;
	pd->cross[crA].sign = PD_NEG_ORIENTATION;

	pd_idx_t e1, e2;

	e1 = pd->cross[crB].edge[2];
	e2 = pd->cross[crB].edge[0];

	if (pd->edge[e1].head != crA) { pd_reorient_edge(pd,e1,PD_NEG_ORIENTATION); }
	if (pd->edge[e2].head != crB) { pd_reorient_edge(pd,e2,PD_NEG_ORIENTATION); }

      }


      /*       +------<-e2------+
	       |                |
	       |0               |3
	   1   |    3        0  |   2
	+-------------->--------|------+
	 crA   |                |   crB
	       |2               |1
	       |                |
	       +------->-e1-----+

      */

      if (pd->cross[crB].edge[3] == pd->cross[crA].edge[0] &&
	  pd->cross[crB].edge[1] == pd->cross[crA].edge[2] &&
	  pd->cross[crB].edge[0] == pd->cross[crA].edge[3]) {

	pd->cross[crB].sign = PD_NEG_ORIENTATION;
	pd->cross[crA].sign = PD_NEG_ORIENTATION;

	pd_idx_t e1, e2;

	e1 = pd->cross[crB].edge[1];
	e2 = pd->cross[crB].edge[3];

	if (pd->edge[e1].head != crB) { pd_reorient_edge(pd,e1,PD_NEG_ORIENTATION); }
	if (pd->edge[e2].head != crA) { pd_reorient_edge(pd,e2,PD_NEG_ORIENTATION); }

      }

      /*       +-------<--e2----+
	       |                |
	       |0               |1
	   1   |    3        2  |   0
	+--------------<--------|------+
	 crA   |                |   crB
	       |2               |3
	       |                |
	       +------->--e1----+

      */

      if (pd->cross[crB].edge[1] == pd->cross[crA].edge[0] &&
	  pd->cross[crB].edge[3] == pd->cross[crA].edge[2] &&
	  pd->cross[crB].edge[2] == pd->cross[crA].edge[3]) {

	pd->cross[crB].sign = PD_POS_ORIENTATION;
	pd->cross[crA].sign = PD_POS_ORIENTATION;

	pd_idx_t e1, e2;

	e1 = pd->cross[crB].edge[3];
	e2 = pd->cross[crB].edge[1];

	if (pd->edge[e1].head != crB) { pd_reorient_edge(pd,e1,PD_NEG_ORIENTATION); }
	if (pd->edge[e2].head != crA) { pd_reorient_edge(pd,e2,PD_NEG_ORIENTATION); }


      }

      /*       +------->--e2----+
	       |                |
	       |0               |1
	   1   |    3        2  |   0
	+--------------<--------|------+
	 crA   |                |   crB
	       |2               |3
	       |                |
	       +-------<--e1----+

      */

      if (pd->cross[crB].edge[1] == pd->cross[crA].edge[0] &&
	  pd->cross[crB].edge[3] == pd->cross[crA].edge[2] &&
	  pd->cross[crB].edge[2] == pd->cross[crA].edge[3]) {

	pd->cross[crB].sign = PD_NEG_ORIENTATION;
	pd->cross[crA].sign = PD_NEG_ORIENTATION;

	pd_idx_t e1, e2;

	e1 = pd->cross[crB].edge[3];
	e2 = pd->cross[crB].edge[1];

	if (pd->edge[e1].head != crA) { pd_reorient_edge(pd,e1,PD_NEG_ORIENTATION); }
	if (pd->edge[e2].head != crB) { pd_reorient_edge(pd,e2,PD_NEG_ORIENTATION); }


      }

      /*       +------->--e2----+
	       |                |
	       |0               |3
	   1   |    3        0  |   2
	+-------------->--------|------+
	 crA   |                |   crB
	       |2               |1
	       |                |
	       +-------<--e1----+

      */

      if (pd->cross[crB].edge[3] == pd->cross[crA].edge[0] &&
	  pd->cross[crB].edge[1] == pd->cross[crA].edge[2] &&
	  pd->cross[crB].edge[0] == pd->cross[crA].edge[3]) {

	pd->cross[crB].sign = PD_POS_ORIENTATION;
	pd->cross[crA].sign = PD_POS_ORIENTATION;

	pd_idx_t e1, e2;

	e1 = pd->cross[crB].edge[1];
	e2 = pd->cross[crB].edge[3];

	if (pd->edge[e1].head != crA) { pd_reorient_edge(pd,e1,PD_NEG_ORIENTATION); }
	if (pd->edge[e2].head != crB) { pd_reorient_edge(pd,e2,PD_NEG_ORIENTATION); }

      }

    }

  }


  /* We can now detect everything else by regenerating... */

  pd_regenerate(pd);

  /* We double check that we got crossing signs... */

  for(cr=0;cr<pd->ncross;cr++) {

    if (pd->cross[cr].sign == PD_UNSET_ORIENTATION) {

      pd_error(SRCLOC,
	       "read_KnotTheory unable to maintain crossing signs across pd_regenerate call\n"
	       "there must be a problem with edges generated by pd_read_KnotTheory.",pd);
      exit(1);

    }

  }

  /* We're done with the crossing buffer. */
  free(crbuf);

  /* We now make a final check, then terminate... */

  if (!pd_ok(pd)) {

    printf("pd_read_KnotTheory: Regenerated PD does not pass pd_ok. Suspect data is corrupt.\n");
    return NULL;

  }

  return pd;

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

char *pd_print_idx(pd_idx_t idx)
/* Returns the index, either sprintf'd to an unsigned type,
   or the string "PD_UNSET_IDX" if this is equal to PD_UNSET_IDX.
   It's the user's responsibility to dispose of the char buffer,
   unfortunately. */
{
  char *out;
  out = calloc(64,sizeof(char));
  assert(out != NULL);

  if (idx != PD_UNSET_IDX) {
    sprintf(out,"%u",(unsigned int)(idx));
  } else {
    sprintf(out,"PD_UNSET_IDX");
  }
  return out;
}

char *pd_print_boundary_or(pd_boundary_or_t or)
/* Returns a character string containing "in", "out", or "?".
   It's the user's responsibility to dispose of the buffer. */
{
  char *outbuf;
  outbuf = calloc(32,sizeof(char));
  assert(outbuf != NULL);

  if (or == in) {
    sprintf(outbuf," in");
  } else if (or == out) {
    sprintf(outbuf,"out");
  } else {
    sprintf(outbuf,"  ?");
  }

  return outbuf;
}


char *pd_print_pos(pd_idx_t pos)
/* Returns the index, either sprintf'd to an unsigned type,
   or the string "PD_UNSET_POS" if this is equal to PD_UNSET_POS.
   It's the user's responsibility to dispose of the char buffer,
   unfortunately. */
{
  char *out;
  out = calloc(32,sizeof(char));
  assert(out != NULL);

  if (pos != PD_UNSET_POS) {
    sprintf(out,"%1u",(unsigned int)(pos));
  } else {
    sprintf(out,"PD_UNSET_POS");
  }
  return out;
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
   %COMP      comp number        compnum (e1 -> e2 -> e3 -> ..... -> e1 (or1)) tag (tag)
   %PD        (no argument)      (\n\n output of pd_write \n\n)

   We also have conversions for pointers.

   %OR        *pd_or_t           +, -, U (unset), or ? (anything else)
   %MULTIDX   *pd_multidx_t      multidx (i[0] i[1] ... i[n-1])
   %COMPGRP   *pd_compgrp_t      compgrp (comp[0] .. comp[n-1])
   %ISO       *pd_iso_t          <a representation of the isomorphism type>
   %PERM      *pd_perm_t         <a representation of the permutation type>
   %DIHEDRAL  *pd_dihedral_t     <a representation of an element of the dihedral group>
   %CYCLIC    *pd_cyclic_t       <a representation of an element of the cyclic group>

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

    if (!strncmp(nxtconv,"%FACEMAP",8)) { /* Face map */

      pd_facemap_t *facemap = (pd_facemap_t *) va_arg(ap,void *);
      char *printed;

      printed = pd_print_facemap(facemap);
      fprintf(stream,"%s",printed);
      free(printed);

      nxtconv += 8;

    } else if (!strncmp(nxtconv,"%FACE",5)) { /* %FACE conversion */

      if (pd != NULL) {

	pd_idx_t face = (pd_idx_t) va_arg(ap,int);
	pd_idx_t edge;

	fprintf(stream,"face %d (",face);
	for(edge=0;edge<pd->face[face].nedges-1 && edge<pd->MAXEDGES;edge++) {

	  char *edge_idx = pd_print_idx(pd->face[face].edge[edge]);
	  fprintf(stream," (%c) %s ->",
		  pd->face[face].or[edge] == PD_POS_ORIENTATION ? '+':'-',
		  edge_idx);
	  free(edge_idx);

	}

	char *edge_idx = pd_print_idx(pd->face[face].edge[edge]);
	fprintf(stream,"(%c) %s) ",
		pd->face[face].or[edge] == PD_POS_ORIENTATION ? '+':'-',
		edge_idx);
	free(edge_idx);

	nxtconv += 5;

      } else {

	fprintf(stderr,"pd_printf: Can't print %%FACE primitive without a pd.\n");
	exit(1);

      }

    } else if (!strncmp(nxtconv,"%COMPGRP",8)) { /* Component group */

      pd_idx_t comp;
      pd_compgrp_t *grp = (pd_compgrp_t *) va_arg(ap,void *);

      fprintf(stream,"compgrp ( ");
      for(comp=0;comp<grp->ncomps;comp++) { fprintf(stream,"%d ",grp->comp[comp]); }
      fprintf(stream,")");

      nxtconv += 8;

    } else if (!strncmp(nxtconv,"%EDGEMAP",8)) { /* Edge map */

      pd_edgemap_t *edgemap = (pd_edgemap_t *) va_arg(ap,void *);
      char *printed_form;

      printed_form = pd_print_edgemap(edgemap);
      fprintf(stream,"%s",printed_form);
      free(printed_form);

      nxtconv += 8;

    } else if (!strncmp(nxtconv,"%EDGE",5)) { /* %EDGE conversion */

      if (pd != NULL) {

	pd_idx_t edge = (pd_idx_t) va_arg(ap,int);

	char *ed = pd_print_idx(edge);
	char *head, *headpos, *tail, *tailpos;
	tail = pd_print_idx(pd->edge[edge].tail);
	tailpos = pd_print_pos(pd->edge[edge].tailpos);
	head = pd_print_idx(pd->edge[edge].head);
	headpos = pd_print_pos(pd->edge[edge].headpos);

	fprintf(stream,"edge %s (%s,%s -> %s,%s)",ed,
		tail,tailpos,head,headpos);

	free(ed);
	free(head); free(headpos); free(tail); free(tailpos);
	nxtconv += 5;

      } else {

	fprintf(stderr,"pd_printf: Can't print %%EDGE primitive without a pd.\n");
	exit(1);

      }

    } else if (!strncmp(nxtconv,"%COMP",5)) { /* %COMP conversion */

      if (pd != NULL) {

	pd_idx_t comp = (pd_idx_t) va_arg(ap,int);
	pd_idx_t edge;

	fprintf(stream,"comp %d (",comp);
	for(edge=0;edge<pd->comp[comp].nedges-1 && edge<pd->MAXEDGES;edge++) {

	  char *ed;
	  ed = pd_print_idx(pd->comp[comp].edge[edge]);
	  fprintf(stream," %s ->",ed);
	  free(ed);

	}

	char *ed;
	ed = pd_print_idx(pd->comp[comp].edge[edge]);
	fprintf(stream," %s ) ",ed);
	free(ed);

	fprintf(stream," tag %c",pd->comp[comp].tag);

	nxtconv += 5;

      } else {

	fprintf(stderr,"pd_printf: Can't print %%COMP primitive without a pd.\n");
	exit(1);

      }

    } else if (!strncmp(nxtconv,"%TANGLE",7)) { /* %TANGLE conversion */

      pd_tangle_t *t = (pd_tangle_t *) va_arg(ap,void *);
      pd_idx_t edge,face;

      fprintf(stream,"tangle (\n\tedge loop: ");
      for(edge=0;edge<t->nedges;edge++) {

	char *ed,*bdo;
	ed = pd_print_idx(t->edge[edge]);
	bdo = pd_print_boundary_or(t->edge_bdy_or[edge]);
	fprintf(stream," %s (%s) ->",ed,bdo);
	free(ed);
	free(bdo);

      }

      fprintf(stream,"\n\tface loop: ");
      for(face=0;face<t->nedges;face++) {

	char *ed;
	ed = pd_print_idx(t->face[face]);
	fprintf(stream," %s ->",ed);
	free(ed);

      }

      pd_idx_t strand;
      for(strand=0;strand<t->nstrands;strand++) {

	char *se,*ee,*ne,*cp,*spe,*epe;

	se = pd_print_idx(t->strand[strand].start_edge);
	ee = pd_print_idx(t->strand[strand].end_edge);
	ne = pd_print_idx(t->strand[strand].nedges);
	cp = pd_print_idx(t->strand[strand].comp);
	spe = pd_print_idx(t->edge[t->strand[strand].start_edge]);
	epe = pd_print_idx(t->edge[t->strand[strand].end_edge]);

	fprintf(stream,"\n\tstrand: tangle edge %s (pd edge %s) -> tangle edge %s (pd edge %s) (%s edges, comp %s)",se,spe,ee,epe,ne,cp);

	free(se); free(ee); free(ne); free(cp); free(spe); free(epe);

      }

      char *nic;
      nic = pd_print_idx(t->ninterior_cross);
      fprintf(stream,"\n\t%s interior cross: ",nic);
      free(nic);

      pd_idx_t cr;
      for(cr=0;cr<t->ninterior_cross;cr++) {

	char *ic;
	ic = pd_print_idx(t->interior_cross[cr]);

	fprintf(stream,"%s ",ic);
	free(ic);

      }

      char *nie;
      nie = pd_print_idx(t->ninterior_edges);
      fprintf(stream,"\n\t%s interior edges: ",nie);
      free(nie);

      for(edge=0;edge<t->ninterior_edges;edge++) {

	char *ie;
	ie = pd_print_idx(t->interior_edge[edge]);

	fprintf(stream,"%s ",ie);
	free(ie);

      }

      fprintf(stream,")\n");
      nxtconv += 7;

    } else if (!strncmp(nxtconv,"%CROSSPTR",9)) { /* %CROSSPTR conversion */

      pd_crossing_t *cross = (pd_crossing_t *) va_arg(ap,void *);

      fprintf(stream,"cross (%d %d %d %d) %c",
	      cross->edge[0],
	      cross->edge[1],
	      cross->edge[2],
	      cross->edge[3],
	      pd_print_or(cross->sign));

	nxtconv += 8;

    }  else if (!strncmp(nxtconv,"%CROSSMAP",9)) { /* Crossing map */

      pd_crossmap_t *crmap = (pd_crossmap_t *) va_arg(ap,void *);
      char *printed;

      printed = pd_print_crossmap(crmap);
      fprintf(stream,"%s",printed);
      free(printed);

      nxtconv += 9;

    } else if (!strncmp(nxtconv,"%CROSS",6)) { /* %CROSS conversion */

      if (pd != NULL) {

	pd_idx_t cross = (pd_idx_t) va_arg(ap,int);

	char *cr,*e[4];
	cr = pd_print_idx(cross);
	int i;
	for(i=0;i<4;i++) { e[i] = pd_print_idx(pd->cross[cross].edge[i]); }

	fprintf(stream,"cross %s (%s %s %s %s) %c",cr,e[0],e[1],e[2],e[3],
		pd_print_or(pd->cross[cross].sign));

	for(i=0;i<4;i++) { free(e[i]); }
	free(cr);

	nxtconv += 6;

      } else {

	fprintf(stderr,"pd_printf: Can't print %%CROSS primitive without a pd.\n");
	exit(1);

      }

    } else if (!strncmp(nxtconv,"%PD",3)) { /* %PD conversion */

      if (pd != NULL) {

	fprintf(stream,"\n");
	pd_write(stream,pd);
	fprintf(stream,"\n");

	nxtconv += 3;

      } else {

	fprintf(stderr,"pd_printf: Can't print %%PD primitive without a pd.\n");
	exit(1);

      }

    } else if (!strncmp(nxtconv,"%FEDGE",6)) { /* FACE/EDGE conversion */

      if (pd != NULL) {

	pd_idx_t face = (pd_idx_t) va_arg(ap,int);
	pd_idx_t edge = (pd_idx_t) va_arg(ap,int);

	fprintf(stream,"%d (%c)",pd->face[face].edge[edge],
		pd->face[face].or[edge] == PD_POS_ORIENTATION ? '+':'-');

	nxtconv += 6;

      } else {

	fprintf(stderr,"pd_printf: Can't print %%FEDGE primitive without a pd.\n");
	exit(1);

      }

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

    } else if (!strncmp(nxtconv,"%ORIENTATION",12)) { /* Orientation group element */

      pd_orientation_t *d = (pd_orientation_t *) va_arg(ap,void *);

      fprintf(stream,"multiorientation ");

      char *dprint = pd_print_orientation(d);
      fprintf(stream,"%s",dprint);
      free(dprint);

      nxtconv += 12;

    } else if (!strncmp(nxtconv,"%OR ",3)) { /* Edge map */

      pd_or_t *or = (pd_or_t *) va_arg(ap,void *);
      fprintf(stream,"or (%c)",pd_print_or(*or));

      nxtconv += 3;

    } else if (!strncmp(nxtconv,"%BDY_OR",7)) { /* Boundary orientation */

      pd_boundary_or_t *bdyor = (pd_boundary_or_t *) va_arg(ap,void *);
      char *printed;

      printed = pd_print_boundary_or(*bdyor);
      fprintf(stream,"%s",printed);
      free(printed);

      nxtconv += 7;

    } else if (!strncmp(nxtconv,"%ISO ",5)) { /* pd_code -> pd_code isomorphism */

      pd_iso_t *iso = (pd_iso_t *) va_arg(ap,void *);
      char *printed;

      printed = pd_print_iso(iso);
      fprintf(stream,"%s",printed);
      free(printed);

      nxtconv += 4;

    } else if (!strncmp(nxtconv,"%DIHEDRAL",9)) { /* Dihedral group element */

      pd_dihedral_t *d = (pd_dihedral_t *) va_arg(ap,void *);

      fprintf(stream,"dihedral ");

      char *dprint = pd_print_dihedral(d);
      fprintf(stream,"%s",dprint);
      free(dprint);

      nxtconv += 9;

    } else if (!strncmp(nxtconv,"%CYCLIC",7)) { /* Cyclic group element */

      pd_cyclic_t *c = (pd_cyclic_t *) va_arg(ap,void *);

      fprintf(stream,"cyclic ");

      char *cprint = pd_print_cyclic(c);
      fprintf(stream,"%s",cprint);
      free(cprint);

      nxtconv += 7;

    } else if (!strncmp(nxtconv,"%PERM",5)) { /* Permutation group elt */

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

    } else if (!strncmp(nxtconv,"%s",2)) { /* Standard string conversion */

      char *str = (char *) va_arg(ap,char *);
      fprintf(stream,"%s",str);
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

  if (PD_VERBOSE < 20) { return false; }

  exit(1);
}

void pd_check_cr(char *file, int line, pd_code_t *pd, pd_idx_t cr)

/* Checks to see if cr is a legitimate crossing number for pd. */

{

  if (cr >= pd->ncross) {

    pd_error(file,line,"pd_code %PD has %d crossings, and variable attempted to reference crossing %d.",pd,pd->ncross,cr);

  }

}

void pd_check_edge(char *file, int line, pd_code_t *pd, pd_idx_t edge)

/* Checks to see if cr is a legitimate crossing number for pd. */

{

  if (edge >= pd->nedges) {

    pd_error(file,line,"pd_code %PD has %d edges, and variable attempted to reference edge %d.",pd,pd->nedges,edge);

  }

}

void pd_check_cmp(char *file, int line, pd_code_t *pd, pd_idx_t cmp)

/* Checks to see if cmp is a legitimate component number for pd. */

{

  if (cmp >= pd->ncomps) {

    pd_error(file,line,"pd_code %PD has %d components, and variable attempted to reference component %d.",pd,pd->ncomps,cmp);

  }

}

void pd_check_face(char *file, int line, pd_code_t *pd, pd_idx_t face)

/* Checks to see if face is a legitimate face number for pd. */

{

  if (face >= pd->nfaces) {

    pd_error(file,line,"pd_code %PD has %d faces, and variable attempted to reference face %d.",pd,pd->nfaces,face);

  }

}

void pd_check_notnull(char *file, int line, char *varname, void *ptr)
{

  if (ptr == NULL) {

    pd_error(file,line,"output variable %s is a NULL pointer",NULL,varname);

  }
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
  pd_code_t *pdA = pd_code_new(nA);
  assert(pdA != NULL);

  pdA->ncross = nA;

  pd_code_t *pdB = pd_code_new(nB);
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

  // Store the call to pd_isomorphic
  bool is_iso = pd_isomorphic(pdA,pdB);

  pd_code_free(&pdA); pd_code_free(&pdB);
  return is_iso;
}


pd_idx_t pd_previous_edge(pd_code_t *pd, pd_idx_t edge)

/* Find the previous edge number */

{
  pd_check_edge(SRCLOC,pd,edge);
  pd_idx_t candidate_edge = pd->cross[pd->edge[edge].tail].edge[(pd->edge[edge].tailpos+2)%4];
  assert(pd->edge[candidate_edge].head == pd->edge[edge].tail);
  return candidate_edge;
}

pd_idx_t pd_next_edge(pd_code_t *pd,pd_idx_t edge)

/* Find the next edge number */

{
  pd_check_edge(SRCLOC,pd,edge);
  pd_idx_t candidate_edge = pd->cross[pd->edge[edge].head].edge[(pd->edge[edge].headpos+2)%4];
  assert(pd->edge[candidate_edge].tail == pd->edge[edge].head);
  return candidate_edge;
}
