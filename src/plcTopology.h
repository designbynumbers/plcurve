/* plcTopology.h.  Generated by configure.  */

/* This header file contains definitions and an API to deal with piecewise
   linear curves as topological knots, including functionality to deal with
   knot diagrams and compute HOMFLY polynomials. */

/* Copyright 2004 The University of Georgia */

/* This file is part of plCurve.

   plCurve is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   plCurve is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with plCurve; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

/*@-exportlocal@*/
#ifndef PLCTOPOLOGY_H
#define PLCTOPOLOGY_H

#if (__cplusplus || c_plusplus)
extern "C" {
#endif

  /*
   * Take our chances that stdio.h and stdbool.h exist.  We need stdio to define
   * the FILE type and stdbool to define bool.  If stdbool doesn't exist, one
   * option is to just
   *   typedef int bool;
   * and
   *   #define true 1
   *   #define false 0
   */
#include <stdio.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>  /* We are going to need the gsl_rng type to be defined below. */
#include <plCurve.h>
#include <limits.h>

  /* This is going to take an actual fork to deal with the fact that
     we can't keep the data in static memory anymore. */

/* #define PD_MAXVERTS         1024   /\* We are only going to deal with diagrams with <= 1024 crossings. *\/ */
/* #define PD_MAXEDGES         (int)(PD_MAXVERTS*2 + 1) */
/* #define PD_MAXCOMPONENTS    (int)(PD_MAXVERTS/2) */
/* #define PD_MAXFACES         PD_MAXVERTS+2 */

  #define PD_HASHSIZE 32

  extern int PD_VERBOSE;

  /* For the sake of compactness and speed, everything
     is written on top of three "abstract" integer types,

     1. the pd "index" type, which is used to index faces,
     edges, or crossings in a pd code,

     2. the pd "orientation" type, which is used to orient
     faces or edges as needed, and only takes the values
     PD_POS_ORIENTATION or PD_NEG_ORIENTATION, and

     3. the pd "position" type, which should contain only
     values 0..3 corresponding to positions around a
     crossing.

     We should use these by default in all the pdcode.c
     functions. */

#define PD_POS_ORIENTATION 1
#define PD_NEG_ORIENTATION 0
#define PD_UNSET_ORIENTATION 2

  typedef unsigned int      pd_idx_t ;  /* pd "index" type */
  typedef unsigned char     pd_or_t;    /* pd "orientation" type */
  typedef unsigned int      pd_pos_t;   /* pd "position" type */
  typedef unsigned long int pd_uid_t;   /* pd "uid" type */
  typedef char              pd_tag_t;   /* pd "component tag" type */

  /* We are going to need special values of these to represent
     UNSET indices at times. */

#define PD_UNSET_IDX UINT_MAX
#define PD_UNSET_POS 4
#define PD_UNSET_TAG 0
#define PD_UNSET_UID ULONG_MAX

  /* The basic architecture of the pd_code is kind of intricate.  The
     problem is that we have to keep track of labelled diagram data in
     such a way that we can access things relatively quickly as we
     work and make an efficient search for isomorphisms, but also not
     lose track of component identities as we go. */


/* Error codes for dynamic error handling */
#define pd_err_check(CHECK,ERR,CODE,RET) { \
    if(ERR) { if(!CHECK) {*ERR = CODE; return RET;} } \
    else { assert(CHECK); }       \
    }
#define pd_err_set(ERR,CODE) {                  \
    if(ERR) {*ERR = CODE;} \
    }
#define PD_NO_ERROR 0
#define PD_NOT_OK 1
#define PD_BAD_FORMAT 2
#define PD_EOF 3
#define PD_TANGLE_ERROR 4


  typedef struct pd_edge_struct {

    /* An oriented edge, joining two verts tail -> head */

    pd_idx_t head;
    pd_pos_t headpos;
    /* Pos [0..3] in crossing record of head vertex. */

    pd_idx_t tail;
    pd_pos_t tailpos;
    /* Pos [0..3] in crossing record of tail vertex. */

  } pd_edge_t;

  /* Since we may have loop edges, we need to record
     the position on the crossing where the head and
     tail of the edge come in.

     For a loop edge, the edge index occurs twice in
     the crossing record, so we can't determine
     orientation otherwise. */

  typedef struct pd_component_struct {

    pd_idx_t nedges;
    pd_idx_t *edge;   // Should be allocated/deallocated as needed.

    /* Edge indices in orientation order
       around a component. These are expected
       to be consecutive. */

    pd_tag_t tag;    /* This tag keeps track of the identity
			of a component as we do crossing moves.
			It is a character from the set "A..Z"
			and proceeding in ASCII order from there.

			The set of tags for an n component link
			must be unique elements with ASCII codes
			>= 'A'.*/

  } pd_component_t;

  typedef struct pd_face_struct {

    pd_idx_t    nedges;
    pd_idx_t    *edge;
    pd_or_t     *orient;

    /* Edge indices around the face in counterclockwise order. These
       are NOT consecutive, nor are they always positively oriented
       (according to their intrinsic edge orientation), so we store
       their orientation as well as their edge number. */

  } pd_face_t;

  typedef struct pd_crossing_struct {

    pd_idx_t edge[4]; /* Edge indices around crossing in counterclockwise order */
    pd_or_t  sign;    /* Whether the crossing is positively or negatively oriented */

    /* The convention used to determine sign is this:

            ^
            |
       ----------->
            |
            |

      positive crossing
      (upper tangent vector) x (lower tangent vector) points OUT of screen.

            ^
            |
       -----|----->
            |
            |

      negative crossing
      (upper tangent vector) x (lower tangent vector) points INTO screen.


      You often simply want to know which of the strands (0-2) or (1-3) is
      on top. There are several cases, because it depends BOTH on the sign of the crossing
      and the orientation of the strands. It's recommended to use the function

      pd_overstrand(pd_code_t *pd,pd_idx_t cr,pd_idx_t *incoming_edgenum, pd_idx_t *outgoing_edgenum)

      to determine which is which. This returns the edge number  (that is, a number in 0..pd->nedges)
      of the incoming and outgoing edges of the strand going over at this crossing.

      pd_overstrand_pos(pd_code_t *pd,pd_idx_t cr,pd_idx_t *incoming_edgepos, pd_idx_t *outgoing_edgepos)

      returns the position in this crossing (that is, a number in 0..3) of the incoming and outgoing
      edges of the strand going over at this crossing.

    */

  } pd_crossing_t;

  typedef struct pd_code_struct {

    pd_uid_t uid;
    /* Unique diagram id number (among nverts diagrams WITH THIS HASH). */

    pd_idx_t MAXVERTS;
    /* Space for this number of vertices has been allocated. */

    pd_idx_t MAXEDGES;
    /* Space for this number of edges has been allocated. (Usually pd->MAXVERTS*2 + 1). */

    pd_idx_t MAXCOMPONENTS;
    /* Space for this number of components has been allocated. (Usually pd->MAXVERTS/2.) */

    pd_idx_t MAXFACES;
    /* Space for this number of faces has been allocated. (Usually pd->MAXVERTS+2). */

    pd_idx_t  ncross;
    /* The total number of crossings in the PD-Code */

    pd_idx_t  nedges;
    /* The total number of edges in the PD-Code */

    pd_idx_t  ncomps;
    /* The number of components in the PD-Code */

    pd_idx_t  nfaces;
    /* The number of faces in the PD-Code */

    char hash[PD_HASHSIZE];
    /* A printable 32 char hash string (incl trailing 0) from generate_pd_hash */

    pd_edge_t      *edge;
    /* nedges entries: tail and head vertices. */

    pd_component_t *comp;
    /* ncomps entries: edge indices/orientations in order around comp */

    pd_crossing_t  *cross;
    /* nverts entries: 4 edge indices (ccw around cross)*/

    pd_face_t      *face;
    /* nfaces entries: edge indices/orientations ccw around face */

  } pd_code_t;

  /*
     A pd_code is a representation of a link diagram with
     assignments of orientations and numberings to the
     components of the link, and numberings for the
     crossings, faces, and edges.

     The edges specified in the crossing and face data are
     always indices into the master list pd->edge[X], which
     stores the "complementary" information for each edge
     about which pair of crossings it connects.

     There are some special case pd codes.

     1) A diagram of an unknot with no crossings is denoted by
        a pd_code with ncross = 0, nedges = 1, 1 edge with
	everything (head, tail, headpos, tailpos) set to
	PD_UNSET_IDX or PD_UNSET_POS, 2 faces (each with
	1 edge), and 1 component.

     2) Split diagrams cannot be encoded by pd_code_t. The
        recommended way to handle split diagrams is to have
	an array of pd_code_t to cover the various connected
	components of the diagram.
  */

  /* A pdCode is like a plCurve... we always deal with a pointer
     which is supposed to be created by the constructor function
     pd_code_new.

     The tricky bit is that the actual list of edges along each face
     is not known in advance (and we don't want to allocate the memory
     for each face data type in advance, either). So we DO have an
     array of edges and faces (though we might not use them all), but
     we DON'T allocate the face records themselves.
  */

  pd_code_t *pd_code_new(pd_idx_t MAXVERTS);
  void       pd_code_free(pd_code_t **pd);
  void       pd_code_eltfree(void **PD);  /* Used when we make a pd_container of pd codes */

  /* Utility Functions For Dealing With PD-code primitives */

  int  pd_idx_cmp(const void *A, const void *B);
  /* The usual comparison function for sorting and searching */

  char pd_print_or(pd_or_t orient);
  /* Returns a one-character printed form for "or": +, -, U (unset), or ? (anything else) */

  pd_or_t pd_compose_or(pd_or_t a,pd_or_t b);
  /* Returns the composition of the two orientation changes: ++ = -- = +, +- = -+ = - */

  bool pd_or_ok(pd_or_t orient); /* Check whether or has a legal value. */

  int pd_or_cmp(const void *A,const void *B);
  /* Compare two *pd_or_t. */

  pd_crossing_t pd_build_cross(pd_idx_t e0,pd_idx_t e1,pd_idx_t e2,pd_idx_t e3);
  /* Builds a crossing from the given edge indices */

  void pd_canonorder_cross(pd_crossing_t *cr, pd_or_t orient);
  /* Reverses (if or == PD_NEG_ORIENTATION) and then rotates
     a crossing into canonical order: cr->edge[0] is the
     lowest edge # */

  void pd_canonorder_face(pd_face_t *face, pd_or_t orient);
  /* Reverses (if or == PD_NEG_ORIENTATION) and then rotates
     a face into canonical order: face->edge[0] is the
     lowest edge # */

  int  pd_cross_cmp(const void *A, const void *B);
  /* When passed two pd_crossing_t, compare order.
     (used for searching, sorting) */

  int  pd_face_cmp(const void *A, const void *B);
  /* When passed two pd_face_t, sort by number of components, then dictionary order
     on edge numbers. Used for searching, sorting. */

  int  pd_component_cmp(const void *A, const void *B);
  /* Sort by number of edges, then dictionary order on the edge numbers.
     Used for searching, sorting. */

  void pd_component_and_pos(pd_code_t *pd,pd_idx_t edge,
			    pd_idx_t *comp,pd_idx_t *comp_pos);
  /* Returns component number and position on component of a
     given edge number (or die) */

  pd_idx_t pd_previous_edge(pd_code_t *pd, pd_idx_t edge);
  pd_idx_t pd_next_edge(pd_code_t *pd,pd_idx_t edge);

  /* Finds the edge number of the previous or next edge from edge.
     according to the orientation of edge. */

  void pd_face_and_pos(pd_code_t *pd, pd_idx_t edge,
		       pd_idx_t *posface, pd_idx_t *posface_pos,
		       pd_idx_t *negface, pd_idx_t *negface_pos);

  /* Finds the two faces which edge occurs on, which should include
     one face where edge appears in positive orientation (posface)
     and one where edge appears with negative orientation (negface).
     If we don't find two with this description, die.

     Return the position (index) of the edge on each face. */

  bool pd_edge_on_face(pd_code_t *pd, pd_idx_t edge, pd_idx_t face);
  /* Returns true if the edge is on the face (with either sign). */

  pd_edge_t pd_oriented_edge(pd_edge_t e,pd_or_t orient);
  /* Returns original edge if orient = PD_POS_ORIENTATION,
     reversed edge if orient = PD_NEG_ORIENTATION */

  void pd_reorient_edge(pd_code_t *pd,pd_idx_t edge,pd_or_t orient);
  /* Flips the edge in pd->edges[] if orient == PD_NEG_ORIENTATION */

  /* Over and under information at a crossing. */

  void pd_overstrand(pd_code_t *pd,pd_idx_t cr,pd_idx_t *incoming_edgenum, pd_idx_t *outgoing_edgenum);

  /*  Returns the edge number  (that is, a number in 0..pd->nedges)
      of the incoming and outgoing edges of the strand going OVER at crossing cr of pd,
      using the sign of the crossing to determine. */

  void pd_understrand(pd_code_t *pd,pd_idx_t cr,pd_idx_t *incoming_edgenum, pd_idx_t *outgoing_edgenum);

  /*  Returns the edge number (that is, a number in 0..pd->nedges) of
      the incoming and outgoing edges of the strand going UNDER at
      crossing cr of pd, using the sign of the crossing to
      determine. */

  void pd_overstrand_pos(pd_code_t *pd,pd_idx_t cr,pd_idx_t *incoming_edgepos, pd_idx_t *outgoing_edgepos);

  /* Returns the position in crossing cr of pd (that is, a number in
     0..3) of the incoming and outgoing edges of the strand going over
     at this crossing, using the crossing sign data and edge
     orientations in order to compute the answer. */

  void pd_understrand_pos(pd_code_t *pd,pd_idx_t cr,pd_idx_t *incoming_edgepos, pd_idx_t *outgoing_edgepos);

  /* Returns the position in crossing cr of pd (that is, a number in
     0..3) of the incoming and outgoing edges of the strand going
     under at this crossing, using the crossing sign data and edge
     orientations in order to compute the answer. */

  /* Functions to compute data from crossing (and other) information */

  void pd_regenerate_crossings(pd_code_t *pd);
  /* Reorders the crossings cyclically to put the
     lowest index edge first and sorts crossings
     in dictionary order based on this reordering. */

  void pd_regenerate_edges(pd_code_t *pd);
  /* Uses crossing information to create edge records
     with random numberings and orientations, then walks
     along components to make sure the edges have consistent
     orientations. */

  void pd_regenerate_comps(pd_code_t *pd);
  /* Strings edges into components, which are then given
     new component tags. Renumbers the edges so that they
     appear consecutively along components, fixing any
     orientations needed along the way. Updates the edge
     references in the crossing data to the new numbering
     scheme. */

  void pd_regenerate_faces(pd_code_t *pd);
  /* Fills in faces from crossing, edge
     information, including orientation of
     each edge along the face. */

  void pd_regenerate_hash(pd_code_t *pd);
  /* Fill in (printable) hash value from
     comps, faces, and crossings */

  void pd_regenerate(pd_code_t *pd);
  /* Regenerates everything from cross data. This needs to be
     used with some care if the crossings have signs. The crossing
     signs depend on edge orientations, which are not uniquely specified
     by crossing data.

     Therefore, pd_regenerate will try to preserve an existing edge set
     if it passes pd_edges_ok (and if it can, will preserve crossing sign
     data attached to the crossings).

     If there is no edge set, or the edge set is inconsistent, pd_regenerate
     will rebuild the edge set (if possible), but will discard the crossing
     sign information as it has no way to transfer crossing sign information
     without knowing the intended orientation of the edges.

     This means that if you're writing code to generate diagrams, you either
     need to

     a) build the edges, too (but not the faces/components/hash)

     or

     b) assign crossing SIGNS after the pd_regenerate call
  */

  /* pd sanity checking */

  /* The ``cross'' data in a pd-code is the
     fundamental object. Everything else is
     computable from this. So these funcs
     check the derived data against
     the cross data.

     If PD_VERBOSE > 10 they terminate on
     error with a helpful message. Otherwise
     they just return true/false.
  */

  bool pd_cross_ok(pd_code_t *pd);
  bool pd_edges_ok(pd_code_t *pd);
  bool pd_faces_ok(pd_code_t *pd);
  bool pd_comps_ok(pd_code_t *pd);
  bool pd_ok(pd_code_t *pd);

  /* entire pd operations */

  /* These functions read and write an internal text format for pd codes: */

  void       pd_write_KnotTheory(FILE *of, pd_code_t *pd);
  void       pd_write(FILE *outfile,pd_code_t *pd);
  void       pd_write_c(FILE *outfile, pd_code_t *pd, char *name);
  /* Writes a c procedure which recreates the pd code pd.
     The procedure will be called pd_create_name() and take
     no arguments. */

  pd_code_t *pd_read(FILE *infile); /* Returns NULL if the file is corrupt */
  pd_code_t *pd_read_err(FILE *infile, int *err); // Same as above, but with err checking

  /* This function reads a pdcode from the Mathematica package KnotTheory,
     exported as text with something like:

     Export["7_2.txt",PD[Knot[7,2]]]

     These PD codes don't have component or face information, so that is
     all regenerated once the crossings have been loaded from the file.
     This will only read one PD code per file.
  */

  pd_code_t *pd_read_KnotTheory(FILE *infile); /* Returns NULL if the file is corrupt */


  bool pd_diagram_isotopic(pd_code_t *A, pd_code_t *B);

  /* Detect whether two pd_codes correspond to diagrams of labelled,
     oriented components related by an isotopy of the 2-sphere or an
     isotopy of the 2-sphere composed with a reflection ("reshaping
     the diagram" or "turning the diagram inside out").

     1. Corresponding crossings are required to have the same sign
     (positive, negative, or unset).
     2. Corresponding components are required to have the same tag.
     3. Corresponding edges are required to have the same orientation.

     In particular, this means that the orientations of each component
     of the pd-codes have to be the same.

     This is the strongest kind of diagram equivalence.
  */

    bool pd_map_isomorphic(pd_code_t *pdA,pd_code_t *pdB);

  bool pd_isomorphic(pd_code_t *pdA,pd_code_t *pdB);
  /* Detect whether two pd codes are correspond to the same polyhedral
     decomposition of the 3-sphere (or are mirror images of each
     other). This is much weaker than being diagram-isotopic;

     1. Crossing signs are ignored.
     2. Component tags are ignored.
     3. The orientations of all edges in each component are either
     (all) preserved or (all) reversed.  However, some component
     orientations may be preserved while others are reversed.

     This is the weakest kind of diagram equivalence.
  */

  bool pd_isomorphic_strings(char *pdcodeA, int nA, char*pdcodeB, int nB);
  /* This method will generate pd code objects and return the call
     to pd_isomorphic. The purpose is to call this from Python via
     SWIG. */

  /* PD storage. We define a flexible container type. It can be
     queried by isomorphism type and inserted in almost linear time,
     and comfortably holds millions of pd_code_t objects in memory.

     The container holds pd_code_t's up to an equivalence relation,
     and can be searched with respect to this relation. */

  typedef struct pdstorage_struct pd_stor_t;
  typedef enum {NONE, ISOMORPHISM, DIAGRAM_ISOTOPY} pd_equivalence_t;

  pd_stor_t *pd_new_pdstor(void);
  /* Create new (empty) pdstor */

  void pd_free_pdstor(pd_stor_t **pdstor);
  /* Delete pdstor (and all associated memory) */

  void pd_addto_pdstor(pd_stor_t *pdstor, pd_code_t *pd,pd_equivalence_t eq);
  /* Add a new-memory copy of pd to pdstor, unless an equivalent pd_code_t
     is already stored. Note: If eq == NONE, then we always store, even if
     this is an exact duplicate of a previous entry. 

     The uid of pd is

     1) left alone if the pd code was inserted successfully
     2) reset to the UID of the equivalent diagram in the pdstor if a 
        duplicate was found in storage. 

  */

  pd_code_t *pd_stor_firstelt(pd_stor_t *pdstor);

  /* Finds the first element of pdstor, and initializes the internal
     state of pdstor. Note that this internal state will go stale if an
     insert or delete operation is performed, and so it's reset if we do
     an insert. This returns a new-memory copy of the element in the
     pdstor, which is the caller's responsibility to free.*/

  pd_code_t *pd_stor_nextelt(pd_stor_t *pdstor);
  /* Can be repeatedly called until it returns NULL to iterate over a
     pdstor. As with pdstor_firstelt, this depends on an internal state
     which goes stale and is reset whenever an insert or delete
     operation is performed. This returns a new-memory copy of the element
     in the pdstor, which is the caller's responsibility to free.*/

  void pd_write_pdstor(FILE *stream,pd_stor_t *pdstor);

  pd_stor_t *pd_read_pdstor(FILE *stream, pd_equivalence_t eq);
  /* Reads pdstor from file, storing elements up to equivalence relation eq.*/


  /* TOPOLOGICAL OPERATIONS.

     Perform various knot-theoretic operations on a pd_code_t.

  */

  bool pd_is_alternating(pd_code_t *pd);
  /* Tests whether the pd code is alternating. Assumes that all the crossing
     information is set. If some crossings aren't set, it will return a pd_error. */

  pd_code_t *pd_copy_newsize(pd_code_t *pd, pd_idx_t MAXVERTS);

  pd_code_t *pd_copy(pd_code_t *pd);
  /* Make a new-memory copy of pd */

  void pd_reorient_component(pd_code_t *pd, pd_idx_t cmp, pd_or_t orient);
  /* Reverse the orientation of component cmp iff or == PD_NEG_ORIENTATION */

  void pd_renumber_component(pd_code_t *pd, pd_idx_t cmp,pd_idx_t ofs);
  /* Changes the numbering of edges in a component by adding "ofs" cyclically
     to each edge number. */

  pd_code_t *pd_R1_loopdeletion(pd_code_t *pd,pd_idx_t cr);
  /* Performs a loop deletion Reidemeister 1 move at the crossing cr.

            +                   +
            |                   |
            |                   |
            |                   |
      +-------------+    ->     +-------+
      |     | cr
      |     |
      +-----+

   (A loop addition is a really different move, computationally speaking.) */

  pd_code_t *pd_R1_loop_addition(pd_code_t *pd,
                                   pd_idx_t f, pd_idx_t e_on_f);

  void pd_R2_bigon_elimination(pd_code_t *pd,pd_idx_t cr[2],
			       pd_idx_t  *noutpd,
			       pd_code_t ***outpd);

  /* Performs a bigon-elimination Reidemeister 2 move.

     |                    |     |               |
     |                    |     |               |
     |   +-----------+    |     +---------------+
     |   |           |    | ->
     +-cr[0]-------cr[1]--+
         |           |          +---------------+
         |           |          |               |


  Input is a pd code and two crossings defining the bigon. Output is a pair of
  pointers to child pd codes. There are several possibilities, because a bigon
  elimination may split the pd code into a pair of disconnected pd codes.

  noutpd counts the number of connected components of the output pd.
  outpd is set to a buffer of pd codes of size noutpd;

  outpd[0] is the pd code containing the component "on top" in the bigon. It may
  be a 0-crossing diagram.

  outpd[1] isn't allocated if noutpd == 1 (in this case, both components are part of the same code).

  if noutpd == 2, then

    outpd[1] contains the pd code of the component "on the bottom"
    in the bigon. Again, it might be a 0-crossing diagram.

  */

  pd_code_t *pd_R2_bigon_addition(pd_code_t *pd, pd_idx_t f,
                                    pd_idx_t e1_on_f, pd_idx_t e2_on_f,
                                    pd_or_t e1_over_e2_or);


  pd_code_t *pd_R3_triangle_flip(pd_code_t *pd, pd_idx_t f);

  pd_code_t *pd_connect_sum(pd_code_t *pdA, pd_idx_t edgeA,
			    pd_code_t *pdB, pd_idx_t edgeB);

  /*
  	  +----+	       +-------+
       	  |    |       	   +-<-|---+   |
       	  |    |  edgeA	   |   |   |   |
       	+-|------<---+ 	   |   |   |   |
      	| |    |     | 	   v   +-------+
      	| +------>---+ 	   +-->----+
       	+------+       	   edgeB

            pdA	       	       pdB

                      ||
                      vv

  	  +----+      	       +-------+
       	  |    |       	   +-<-|---+   |
       	  |    |       	   |   |   |   |
       	+-|------<---------+   |   |   |
      	| |    |       	       +-------+
      	| +------>------------>----+
       	+------+

                  (output pd)

  */

  pd_code_t *pd_simplify(pd_code_t *pd);

  /* Simplify the pd code using combinations of the moves above to
     reduce crossing number as far as possible. */

  /******************** pd human output ***************/

  /* These are functions which implement a superset of printf in order
     to display human-readable output and error messages.

     Tag        pd_idx_t           output

     %FACE      face number        face fnum (e1 (or1) -> e2 (or2) -> ... -> e1 (or1))
     %EDGE      edge number        edge enum (tail (tailpos) -> head (headpos) )
     %CROSS     cross number       cross cnum (e0 e1 e2 e3) +/-/U
     %COMP      comp number        compnum (e1 -> e2 -> e3 -> ..... -> e1 (or1))
     %PD        (no argument)      (\n\n output of pd_write \n\n)
     %FEDGE     face, edge         edge number (orientation on face)
                (2 pd_idx_ts)

     We also have conversions for pointers.

     %ORIENTATION *pd_orientation_t  multiorientation
     %OR          *pd_or_t           +, -, or U (unset)
     %MULTIDX     *pd_multidx_t      multidx (i[0] i[1] ... i[n-1])
     %COMPGRP     *pd_compgrp_t      compgrp (comp[0] .. comp[n-1])
     %COMPMAP     *pd_compmap_t      compmap (0 -> (+/-) comp[0], ..., n-1 -> (+/-) comp[n-1])
     %EDGEMAP     *pd_edgemap_t      edgemap (+/- e1 +/- e2 ... +/- en)
     %TANGLE      *pd_tangle_t       tangle display format
     %CROSSPTR    *pd_crossing_t     cross (e0 e1 e2 e3) (+/-/U)
     %ISO         *pd_iso_t          iso cr (# cross) e (# edges) f (# faces) cmps (# comps)
                                        compperm (permutation of components)
				        (edgemap, as above)
                                        (crossmap, as above)
					(facemap, as above)
     %DIHEDRAL    *pd_dihedral_t     rot (target of element 0) or ref (target of element 0)
     %CYCLIC      *pd_cyclic_t       rot (target of element 0)
     %PERM        *pd_perm_t         perm (map[0] ... map[n-1]) idx (precomputed perm index)
     %BDY_OR      *pd_boundary_or_t  in/out/?

     The function also converts %d and %s specifications in the usual way.

     We ignore any other format conversions present in fmt,
     passing them unchanged into the output stream. The intention
     is that these would have already been processed in fmt (using
     sprintf) before calling pd_vfprintf.

  */

  void pd_printf(char *fmt,pd_code_t *pd, ... );

#define SRCLOC __FILE__, __LINE__   /* A convenience tag to help us
				       specify a location for the error. */

  bool pd_error(char *file, int line, char *fmt, pd_code_t *pd, ...);

  /* If PD_VERBOSE > 10, outputs the error string in fmt,
     converted with the special format conversions
     to stderr, and then exits with error code 1.

     Otherwise, simply returns false. */

  /* There are some standard error checks which we give convenience functions for: */

  void pd_check_cr(char *file, int line, pd_code_t *pd, pd_idx_t cr);
  /* Checks if crossing number is legal, dies with error if not.
     Expected to be called pd_check_cr(SRCLOC,pd,cr). */

  void pd_check_edge(char *file, int line, pd_code_t *pd, pd_idx_t edge);
  /* Checks if edge number is legal, dies with error if not.
     Expected to be called pd_check_edge(SRCLOC,pd,edge). */

  void pd_check_cmp(char *file, int line, pd_code_t *pd, pd_idx_t cmp);
  /* Checks if component number is legal, dies with error if not.
     Expected to be called pd_check_cmp(SRCLOC,pd,cmp); */

  void pd_check_face(char *file, int line, pd_code_t *pd, pd_idx_t face);
  /* Checks if face number is legal, dies with error if not.
     Expected to be called pd_check_cmp(SRCLOC,pd,face); */

  void pd_check_notnull(char *file, int line, char *varname, void *ptr);
  /* Checks if pointer is null, dies with error if so. The field "varname" is a string giving
     the name of the pointer. Should be called like pd_check_notnull(SRCLOC,"invar",invar); */

  /********** Tangle Operations **************/

  /* It's often useful to cut "tangles" out of the center of pd codes.
     A tangle is determined by a cycle of edges e[0] -> e[nedges-1] and a
     corresponding cycle of faces f[0] -> f[nfaces-1], as below:

               f[0]

          +-------------+
          |             |
    e[0]--+             +---e[nedges-1]
          |             |
     f[1] |     T       |  f[nfaces-1]
          |             |
    e[1]--+             +---e[nedges-2]
          |             |
          +-------------+  f[nfaces-2]
     f[2]   |  ....   |


     Tangles must have unique edge sets (that is, the same edge is not
     to occur twice in the boundary of a tangle). The face set need
     not be unique: equivalently, a tangle a cycle in the dual graph
     of the diagram which does not revisit the same edge, but may
     pass through a vertex more than once.

       	....................   	       	..................
       	.      	/----\	   .		.  		 .
       	. 	|    |	   .		.	 /--------------
      --.-------|----/	   .	    ---------\	 |	 .
       	.       |          .   	       	.    |	 |	 .
       	.      	\----------.----	.    |	 |	 .
       	. 		   .		.    |	 |	 .
       	.   /--------------.----    ---------/ 	 |   	 .
     ---.---|---\    	   .		.      	 |   	 .
       	.   |  	|	   .	   --------------/	 .
	.   \---/	   .		..................
	....................

        ok (even through face            not ok (same edge occurs
	  occurs twice on bdy)             twice on the boundary)

     The edges which enter or leave the tangle can be joined
     (pairwise) by chains of edges inside the tangle. Each such chain
     of edges is called a "strand" of the tangle. Note that the union
     of the edges in the strands is NOT always the complete collection
     of edges inside the tangle because there may be entire components
     contained within the tangle which don't cross the boundary.

     Each edge which crosses the boundary of the tangle has a boundary
     orientation on the tangle, which records whether it is heading
     into or out of the tangle.

  */

  typedef enum {in,out,unset} pd_boundary_or_t;

  typedef struct tangle_strand_struct {

    pd_idx_t start_edge;  /* Edge number (in tangle) where strand enters tangle. */
    pd_idx_t end_edge;    /* Edge number (in tangle) where strand exits tangle. */
    pd_idx_t nedges;      /* Number of edges in tangle (counting start, end) */
    pd_idx_t comp;        /* Component (in pd) containing this strand. */

    /* The array of strands is stored in canonical order, sorted by
       the edge number (in the tangle) of the start edge. */

  } pd_tangle_strand_t;

  typedef struct pd_tangle_struct {

    pd_idx_t nedges;   /* Should always be even. */
    pd_idx_t nstrands; /* Always equal to nedges/2 */

    pd_idx_t *edge;
    pd_idx_t *face;

    pd_boundary_or_t *edge_bdy_or; /* in/out orientations for the edges e */
    pd_tangle_strand_t *strand;

    pd_idx_t  ninterior_cross;   /* Crossings in the interior of the tangle. */
    pd_idx_t *interior_cross;    /* List of interior crossing indices */
                                 /* NULL if there are no interior crossings */

    pd_idx_t  ninterior_edges;   /* Edges in the interior of the tangle. */
    pd_idx_t *interior_edge;     /* List of interior edge indices. */
                                 /* NULL if there are no interior edges */

  } pd_tangle_t;

  bool pd_tangle_ok(pd_code_t *pd,pd_tangle_t *t);

  pd_tangle_t *pd_tangle_new(pd_idx_t nedges);
  void pd_tangle_free(pd_tangle_t **t);

  void pd_regenerate_tangle(pd_code_t *pd,pd_tangle_t *t);
  void pd_regenerate_tangle_err(pd_code_t *pd, pd_tangle_t *t, int *err);
  /* The usual procedure for generating a tangle is to specify the
     loop of edges and faces and call "pd_regenerate_tangle" in order
     to reconstruct the remaining data. */

  pd_boundary_or_t pd_tangle_bdy_or(pd_code_t *pd,pd_tangle_t *t, pd_idx_t pd_edge_num);
  /* Finds the boundary orientation of an edge on the tangle given the
     index of the edge *in the pd*. (If you have the index of the edge
     in the t->edge array, you can just look up the orientation in the
     corresponding t->edge_bdy_or array directly.) */

  void pd_tangle_slide(pd_code_t *pd,pd_tangle_t *t,
		       pd_idx_t n,
		       pd_idx_t *overstrand_edges,
		       pd_idx_t *border_faces,
		       pd_idx_t *npieces,
		       pd_code_t ***pd_pieces);
  void pd_tangle_slide_err(pd_code_t *pd,pd_tangle_t *t,
                           pd_idx_t n,
                           pd_idx_t *overstrand_edges,
                           pd_idx_t *border_faces,
                           pd_idx_t *npieces,
                           pd_code_t ***pd_pieces,
                           int *err);


 /* Given a list of edges overstrand_edges (e[0]...e[n-1], below) and
    corresponding faces bordering the tangle (f[0]...f[n-1], below),
    slide the strand over the tangle to cross the remaining edges
    of the tangle, as below. The edges e[0]..e[n-1] are supposed to
    occur in orientation order along their component.


		  |  	   |
		  |  	   |
       	     +-------------------+
	     |		    	 |
             | 	    Tangle     	 |
    ---+     | 	       	       	 |    +---
       |     | 	    	    	 |    |
       |     +-------------------+    |
       | f[n-1]  |   |  f[1] | 	 f[0] |
       +--e[n-1]---<----e[1]-----e[0]-+
       	       	 |   | 	     |


    becomes

		  |    	   |
       +------------------------------+
       |	  |  	   |	      |
       |     +-------------------+    |
       |     |		    	 |    |
       |     | 	    Tangle     	 |    |
    ---+     | 	       	       	 |    +---
             | 	    	       	 |
             +-------------------+
                 |   |       |
                 |   |       |
       	       	 |   | 	     |

    We handle correctly the case where the initial and/or final
    edges of the strand are tangle edges themselves. We also note
    that while we call the strand the "overstrand", we also handle
    the case where the strand goes UNDER all of the tangle strands.

    This operation can disconnect the diagram, potentially into many
    pieces. We return the number of connected components of the
    diagram in "npieces" and the components themselves in
    "pd_pieces". The buffer of pd_code_t pointers pd_pieces is
    allocated internally and is the caller's responsibility to
    dispose of.

   */



  pd_code_t *pd_tangle_flype(pd_code_t *pd,pd_tangle_t *t);

  /*  This is a classical "flype" move, which we define by giving input
      and output edges. These must be connected in the (topological)
      situation below.
                           +---------+
      +in[0]---+    +------+         +----out[0]--+
               |    |      |         |
               |    |      |         |
            +--------+     |  VVVVV  |
            |  | cr        |         |
            |  |           |         |
      +in[1]+  +-----------+         +----out[1]--+
                           +---------+
                        |
                        v
         +---------+
         |         |
      +--+         +------+  +-----------------+
         |         |      |  |
         |  ^^^^^  |      |  | cr
         |         |      +------+
         |         |         |   |
      +--+         +---------+   +-------------+
         +---------+

  */

  /* Standard, valid pd codes (for test purposes). */

  pd_code_t *pd_build_twist_knot(pd_idx_t n);            /* Twist knot with n twists */
  pd_code_t *pd_build_torus_knot(pd_idx_t p,pd_idx_t q); /* Only 2, q implemented now */

  pd_code_t *pd_build_simple_chain(pd_idx_t n);          /* An n-link chain. */

  pd_code_t *pd_build_unknot(pd_idx_t n);                /* An n-crossing unknot diagram, where n >= 0 */

  /* An n-crossing unknot diagram in the form


         +--<--+       +--<--+            +-<--+
         |     |       |     |                 |
         |     |       |     |                 |
  +------|->-----------|->-------+  ...    +----->---+
  |      |     |       |     |                 |     |
  |      |     |       |     |                 ^     |
  +--<---+     +---<---+     +---+             +-----+

  with all crossings positive. This generates the
  case with n == 0 and n == 1 correctly. */

  pd_code_t *pd_build_unknot_wye(pd_idx_t a,pd_idx_t b,pd_idx_t c); /* An unknot diagram designed for hash collisions */

  /******* Interface with traditional plCurve types ******/

  pd_code_t *pd_code_from_plCurve(gsl_rng *rng, plCurve *L);

  /************************ Topological Invariants ********************/

  /* Compute the HOMFLY polynomial of a pd_code (returned as string) */
  /* Returns NULL if the timeout limit is reached. */

  char *pd_homfly( pd_code_t *pdC);
  char *pd_homfly_timeout( pd_code_t *pdC, int timeout);

  /* Compute the linking number of two components of a plcurve. */
  /* (Requires that crossing signs be set; otherwise, fails out.) */

    int pd_linking_number(pd_code_t *L,pd_idx_t c1,pd_idx_t c2);
    unsigned int pd_unsigned_linking_number(pd_code_t *L,pd_idx_t c1,pd_idx_t c2);

  /* Compute the HOMFLY polynomial of a plCurve (returned as string) */

  char *plc_homfly( gsl_rng *rng, plCurve *L);


#define MAXPRIMEFACTORS 10
#define MAXHOMFLY       1024

  typedef struct knottypestruct {

    int  nf;                            /* Number of prime factors */
    int  cr[MAXPRIMEFACTORS];           /* Crossing number of each prime factor */
    int  ind[MAXPRIMEFACTORS];          /* Index (in Rolfsen or Cerf) of each prime factor */
    char sym[MAXPRIMEFACTORS][128];     /* Symmetry tag (Whitten group element) for each prime factor */
    char homfly[MAXHOMFLY];             /* Homfly polynomial (as plc_lmpoly output) */

  } plc_knottype;

  /* Prints the knot type in a neatly formatted human-readable version */
  void plc_write_knottype(FILE *out,plc_knottype kt);

  /* Reads a knot-type as written by a human. Returns NULL if it can't parse the input string.*/
  plc_knottype *plc_read_knottype(const char *kt);


  /* Find the knot type of a single component plCurve */
  /* Sets nposs to the number of possible knottypes found for the curve. If we cannot
     classify the knot, return 0 for nposs and NULL for the buffer of knot types. */
  plc_knottype *plc_classify( gsl_rng *rng, plCurve *L, int *nposs);

  /* Find the knot type of a single component pdcode */
  /* Same return data as from plc_classify */
  plc_knottype *pd_classify(pd_code_t *pdC, int *nposs);

#if (__cplusplus || c_plusplus)
   }
#endif

#endif

