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

  /* The basic architecture of the pd_code is kind of intricate.  The
     problem is that we have to keep track of labelled diagram data in
     such a way that we can access things relatively quickly as we
     work and make an efficient search for isomorphisms, but also not
     lose track of component identities as we go. */

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
			It is a character, usually "A..Z" followed
			by lower case "a..z". It is independent
			of the position of the component in the component
			array because the component array gets resorted
			to be in canonical order. */
  } pd_component_t;

  typedef struct pd_face_struct {

    pd_idx_t    nedges;
    pd_idx_t    *edge;
    pd_or_t     *or;

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

  char pd_print_or(pd_or_t or);
  /* Returns a one-character printed form for "or": +, -, U (unset), or ? (anything else) */

  pd_or_t pd_compose_or(pd_or_t a,pd_or_t b);
  /* Returns the composition of the two orientation changes: ++ = -- = +, +- = -+ = - */

  bool pd_or_ok(pd_or_t or); /* Check whether or has a legal value. */

  int pd_or_cmp(const void *A,const void *B);
  /* Compare two *pd_or_t. */

  pd_crossing_t pd_build_cross(pd_idx_t e0,pd_idx_t e1,pd_idx_t e2,pd_idx_t e3);
  /* Builds a crossing from the given edge indices */

  void pd_canonorder_cross(pd_crossing_t *cr, pd_or_t or);
  /* Reverses (if or == PD_NEG_ORIENTATION) and then rotates
     a crossing into canonical order: cr->edge[0] is the
     lowest edge # */

  void pd_canonorder_face(pd_face_t *face, pd_or_t or);
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
  /* Sort by number of components, then dictionary order on the edge numbers.
     Used for searching, sorting. */

  void pd_component_and_pos(pd_code_t *pd,pd_idx_t edge,
			    pd_idx_t *comp,pd_idx_t *comp_pos);
  /* Returns component number and position on component of a
     given edge number (or die) */

  pd_idx_t pd_previous_edge(pd_code_t *pd, pd_idx_t edge);
  pd_idx_t pd_next_edge(pd_code_t *pd,pd_idx_t edge);

  /* Finds the number of the previous or next edge along the component containing edge. */

  void pd_face_and_pos(pd_code_t *pd, pd_idx_t edge,
		       pd_idx_t *posface, pd_idx_t *posface_pos,
		       pd_idx_t *negface, pd_idx_t *negface_pos);

  /* Finds the two faces which edge occurs on, which should include
     one face where edge appears in positive orientation (posface)
     and one where edge appears with negative orientation (negface).
     If we don't find two with this description, die.

     Return the position (index) of the edge on each face. */

  pd_edge_t pd_oriented_edge(pd_edge_t e,pd_or_t or);
  /* Returns original edge if or = PD_POS_EDGE_ORIENTATION,
     reversed edge if or = PD_NEG_EDGE_ORIENTATION */

  void pd_reorient_edge(pd_code_t *pd,pd_idx_t edge,pd_or_t or);
  /* Flips the edge in pd->edges[] if or == PD_NEG_ORIENTATION */

  /* Over and under information at a crossing. */

  void pd_overstrand(pd_code_t *pd,pd_idx_t cr,pd_idx_t *incoming_edgenum, pd_idx_t *outgoing_edgenum);

  /*  Returns the edge number  (that is, a number in 0..pd->nedges)
      of the incoming and outgoing edges of the strand going OVER at crossing cr of pd,
      using the sign of the crossing to determine. */

  void pd_understrand(pd_code_t *pd,pd_idx_t cr,pd_idx_t *incoming_edgenum, pd_idx_t *outgoing_edgenum);

  /*  Returns the edge number  (that is, a number in 0..pd->nedges)
      of the incoming and outgoing edges of the strand going UNDER at crossing cr of pd,
      using the sign of the crossing to determine. */

  void pd_overstrand_pos(pd_code_t *pd,pd_idx_t cr,pd_idx_t *incoming_edgepos, pd_idx_t *outgoing_edgepos);

  /* Returns the position in crossing cr of pd (that is, a number in 0..3) of the incoming and outgoing
     edges of the strand going over at this crossing, using the crossing sign data
     and edge orientations in order to compute the answer. */

  void pd_understrand_pos(pd_code_t *pd,pd_idx_t cr,pd_idx_t *incoming_edgepos, pd_idx_t *outgoing_edgepos);

  /* Returns the position in crossing cr of pd (that is, a number in 0..3) of the incoming and outgoing
     edges of the strand going under at this crossing, using the crossing sign data
     and edge orientations in order to compute the answer. */

  /* Functions to compute data from crossing (and other) information */

  void pd_regenerate_crossings(pd_code_t *pd);
  /* Reorders the crossings cyclically to put the
     lowest index edge first and sorts crossings
     in dictionary order based on this reordering. */

  void pd_regenerate_comps(pd_code_t *pd);
  /* Generates randomly oriented and numbered
     edges from crossings, then strings them
     into components, sorts components by
     size, renumbers edges and orients them along
     components. Updates and regenerates crossings. */

  void pd_regenerate_faces(pd_code_t *pd);
  /* Fills in faces from crossing, edge
     information, including orientation of
     each edge along the face. */

  void pd_regenerate_hash(pd_code_t *pd);
  /* Fill in (printable) hash value from
     comps, faces, and crossings */

  void pd_regenerate(pd_code_t *pd);
  /* Regenerates everything from cross data. */

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

  void       pd_write(FILE *outfile,pd_code_t *pd);
  void       pd_write_c(FILE *outfile, pd_code_t *pd, char *name);
  /* Writes a c procedure which recreates the pd code pd.
     The procedure will be called pd_create_name() and take 
     no arguments. */  

  pd_code_t *pd_read(FILE *infile); /* Returns NULL if the file is corrupt */

  /* This function reads a pdcode from the Mathematica package KnotTheory,
     exported as text with something like:

     Export["7_2.txt",PD[Knot[7,2]]]

     These PD codes don't have component or face information, so that is
     all regenerated once the crossings have been loaded from the file.
     This will only read one PD code per file.
  */

  pd_code_t *pd_read_KnotTheory(FILE *infile); /* Returns NULL if the file is corrupt */

  bool pd_isomorphic(pd_code_t *pdA,pd_code_t *pdB);
  /* Test whether two pd codes are isomorphic. */

  bool pd_isomorphic_strings(char *pdcodeA, int nA, char*pdcodeB, int nB);
  /* This method will generate pd code objects and return the call
     to pd_isomorphic. The purpose is to call this from Python via
     SWIG. */

  pd_code_t *pd_copy(pd_code_t *pd);
  /* Make a new-memory copy of pd */

  void pd_reorient_component(pd_code_t *pd, pd_idx_t cmp, pd_or_t or);
  /* Reverse the orientation of component cmp iff or == PD_NEG_ORIENTATION */

  pd_code_t *pd_R1_loopdeletion(pd_code_t *pd,pd_idx_t cr);
  /* Performs a loop deletion Reidemeister 1 move at the crossing cr. */
  /* (A loop addition is a really different move, computationally speaking.) */

  pd_code_t *pd_simplify(pd_code_t *pd);
  /* Simplify the pd code by eliminating loops and ``generalized loops''. */

  /* pd human output */

  /* These are functions which implement a superset of printf in order
     to display human-readable output and error messages.

     Tag        pd_idx_t           output

     %FACE      face number        fnum (e1 (or1) -> e2 (or2) -> ... -> e1 (or1))
     %EDGE      edge number        enum (tail (tailpos) -> head (headpos) )
     %CROSS     cross number       cnum (e0 e1 e2 e3)
     %COMP      comp number        compnum (e1 -> e2 -> e3 -> ..... -> e1 (or1))
     %PD        (no argument)      (\n\n output of pd_write \n\n)

     We also have conversions for pointers.

     %MULTIDX   *pd_multidx_t      multidx (i[0] i[1] ... i[n-1])
     %COMPGRP   *pd_compgrp_t      compgrp (comp[0] .. comp[n-1])
     %COMPMAP   *pd_compmap_t      compmap (0 -> (+/-) comp[0], ..., n-1 -> (+/-) comp[n-1])

     The functions also convert %d specifications in the usual way.

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



  /* Standard, valid pd codes (for test purposes). */

  pd_code_t *pd_build_twist_knot(pd_idx_t n);            /* Twist knot with n twists */
  pd_code_t *pd_build_torus_knot(pd_idx_t p,pd_idx_t q); /* Only 2, q implemented now */
  pd_code_t *pd_build_simple_chain(pd_idx_t n);          /* An n-link chain. */
  pd_code_t *pd_build_unknot(pd_idx_t n);                /* An n-crossing unknot diagram */
  pd_code_t *pd_build_unknot_wye(pd_idx_t a,pd_idx_t b,pd_idx_t c); /* An unknot diagram designed for hash collisions */

  /**************************** Interface with traditional plCurve types **********************************/

  pd_code_t *pd_code_from_plCurve(gsl_rng *rng, plCurve *L);

  /* Compute the HOMFLY polynomial of a pd_code (returned as string) */

  char *pd_homfly( pd_code_t *pdC);

  /* Compute the HOMFLY polynomial of a plCurve (returned as string) */

  char *plc_homfly( gsl_rng *rng, plCurve *L);

  /************************ plCurve Topology Library ********************/

  /* This contains some functionality designed to work with plCurves as knots,
     including converting them to an abstract ``crossing'' representation,
     computing their HOMFLY polynomials (using lmpoly) and identifying their
     knot types (by HOMFLY). */

#define MAXPRIMEFACTORS 10
#define MAXHOMFLY       1024

  typedef struct knottypestruct {

    int  nf;                            /* Number of prime factors */
    int  cr[MAXPRIMEFACTORS];           /* Crossing number of each prime factor */
    int  ind[MAXPRIMEFACTORS];           /* Index (in Rolfsen or Cerf) of each prime factor */
    char sym[MAXPRIMEFACTORS][128];     /* Symmetry tag (Whitten group element) for each prime factor */
    char homfly[MAXHOMFLY];             /* Homfly polynomial (as plc_lmpoly output) */

  } plc_knottype;


  /* Find the knot type of a single component plCurve */
  /* Sets nposs to the number of possible knottypes found for the curve. If we cannot
     classify the knot, return 0 for nposs and NULL for the buffer of knot types. */
  plc_knottype *plc_classify( gsl_rng *rng, plCurve *L, int *nposs);

#endif
