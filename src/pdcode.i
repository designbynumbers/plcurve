%module pdcode
%{
#include "stdbool.h"
#include "pdcode.h"
%}
%include "typemaps.i"

typedef int pd_idx_t ;  /* pd "index" type */
typedef int pd_or_t;    /* pd "orientation" type */
typedef int pd_pos_t;   /* pd "position" type */
typedef int pd_uid_t;   /* pd "uid" type */

%rename(ShadowEdge) pd_edge_struct;
typedef struct pd_edge_struct {

  /* An oriented edge, joining two verts tail -> head */

  pd_idx_t head;
  pd_pos_t  headpos;
  /* Pos [0..3] in crossing record of head vertex. */

  pd_idx_t tail;
  pd_pos_t  tailpos;
  /* Pos [0..3] in crossing record of tail vertex. */

} pd_edge_t;

/* Since we may have loop edges, we need to record
   the position on the crossing where the head and
   tail of the edge come in.

   For a loop edge, the edge index occurs twice in
   the crossing record, so we can't determine
   orientation otherwise. */

%rename(ShadowComponent) pd_component_struct;
typedef struct pd_component_struct {

  pd_idx_t nedges;
  pd_idx_t edge[PD_MAXEDGES];

  /* Edge indices in orientation order
     around a component. These are expected
     to be consecutive. */

} pd_component_t;

%rename(ori) or;
%rename(ShadowFace) pd_face_struct;
typedef struct pd_face_struct {

  pd_idx_t    nedges;
  pd_idx_t    edge[PD_MAXEDGES];
  pd_or_t     or[PD_MAXEDGES];

  /* Edge indices around the face in counterclockwise
     order. These are NOT consecutive, nor are they always
     positively oriented (according to their intrinsic
     edge orientation), so we store their orientation
     as well as their edge number. */

} pd_face_t;

%rename(ShadowCrossing) pd_crossing_struct;
typedef struct pd_crossing_struct {

  pd_idx_t edge[4];
  /* Edge indices around crossing
     in counterclockwise order */

} pd_crossing_t;

// Shadow is the python wrapper for pd_code_t
// SWIG argument output typemaps
%apply int *OUTPUT { pd_idx_t *comp, pd_idx_t *comp_pos };

// SWIG wrapper declaration
%rename(Shadow) pd_code_struct;
typedef struct pd_code_struct {
  pd_uid_t uid;
  pd_idx_t  ncross;
  pd_idx_t  nedges;
  pd_idx_t  ncomps;
  pd_idx_t  nfaces;
  char hash[32];
  pd_edge_t      edge[PD_MAXEDGES];
  pd_component_t comp[PD_MAXCOMPONENTS];
  pd_crossing_t  cross[PD_MAXVERTS];
  pd_face_t      face[PD_MAXFACES];

} pd_code_t;

// Extend Shadow with methods
%extend pd_code_t {
    // Standard, valid pd code builders
    inline static pd_code_t *from_twist_knot(pd_idx_t n) {
	return pd_build_twist_knot(n);
    }
    inline static pd_code_t *from_torus_knot(pd_idx_t p,pd_idx_t q) {
	return pd_build_torus_knot(p, q);
    }
    inline static pd_code_t *from_simple_chain(pd_idx_t n) {
	return pd_build_simple_chain(n);
    }
    inline static pd_code_t *from_unknot(pd_idx_t n) {
	return pd_build_unknot(n);
    }
    inline static pd_code_t *from_unknot_wye(pd_idx_t a,pd_idx_t b,pd_idx_t c) {
	return pd_build_unknot_wye(a, b, c);
    }

    // Bind functions taking a pd_code_t* as members of python objects
    //  Miscellaneous methods
    inline void component_and_pos(pd_idx_t edge,
				  pd_idx_t *comp, pd_idx_t *comp_pos) {
	pd_component_and_pos($self, edge, comp, comp_pos);
    }
    inline void reorient_edge(pd_idx_t edge, pd_or_t or) {
	pd_reorient_edge($self, edge, or);
    }

    //  Recomputing methods
    inline void regenerate_crossings() {
	pd_regenerate_crossings($self);
    }
    inline void regenerate_comps() {
	pd_regenerate_comps($self);
    }
    inline void regenerate_faces() {
	pd_regenerate_faces($self);
    }
    inline void regenerate_hash() {
	pd_regenerate_hash($self);
    }
    inline void regenerate() {
	pd_regenerate($self);
    }

    //  Sanity checking methods
    inline bool cross_ok() {
	return pd_cross_ok($self);
    }
    inline bool edges_ok() {
	return pd_edges_ok($self);
    }
    inline bool faces_ok() {
	return pd_faces_ok($self);
    }
    inline bool comps_ok() {
	return pd_comps_ok($self);
    }
    inline bool is_ok() {
	return pd_ok($self);
    }
    
    //  pd_isomorphic
    //inline bool isomorphic(char *pdcodeA, int nA, char *pdcodeB, int nB) {
    //  return pd_isomorphic_strings(pdcodeA, nA, pdcodeB, nB);
    //}
};

//  pd_isomorphic
%inline %{ 
  bool isomorphic(char *pdcodeA, int nA, char *pdcodeB, int nB) {
    return pd_isomorphic_strings(pdcodeA, nA, pdcodeB, nB);
  }
  %}
