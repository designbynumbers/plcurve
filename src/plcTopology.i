%module(package="plcurve.pd") pd
%feature("autodoc", "1");
%{
#include "stdbool.h"
#include "plcTopology.h"
#include "pd_multidx.h"
#include "pd_perm.h"
#include "pd_isomorphisms.h"
#include "pd_storage.h"
%}
%include "typemaps.i"

typedef int pd_idx_t ;  /* pd "index" type */
typedef int pd_or_t;    /* pd "orientation" type */
typedef int pd_pos_t;   /* pd "position" type */
typedef int pd_uid_t;   /* pd "uid" type */

%rename(Edge) pd_edge_struct;
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

%rename(Component) pd_component_struct;
typedef struct pd_component_struct {

  pd_idx_t nedges;
  pd_idx_t *edge;

  /* Edge indices in orientation order
     around a component. These are expected
     to be consecutive. */

} pd_component_t;

%rename(ori) or;
%rename(Face) pd_face_struct;
typedef struct pd_face_struct {

  pd_idx_t    nedges;
  pd_idx_t    *edge;
  pd_or_t     *or;

  /* Edge indices around the face in counterclockwise
     order. These are NOT consecutive, nor are they always
     positively oriented (according to their intrinsic
     edge orientation), so we store their orientation
     as well as their edge number. */

} pd_face_t;

// Crossing structs are really just 4-tuples, so let's treat it as such.
#define _PY_PD_CROSSING_LEN 4
%{
#define _PY_PD_CROSSING_LEN 4
  inline PyObject *pd_crossing_as_pytuple(pd_crossing_t *c) {
    PyObject *result;
    int i;
    result = PyTuple_New(_PY_PD_CROSSING_LEN);
    for (i = 0; i < _PY_PD_CROSSING_LEN; i++) {
      PyTuple_SetItem(result, i, PyInt_FromLong((long)(c->edge[i])));
    }
    return result;
  }
%}

%typemap(in) pd_crossing_t (pd_crossing_t v) {
  int i;
  if (!PySequence_Check($input)) {
    PyErr_SetString(PyExc_TypeError,
                    "Expecting a sequence with _PY_PD_CROSSING_LEN elements");
    return NULL;
  }
  if (PyObject_Length($input) != _PY_PD_CROSSING_LEN) {
    PyErr_SetString(PyExc_ValueError,
                    "Expecting a sequence with _PY_PD_CROSSING_LEN elements");
    return NULL;
  }
  for (i =0; i < _PY_PD_CROSSING_LEN; i++) {
    PyObject *o = PySequence_GetItem($input,i);
    // TODO: Better handle pd_idx_t
    if (PyInt_Check(o)) {
      v.edge[i] = (pd_idx_t)PyInt_AsLong(o);
      Py_DECREF(o);
    } else if (PyLong_Check(o)) {
      v.edge[i] = (pd_idx_t)PyLong_AsLong(o);
      Py_DECREF(o);
    } else {
      Py_XDECREF(o);
      PyErr_SetString(PyExc_ValueError,"Expecting a sequence of ints or longs");
      return NULL;
    }
  }
  $1 = v;
 }
%typemap(out) pd_crossing_t {
  int i;
  $result = pd_crossing_as_pytuple(&$1);
 }
%typemap(out) pd_crossing_t *{
  int i;
  if ($1 == NULL) {
    $result = Py_None;
    return;
  }
  $result = pd_crossing_as_pytuple($1);
 }
/* %rename(Crossing) pd_crossing_struct; */
/* typedef struct pd_crossing_struct { */

/*   pd_idx_t edge[4]; */
/*   /\* Edge indices around crossing */
/*      in counterclockwise order *\/ */

/* } pd_crossing_t; */

// Out of bounds exception typemap pair
%typemap(in, numinputs=0) bool *OOB (bool e) {
  $1 = &e;
 }
%typemap(argout) bool *OOB {
  if ((*$1)) {
    PyErr_SetString(PyExc_IndexError, "Index is out of bounds");
    return NULL;
  }
 }

//  is the python wrapper for pd_code_t
// SWIG argument output typemaps
%apply int *OUTPUT { pd_idx_t *comp, pd_idx_t *comp_pos };

%typemap(in) FILE * {
  if (!PyFile_Check($input)) {
    PyErr_SetString(PyExc_Exception, "Expecting a file object");
    return NULL;
  }
  $1 = PyFile_AsFile($input);
 }

// SWIG wrapper declaration
%rename(PlanarDiagram) pd_code_struct;
typedef struct pd_code_struct {
  pd_uid_t uid;
  pd_idx_t  ncross;
  pd_idx_t  nedges;
  pd_idx_t  ncomps;
  pd_idx_t  nfaces;
  char hash[32];

  %typemap(out) pd_code_t *crossings {
    int i;
    PyObject *cr_tuple;
    $result = PyTuple_New($1->ncross);
    for (i = 0; i < $1->ncross; i++) {
      cr_tuple = pd_crossing_as_pytuple(&($1->cross[i]));
      PyTuple_SetItem($result, i, cr_tuple);
    }
  }
  %typemap(out) pd_code_t *edges {
    int i;
    $result = PyTuple_New($1->nedges);
    for (i = 0; i < $1->nedges; i++) {
      PyTuple_SetItem($result, i,
                      SWIG_NewPointerObj(((pd_edge_t*)$1->edge)+i,
                                         SWIGTYPE_p_pd_edge_struct, 0));
    }
  }
  %typemap(out) pd_code_t *components {
    int i;
    $result = PyTuple_New($1->ncomps);
    for (i = 0; i < $1->ncomps; i++) {
      PyTuple_SetItem($result, i,
                      SWIG_NewPointerObj(((pd_component_t*)$1->comp)+i,
                                         SWIGTYPE_p_pd_component_struct, 0));
    }
  }
  %typemap(out) pd_code_t *faces {
    int i;
    $result = PyTuple_New($1->nfaces);
    for (i = 0; i < $1->nfaces; i++) {
      PyTuple_SetItem($result, i,
                      SWIG_NewPointerObj(((pd_face_t*)$1->face)+i,
                                         SWIGTYPE_p_pd_face_struct, 0));
    }
  }

  // Extend Shadow with methods
  %extend {
    // Constructor
    pd_code_struct() { return (pd_code_t *)calloc(1, sizeof(pd_code_t)); }
    // Copy constructor
    pd_code_struct(pd_code_t *to_copy) { return pd_copy(to_copy); }
    // From file pointer constructor
    pd_code_struct(FILE *f) {
      return pd_read(f);
    }
    // Destructor
    ~pd_code_struct() { free($self); }

    // Copy as a method
    %newobject copy;
    inline pd_code_t *copy() { return pd_copy($self); }
    // Simplify a pd_code, returning a copy
    //  TODO? use name which better describes that this is not in place
    //%newobject simplify;
    // Not actually implemented...
    //inline pd_code_t *simplify() { return pd_simplify($self); }

    // Standard, valid pd code builders
    // Be sure to add %newobject directives
    %newobject from_twist_knot;
    inline static pd_code_t *from_twist_knot(pd_idx_t n) {
      return pd_build_twist_knot(n);
    }
    %newobject from_torus_knot;
    inline static pd_code_t *from_torus_knot(pd_idx_t p,pd_idx_t q) {
      return pd_build_torus_knot(p, q);
    }
    %newobject from_simple_chain;
    inline static pd_code_t *from_simple_chain(pd_idx_t n) {
      return pd_build_simple_chain(n);
    }
    %newobject_from_unknot;
    inline static pd_code_t *from_unknot(pd_idx_t n) {
      return pd_build_unknot(n);
    }
    %newobject from_unknot_wye;
    inline static pd_code_t *from_unknot_wye(pd_idx_t a,pd_idx_t b,pd_idx_t c) {
      return pd_build_unknot_wye(a, b, c);
    }

    // Fake variables which return member arrays
    const pd_code_t *const edges;
    const pd_code_t *const crossings;
    const pd_code_t *const components;
    const pd_code_t *const faces;

    // Get/set/add crossings
    inline pd_crossing_t get_crossing(pd_idx_t i, bool *OOB) {
      if(i < 0 || i >= $self->ncross) { *OOB = true; return; }
      *OOB = false;
      return ($self->cross[i]);
    }
    inline void set_crossing(pd_idx_t i, pd_crossing_t cross, bool *OOB) {
      if(i < 0 || i >= $self->ncross) { *OOB = true; return; }
      *OOB = false;
      $self->cross[i] = cross;
    }
    /* inline void add_crossing(pd_crossing_t cross, bool *OOB) { */
    /*   if($self->ncross >= PD_MAXVERTS) { *OOB = true; return; } */
    /*   *OOB = false; */
    /*   $self->cross[$self->ncross] = cross; */
    /*   $self->ncross++; */
    /* } */


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

    // Entire pd operation methods
    inline void write(FILE *f) { return pd_write(f, $self); }
    %newobject read;
    inline static pd_code_t *read(FILE *f) {
      return pd_read(f);
    }
    inline bool isomorphic(pd_code_t *to) {
      return pd_isomorphic($self, to);
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
  };


} pd_code_t;

%{
  inline pd_code_t *pd_code_struct_crossings_get(pd_code_t *pd) {
    return pd;
  }
  inline pd_code_t *pd_code_struct_edges_get(pd_code_t *pd) {
    return pd;
  }
  inline pd_code_t *pd_code_struct_components_get(pd_code_t *pd) {
    return pd;
  }
  inline pd_code_t *pd_code_struct_faces_get(pd_code_t *pd) {
    return pd;
  }

  %}
