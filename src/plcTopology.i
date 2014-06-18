%module(package="libplcurve") pd
%feature("autodoc", "1");
%{
#include "stdbool.h"
#include "plcTopology.h"
#include "pd_multidx.h"
#include "pd_perm.h"
#include "pd_isomorphisms.h"
#include <stddef.h> // SWIG should include this itself but Debian version does not
  static char eoi = 0; // end of iteration exception
%}
%include "typemaps.i"

%{
// Secret function which is used in test_homfly
char *pdcode_to_ccode(pd_code_t *pd);
int __set_gc(PyObject *self) {
  return 1;
}
%}

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

%feature("python:slot", "sq_length", functype="lenfunc") pd_component_struct::__len__;
%feature("python:slot", "sq_item", functype="ssizeargfunc") pd_component_struct::__getitem__;
%rename(Component) pd_component_struct;
typedef struct pd_component_struct {

  pd_idx_t nedges;
  pd_idx_t *edge;

  /* Edge indices in orientation order
     around a component. These are expected
     to be consecutive. */

  %exception __getitem__ {
    assert(!eoi);
    $action
      if (eoi) {
        eoi = 0; // clear flag for next time
        PyErr_SetString(PyExc_StopIteration, "End of iterator");
        return NULL;
      }
  }

  %extend{
    const pd_idx_t __len__() const { return (Py_ssize_t)($self->nedges); }
    const pd_idx_t __getitem__(size_t j) const {
      if (j >= $self->nedges) {
        eoi=1;
        return -1;
      } else {
        return $self->edge[j];
      }
    }
  }

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

%feature("python:slot", "mp_subscript", functype="binaryfunc") pd_crossing_struct::__getitem__;
%rename(Crossing) pd_crossing_struct;
typedef struct pd_crossing_struct {
  pd_idx_t edge[4]; /* Edge indices around crossing in counterclockwise order */
  pd_or_t sign; /* Whether the crossing is positively or negatively oriented */

  %extend{
    const int __getitem__(int pos) const { return $self->edge[pos]; }
  }
} pd_crossing_t;

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
    $result = PyTuple_New($1->ncross);
    for (i = 0; i < $1->ncross; i++) {
      PyTuple_SetItem($result, i,
                      SWIG_NewPointerObj(((pd_crossing_t*)$1->cross)+i,
                                         SWIGTYPE_p_pd_crossing_struct, 0));
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

  %typemap(in) (PyObject *_secret_to_replace, int success) {
    SwigPyObject *sobj = (SwigPyObject *)$input;
    SwigPyObject *swig_self = (SwigPyObject *)self;
    sobj->ptr  = swig_self->ptr;
    sobj->ty   = swig_self->ty;
    sobj->own  = swig_self->own;
    sobj->next = 0;
    swig_self->own = false;
    $1 = NULL;
    $2 = 1;
  }

  // Extend Shadow with methods
  %extend {
    %feature("python:slot", "tp_is_gc", functype="inquiry") __set_gc;

    // Constructor
    pd_code_struct() { return (pd_code_t *)calloc(1, sizeof(pd_code_t)); }
    // Copy constructor
    pd_code_struct(pd_code_t *to_copy) { return pd_copy(to_copy); }

    int _secret_swig_usurp(PyObject *_secret_to_replace, int success) {
      return success;
    }

    // Destructor
    ~pd_code_struct() {
        pd_code_free(&$self);
    }

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
    %feature("docstring") from_twist_knot
         "Twist knot with ``n`` twists";
    %newobject from_twist_knot;
    inline static pd_code_t *from_twist_knot(pd_idx_t n) {
      return pd_build_twist_knot(n);
    }
    %feature("docstring") from_torus_knot
        "Torus knot. Only (2,q) is implemented for now.";
    %newobject from_torus_knot;
    inline static pd_code_t *from_torus_knot(pd_idx_t p,pd_idx_t q) {
      return pd_build_torus_knot(p, q);
    }
    %feature("docstring") from_twist_knot
        "An ``n``-link chain";
    %newobject from_simple_chain;
    inline static pd_code_t *from_simple_chain(pd_idx_t n) {
      return pd_build_simple_chain(n);
    }
    %feature("docstring") from_unknot
        "An ``n``-crossing unknot diagram";
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
    inline pd_crossing_t *get_crossing(pd_idx_t i, bool *OOB) {
      if(i >= $self->ncross) { *OOB = true; return NULL; }
      *OOB = false;
      return &($self->cross[i]);
    }
    inline void set_crossing(pd_idx_t i, pd_crossing_t cross, bool *OOB) {
      if(i >= $self->ncross) { *OOB = true; return; }
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
    inline void reorient_component(pd_idx_t cmp, pd_or_t or) {
      pd_reorient_component($self, cmp, or);
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

    // Topology methods
    %newobject homfly;
    inline char *homfly() {
        return pd_homfly($self);
    }
    %newobject ccode;
    char *ccode() {
        return pdcode_to_ccode($self);
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

    // String & display methods
    // TODO: Finish this implementation. May need to be done with Python? instead of SWIG
    //  or require different implementation of pd_printf.
    void printf(char *fmt, ...) {
        pd_printf(fmt,$self);
    }
    %feature("python:slot", "tp_str", functype="reprfunc") __str__;
    %feature("python:slot", "tp_repr", functype="reprfunc") __repr__;
    %newobject __str__;
    char *__str__() {
        char *buf;
        buf = malloc(255*sizeof(char)); // TODO: make more generic

        snprintf(buf, 255, "PlanarDiagram of %d crossings made up of %d components",
                $self->ncross,
                $self->ncomps);
        return buf;
    }
    %newobject __repr__;
    char *__repr__() {
        const int BUFLEN=500;
        char *buf, *nbuf;
        int edge, cross;
        size_t pos;
        buf = malloc(BUFLEN*sizeof(char)); // TODO: make more generic

        pos = snprintf(buf, BUFLEN, "PlanarDiagram(edges=[");
        for (edge = 0; edge<$self->nedges; edge++) {
            nbuf = buf + pos*sizeof(char);
            pos += snprintf(nbuf, BUFLEN-pos, "%u_%u->%u_%u, ",
                           (unsigned int)($self->edge[edge].tail),
                           (unsigned int)($self->edge[edge].tailpos),
                           (unsigned int)($self->edge[edge].head),
                           (unsigned int)($self->edge[edge].headpos));
        }
        pos -= 2;
        nbuf = buf + pos*sizeof(char);
        pos += snprintf(nbuf, BUFLEN-pos, "], cross=[");
        for (cross = 0; cross<$self->ncross; cross++) {
            nbuf = buf + pos*sizeof(char);
            pos += snprintf(nbuf, BUFLEN-pos, "%u.%u.%u.%u%s, ",
                            (unsigned int)($self->cross[cross].edge[0]),
                            (unsigned int)($self->cross[cross].edge[1]),
                            (unsigned int)($self->cross[cross].edge[2]),
                            (unsigned int)($self->cross[cross].edge[3]),
                            $self->cross[cross].sign != PD_UNSET_ORIENTATION ?
                            ($self->cross[cross].sign == PD_POS_ORIENTATION ? "+" : "-") : "x");
        }
        pos -= 2;
        nbuf = buf + pos*sizeof(char);
        pos += snprintf(nbuf, BUFLEN-pos, "])");

        return buf;
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
