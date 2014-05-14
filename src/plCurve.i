%module plcurve
%feature("autodoc", "1");
%{
#include "plCurve.h"
#include "matrix.h"
#include "plcTopology.h"
    //#include "homfly.h"
#include <gsl/gsl_rng.h>
#include <stddef.h> // SWIG should include this itself but Debian version does not

  // varray is a hack-in to support passing variable length arrays to Python
  typedef struct variable_array_type {
    int len;
    void *array;
  } varray;
  %}

%include "carrays.i"
%include "typemaps.i"
%import "plcTopology.i"

typedef double plc_matrix[3][3];

struct plc_vertex_loc {

  int cp;
  int vt;

};

typedef struct plc_type plCurve; /* We need to forward declare the plCurve type. */

%rename(Symmetry) plc_symmetry_type;
typedef struct plc_symmetry_type {

  plc_matrix              *transform;
  plCurve                 *curve;
  struct plc_vertex_loc  **target; /* Array of curve->nc arrays of curve->cp[cp].nv arrays of plc_vertex_loc */

} plc_symmetry;

%rename(SymmetryGroup) plc_symmetry_group_type;
typedef struct plc_symmetry_group_type { /* This requires a little bit of the group structure to be specified. */

  int n;
  plc_symmetry **sym; /* This is just a list of the group elements. */
  int *inverse; /* the (index of) the inverse of symmetry i is given by G->inverse[i] */

} plc_symmetry_group;

// It is sufficient to pretend that plc_vectors are just python sequences of length 3
%typemap(out) plc_vector {
  int i;
  $result = PyTuple_New(3);
  for (i = 0; i < 3; i++) {
    PyTuple_SetItem($result, i, PyFloat_FromDouble($1.c[i]));
  }
 }
%typemap(in) plc_vector(plc_vector v) {
  int i;
  if (!PySequence_Check($input)) {
    PyErr_SetString(PyExc_TypeError,"Expecting a sequence");
    return NULL;
  }
  if (PyObject_Length($input) != 3) {
    PyErr_SetString(PyExc_ValueError,"Expecting a sequence with 3 elements");
    return NULL;
  }
  for (i =0; i < 3; i++) {
    PyObject *o = PySequence_GetItem($input,i);
    if (!PyFloat_Check(o) && !PyInt_Check(o) && !PyLong_Check(o)) {
      Py_XDECREF(o);
      PyErr_SetString(PyExc_ValueError,"Expecting a sequence of floats or ints");
      return NULL;
    }
    v.c[i] = PyFloat_AsDouble(o);
    Py_DECREF(o);
  }
  $1 = v;
 }

%{
  static int fail_if_non_py_numeric(PyObject *o, char *fail_msg) {
    if (!PyFloat_Check(o) && !PyInt_Check(o) && !PyLong_Check(o)) {
      PyErr_SetString(PyExc_TypeError, fail_msg);
      return 0;
    }
    return 1;
  }
  static int convert_tuple_to_color(PyObject *py_color, plc_color *clr) {
    PyObject *n;
    if (!PySequence_Check(py_color) || PyObject_Length(py_color) < 4) {
      PyErr_SetString(PyExc_TypeError,
                      "Vertex list must contain ordered quadruple");
      return 0;
    }

    // Red
    n = PySequence_GetItem(py_color, 0);
    if (fail_if_non_py_numeric(n, "R Color components must be numerical")) {
      clr->r = PyFloat_AsDouble(n);
    } else { return 0; }
    // Green
    n = PySequence_GetItem(py_color, 1);
    if (fail_if_non_py_numeric(n, "G Color components must be numerical")) {
      clr->g = PyFloat_AsDouble(n);
    } else { return 0; }
    // Blue
    n = PySequence_GetItem(py_color, 2);
    if (fail_if_non_py_numeric(n, "B Color components must be numerical")) {
      clr->b = PyFloat_AsDouble(n);
    } else { return 0; }
    // Alpha
    n = PySequence_GetItem(py_color, 3);
    if (fail_if_non_py_numeric(n, "A Color components must be numerical")) {
      clr->alpha = PyFloat_AsDouble(n);
    } else { return 0; }

    return 1; // Success!
  }

  %}

// Example functions to test plc_vector typemaps (TODO: remove later)
plc_vector plc_build_vect(const double x, const double y, const double z);
plc_vector plc_vect_sum(plc_vector A, plc_vector B);

%feature("python:slot", "tp_str", functype="reprfunc") plc_strand_type::__str__;
%feature("python:slot", "tp_repr", functype="reprfunc") plc_strand_type::__repr__;


%rename(Strand) plc_strand_type;
typedef struct plc_strand_type {
  %rename(num_vertices) nv;
  %rename(is_open) open;
  %rename(num_colors) cc;
  int            nv;     /* Number of vertices */
  const bool     open;   /* This is an "open" strand (with distinct ends) */
  int            cc;     /* Color count (number of colors) */

  // Obscure direct pointer access to these variable-length arrays
  //plc_vector   *vt;     /* Actual vertices */
  //plc_color    *clr;    /* Colors */

  %typemap(out) varray vertices {
    PyObject *temp_entry;
    int i, j;
    $result = PyList_New(0);
    for (i = 0; i < $1.len; i++) {
      temp_entry = PyTuple_New(3);
      for (j = 0; j < 3; j++) {
        PyTuple_SetItem(temp_entry, j,
                        PyFloat_FromDouble(((plc_vector*)$1.array)[i].c[j]));
      }
      PyList_Append($result, temp_entry);
    }
  }
  %typemap(out) varray colors {
    PyObject *temp_entry;
    int i, j;
    $result = PyList_New(0);
    for (i = 0; i < $1.len; i++) {
      temp_entry = PyTuple_New(4);
      PyTuple_SetItem(temp_entry, 0,
                      PyFloat_FromDouble(((plc_color*)$1.array)[i].r));
      PyTuple_SetItem(temp_entry, 1,
                      PyFloat_FromDouble(((plc_color*)$1.array)[i].g));
      PyTuple_SetItem(temp_entry, 2,
                      PyFloat_FromDouble(((plc_color*)$1.array)[i].b));
      PyTuple_SetItem(temp_entry, 3,
                      PyFloat_FromDouble(((plc_color*)$1.array)[i].alpha));
      PyList_Append($result, temp_entry);
    }
  }

  %extend {
    const varray vertices;
    const varray colors;

    // Python special methods
    //
    const char *__str__() {
      char buf[255];

      sprintf(buf, "Strand (%s) with %d vertices",
              $self->open ? "open" : "closed",
              $self->nv);
      return buf;
    }
    const char *__repr__() {
		char buf[255];

		sprintf(buf, "Strand (%s) with %d vertices",
				$self->open ? "open" : "closed",
				$self->nv);
		return buf;
    }
  }
} plc_strand;

%{
  const varray plc_strand_type_vertices_get(plc_strand *st) {
    varray arr;
    arr.len = st->nv;
    arr.array = st->vt;
    return arr;
  }
  const varray plc_strand_type_colors_get(plc_strand *st) {
    varray arr;
    arr.len = st->cc;
    arr.array = st->clr;
    return arr;
  }
  %}

/* Curve constraint kind */
%rename(ConstantKind) plc_cst_kind_type;
typedef enum plc_cst_kind_type {
  unconstrained = 0,
  fixed,
  line,
  plane,
} plc_cst_kind;

%rename(Constraint) plc_constraint_type;
typedef struct plc_constraint_type {
  plc_cst_kind kind;    /* What kind of constraint */
  plc_vector   vect[2]; /* Vectors to define plane, line or fixed point */
  int           cmp;     /* Component */
  int           vert;    /* Starting vertex */
  int           num_verts; /* Length of run */
  /*@only@*/ /*@null@*/ struct plc_constraint_type *next;
} plc_constraint;

/* The vectors in vect change meaning depending on the type of constraint:

   fixed : vect[0] is the fixed point, vect[1] is ignored
   line:   vect[0] is the tangent vector, vect[1] is a point on the line
   plane:  vect[0] is the normal vector, vect[1].c[0] is the distance from origin.

*/

%rename(VertexQuantifier) plc_vert_quant_type;
typedef struct plc_vert_quant_type { /* Vertex quantifiers */
  int              cmp;    /* Component */
  int              vert;   /* Vertex */
  char             tag[4]; /* 3-character tag */
  double           quant;  /* Quantifier */
  /*@only@*/ /*@null@*/ struct plc_vert_quant_type *next;
} plc_vert_quant;

%{
  //char *plc_serialize(plCurve *const L, int *ret_len);
  %}

/*
  plCurve wrapper SWIG declaration
  ================================================
  Wraps most functions which take (at least one) plCurve* argument
*/
%feature("python:slot", "tp_str", functype="reprfunc") plc_type::__str__;
%feature("python:slot", "tp_repr", functype="reprfunc") plc_type::__repr__;
%rename(PlCurve) plc_type;
struct plc_type {
  %rename(num_components) nc;
  int         nc;                              /* Number of components */

  // Hide raw components pointer
  // plc_strand *cp;                              /* Components */
  /*@only@*/ /*@null@*/ plc_constraint *cst;   /* Constraints */
  /*@only@*/ /*@null@*/ plc_vert_quant *quant; /* per-vertex quantifiers */

  plc_symmetry_group *G;                       /* Symmetry group (may be null) */

  %typemap(in, numinputs=0) (char **var_buf, int *len) (char *v, int l){
    $1 = &v; $2 = &l;
  }
  %typemap(argout) (char **var_buf, int *len) {
    $result = PyString_FromStringAndSize(*$1, *$2);
  }
  %typemap(out) varray components {
    int i, j;
    $result = PyList_New(0);
    for (i = 0; i < $1.len; i++) {
      PyList_Append($result,
                    SWIG_NewPointerObj(((plc_strand*)$1.array)+i,
                                       SWIGTYPE_p_plc_strand_type, 0));
    }
  }
  %typemap(in) (const int nv, const plc_vector *const vt) {
    PyObject *py_vtx, *n;
    // Assert that input is a sequence
    if (PySequence_Check($input)) {
      int i, j;
      $1 = PyObject_Length($input); // Set length argument
      $2 = malloc(($1)*sizeof(plc_vector)); // Malloc array
      for (i = 0; i < $1; i++) {
        py_vtx = PySequence_GetItem($input,i);
        // Assert that input sequence contains a triple
        if (!PySequence_Check(py_vtx) || PyObject_Length(py_vtx) < 3) {
          PyErr_SetString(PyExc_TypeError,
                          "Vertex list must contain ordered triples");
          free($2);
          return NULL;
        }

        for (j = 0; j < 3; j++) {
          n = PySequence_GetItem(py_vtx, j);
          if (PyFloat_Check(n) || PyInt_Check(n)) {
            $2[i].c[j] = PyFloat_AsDouble(n);
          } else {
            PyErr_SetString(PyExc_TypeError,
                            "Vector coordinates must be numerical");
            free($2);
            return NULL;
          }
        }
      }
    } else {
      PyErr_SetString(PyExc_TypeError,"nv must be a sequence");
      return NULL;
    }
  }
  %typemap(in) (const int cc, const plc_color *const clr) {
    PyObject *py_color;
    // Assert that input is a sequence
    if (PySequence_Check($input)) {
      int i, j;
      $1 = PyObject_Length($input); // Set length argument
      $2 = malloc(($1)*sizeof(plc_color)); // Malloc array
      for (i = 0; i < $1; i++) {
        py_color = PySequence_GetItem($input,i);
        if(!convert_tuple_to_color(py_color, &($2[i]))) {
          free($2);
          PyErr_SetString(PyExc_TypeError,"I'm confused.");
          return NULL;
        }

      }
    } else if(!$input || $input == Py_None) {
      $1 = 0; $2 = NULL;
    } else {
      PyErr_SetString(PyExc_TypeError,"clr must be a sequence, or none");
      return NULL;
    }
  }

  %typemap(in, numinputs=0) int *nposs (int temp) {
    $1 = &temp;
  }
  %typemap(argout) int *nposs {
    PyObject *knottype, *np, *o3;

    np = PyInt_FromLong(*$1);
    if(!PyTuple_Check($result)) {
      knottype = $result;
      $result = PyTuple_New(1);
      PyTuple_SetItem($result,0,knottype);
    }
    o3 = PyTuple_New(1);
    PyTuple_SetItem(o3,0,np);
    knottype = $result;
    $result = PySequence_Concat(knottype,o3);
    Py_DECREF(knottype);
    Py_DECREF(o3);
  }

  %typemap(in, numinputs=0) bool *ok (bool success) {
    $1 = &success;
  }
  %typemap(argout) bool *ok {
    if (!(*$1)) {
      PyErr_SetString(PyExc_Exception,"Error in input to function");
      return NULL;
    }

  }

  %typemap(in, numinputs=0) (int *error_num, char error_str[], size_t error_str_len)
    (int i, char c[256]) {
    $1 = &i; $2 = c; $3 = 256;
  }
  %typemap(out) (int *error_num, char error_str[], size_t error_str_len) {
    if ((*$1)) { // If there is an error
      // TODO: return better exceptions depending on error number
      // TODO: Only return string up until end of error
      PyErr_SetString(PyExc_Exception,PyString_FromStringAndSize($2, 256));
      return NULL;
    }
  }

  %typemap(in) FILE * {
    if (!PyFile_Check($input)) {
      PyErr_SetString(PyExc_Exception, "Expecting a file object");
      return NULL;
    }
    $1 = PyFile_AsFile($input);
  }

  // SWIG extensions
  %extend {
    // Virtual class members
    //
    const varray components;

    // Constructors and destructor
    //
    plc_type(const int components,
             const int * const nv,
             const bool * const open,
             const int * const cc) { return plc_new(components, nv, open, cc); }
    plc_type(const plCurve * const L) { return plc_copy(L); }
    ~plc_type() {
      plc_free($self);
    }

    // File i/o
    %newobject from_file;
    static plCurve *from_file(FILE *file, int *error_num,
                              char error_str[], size_t error_str_len) {
      return plc_read(file, error_num, error_str, error_str_len);
    }
    void write(FILE *outfile) {
      plc_write(outfile, $self);
    }

    // Random PlCurve creators
    //
    %feature("docstring") random_closed_polygon
         "Generate random length 2 space polygons of nEdges edges using
the symmetric measure of Cantarella, Deguchi, Shonkwiler";
    %newobject random_closed_polygon;
    static plCurve *random_closed_polygon(gsl_rng *r, int nEdges)
    {
        return plc_random_closed_polygon(r, nEdges); }

    %feature("docstring") random_open_polygon
         "Generate random length 2 space polygons of nEdges edges using
the symmetric measure of Cantarella, Deguchi, Shonkwiler";
    %newobject random_open_polygon;
    static plCurve *random_open_polygon(gsl_rng *r,int nEdges)
    { return plc_random_open_polygon(r, nEdges); }

    %feature("docstring") random_closed_plane_polygon
         "Generate random length 2 planar polygons of nEdges edges using
the symmetric measure of Cantarella, Deguchi, Shonkwiler";
    %newobject random_closed_plane_polygon;
    static plCurve *random_closed_plane_polygon(gsl_rng *r,int nEdges)
    { return plc_random_closed_plane_polygon(r, nEdges); }

    %feature("docstring") random_open_plane_polygon
         "Generate random length 2 planar polygons of nEdges edges using
the symmetric measure of Cantarella, Deguchi, Shonkwiler";
    %newobject random_open_plane_polygon;
    static plCurve *random_open_plane_polygon(gsl_rng *r,int nEdges)
    { return plc_random_open_plane_polygon(r, nEdges); }

    %feature("docstring") random_equilateral_closed_polygon
         "Random equilateral polygons (this uses a hedgehog/fold
combination method).  They are scaled to total length 2 in order
to be comparable to the polygons generated by the other methods.";
    %newobject random_equilateral_closed_polygon;
    static plCurve *random_equilateral_closed_polygon(gsl_rng *r,int nEdges)
    { return plc_random_equilateral_closed_polygon(r, nEdges); }

    %feature("docstring") random_equilateral_open_polygon
         "Random equilateral polygons (this uses a hedgehog/fold
combination method).  They are scaled to total length 2 in order
to be comparable to the polygons generated by the other methods.";
    %newobject random_equilateral_open_polygon;
    static plCurve *random_equilateral_open_polygon(gsl_rng *r,int nEdges)
    { return plc_random_equilateral_open_polygon(r, nEdges); }

    %newobject loop_closure;
    static plCurve *loop_closure(gsl_rng *r,int cp,plCurve *openL,int nEdges)
    { return plc_loop_closure(r, cp, openL, nEdges); }

    // Data operation methods
    //
    void add_component(const int add_as, const int nv, const plc_vector *const vt,
                       const bool open, const int cc, const plc_color *const clr) {
      plc_add_component($self, add_as, nv, open, cc, vt, clr);
    }
    inline void drop_component(const int cmp) {
      plc_drop_component($self, cmp);
    }
    inline void resize_colorbuf(const int cp, const int cc) {
      plc_resize_colorbuf($self, cp, cc);
    }
    inline void set_fixed(const int cmp, const int vert, const plc_vector point) {
      plc_set_fixed($self, cmp, vert, point);
    }
    inline void constrain_to_line(const int cmp,
                                  const int vert,
                                  const int num_verts,
                                  const plc_vector tangent,
                                  const plc_vector point_on_line) {
      plc_constrain_to_line($self, cmp, vert,
                            num_verts, tangent, point_on_line);
    }
    inline void constrain_to_plane(const int cmp,
                                   const int vert,
                                   const int num_verts,
                                   const plc_vector normal,
                                   const double dist_from_origin) {
      plc_constrain_to_plane($self, cmp, vert,
                             num_verts, normal, dist_from_origin);
    }
    inline void unconstrain(const int cmp, const int vert,
                            const int num_verts) {
      plc_unconstrain($self, cmp, vert, num_verts);
    }
    inline void remove_constraint(const plc_cst_kind kind,
                                  const plc_vector vect[]) {
      plc_remove_constraint($self, kind, vect);
    }
    inline void remove_all_constraints() {
      plc_remove_all_constraints($self);
    }
    // TODO: make this just return the constraints? (empty list means false)
    // inline bool is_constrained(int cp, int vt, plc_constraint **constraint)
    inline void fix_wrap() { plc_fix_wrap($self); }
    inline void fix_constraints() { plc_fix_cst($self); }
    %newobject double_vertices;
    inline plCurve *double_vertices() { return plc_double_verts($self); }

    // Spline methods
    %newobject convert_to_spline;
    plc_spline *convert_to_spline(bool *ok) {
      return plc_convert_to_spline($self, ok);
    }

    // Geometric information methods
    inline int num_edges() { return plc_edges($self, NULL); }
    inline int num_vertices() { return plc_num_verts($self); }
    inline int vertex_num(const int cp, const int vt) {
      return plc_vertex_num($self, cp, vt);
    }
    inline int cp_num(int wrapVt) { return plc_cp_num($self, wrapVt); }
    inline int vt_num(int wrapVt) { return plc_vt_num($self, wrapVt); }
    inline double turning_angle(const int component, const int vertex, bool *ok) {
      return plc_turning_angle($self, component, vertex, ok);
    }
    inline double MR_curvature(const int component, const int vertex) {
      return plc_MR_curvature($self, component, vertex);
    }
    // TODO: Overload these 2 for possible component-wise return values
    double total_curvature() {
      return plc_totalcurvature($self, NULL);
    }
    double total_torsion() {
      return plc_totaltorsion($self, NULL);
    }
    //float *plc_dihedral_angles() needs out array typemap...
    inline plc_vector mean_tangent(const int component, const int vertex, bool *ok) {
      return plc_mean_tangent($self, component, vertex, ok);
    }
    // arclength
    inline double subarc_length(const int cmp, const int v1, const int v2) {
      return plc_subarc_length($self, cmp, v1, v2);
    }
    inline double s(const int cmp, const int vert) {
      return plc_s($self, cmp, vert);
    }
    // inline void edgelength_stats
    inline double check_constraint() { return plc_check_cst($self); }
    const double pointset_diameter;
    const plc_vector center_of_mass;
    const double gyradius;
    inline double mean_squared_chordlengths(int cp, int *skips, int nskips) {
      plc_mean_squared_chordlengths($self, cp, skips, nskips); // TODO: typemap
    }
    // nearest_neighbor
    // nearest_vertex

    // Geometric Operations
    //
    inline void scale(const double alpha) {
      plc_scale($self, alpha);
    }
    //inline void whitten(int mirror, int *eps, int *perm);
    inline void pfm(int cp, int vt0, int vt1, double angle) {
      plc_pfm($self, cp, vt0, vt1, angle);
    }
    inline void rotate(plc_vector axis, double angle) {
      plc_rotate($self, axis, angle);
    }
    inline void random_rotate(plc_vector axis) {
      plc_random_rotate($self, axis);
    }
    inline void translate(plc_vector translation) {
      plc_translate($self, translation);
    }
    inline void perturb(double radius) {
      plc_perturb($self, radius);
    }
    inline void project(plc_vector normal) {
      plc_project($self, normal);
    }
    inline void delete_arc(int cp, int vt1, int vt2) {
      plc_delete_arc($self, cp, vt1, vt2);
    }

    // Symmetry methods
    // presently omitted
    //

    // Topology methods
    //
    /* const char *const ccode; */
    /* const char *const homfly; */
    %newobject classify;
    plc_knottype *classify(gsl_rng *r, int *nposs)
    { return plc_classify(r, $self, nposs); }

    // Python special methods
    //
    const char *__str__() {
      char buf[255];

      sprintf(buf, "PlCurve with %d components", $self->nc);
      return buf;
    }
    const char *__repr__() {
		char buf[255];

		sprintf(buf, "PlCurve with %d components", $self->nc);
		return buf;
    }

    //inline void serialize(char **var_buf, int *len) {
    //  *var_buf = (char *)(plc_serialize($self, len));
    //}
  }
};

%{
  const varray plc_type_components_get(plCurve *curve) {
    varray arr;
    arr.len = curve->nc;
    arr.array = curve->cp;
    return arr;
  }
  const double plc_type_pointset_diameter_get(plCurve *c) {
    return plc_pointset_diameter(c);
  }
  const plc_vector plc_type_center_of_mass_get(plCurve *c) {
    return plc_center_of_mass(c);
  }
  const double plc_type_gyradius_get(plCurve *c) {
    return plc_gyradius(c);
  }
  inline char *plc_type_ccode_get(plCurve *c) {
    return old_plc_ccode(c);
  }
  inline char *plc_type_homfly_get(plCurve *c) {
      //return plc_homfly(c);
  }
  %}

/* PlCurve_spline types */
%rename(SplineStrand) plc_spline_strand_type;
typedef struct plc_spline_strand_type {
  bool         open;     /* This is an "open" strand (with distinct ends) */
  int          ns;       /* Number of samples used to build spline. */
  double      *svals;    /* s values at samples */
  plc_vector *vt;       /* positions at these s values */
  plc_vector *vt2;      /* _second_ derivatives at these s values */
  int          cc;       /* Number of colors */
  plc_color  *clr;      /* color values at samples */
} plc_spline_strand;

%rename(Spline) plc_spline_type;
typedef struct plc_spline_type {
  int                 nc;     /* Number of components */
  plc_spline_strand  *cp;     /* Components */
  plc_constraint     *cst;    /* Constraints */
} plc_spline;

#define MAXPRIMEFACTORS 10
#define MAXHOMFLY       1024

%{
  typedef struct {
    int cr;
    int ind;
    char sym[128];
  } prime_factor;
  %}
typedef struct {
  int cr;
  int ind;
  char sym[128];
	%extend {
		~prime_factor() {
			free($self);
		}
	}
} prime_factor;

%rename(KnotType) knottypestruct;
typedef struct knottypestruct {
  %rename(num_factors) nf;

  int  nf;                            /* Number of prime factors */
  int  cr[MAXPRIMEFACTORS];           /* Crossing number of each prime factor */
  int  ind[MAXPRIMEFACTORS];           /* Index (in Rolfsen or Cerf) of each prime factor */
  char sym[MAXPRIMEFACTORS][128];     /* Symmetry tag (Whitten group element) for each prime factor */
  char homfly[MAXHOMFLY];             /* Homfly polynomial (as plc_lmpoly output) */

	%typemap(out) plc_knottype *factors {
		int i;
		prime_factor *f;
		$result = PyTuple_New($1->nf);
		for (i = 0; i < $1->nf; i++) {
			f = calloc(1, sizeof(prime_factor));
			f->cr = $1->cr[i];
			f->ind = $1->ind[i];
			memcpy(f->sym, $1->sym[i], 128*sizeof(char));
			PyTuple_SetItem($result, i,
							SWIG_NewPointerObj((prime_factor*)f,
											   SWIGTYPE_p_prime_factor, 1));
		}
	}

  %extend {
    ~knottypestruct() {
      free($self);
    }

    const plc_knottype *const factors;
  }

} plc_knottype;

%{
    const plc_knottype *knottypestruct_factors_get(plc_knottype *kt) {
        return kt;
    }
%}

%rename(RandomGenerator) gsl_rng;
typedef struct {
  %typemap(out) gsl_rng *get_state {
    $result = PyString_FromStringAndSize(gsl_rng_state($1),
                                         gsl_rng_size($1));
  }
  %apply gsl_rng *get_state { gsl_rng *__getstate__ }
  %typemap(in) (char *stream, size_t items) {
    // TODO: Error checking
    if (!PyString_Check($input)) {
      PyErr_SetString(PyExc_TypeError,"Expecting a string");
      return NULL;
    }

    $1 = PyString_AsString($input);
    $2 = PyString_Size($input);
  }
  %exception from_state {
    $action
      if (!result) {
        PyErr_SetString(PyExc_Exception,"Malformed rng state.");
        return NULL;
      }
  }
  %exception __setstate__ {
    $action
      if (!result) {
        PyErr_SetString(PyExc_Exception,"Malformed rng state.");
        return NULL;
      }
  }

  %extend {
    gsl_rng() {
      gsl_rng *r;
      const gsl_rng_type *T;

      gsl_rng_env_setup();
      T = gsl_rng_default;
      r = gsl_rng_alloc(T);

      return r;
    }
    ~gsl_rng() {
      gsl_rng_free($self);
    }

    void set(unsigned long int s);

    gsl_rng *__getstate__() {
      return $self;
    }
    gsl_rng *get_state() {
      return $self;
    }
    int __setstate__(char *stream, size_t items) {
      size_t n = $self->type->size;
      char *state = (char *)$self->state;

      if (items != n) {
        // Error parsing state. Cleanup. (Throws an exception)
        return 0;
      }

      memcpy(state, stream, items);
      return 1;
    }
    static gsl_rng *from_state(char *stream, size_t items) {
      gsl_rng *r = new_gsl_rng();
      if (!gsl_rng___setstate__(r, stream, items)) {
        // Error parsing state. Cleanup.
        gsl_rng_free(r);
        return NULL;
      }
      return r;
    }
  }

} gsl_rng;
