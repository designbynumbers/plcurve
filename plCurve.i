%module plcurve
%{
#include "plCurve.h"
#include "matrix.h"
#include <gsl/gsl_rng.h>

    // varray is a hack-in to support passing variable length arrays to Python
    typedef struct variable_array_type {
	int len;
	void *array;
    } varray;
%}
%include "carrays.i"
%include "typemaps.i"

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
	if (!PyFloat_Check(o) && !PyInt_Check(o)) {
	    Py_XDECREF(o);
	    PyErr_SetString(PyExc_ValueError,"Expecting a sequence of floats or ints");
	    return NULL;
	}
	v.c[i] = PyFloat_AsDouble(o);
	Py_DECREF(o);
    }
    $1 = v;
}

// Example functions to test plc_vector typemaps (TODO: remove later)
plc_vector plc_build_vect(const double x, const double y, const double z);
plc_vector plc_vect_sum(plc_vector A, plc_vector B);

%rename(Color) plc_color_type;
typedef struct plc_color_type {
    double r;
    double g;
    double b;
    double alpha;
} plc_color;

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

%rename(PlCurve) plc_type;
struct plc_type {
    %rename(num_components) nc;
    int         nc;                              /* Number of components */

    // Hide raw components pointer
    // plc_strand *cp;                              /* Components */
    /*@only@*/ /*@null@*/ plc_constraint *cst;   /* Constraints */
    /*@only@*/ /*@null@*/ plc_vert_quant *quant; /* per-vertex quantifiers */

    plc_symmetry_group *G;                       /* Symmetry group (may be null) */

    %typemap(out) varray components {
	int i, j;
	$result = PyList_New(0);
	for (i = 0; i < $1.len; i++) {

	    PyList_Append($result,
			  SWIG_NewPointerObj(((plc_strand*)$1.array)+i,
					     SWIGTYPE_p_plc_strand_type, 0));
	}
    }

    // SWIG extensions
    %extend {
	// Virtual class members
	const varray components;

	// Constructors and destructor
	//
	plc_type(const int components,
		 const int * const nv,
		 const bool * const open,
		 const int * const cc) { return plc_new(components, nv, open, cc); }
	plc_type(const plCurve * const L) { return plc_copy(L); }
	~plc_type() { plc_free($self); }

	// Random PlCurve creators
	//
	static plCurve *random_closed_polygon(gsl_rng *r, int nEdges)
	{ return plc_random_closed_polygon(r, nEdges); }

	static plCurve *random_open_polygon(gsl_rng *r,int nEdges)
	{ return plc_random_open_polygon(r, nEdges); }

	static plCurve *random_closed_plane_polygon(gsl_rng *r,int nEdges)
	{ return plc_random_closed_plane_polygon(r, nEdges); }

	static plCurve *random_open_plane_polygon(gsl_rng *r,int nEdges)
	{ return plc_random_open_plane_polygon(r, nEdges); }

	static plCurve *random_equilateral_closed_polygon(gsl_rng *r,int nEdges)
	{ return plc_random_equilateral_closed_polygon(r, nEdges); }

	static plCurve *random_equilateral_open_polygon(gsl_rng *r,int nEdges)
	{ return plc_random_equilateral_open_polygon(r, nEdges); }

	static plCurve *loop_closure(gsl_rng *r,int cp,plCurve *openL,int nEdges)
	{ return plc_loop_closure(r, cp, openL, nEdges); }

	// Python special methods
	//
	const char *__str__() {
	    char buf[255];

	    sprintf(buf, "PlCurve with %d components", $self->nc);
	    return buf;
	}
    }
};

%{
    const varray plc_type_components_get(plCurve *curve) {
	varray arr;
	arr.len = curve->nc;
	arr.array = curve->cp;
	return arr;
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

/* %inline %{ */
/*     double darray_get(double *a, int idx) { */
/* 	return a[idx]; */
/*     } */
/*     int iarray_get(int *a, int idx) { */
/* 	return a[idx]; */
/*     } */

%inline %{
    gsl_rng *make_gsl_rng() {
	gsl_rng *r;
	const gsl_rng_type *T;

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	return r;
    }
    %}

/*     void free_knottype_struct(plc_knottype *kt) { */
/* 	free(kt); */
/*     } */
/*     %} */

void gsl_rng_set(const gsl_rng *r, unsigned long int s);
void gsl_rng_free(gsl_rng *r);

/* %typemap(in, numinputs=0) int *nposs (int temp) { */
/*     $1 = &temp; */
/* } */

/* %typemap(argout) int *nposs { */
/*     PyObject *knottype, *np, *o3; */

/*     np = PyInt_FromLong(*$1); */
/*     if(!PyTuple_Check($result)) { */
/* 	knottype = $result; */
/* 	$result = PyTuple_New(1); */
/* 	PyTuple_SetItem($result,0,knottype); */
/*     } */
/*     o3 = PyTuple_New(1); */
/*     PyTuple_SetItem(o3,0,np); */
/*     knottype = $result; */
/*     $result = PySequence_Concat(knottype,o3); */
/*     Py_DECREF(knottype); */
/*     Py_DECREF(o3); */
/*  } */



/* %extend plc_strand_type { */
/*     plc_vector *get_edge(int n) { assert(n < $self->nv); return &($self->vt[n]); } */
/* }; */
