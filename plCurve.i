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

%{
    static int fail_if_non_py_numeric(PyObject *o, char *fail_msg) {
	if (!PyFloat_Check(o) || !PyInt_Check(o)) {
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
	if (!fail_if_non_py_numeric(n, "Color components must be numerical")) {
	    clr->r = PyFloat_AsDouble(n);
	} else { return 0; }
	// Green
	n = PySequence_GetItem(py_color, 1);
	if (!fail_if_non_py_numeric(n, "Color components must be numerical")) {
	    clr->g = PyFloat_AsDouble(n);
	} else { return 0; }
	// Blue
	n = PySequence_GetItem(py_color, 2);
	if (!fail_if_non_py_numeric(n, "Color components must be numerical")) {
	    clr->b = PyFloat_AsDouble(n);
	} else { return 0; }
	// Alpha
	n = PySequence_GetItem(py_color, 3);
	if (!fail_if_non_py_numeric(n, "Color components must be numerical")) {
	    clr->alpha = PyFloat_AsDouble(n);
	} else { return 0; }

	return 1; // Success!
    }

%}

// Example functions to test plc_vector typemaps (TODO: remove later)
plc_vector plc_build_vect(const double x, const double y, const double z);
plc_vector plc_vect_sum(plc_vector A, plc_vector B);

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

	// Data operation methods
	//
	void add_component(const int add_as, const int nv, const plc_vector *const vt,
			   const bool open, const int cc, const plc_color *const clr) {
	    plc_add_component($self, add_as, nv, open, cc, vt, clr);
	}

	// Topology methods
	//
	plc_knottype *classify(int *nposs) { return plc_classify($self, nposs); }

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
} prime_factor;

%rename(KnotType) knottypestruct;
typedef struct knottypestruct {
    %rename(num_factors) nf;

    int  nf;                            /* Number of prime factors */
    int  cr[MAXPRIMEFACTORS];           /* Crossing number of each prime factor */
    int  ind[MAXPRIMEFACTORS];           /* Index (in Rolfsen or Cerf) of each prime factor */
    char sym[MAXPRIMEFACTORS][128];     /* Symmetry tag (Whitten group element) for each prime factor */
    char homfly[MAXHOMFLY];             /* Homfly polynomial (as plc_lmpoly output) */

    %extend {
	~knottypestruct() {
	    free($self);
	}

	const prime_factor *const factors;
    }

} plc_knottype;

%{
    const prime_factor *knottypestruct_factors_get(plc_knottype *kt) {
	int i;
	prime_factor *pfs = malloc(sizeof(prime_factor)*kt->nf);
	for (i = 0; i < kt->nf; i++) {
	    pfs[i].cr = kt->cr[i];
	    pfs[i].ind = kt->ind[i];
	    memcpy(pfs[i].sym, kt->sym[i], sizeof(char)*128);
	}
	return pfs;
    }
%}

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
