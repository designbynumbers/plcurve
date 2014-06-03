%module(package="libplcurve") tsmcmc
%feature("autodoc", "1");
%{
#include "stdbool.h"
#include "plcTopology.h"
#include "plcRandomPolygon.h"
#include "pd_multidx.h"
#include "pd_perm.h"
#include "pd_isomorphisms.h"
#include <stddef.h> // SWIG should include this itself but Debian version does not
//static char eoi = 0; // end of iteration exception
static char _exception =0;
%}
%include "typemaps.i"
%include "exception.i"
%import "plCurve.i"

// typedef enum chord_enum {diagonal,edge} chordtype_t;

%rename(Chord) chord_struct;
typedef struct chord_struct {
  int vt[2];    /* An chord joins two vertices. Edges ARE chords. */

} tsmcmc_chord_t;

%rename(Triangle) triangle_struct;
typedef struct triangle_struct {

  int         parent_chord;        /* The parent chord index.                              */
  chordtype_t parent_type;         /* Every chord is either a "diagonal" or "edge" chord.  */

  int         daughter_chord[2];   /* The other two chords in the triangle are daughters.  */
  chordtype_t daughter_type[2];    /* The daughters are either "diagonal" or "edge" chords */
  int         daughter_tri[2];     /* These are indices into the
                                      "triangles" array for each
                                      daughter chord of type
                                      "diagonal", or (-1) for each daughter chord of type "edge". */

  /* The daughter chords are stored in order by the orientation of the triangle so that
     the order goes "parent_diag->daughter_chord[0]->daughter_chord[1]". The chords are
     stored oriented so that the daughter chords are oriented positively on this triangle,
     so the "free" vertex is the second vertex of daughter_chord[0] or the first vertex
     of daughter_chord[1]. */

} tsmcmc_triangle_t;

%rename(Triangulation) triangulation_struct;
typedef struct triangulation_struct {

  int                ntri;      /* This is a master list of (n-2) triangles. */
  tsmcmc_triangle_t *triangles; /* They have an internal tree structure,
                                   but keeping them together in memory
                                   improves performance slightly and
                                   simplifies memory management
                                   greatly. */

  int                ndiags;    /* This is a master list of (n-3) "internal" diagonals. */
  tsmcmc_chord_t    *diags;     /* We use this to index into the array diagonal_lengths. */

  int                nedges;    /* This is a master list of polygon edges. */
  tsmcmc_chord_t    *edges;     /* We use this to index into the array edge_lengths. */

  %extend{
    ~triangulation_struct() {
      tsmcmc_triangulation_free(*$self);
    }

    /* Constructs the fan triangulation for an n-edge polygon. */
    static tsmcmc_triangulation_t new_fan(int nedges)
    { return tsmcmc_fan_triangulation(nedges); }
  }

} tsmcmc_triangulation_t;

  tsmcmc_triangulation_t tsmcmc_spiral_triangulation(int nedges);
  /* Constructs the spiral triangulation which recursively joins adjacent edges with a diagonal. */

  tsmcmc_triangulation_t tsmcmc_random_triangulation(gsl_rng *rng,int nedges);
  /* Constructs a random triangulation of an n-edge polygon. */

  tsmcmc_triangulation_t tsmcmc_teeth_triangulation(int nedges);
  /* Constructs a triangulation by "walking" across the polygon. */

  bool tsmcmc_triangulation_ok(tsmcmc_triangulation_t T);
  /* Self-tests on a triangulation. */

  void tsmcmc_triangulation_print(FILE *outfile,tsmcmc_triangulation_t T);
  /* Prints to a text file to be picked up by Mathematica */

  char *tsmcmc_triangulation_MathematicaForm(tsmcmc_triangulation_t T);
  /* Allocates a string for a Mathematica representation of the triangulation */

  void tsmcmc_triangulation_polymake(FILE *outfile,tsmcmc_triangulation_t T,double *edge_lengths,bool integral);
  /* Outputs the moment polytope of the triangulation in polymake form */

  /*-

    The equilateral polygon sampling algorithms are not direct sampling
    algorithms; but rather produce a Markov chain of polygons which is
    guaranteed to converge.  As a result, there's no direct API for
    accessing the individual polygons in the chain directly- you just
    compute a single (real-valued) function on each polygon and
    integrate it. The code below computes 95% confidence error bars
    for the result of the integral using the "Geyer IPS estimators".

    Notice that this explicitly works (including convergence and
    error bars) for functions like:

    double f(plCurve *L) {

    if (is_trefoil_knot(L)) { return 1.0; }
    else {return 0.0;}

    }

    Note: you can get to the chain itself using plcrwalk from the command line;
    since the chains are usually very large, they need to be written to
    disk and then post-processed from there.

    -*/


  /* These structures deal with internals of the algorithm. It's ok to
     use the default values. */

%rename(RunStats) tsmcmc_run_stats_struct;
typedef struct tsmcmc_run_stats_struct {

  int max_lagged_covariance_used;
  int lagged_covariances_available;

  int dihedral_steps;
  int mp_steps;
  int permute_steps;

  double total_seconds;
  double geyer_ips_seconds;

} tsmcmc_run_stats;

%rename(RunParams) tsmcmc_run_parameters_struct;
typedef struct tsmcmc_run_parameters_struct {

  int    burn_in;
  double delta;                  /* The fraction of permute steps. */
  double beta;                   /* The fraction of moment polytope steps among non-permute steps. */
  int    moment_polytope_repeat; /* Repeat moment polytope steps this many times. */

  int    log_interval;           /* Log data every this many steps */
  char   logfile_name[4096];     /* Filename to use for logging data (if any) */
  FILE   *logfile;               /* Filehandle for log. Set to NULL to turn logging off. */

  %extend{
    static tsmcmc_run_parameters default_unconfined() {
      return tsmcmc_default_unconfined_parameters();
    }
    static tsmcmc_run_parameters default_confined() {
      return tsmcmc_default_confined_parameters();
    }
  }

} tsmcmc_run_parameters;


%{
  typedef struct {
    PyObject *py_integrand;
    PyObject *py_args;
    PyObject *py_kwargs;
  } py_integrand_args;
%}

%typecheck(SWIG_TYPECHECK_POINTER) (double (*integrand)(plCurve *L, void *args), void *args) (py_integrand_args arg_struct) {
  $1 = ((PySequence_Check($input) &&
         (0 < (len = PySequence_Length($input))) && (len) <= 3) ||
        (PyCallable_Check($input))) ? 1 : 0;
}
%typemap(in, doc="Callable integrand", numargs=1) (double (*integrand)(plCurve *L, void *args), void *args) (py_integrand_args arg_struct) {
  Py_ssize_t len;
  PyObject *o;
  arg_struct.py_integrand = NULL;
  arg_struct.py_args = NULL;
  arg_struct.py_kwargs = NULL;
  $2 = &arg_struct;
  // The $input must be a callable or a sequence (which contains a callable!)
  if (PySequence_Check($input) && (0 < (len = PySequence_Length($input))) && (len) <= 3) {
    // A sequence passed as an integrand is of the form
    //  (integrand, [args[, kwargs]])
    // Read the callable
    if (!PyCallable_Check(o = PySequence_GetItem($input, 0))) {
      SWIG_exception(SWIG_TypeError,
                     "Argument `integrand` must be callable or a tuple of the form (callable[, args[, kwargs]]).");
    } else {
      arg_struct.py_integrand = o;
    }
    if (len > 1) {
      // If there's a list of args, take note of it; makes a new ref
      if (PySequence_Check(o = PySequence_GetItem($input, 1))) {
      arg_struct.py_args = PyList_New(len+1);
      // We have to do a little magic since we pass the plCurve as args[0]
      if(-1 == PyList_SetSlice(arg_struct.py_args, 1, len+1, o)) {
        SWIG_exception(SWIG_TypeError,
                   "Something went terribly awry setting up the callback default arguments");
      }
      Py_XDECREF(o); // We don't need this new ref anymore since we copied it all
      } else {
        SWIG_exception(SWIG_TypeError,
                       "Argument `integrand` must be callable or a tuple of the form (callable[, args[, kwargs]]).");
      }
      if (len > 2) {
        // If there's a dict of kwargs, hold on to it
        if (PyDict_Check(o = PySequence_GetItem($input, 2))) {
          arg_struct.py_kwargs = o;
        } else {
          SWIG_exception(SWIG_TypeError,
                         "Argument `integrand` must be callable or a tuple of the form (callable[, args[, kwargs]]).");
        }
      } else {
        arg_struct.py_kwargs = NULL;
      }
    } else {
      arg_struct.py_args = NULL;
      arg_struct.py_kwargs = NULL;
    }
  } else if (PyCallable_Check($input)) {
    // Integrand is a callable
    arg_struct.py_integrand = $input;
    Py_INCREF($input); // kind-of a hack; we promise to decref in freearg
    arg_struct.py_args = NULL;
    arg_struct.py_kwargs = NULL;
  } else {
    SWIG_exception(SWIG_TypeError,
                   "Argument `integrand` must be callable or a tuple of the form (callable[ args[, kwargs]]).");
  }
  $1 = py_integrand_helper;
}
%typemap(freearg) (double (*integrand)(plCurve *L, void *args), void *args) {
  // Free all the python refs we made in the typemap
  if($2) {
    Py_XDECREF(((py_integrand_args *)$2)->py_integrand);
    Py_XDECREF(((py_integrand_args *)$2)->py_args);
    Py_XDECREF(((py_integrand_args *)$2)->py_kwargs);
  }
}
%typemap(in, numinputs=0) double *error_OUT (double temp) {
  $1 = &temp;
}
%typemap(in, numinputs=0) tsmcmc_run_stats *run_stats_OUT (tsmcmc_run_stats temp) {
  $1 = &temp;
}
%typemap(argout) tsmcmc_run_stats *run_stats_OUT {
// It's important that any functions to which this typemap applies NOT return tuples.
// If they do, you will need to rework this (you probably don't want to append to it...)
  PyObject *o, *o2, *o3;
  $1_type real_out = malloc(sizeof($*1_type));
  if ($1) {
    memcpy(real_out, $1, sizeof($*1_type));
    o = SWIG_NewPointerObj(real_out, $1_descriptor, 1);
  } else {
    o = Py_None;
  }
  // If the result is not a PyTuple (this is the FIRST argout; see comment above)
  if (!PyTuple_Check($result)) {
    PyObject *o2 = $result;
    $result = PyTuple_New(1);
    PyTuple_SetItem($result,0,o2);
  }
  o3 = PyTuple_New(1);
  PyTuple_SetItem(o3,0,o);
  o2 = $result;
  $result = PySequence_Concat(o2,o3);
  Py_DECREF(o2);
  Py_DECREF(o3);
}
%typemap(argout) double *error_OUT {
  // It's important that any functions to which this typemap applies NOT return tuples.
  // If they do, you will need to rework this (you probably don't want to append to it...)
  PyObject *o, *o2, *o3;
  o = PyFloat_FromDouble(*$1);

  // If the result is not a PyTuple (this is the FIRST argout; see comment above)
  if (!PyTuple_Check($result)) {
    PyObject *o2 = $result;
    $result = PyTuple_New(1);
    PyTuple_SetItem($result,0,o2);
  }
  o3 = PyTuple_New(1);
  PyTuple_SetItem(o3,0,o);
  o2 = $result;
  $result = PySequence_Concat(o2,o3);
  Py_DECREF(o2);
  Py_DECREF(o3);
 }

%exception py_integrand_helper {
  assert(!_exception);
  $action
    if (_exception) {
      _exception = 0;
      //SWIG_fail;
      SWIG_exception(SWIG_RuntimeError, "There was an error calling the Python callback.");
      SWIG_fail;
    }
}

%inline %{
  double py_integrand_helper(plCurve *L, void *argptr) {
    PyObject *py_plc;
    PyObject *result;
    double ret;
    PyObject *py_integrand;
    PyObject *py_args;
    py_integrand_args *args = (py_integrand_args *)argptr;
    PyObject *err;

    py_integrand = args->py_integrand;
    py_args = args->py_args;
    py_plc = SWIG_InternalNewPointerObj(L, SWIGTYPE_p_plc_type, 0);
    if (py_args) {
      PyList_SetItem(py_args, 0, py_plc);
      result = PyObject_Call(py_integrand, py_args=PyList_AsTuple(py_args), args->py_kwargs);
      Py_DECREF(py_args);
    } else {
      result = PyObject_CallFunctionObjArgs(py_integrand, py_plc, NULL);
      Py_DECREF(py_plc);
    }
    err = PyErr_Occurred();
    if (err != NULL) {
        PyErr_PrintEx(0);
    }
    if (result != NULL) {
      if (!PyFloat_Check(result) && !PyLong_Check(result) && !PyInt_Check(result)) {
        _exception = 1;
        return NAN;
      }
      ret = PyFloat_AsDouble(result);
      Py_DECREF(result);
      return ret;
    } else {
      _exception = 1;
      return NAN;
    }
  }
%}


  char    *tsmcmc_run_stats_MathematicaForm(tsmcmc_run_stats run_stats);
  char    *tsmcmc_run_params_MathematicaForm(tsmcmc_run_parameters run_params);

  /*-
    Master functions for integrating over polygon space.
    -*/

%typemap(check) PyObject *integrand
{
  if (!PyCallable_Check($1)) {
    SWIG_exception(SWIG_TypeError, "Argument `integrand` must be callable");
  }
}

%rename(equilateral_expectation) tsmcmc_equilateral_expectation;
%feature("autodoc", "equilateral_expectation(rng, df, max_steps, max_secs, T, run_params) \
-> result, stats, error")
tsmcmc_equilateral_expectation;
%feature("docstring") tsmcmc_equilateral_expectation
"This is the 'master' driver function for computing an expectation over
equilateral unconfined polygons. It uses the Geyer ips estimator to
compute error bars.

We set parameters for the algorithm with run_params (intended to be
one of the predefined profiles for a run), and return a lot of
detailed information about the run in run_stats (optional, set to NULL
if you don't care).

Note that the number of edges is set (implicitly) by the triangulation.

``df`` is either a callable object (i.e. a [lambda] function) or a tuple of the form
``(callable[, arg_list[, kwarg_dict]])``";
double tsmcmc_equilateral_expectation(gsl_rng *rng,
                                      double (*integrand)(plCurve *L,void *a),
                                      void *args,
                                      int max_steps,
                                      int max_seconds,
                                      tsmcmc_triangulation_t T,
                                      tsmcmc_run_parameters rp,
                                      tsmcmc_run_stats *run_stats_OUT,
                                      double *error_OUT); %{%}

%rename(confined_equilateral_expectation) tsmcmc_confined_equilateral_expectation;
%feature("autodoc", "confined_equilateral_expectation(rng, df, n_edges, max_steps, max_secs, run_params) \
-> result, stats, error")
tsmcmc_confined_equilateral_expectation;
%feature("docstring") tsmcmc_confined_equilateral_expectation
"This is the 'master' driver function for computing an expectation over
equilateral polygons in 'rooted' spherical confinement (the
confinement is 'rooted' when the first vertex of the polygon is at the
center of the sphere). Since only one triangulation is possible, we
don't pass in a triangulation.  It uses the Geyer ips estimator to
compute error bars.

We set the run parameters with the usual run_params struct, but notice
that permutation steps aren't possible, so we ignore delta (if
set). Again, this is usually intended to be one of the predefined
defaults..";
double tsmcmc_confined_equilateral_expectation(gsl_rng *rng,
                                               double (*integrand)(plCurve*,void*),
                                               void *args,
                                               double confinement_radius, int nedges,
                                               int max_steps,int max_seconds,
                                               tsmcmc_run_parameters run_params,
                                               tsmcmc_run_stats *run_stats_OUT,
                                               double *error_OUT); %{%}

%rename(fixed_ftc_expectation) tsmcmc_fixed_ftc_expectation;
%feature("autodoc", "fixed_ftc_expectation(rng, df, close_failure, max_steps, max_secs, T, run_params) \
-> result, stats, error")
tsmcmc_fixed_ftc_expectation;
%feature("docstring") tsmcmc_fixed_ftc_expectation
"This is the 'master' driver function for computing an expectation over
equilateral polygons with a fixed failure to close (one long edge of
length failure_to_close). These can be triangulated any way that is
desired, but the edge from vertex 0 to vertex 1 is still the long
edge. It uses the Geyer ips estimator to compute error bars.

We set the run parameters with the usual ``run_params``
argument. Permutations only permute the equal edge-length steps.

Note that the number of edges is set (implicitly) by the
triangulation.";
double tsmcmc_fixed_ftc_expectation(gsl_rng *rng,double integrand(plCurve *L, void* args),
                                    void *args,
                                    double failure_to_close,
                                    int max_steps,int max_seconds,
                                    tsmcmc_triangulation_t T,
                                    tsmcmc_run_parameters run_params,
                                    tsmcmc_run_stats *run_stats_OUT,
                                    double *error_OUT); %{%}

  /*---*/
