%module plCurve
%{
#include "plCurve.h"
#include "matrix.h"
#include <gsl/gsl_rng.h>
%}


%inline %{
    double darray_get(double *a, int idx) {
	return a[idx];
    }
    int iarray_get(int *a, int idx) {
	return a[idx];
    }

    gsl_rng *make_gsl_rng() {
	gsl_rng *r;
	const gsl_rng_type *T;

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	return r;
    }

    void free_knottype_struct(plc_knottype *kt) {
	free(kt);
    }
    %}

void gsl_rng_set(const gsl_rng *r, unsigned long int s);
void gsl_rng_free(gsl_rng *r);

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

%include "plCurve.h"
