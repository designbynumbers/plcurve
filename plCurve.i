%module plCurve
%{
#include "plCurve.h"
#include "matrix.h"
#include <gsl/gsl_rng.h>
%}

%include "plCurve.h"

%inline %{
    double darray_get(double *a, int idx) {
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
    %}
