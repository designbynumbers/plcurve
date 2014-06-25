cdef extern from "gsl/gsl_rng.h":
    ctypedef struct gsl_rng:
        pass
    ctypedef struct gsl_rng_type:
        pass

    gsl_rng_type* gsl_rng_default
    void gsl_rng_env_setup()
    gsl_rng* gsl_rng_alloc(gsl_rng_type* T)
    void gsl_rng_set(const gsl_rng *r, unsigned long int seed)
