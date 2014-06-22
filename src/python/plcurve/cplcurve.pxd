cdef extern from "gsl/gsl_rng.h":
    ctypedef struct gsl_rng:
        pass
    ctypedef struct gsl_rng_type:
        pass

    gsl_rng_type* gsl_rng_default
    void gsl_rng_env_setup()
    gsl_rng* gsl_rng_alloc(gsl_rng_type* T)

cdef extern from "plCurve.h":
    ctypedef struct plc_strand:
        pass
    ctypedef struct plc_constraint:
        pass
    ctypedef struct plc_vert_quant:
        pass
    ctypedef struct plc_symmetry_group:
        pass
    ctypedef struct plCurve:
        int nc
        plc_strand* cp
        plc_constraint* cst
        plc_vert_quant* quant
        plc_symmetry_group* G
    plCurve* plc_random_closed_polygon(gsl_rng* r, int n_edges)
    void plc_free(plCurve* L)
