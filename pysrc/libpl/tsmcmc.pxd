from plcurve cimport *

cdef extern from "plcRandomPolygon.h":
    ctypedef struct tsmcmc_run_parameters:
        pass
    ctypedef struct tsmcmc_run_stats:
        pass
    ctypedef struct tsmcmc_triangulation_t:
        pass

    cdef tsmcmc_triangulation_t tsmcmc_fan_triangulation(int nedges)
    cdef void tsmcmc_triangulation_free(tsmcmc_triangulation_t)

    cdef tsmcmc_run_parameters tsmcmc_default_unconfined_parameters()
    cdef tsmcmc_run_parameters tsmcmc_default_confined_parameters()

    cdef double tsmcmc_equilateral_expectation(
        gsl_rng *rng,
        double integrand(plCurve *L, void* userdata), void* args,
        int max_steps, int max_seconds,
        tsmcmc_triangulation_t T,
        tsmcmc_run_parameters run_params,
        tsmcmc_run_stats *run_stats,
        double *error)

    cdef double tsmcmc_confined_equilateral_expectation(
                 gsl_rng *rng,
                 double integrand(plCurve *L, void *args),
                 void *args,
    						 double confinement_radius, int nedges,
    						 int max_steps,int max_seconds,
    						 tsmcmc_run_parameters run_params,
    						 tsmcmc_run_stats *run_stats,
    						 double *error)

cdef class Triangulation:
    cdef tsmcmc_triangulation_t c
cdef class RunParameters:
    cdef tsmcmc_run_parameters c
