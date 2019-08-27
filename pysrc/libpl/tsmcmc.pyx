from tsmcmc cimport *
from cython.operator cimport dereference as deref

cdef double _py_integrand(plCurve* L, void* userdata):
    cdef tuple pycb_data = <object>userdata
    cdef object cb
    cdef object args
    cdef object kwargs
    cb, args, kwargs = pycb_data
    cdef PlCurve pywrap = PlCurve.__new__(PlCurve)
    pywrap.p = L
    return cb(pywrap, *args, **kwargs)

cdef class Triangulation:
    def __dealloc__(self):
        tsmcmc_triangulation_free(self.c)

    @classmethod
    def fan(cls, int num_edges):
        cdef Triangulation ret = Triangulation.__new__(cls)
        ret.c = tsmcmc_fan_triangulation(num_edges)
        return ret

cdef class RunParameters:
    @classmethod
    def default_unconfined(cls):
        cdef RunParameters ret = RunParameters.__new__(cls)
        ret.c = tsmcmc_default_unconfined_parameters()
        return ret
    @classmethod
    def default_confined(cls):
        cdef RunParameters ret = RunParameters.__new__(cls)
        ret.c = tsmcmc_default_confined_parameters()
        return ret
    

def equilateral_expectation(
        RandomGenerator rng, integrand,
        int max_steps, int max_seconds,
        Triangulation triangulation,
        RunParameters run_params, get_stats=False,
        cb_args=(), cb_kwargs={},):
    cdef double error = 0
    cdef object args = (integrand, cb_args, cb_kwargs)
    return tsmcmc_equilateral_expectation(
        rng.p, _py_integrand, <void*>args,
        max_steps, max_seconds,
        triangulation.c,
        run_params.c,
        NULL, &error
    ), error

def confined_equilateral_expectation(
            RandomGenerator rng, integrand,
            float confinement_radius,
            int nedges,
            int max_steps, int max_seconds,
            RunParameters run_params, get_stats=False,
            cb_args=(), cb_kwargs={},):
        cdef double error = 0
        cdef object args = (integrand, cb_args, cb_kwargs)
        return tsmcmc_confined_equilateral_expectation(
            rng.p, _py_integrand, <void*>args,
            confinement_radius, nedges,
            max_steps, max_seconds,
            run_params.c,
            NULL, &error
        ), error
