cimport cplcurve

cpdef long new_rng():
    cplcurve.gsl_rng_env_setup()
    return <long>cplcurve.gsl_rng_alloc(cplcurve.gsl_rng_default)

cdef class PlCurve:
    """
    It's just a PlCurve.
    """
    cdef cplcurve.plCurve* _cptr
    def __cinit__(self):
        pass
    def __dealloc__(self):
        if self._cptr is not NULL:
            # really augmented free
            cplcurve.plc_free(self._cptr)
    cdef _setptr(self, cplcurve.plCurve* ptr):
        self._cptr = ptr

    @classmethod
    def random_closed_polygon(cls, r, N):
        """
        Create a random closed polygon.
        """
        plc = PlCurve.__new__(cls)
        plc._init_random_closed_polygon(r, N)
        return plc

    cpdef _init_random_closed_polygon(self, long r, int N):
        self._setptr(cplcurve.plc_random_closed_polygon(<cplcurve.gsl_rng*>r, N))
