ctypedef bint bool

cdef class HOMFLYTerm:
    cdef int C, alpha, zeta
    cpdef bool equals(HOMFLYTerm x, HOMFLYTerm y)

cdef class HOMFLYPolynomial:
    """An [immutable] HOMFLY polynomial object."""
    cdef readonly tuple terms

