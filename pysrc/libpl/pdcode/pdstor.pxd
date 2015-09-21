from .plctopology cimport pd_equivalence_t
from .pd_storage cimport pd_stor_t
from .diagram cimport PlanarDiagram

cdef class PDStorage:
    cdef pd_stor_t *p

cdef PDStorage PDStorage_wrap(pd_stor_t *stor)
