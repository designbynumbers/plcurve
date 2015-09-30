from .plctopology cimport pd_code_t

cdef extern from "pd_invariants.h":
    int *pd_interlaced_crossings(pd_code_t *pd)
    unsigned int *pd_interlaced_crossings_unsigned(pd_code_t *pd)
