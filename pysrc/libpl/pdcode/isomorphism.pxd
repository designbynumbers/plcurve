from pdisomorphism cimport pd_iso_t

cdef class PlanarIsomorphism:
    cdef pd_iso_t *p
