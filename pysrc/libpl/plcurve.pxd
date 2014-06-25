from plCurve_h cimport *
from numpy cimport ndarray

cdef plc_vector new_plc_vector(list l)
cdef ndarray plc_vector_as_pyo(plc_vector *v)
cdef ndarray plc_vector_as_pyo_copy(plc_vector *v)

cdef class RandomGenerator:
    cdef gsl_rng* p

cdef class Component:
    cdef plc_strand* p

cdef class PlCurve:
    cdef plCurve* p
    cdef bool own
    cdef list _components
    cdef void generate_py_components(PlCurve self)
