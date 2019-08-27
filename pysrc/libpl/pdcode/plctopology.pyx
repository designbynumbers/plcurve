from .plctopology cimport PD_POS_ORIENTATION, PD_NEG_ORIENTATION, PD_UNSET_ORIENTATION
from .plctopology cimport *
from ..plcurve cimport *
from cython.operator cimport dereference as deref
from libc.stdlib cimport malloc, free

def plc_classify_knot( RandomGenerator rng, PlCurve L):
    cdef int nposs = 1
    cdef plc_knottype * knottypes
    knottypes = plc_classify(rng.p, L.p, &nposs)
    py_knottypes = []
    for i in range(nposs):
        py_knottype = []
        for j in range(knottypes[i].nf):
            py_knottype.append((knottypes[i].cr, knottypes[i].ind))
        py_knottypes.append(py_knottype)

    nf = None
    cr = None
    ind = None
    if knottypes is not NULL:
      nf = int(knottypes.nf)
      cr = list(knottypes.cr)
      ind = list(knottypes.ind)

    #free memory
    free(knottypes)

    #return py_knottypes, int(nposs)
    return nf, cr, ind, int(nposs)

def plc_knot_homfly( RandomGenerator rng, PlCurve L):
    cdef char * ret_string
    ret_string = plc_homfly(rng.p, L.p)
    if ret_string is not NULL:
      return str(ret_string)
    else:
      return None
