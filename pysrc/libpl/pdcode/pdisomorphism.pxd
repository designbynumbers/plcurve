from .plctopology cimport pd_idx_t, pd_or_t, pd_code_t

cdef extern from "pd_multidx.h":
    pass
cdef extern from "pd_perm.h":
    cdef struct pd_perm_struct:
        pd_idx_t n
        pd_idx_t *map
        pd_idx_t pc_idx
    ctypedef pd_perm_struct pd_perm_t

cdef extern from "pd_isomorphisms.h":
    cdef struct pd_edgemap_struct:
        pd_perm_t *perm
        pd_or_t *ori "or"
    ctypedef pd_edgemap_struct pd_edgemap_t
    cdef struct pd_crossmap_struct:
        pd_perm_t *perm
        pd_or_t ori "or"
    ctypedef pd_crossmap_struct pd_crossmap_t
    cdef struct pd_facemap_struct:
        pd_perm_t *perm
        pd_or_t ori "or"
    ctypedef pd_facemap_struct pd_facemap_t

    cdef struct pd_iso_struct:
        pd_perm_t *compperm
        pd_edgemap_t *edgemap
        pd_crossmap_t *crossmap
        pd_facemap_t *facemap
    ctypedef pd_iso_struct pd_iso_t

    void pd_apply_edgemap(pd_code_t *pd, pd_edgemap_t *edgemap)

    void pd_free_iso(pd_iso_t **iso)

    cdef pd_iso_t **pd_build_diagram_isotopies(
        pd_code_t *pdA, pd_code_t *pdB, unsigned int *nisos)
    cdef pd_iso_t **pd_build_isos(
        pd_code_t *pdA, pd_code_t *pdB, unsigned int *nisos)
