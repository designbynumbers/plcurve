from .plctopology cimport pd_equivalence_t, pd_code_t, pd_uid_t
from .pdisomorphism cimport pd_iso_t

cdef extern from "stdio.h":
    ctypedef struct FILE:
        pass

cdef extern from "pd_storage.h":
    ctypedef struct pd_stor_t:
        pass

    pd_stor_t* pd_new_pdstor()
    void pd_free_pdstor(pd_stor_t **pdstor)

    unsigned int pd_stor_nelts(pd_stor_t *pdstor)

    void pd_copyinto_pdstor(pd_stor_t *pdstor, pd_code_t *pd)
    void pd_copyinto_cass(pd_code_t *pd)

    void pd_addto_pdstor(pd_stor_t *pdstor, pd_code_t *pd, pd_equivalence_t eq)

    pd_code_t *pd_stor_firstelt(pd_stor_t *pdstor)
    pd_code_t *pd_stor_nextelt(pd_stor_t *pdstor)
    pd_code_t *pd_search_pdstor_by_isomorphism(
        pd_stor_t *pdstor, pd_code_t *pd, pd_iso_t ***isos, unsigned int *nisos)
    pd_code_t *pd_search_pdstor_by_hash_uid(
        pd_stor_t *pdstor, char *hash, pd_uid_t uid)

    pd_stor_t *pd_read_pdstor(FILE *stream, pd_equivalence_t eq)

    void pd_write_pdstor(FILE *stream, pd_stor_t *pdstor)
    void pd_start_incremental_pdstor(FILE *stream)
    void pd_addto_incremental_pdstor(FILE *stream, pd_stor_t *pdstor,
                                     unsigned int *nhashes,
                                     unsigned int *nelts)
    void pd_finish_incremental_pdstor(FILE *stream,
                                      unsigned int nhashes,
                                      unsigned int nelts)

    void pd_stor_stats(pd_stor_t *pdstor, unsigned int *nhashes,
                      unsigned int *nelts)
    void pd_display_pdstor(FILE *stream, pd_stor_t *pdstor)
