cdef extern from "stdio.h":
    ctypedef struct FILE:
        pass

cdef extern from "plcTopology.h":
    cdef enum:
        PD_HASHSIZE = 32
    cdef enum:
        PD_POS_ORIENTATION = 1
        PD_NEG_ORIENTATION = 0
        PD_UNSET_ORIENTATION = 2
    cdef enum:
        PD_NO_ERROR = 0
        PD_NOT_OK = 1
        PD_BAD_FORMAT = 2
        PD_EOF = 3
    cdef extern int PD_VERBOSE
    cdef extern int PD_LIVE_ON_ERROR

    ctypedef bint bool

    ctypedef unsigned int pd_idx_t
    ctypedef unsigned char pd_or_t
    ctypedef unsigned int pd_pos_t
    ctypedef unsigned long int pd_uid_t

    # Edge
    cdef struct pd_edge_struct:
        pd_idx_t head
        pd_pos_t headpos

        pd_idx_t tail
        pd_pos_t tailpos

    ctypedef pd_edge_struct pd_edge_t

    # Component
    cdef struct pd_component_struct:
        pd_idx_t nedges
        pd_idx_t *edge

    ctypedef pd_component_struct pd_component_t

    # Face
    cdef struct pd_face_struct:
        pd_idx_t nedges
        pd_idx_t *edge
        pd_or_t *ori "or"

    ctypedef pd_face_struct pd_face_t

    # Crossing
    cdef struct pd_crossing_struct:
        pd_idx_t edge[4]
        pd_or_t  sign

    ctypedef pd_crossing_struct pd_crossing_t

    # Tangle
    cdef enum pd_boundary_or_t:
        input "in"
        output "out"
        unset

    cdef struct pd_tangle_strand_struct:
        pd_idx_t start_edge
        pd_idx_t end_edge
        pd_idx_t nedges
        pd_idx_t comp
    ctypedef pd_tangle_strand_struct pd_tangle_strand_t

    cdef struct pd_tangle_struct:
        pd_idx_t nedges
        pd_idx_t nstrands

        pd_idx_t *edge
        pd_idx_t *face

        pd_boundary_or_t *edge_bdy_or
        pd_tangle_strand_t *strand

        pd_idx_t ninterior_cross
        pd_idx_t *interior_cross

        pd_idx_t ninterior_edges
        pd_idx_t *interior_edge
    ctypedef pd_tangle_struct pd_tangle_t

    # Diagram
    cdef struct pd_code_struct:
        pd_uid_t uid

        pd_idx_t MAXVERTS
        pd_idx_t MAXEDGES
        pd_idx_t MAXCOMPONENTS
        pd_idx_t MAXFACES

        pd_idx_t ncross
        pd_idx_t nedges
        pd_idx_t ncomps
        pd_idx_t nfaces

        char hash[PD_HASHSIZE]

        pd_edge_t *edge
        pd_component_t *comp
        pd_crossing_t *cross
        pd_face_t *face

    ctypedef pd_code_struct pd_code_t

    pd_code_t *pd_code_new(pd_idx_t MAXVERTS)
    void       pd_code_free(pd_code_t **pd)
    void       pd_code_eltfree(void **PD)

    char pd_print_or(pd_or_t ori)
    pd_or_t pd_compose_or(pd_or_t a,pd_or_t b)
    bool pd_or_ok(pd_or_t ori)
    int pd_or_cmp(const void *A,const void *B)

    pd_crossing_t pd_build_cross(pd_idx_t e0,pd_idx_t e1,pd_idx_t e2,pd_idx_t e3)
    void pd_canonorder_cross(pd_crossing_t *cr, pd_or_t ori)

    void pd_canonorder_face(pd_face_t *face, pd_or_t ori)

    int  pd_cross_cmp(const void *A, const void *B)
    int  pd_face_cmp(const void *A, const void *B)
    int  pd_component_cmp(const void *A, const void *B)

    void pd_component_and_pos(pd_code_t *pd,pd_idx_t edge,
                              pd_idx_t *comp,pd_idx_t *comp_pos)
    void pd_face_and_pos(pd_code_t *pd, pd_idx_t edge,
                         pd_idx_t *posface, pd_idx_t *posface_pos,
                         pd_idx_t *negface, pd_idx_t *negface_pos)

    pd_edge_t pd_oriented_edge(pd_edge_t e,pd_or_t ori)
    void pd_reorient_edge(pd_code_t *pd,pd_idx_t edge,pd_or_t ori)

    # Over/understrand (actually crossing methods)
    void pd_overstrand(
        pd_code_t *pd,pd_idx_t cr,pd_idx_t *incoming_edgenum, pd_idx_t *outgoing_edgenum)
    void pd_understrand(
        pd_code_t *pd,pd_idx_t cr,pd_idx_t *incoming_edgenum, pd_idx_t *outgoing_edgenum)
    void pd_overstrand_pos(
        pd_code_t *pd,pd_idx_t cr,pd_idx_t *incoming_edgepos, pd_idx_t *outgoing_edgepos)
    void pd_understrand_pos(
        pd_code_t *pd,pd_idx_t cr,pd_idx_t *incoming_edgepos, pd_idx_t *outgoing_edgepos)

    # Regenerate pd code data
    void pd_regenerate_crossings(pd_code_t *pd)
    void pd_regenerate_edges(pd_code_t *pd)
    void pd_regenerate_comps(pd_code_t *pd)
    void pd_regenerate_faces(pd_code_t *pd)
    void pd_regenerate_hash(pd_code_t *pd)
    void pd_regenerate(pd_code_t *pd)

    # Validity checks
    bool pd_cross_ok(pd_code_t *pd)
    bool pd_edges_ok(pd_code_t *pd)
    bool pd_faces_ok(pd_code_t *pd)
    bool pd_comps_ok(pd_code_t *pd)
    bool pd_ok(pd_code_t *pd)

    void       pd_write(FILE *outfile,pd_code_t *pd)
    pd_code_t *pd_read_err(FILE *infile, int *err)

    void       pd_write_KnotTheory(FILE *outfile, pd_code_t *pd)
    pd_code_t *pd_read_KnotTheory(FILE *infile)

    void pd_write_c(FILE *outfile, pd_code_t *pd, char *name)

    bool pd_diagram_isotopic(pd_code_t *pdA,pd_code_t *pdB)
    bool pd_isomorphic(pd_code_t *pdA,pd_code_t *pdB)
    bool pd_isomorphic_strings(char *pdcodeA, int nA, char*pdcodeB, int nB)

    pd_code_t *pd_copy(pd_code_t *pd)

    void pd_reorient_component(pd_code_t *pd, pd_idx_t cmp, pd_or_t ori)

    void pd_printf(char *fmt,pd_code_t *pd, ... )
    bool pd_error(char *file, int line, char *fmt, pd_code_t *pd, ...)

    void pd_check_cr(char *file, int line, pd_code_t *pd, pd_idx_t cr)
    void pd_check_edge(char *file, int line, pd_code_t *pd, pd_idx_t edge)
    void pd_check_cmp(char *file, int line, pd_code_t *pd, pd_idx_t cmp)
    void pd_check_face(char *file, int line, pd_code_t *pd, pd_idx_t face)
    void pd_check_notnull(char *file, int line, char *varname, void *ptr)

    # Tangle utilities
    bool pd_tangle_ok(pd_code_t *pd, pd_tangle_t *tangle)
    pd_tangle_t *pd_tangle_new(pd_idx_t nedges)
    void pd_tangle_free(pd_tangle_t **t)
    void pd_regenerate_tangle(pd_code_t *pd, pd_tangle_t *tangle)
    void pd_regenerate_tangle_err(pd_code_t *pd, pd_tangle_t *tangle, int *err)

    # Knot-isomorphic modifications
    pd_code_t* pd_R1_loopdeletion(pd_code_t *pd, pd_idx_t cr)
    void pd_R2_bigon_elimination(pd_code_t *pd,
                                 pd_idx_t cr[2],
                                 pd_idx_t *outpd,
                                 pd_code_t ***outpd)
    void pd_tangle_slide(pd_code_t *pd, pd_tangle_t *tangle,
                         pd_idx_t strand_n_edges,
                         pd_idx_t *strand_edges,
                         pd_idx_t *border_faces,
                         pd_idx_t *n_pieces,
                         pd_code_t ***pd_pieces)
    void pd_tangle_slide_err(pd_code_t *pd, pd_tangle_t *tangle,
                             pd_idx_t strand_n_edges,
                             pd_idx_t *strand_edges,
                             pd_idx_t *border_faces,
                             pd_idx_t *n_pieces,
                             pd_code_t ***pd_pieces,
                             int *err)


    # Builders for typical pd types
    pd_code_t *pd_build_twist_knot(pd_idx_t n)
    pd_code_t *pd_build_torus_knot(pd_idx_t p,pd_idx_t q)
    pd_code_t *pd_build_simple_chain(pd_idx_t n)
    pd_code_t *pd_build_unknot(pd_idx_t n)
    pd_code_t *pd_build_unknot_wye(pd_idx_t a,pd_idx_t b,pd_idx_t c)

    char *pd_homfly(pd_code_t *pdC)

cdef extern char* pdcode_to_ccode(pd_code_t *pdC)
