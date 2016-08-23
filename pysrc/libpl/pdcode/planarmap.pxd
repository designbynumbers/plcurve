cdef extern from "PMdef.h":
    ctypedef struct pm_vertex:
        pm_edge *root
        pm_vertex *next_v "next"
        long mark
        short type
        long label

    ctypedef struct pm_edge:
        pm_vertex *from_v "from"
        pm_vertex *r_face "face"
        pm_edge *prev_e "prev"
        pm_edge *next_e "next"
        pm_edge *oppo_e "oppo"
        long mark
        short type
        long label

    ctypedef struct pmMap:
        pm_edge *root
        long e, v, f, i
        long geodesic_distance "geodesicDistance"

    ctypedef struct pmMethod:
        char core, pic
        long seed
        char verbose

    ctypedef struct pmMemory:
        pass

    cdef int PM_MAP_TYPE_QUART_2C
    cdef int PM_MAP_TYPE_QUART_2C_2LEG

    cdef int PM_BASIC_TYPE_QUART_2C

    ctypedef struct pmSize:
        char map_type "m"
        char basic_type "b"
        long n_edges "e"
        long n_verts "v"
        long n_faces "f"

        long min_loop_comps "minLoopComps"
        long max_loop_comps "maxLoopComps"

        long r, g, d
        long t

        long *dgArr

    cdef int pmInitRND(pmMethod *method)
    cdef int pmMemoryInit(pmSize *size, pmMethod *method, pmMemory *memory)
    cdef void pmPrintChndVtx(pm_vertex *vtx)

cdef extern from "PMplanmap.h":
    cdef int pmSetParameters(pmSize *size, pmMethod *method)
    cdef int pmExtendMemory(pmSize *size, pmMethod *method, pmMemory *memory, char OtherReason)
    cdef int pmPlanMap(pmSize *size, pmMethod *method, pmMemory *memory, pmMap *output_map)
    cdef int pmFreeMap(pmMap *map_to_free)
