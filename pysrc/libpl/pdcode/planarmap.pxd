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

    ctypedef struct pmMethod:
        char core, pic
        long seed
        char verbose

    ctypedef struct pmMemory:
        pass

    ctypedef struct pmSize:
        char m, b
        long e, v, f, r, g, d, t
        long *dgArr

    cdef int pmInitRND(pmMethod *method)
    cdef int pmMemoryInit(pmSize *size, pmMethod *method, pmMemory *memory)

cdef extern from "PMplanmap.h":
    cdef int pmSetParameters(pmSize *size, pmMethod *method)
    cdef int pmPlanMap(pmSize *size, pmMethod *method, pmMemory *memory, pmMap *output_map)
    cdef int pmFreeMap(pmMap *map_to_free)
