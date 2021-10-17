from planarmap cimport *
from planarmap cimport pmPrintChndVtx
from libc.stdlib cimport malloc, free

cdef class _PlanarMap:
    """Python wrapper for pmMap object of
    Gilles Schaeffer's PlanarMap software"""

    cdef pmMap *p
    """Internal pointer to C object"""

    def __cinit__(self):
        self.p = <pmMap *>malloc(sizeof(pmMap))

    def __dealloc__(self):
        if self.p is not NULL:
            pmFreeMap(self.p)
            self.p = NULL

    def n_loop_components(self):
        return pmStatGauss(self.p)

    def geodesic_distance(self):
        return pmStatGeodesicDistance(self.p)

    def to_combinatorial_map(self):
        """
        Get the \sigma, \tau cycle permutations representing this map
        combinatorially.

        THIS DOES NOT WORK! planarmap allocates some global variables and the actual
        edges and vertices might be dealloc'ed on this thread!!
        """
        cdef pm_edge *Cur1
        cdef pm_vertex *Vtx

        return
        pmPrintChndVtx(self.p.root.from_v)
        sigma = []
        tau = []

        Vtx = self.p.root.from_v

        Cur1 = Vtx.root
        sigma_cycle = []
        while Cur1 != Vtx.root.prev_e:
            sigma_cycle.append(Cur1.label)
            Cur1 = Cur1.next_e
        sigma_cycle.append(Cur1.label)
        sigma.append(sigma_cycle)

        while Vtx.next_v is not NULL:
            Vtx = Vtx.next_v
            Cur1 = Vtx.root
            sigma_cycle = []
            while Cur1 != Vtx.root.prev_e:
                sigma_cycle.append(Cur1.label)
                Cur1 = Cur1.next_e
            sigma_cycle.append(Cur1.label)
            sigma.append(sigma_cycle)

        return sigma, tau

    @classmethod
    def new_random(cls,
                   n_verts=5,
                   map_type=None,
                   seed=None):
        """Generate a new PlanarMap randomly

        THIS DOES NOT WORK! planarmap allocates some global variables and the actual
        edges and vertices might be dealloc'ed on other threads!!

        This is currently only here to serve as an example to how to use the C
        bindings to planarmap to actually generate a random object. Presently,
        I suggest that you just reimplement this in your more specific case...

        See diagram.pyx for a reference example of this!
        """

        import os

        cdef pmSize size
        cdef pmMethod meth
        cdef pmMemory mem
        cdef pm_edge *cur_e
        cdef pm_vertex *cur_v
        cdef _PlanarMap self = _PlanarMap.__new__(cls)

        ## Initialize "size", which determines what type of map to sample
        size.n_edges = 0
        size.n_verts = n_verts
        size.n_faces = 0

        # We don't currently do anything with these parameters
        size.r = 0 # #"red"
        size.g = 0 # #"green"
        size.d = 0 # something to do with degree on other types of maps
        size.t = -1 # allowed error on size

        size.dgArr = NULL

        size.min_loop_comps = 0
        size.max_loop_comps = 0

        if map_type is None or map_type == 'quart_2c':
            # 2-edge-connected quartic map
            size.map_type = PM_MAP_TYPE_QUART_2C
            size.basic_type = PM_BASIC_TYPE_QUART_2C

        elif map_type == 'quart_4c':
            # 4-edge-connected quartic map
            size.map_type = 5
            size.basic_type = 5

        elif map_type == 'quart_6c':
            # 6-edge-connected quartic map
            size.map_type = 6
            size.basic_type = 5

        elif map_type == 'biquart':
            # bi-quartic map
            size.map_type = 9
            size.basic_type = 9

        elif map_type == 'quart_2c_2leg':
            # 2-edge-connected quartic map with 2 legs
            size.map_type = PM_MAP_TYPE_QUART_2C_2LEG
            size.basic_type = PM_BASIC_TYPE_QUART_2C

        else:
            raise Exception("dia_type is not an expected value")

        ## Initialize some parameters of the sampling method (namely, seed)
        meth.core = 0
        meth.pic = 0
        if seed is None:
            meth.seed = int(os.urandom(5).encode('hex'), 16)
        else:
            meth.seed = seed
        meth.verbose = 0

        if not pmInitRND(&meth):
            raise Exception("Failure during init RND")
        if not pmSetParameters(&size, &meth):
            raise Exception("Failure during set size")

        if not pmMemoryInit(&size, &meth, &mem):
            raise Exception("Failure during memory init")
        if not pmExtendMemory(&size, &meth, &mem, 0):
            raise Exception("Failure during memory extend")

        # Generate the map itself
        if not pmPlanMap(&size, &meth, &mem, self.p):
            raise Exception("Failure during map generation")

        #pmPrintChndVtx(self.p.root.from_v)
        #print self.to_combinatorial_map()

        return self

def new_random_combinatorial(
        n_verts=5,
        min_components=0,
        max_components=0,
        map_type=None,
        seed=None):
    """
    new_random_combinatorial(n_verts=5, map_type=None, seed=None) -> \sigma

    Generate a new vertex cycle \sigma randomly. The edge cycle \tau is
    understood to be (1 2)(3 4)...(2n-1 2n).
    """

    import os

    cdef pmSize size
    cdef pmMethod meth
    cdef pmMemory mem
    cdef pm_edge *cur_e
    cdef pm_vertex *cur_v
    cdef pmMap pmap

    ## Initialize "size", which determines what type of map to sample
    size.n_edges = 0
    size.n_verts = n_verts
    size.n_faces = 0

    # We don't currently do anything with these parameters
    size.r = 0 # #"red"
    size.g = 0 # #"green"
    size.d = 0 # something to do with degree on other types of maps
    size.t = -1 # allowed error on size

    size.dgArr = NULL

    size.min_loop_comps = min_components
    size.max_loop_comps = max_components

    if map_type is None or map_type == 'quart_2c':
        # 2-edge-connected quartic map
        size.map_type = PM_MAP_TYPE_QUART_2C
        size.basic_type = PM_BASIC_TYPE_QUART_2C

    elif map_type == 'quart_4c':
        # 4-edge-connected quartic map
        size.map_type = 5
        size.basic_type = 5

    elif map_type == 'quart_6c':
        # 6-edge-connected quartic map
        size.map_type = 6
        size.basic_type = 5

    elif map_type == 'biquart':
        # bi-quartic map
        size.map_type = 9
        size.basic_type = 9

    elif map_type == 'quart_2c_2leg':
        # 2-edge-connected quartic map with 2 legs
        size.map_type = PM_MAP_TYPE_QUART_2C_2LEG
        size.basic_type = PM_BASIC_TYPE_QUART_2C

    else:
        raise Exception("dia_type is not an expected value")

    ## Initialize some parameters of the sampling method (namely, seed)
    meth.core = 0
    meth.pic = 0
    if seed is None:
        meth.seed = int(os.urandom(5).encode('hex'), 16)
    else:
        meth.seed = seed
    meth.verbose = 0

    if not pmInitRND(&meth):
        raise Exception("Failure during init RND")
    if not pmSetParameters(&size, &meth):
        raise Exception("Failure during set size")

    if not pmMemoryInit(&size, &meth, &mem):
        raise Exception("Failure during memory init")
    if not pmExtendMemory(&size, &meth, &mem, 0):
        raise Exception("Failure during memory extend")

    # Generate the map itself
    if not pmPlanMap(&size, &meth, &mem, &pmap):
        raise Exception("Failure during map generation")

    sigma = []

    cur_v = pmap.root.from_v
    cur_e = cur_v.root
    sigma_cycle = []
    while cur_e != cur_v.root.prev_e:
        sigma_cycle.append(cur_e.label)
        cur_e = cur_e.next_e
    sigma_cycle.append(cur_e.label)
    sigma.append(sigma_cycle)

    while cur_v.next_v is not NULL:
        cur_v = cur_v.next_v
        cur_e = cur_v.root
        sigma_cycle = []
        while cur_e != cur_v.root.prev_e:
            sigma_cycle.append(cur_e.label)
            cur_e = cur_e.next_e
        sigma_cycle.append(cur_e.label)
        sigma.append(tuple(sigma_cycle))

    pmFreeMap(&pmap)

    return tuple(tuple(2*(abs(p)-1) + (1 if p<0 else 0) for p in v) for
                 v in sigma)

def generate_random_maps(
        n=1,
        n_verts=5,
        min_components=0,
        max_components=0,
        map_type=None,
        seed=None,):
    """
    generate_random_maps(n_verts=5, min_components=0, max_components=0,
    map_type=None, seed=None) -> \sigma

    Generate new maps, and pass them to cb
    """

    import os

    cdef pmSize size
    cdef pmMethod meth
    cdef pmMemory mem
    cdef pm_edge *cur_e
    cdef pm_vertex *cur_v
    cdef _PlanarMap pm = _PlanarMap()

    ## Initialize "size", which determines what type of map to sample
    size.n_edges = 0
    size.n_verts = n_verts
    size.n_faces = 0

    # We don't currently do anything with these parameters
    size.r = 0 # #"red"
    size.g = 0 # #"green"
    size.d = 0 # something to do with degree on other types of maps
    size.t = -1 # allowed error on size

    size.dgArr = NULL

    size.min_loop_comps = min_components
    size.max_loop_comps = max_components

    if map_type is None or map_type == 'quart_2c':
        # 2-edge-connected quartic map
        size.map_type = PM_MAP_TYPE_QUART_2C
        size.basic_type = PM_BASIC_TYPE_QUART_2C

    elif map_type == 'quart_4c':
        # 4-edge-connected quartic map
        size.map_type = 5
        size.basic_type = 5

    elif map_type == 'quart_6c':
        # 6-edge-connected quartic map
        size.map_type = 6
        size.basic_type = 5

    elif map_type == 'biquart':
        # bi-quartic map
        size.map_type = 9
        size.basic_type = 9

    elif map_type == 'quart_2c_2leg':
        # 2-edge-connected quartic map with 2 legs
        size.map_type = PM_MAP_TYPE_QUART_2C_2LEG
        size.basic_type = PM_BASIC_TYPE_QUART_2C

    else:
        raise Exception("dia_type is not an expected value")

    ## Initialize some parameters of the sampling method (namely, seed)
    meth.core = 0
    meth.pic = 0
    if seed is None:
        meth.seed = int(os.urandom(5).encode('hex'), 16)
    else:
        meth.seed = seed
    meth.verbose = 0

    if not pmInitRND(&meth):
        raise Exception("Failure during init RND")
    if not pmSetParameters(&size, &meth):
        raise Exception("Failure during set size")

    if not pmMemoryInit(&size, &meth, &mem):
        raise Exception("Failure during memory init")
    if not pmExtendMemory(&size, &meth, &mem, 0):
        raise Exception("Failure during memory extend")

    for _ in xrange(n):
        # Generate the map itself
        if not pmPlanMap(&size, &meth, &mem, pm.p):
            raise Exception("Failure during map generation")

        yield pm

        pmFreeMap(pm.p)

    pm.p = NULL
    print "Done..."
