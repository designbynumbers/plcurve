from .pdisomorphism cimport *

cdef class PlanarIsomorphism:
    # cdef pd_iso_t *p ## in .pxd

    property face_maps:
        def __get__(self):
            return tuple(self.p.facemap.perm.map[i] for 
                         i in range(self.p.facemap.perm.n))
    property face_maps_ori:
        def __get__(self):
            return self.p.facemap.ori

    property edge_maps:
        def __get__(self):
            return tuple(self.p.edgemap.perm.map[i] for 
                         i in range(self.p.edgemap.perm.n))
    property edge_maps_ori:
        def __get__(self):
            return tuple(self.p.edgemap.ori[i] for
                         i in range(self.p.edgemap.perm.n))

    property cross_maps:
        def __get__(self):
            return tuple(self.p.crossmap.perm.map[i] for 
                         i in range(self.p.crossmap.perm.n))
    property cross_maps_ori:
        def __get__(self):
            return self.p.crossmap.ori

