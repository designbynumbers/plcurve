from .pdisomorphism cimport *
from .isomorphism cimport *
from .diagram cimport PlanarDiagram

cdef class PlanarIsomorphism:
    # cdef pd_iso_t *p ## in .pxd
    def __dealloc__(self):
        if self.p is not NULL:
            pd_free_iso(&self.p)
            self.p = NULL

    def __call__(self, PlanarDiagram pd):
        cdef PlanarDiagram newobj = pd.copy()
        pd_apply_edgemap(newobj.p, self.p.edgemap)
        newobj.regenerate_py_os()
        return newobj

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

    def __repr__(self):
        return "PlanarIsomorphism(%s)"%(
            ",".join("%s->%s%s"%(i, self.p.edgemap.perm.map[i],
                                 "+" if self.p.edgemap.ori[i] else "-") for
                     i in range(self.p.edgemap.perm.n))
        )
