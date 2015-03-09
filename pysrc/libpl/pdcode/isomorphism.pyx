from .pdisomorphism cimport *
from .isomorphism cimport *
from .diagram cimport PlanarDiagram

cdef PlanarIsomorphism_FromPointer(pd_iso_t *c_iso):
    cdef PlanarIsomorphism new_iso = PlanarIsomorphism.__new__(PlanarIsomorphism)
    new_iso.p = c_iso
    return new_iso

cdef class PlanarIsomorphism:
    # cdef pd_iso_t *p ## in .pxd
    def __dealloc__(self):
        if self.p is not NULL:
            pd_free_iso(&self.p)
            self.p = NULL

    def __call__(self, PlanarDiagram pd):
        """Returns the diagram which results from passing the input diagram
        through this isomorphism."""

        cdef PlanarDiagram newobj = pd.copy()
        self.apply(newobj)
        return newobj

    def __mul__(PlanarIsomorphism self, PlanarIsomorphism other):
        """Composition of isomorphisms."""
        return PlanarIsomorphism_FromPointer(pd_compose_isos(self.p, other.p))

    def __imul__(self, PlanarIsomorphism other):
        """In-place composition of isomorphisms."""
        pd_stareq_iso(self.p, other.p)

    def __cmp__(PlanarIsomorphism self, PlanarIsomorphism other):
        cdef int cmpval = pd_iso_cmp(&self.p, &other.p)
        return (0 if cmpval == 0 else 1 if cmpval > 0 else -1)

    def apply(self, PlanarDiagram pd):
        """apply(PlanarDiagram)

        Applies this isomorphism to the input planar diagram *in-place*."""

        pd_apply_edgemap(pd.p, self.p.edgemap)
        pd.regenerate_py_os()

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
