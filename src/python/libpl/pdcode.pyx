from pdcode cimport *
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.string cimport memcpy
import sys

cdef extern from "stdio.h":
    ctypedef struct FILE:
        pass

cdef extern from "Python.h":
    cdef FILE* PyFile_AsFile(file_obj)

cdef pd_edge_t* Edge_as_pd_edge_t(Edge e):
    cdef pd_edge_t* p = e.p
    return p

def ori_char(pd_or_t ori):
    return chr(pd_print_or(ori))

cdef class PlanarDiagram

class NoParent(Exception):
    pass

cdef class Edge:
    cdef pd_edge_t *p
    cdef readonly PlanarDiagram parent

    property head:
        def __get__(self):
            return (self.p.head, self.p.headpos)
    property tail:
        def __get__(self):
            return (self.p.tail, self.p.tailpos)

    def __cinit__(self):
        pass
    def __init__(self, PlanarDiagram parent):
        self.parent = parent
    def __dealloc__(self):
        pass#if self.parent is None and self.p is not NULL:
        #    PyMem_Free(self.p)

    cpdef disown(self):
        """Disconnect this data from its parent and copy the data."""
        self.p = <pd_edge_t *>memcpy(PyMem_Malloc(sizeof(pd_edge_t)),
                                     self.p, sizeof(pd_edge_t))
        self.parent = None

    def __str__(self):
        return "%d_%d->%d_%d"%(self.p.tail, self.p.tailpos, self.p.head, self.p.headpos)
    def __repr__(self):
        return "Edge(%d_%d->%d_%d)"%(self.p.tail, self.p.tailpos, self.p.head, self.p.headpos)

    def prev_crossing(self):
        if self.parent is None:
            raise NoParent("This Edge is not owned by any PlanarDiagram.")
        return self.parent.crossings[self.p.tail]
    def next_crossing(self):
        if self.parent is None:
            raise NoParent("This Edge is not owned by any PlanarDiagram.")
        return self.parent.crossings[self.p.head]

    def prev_edge(self):
        if self.parent is None:
            raise NoParent("This Edge is not owned by any PlanarDiagram.")
        return self.parent.edges[
            self.parent.p.cross[self.p.tail].edge[(self.p.headpos+2)%4]]
    def next_edge(self):
        if self.parent is None:
            raise NoParent("This Edge is not owned by any PlanarDiagram.")
        return self.parent.edges[
            self.parent.p.cross[self.p.head].edge[(self.p.headpos+2)%4]]

cdef class Component:
    cdef pd_component_t* p
    cdef readonly PlanarDiagram parent

    property nedges:
        def __get__(self):
            return self.p.nedges

    def __init__(self, PlanarDiagram parent):
        self.parent = parent

    def __len__(self):
        return self.p.nedges

    def __getitem__(self, i):
        if self.parent is None:
            raise NoParent("This Component is not owned by any PlanarDiagram.")

        i = -i if i < 0 else i
        if i >= len(self):
            raise IndexError("Out of bounds error on component access")
        return self.parent.edges[self.p.edge[i]]

cdef class Face:
    cdef pd_face_t* p
    cdef readonly PlanarDiagram parent

    property nedges:
        def __get__(self):
            return self.p.nedges

    def __init__(self, PlanarDiagram parent):
        self.parent = parent

cdef class Crossing:
    cdef pd_crossing_t* p
    cdef readonly PlanarDiagram parent

    property sign:
        def __get__(self):
            return self.p.sign
        def __set__(self, pd_or_t sign):
            self.p.sign = sign

    def __cinit__(self):
        pass
    def __init__(self, PlanarDiagram parent):
        self.parent = parent
    def __dealloc__(self):
        pass

    def __str__(self):
        return "%d.%d.%d.%d%s"%(
            tuple(self.p.edge[i] for i in range(4)) +
            (ori_char(self.p.sign),))
    def __repr__(self):
        return "Crossing(%d.%d.%d.%d %s)"%(
            tuple(self.p.edge[i] for i in range(4)) +
            (ori_char(self.p.sign),))

# List subclasses which know how to properly delete/set elements within
cdef class _OwnedObjectList(list):
    def __delitem__(self, i):
        if len(self) == 0:
            return

        if isinstance(i, slice):
            for j in i.indices(len(self)):
                self[j].disown()
            super(_OwnedObjectList, self).__delitem__(i)
        else:
            self[i].disown()
            super(_OwnedObjectList, self).__delitem__(i)
    def __delslice__(self, i,j):
        del self[max(0,i):max(0,j):]

    cdef delitem(self, i):
        _OwnedObjectList.__delitem__(self, i)
    cdef clearout(self):
        _OwnedObjectList.__delitem__(self, slice(None))

cdef class _EdgeList(_OwnedObjectList):
    def __delitem__(self, i):
        raise NotImplementedError("Can't delete edges yet.")
cdef class _CrossingList(_OwnedObjectList):
    def __delitem__(self, i):
        raise NotImplementedError("Can't delete crossings yet.")
cdef class _ComponentList(_OwnedObjectList):
    def __delitem__(self, i):
        raise NotImplementedError("Can't delete components yet.")
cdef class _FaceList(_OwnedObjectList):
    def __delitem__(self, i):
        raise NotImplementedError("Can't delete faces.")

cdef class PlanarDiagram:
    """
    Class which represents a PD code (planar diagram)."""
    cdef pd_code_t *p
    cdef readonly _EdgeList edges
    cdef readonly _CrossingList crossings
    cdef readonly _ComponentList components
    cdef readonly _FaceList faces

    cdef bool hashed

    property ncomps:
        """Number of components in this PlanarDiagram"""
        def __get__(self):
            return self.p.ncomps
    property nedges:
        def __get__(self):
            return self.p.nedges
    property ncross:
        def __get__(self):
            return self.p.ncross
    property nfaces:
        def __get__(self):
            return self.p.nfaces

    property hash:
        def __get__(self):
            if not self.hashed:
                pd_regenerate_hash(self.p)
            return self.p.hash

    def __cinit__(self):
        self.p = NULL
        self.hashed = False
        self.edges = _EdgeList()
        self.crossings = _CrossingList()
        self.components = _ComponentList()
        self.faces = _FaceList()

    def __init__(self, max_verts=15):
        self.p = pd_code_new(max_verts)

    def copy(self):
        """Returns a memory deepcopy of this PlanarDiagram
        """
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(self.__class__)
        newobj.p = pd_copy(self.p)
        newobj.regenerate_py_os()
        return newobj

    def assert_nonnull(self):
        if self.p is NULL:
            raise Exception("Initialization for this PlanarDiagram is incomplete.")

    @classmethod
    def read(cls, f):
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        open_f = None
        to_close = False
        if not isinstance(f, file):
            open_f = open(f)
            to_close = True
        newobj.p = pd_read(PyFile_AsFile(f))
        if to_close:
            open_f.close()
        if newobj.p is NULL:
            return None
        newobj.regenerate_py_os()
        return newobj

    @classmethod
    def read_knot_theory(cls, f):
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        open_f = None
        to_close = False
        if not isinstance(f, file):
            open_f = open(f)
            to_close = True
        newobj.p = pd_read_KnotTheory(PyFile_AsFile(f))
        if to_close:
            open_f.close()
        if newobj.p is NULL:
            return None
        newobj.regenerate_py_os()
        return newobj

    @classmethod
    def twist_knot(cls, n):
        """Create a new `PlanarDiagram` which represents a twist knot with \$n\$ twists."""
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        newobj.p = pd_build_twist_knot(n)
        newobj.regenerate_py_os()
        return newobj
    from_twist_knot = twist_knot # Deprecated

    @classmethod
    def torus_knot(cls, p, q):
        """
        Create a new `PlanarDiagram` which represents a \$(p,q)\$-torus knot.

        *Caveat* Only implemented for \$p=2\$."""
        cdef PlanarDiagram newobj
        if p != 2:
            raise(Exception("torus_knot only implemented for p=2"))

        newobj = PlanarDiagram.__new__(cls)
        newobj.p = pd_build_torus_knot(p,q)
        newobj.regenerate_py_os()
        return newobj
    from_torus_knot = torus_knot # Deprecated

    @classmethod
    def simple_chain(cls, n):
        """Create a new `PlanarDiagram` which represents an \$n\$-link chain."""
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        newobj.p = pd_build_simple_chain(n)
        newobj.regenerate_py_os()
        return newobj
    from_simple_chain = simple_chain # Deprecated

    @classmethod
    def unknot(cls, n):
        """Create a new `PlanarDiagram` which represents an \$n\$-crossing unlink diagram"""
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        newobj.p = pd_build_unknot(n)
        newobj.regenerate_py_os()
        return newobj
    from_unknot = unknot # Deprecated

    @classmethod
    def unknot_wye(cls, a,b,c):
        """Create a new `PlanarDiagram` representing an unknot which is
        designed for hash collisions."""
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        newobj.p = pd_build_unknot_wye(a,b,c)
        newobj.regenerate_py_os()
        return newobj
    from_unknot_wye = unknot_wye # Deprecated

    def reorient_component(self, pd_idx_t component, pd_or_t sign):
        pd_reorient_component(self.p, component, sign)

    def regenerate_crossings(self):
        pd_regenerate_crossings(self.p)

    def regenerate_components(self):
        pd_regenerate_comps(self.p)

    def regenerate_faces(self):
        pd_regenerate_faces(self.p)

    def regenerate_hash(self):
        pd_regenerate_hash(self.p)

    def regenerate(self):
        pd_regenerate(self.p)

    def is_ok(self):
        return pd_ok(self.p)

    def homfly(self):
        return pd_homfly(self.p)

    cdef regenerate_py_os(self):
        self.edges.clearout()
        self.components.clearout()
        self.crossings.clearout()
        self.faces.clearout()
        if self.p is not NULL:
            for i in range(self.p.nedges):
                e = Edge(self)
                e.p = &(self.p.edge[i])
                self.edges.append(e)
            for i in range(self.p.ncross):
                c = Crossing(self)
                c.p = &(self.p.cross[i])
                self.crossings.append(c)
            for i in range(self.p.ncomps):
                m = Component(self)
                m.p = &(self.p.comp[i])
                self.components.append(m)
            for i in range(self.p.nfaces):
                f = Face(self)
                f.p = &(self.p.face[i])
                self.faces.append(f)
        else:
            raise Exception("Pointer for this PD code was not set")

    def __len__(self):
        return self.p.ncomps

    def __str__(self):
        return ("PlanarDiagram with %d crossings made up of %d components"%(
            self.p.ncross, self.p.ncomps))

    def __repr__(self):
        return ("PlanarDiagram(edges=[" +
                ", ".join(str(e) for e in self.edges) +
                "], cross=[" +
                ", ".join(str(c) for c in self.crossings) +
                "])")

    def __dealloc__(self):
        if self.p is not NULL:
            pd_code_free(&self.p)
            self.p = NULL
