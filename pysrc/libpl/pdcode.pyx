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

cdef class _Disownable:
    cdef disown(self):
        raise NotImplementedError("Object must implement a disown method!")

class NoParentException(Exception):
    pass

cdef class Edge(_Disownable):
    """An oriented edge, joining two verts tail -> head"""
    cdef pd_edge_t *p
    cdef readonly PlanarDiagram parent
    """The :py:class:`PlanarDiagram` to which this edge belongs."""

    property head:
        """A tuple (vertex, pos) describing where this edge points.

        The position is an index [0..3] in crossing record of head vertex."""
        def __get__(self):
            return (self.p.head, self.p.headpos)
    property tail:
        """A tuple (vertex, pos) describing where this edge leaves.

        The position is an index [0..3] in crossing record of tail vertex."""
        def __get__(self):
            return (self.p.tail, self.p.tailpos)

    def __cinit__(self):
        pass
    def __init__(self, PlanarDiagram parent):
        self.parent = parent

    def __dealloc__(self):
        pass#if self.parent is None and self.p is not NULL:
        #    PyMem_Free(self.p)

    cdef disown(self):
        """Disconnect this data from its parent and copy the data."""
        self.p = <pd_edge_t *>memcpy(PyMem_Malloc(sizeof(pd_edge_t)),
                                     self.p, sizeof(pd_edge_t))
        self.parent = None

    def __str__(self):
        return "%d_%d->%d_%d"%(self.p.tail, self.p.tailpos, self.p.head, self.p.headpos)
    def __repr__(self):
        return "Edge(%d_%d->%d_%d)"%(self.p.tail, self.p.tailpos, self.p.head, self.p.headpos)

    def prev_crossing(self):
        """prev_crossing() -> Crossing

        Return the crossing from which this edge originates."""
        if self.parent is None:
            raise NoParentException("This Edge is not owned by any PlanarDiagram.")
        return self.parent.crossings[self.p.tail]
    def next_crossing(self):
        """next_crossing() -> Crossing

        Return the crossing towards which this edge points"""
        if self.parent is None:
            raise NoParentException("This Edge is not owned by any PlanarDiagram.")
        return self.parent.crossings[self.p.head]

    def prev_edge(self):
        """prev_edge() -> Edge

        Return the precedent edge along the component"""
        if self.parent is None:
            raise NoParentException("This Edge is not owned by any PlanarDiagram.")
        return self.parent.edges[
            self.parent.p.cross[self.p.tail].edge[(self.p.headpos+2)%4]]
    def next_edge(self):
        """next_edge() -> Edge

        Return the next edge along the component"""
        if self.parent is None:
            raise NoParentException("This Edge is not owned by any PlanarDiagram.")
        return self.parent.edges[
            self.parent.p.cross[self.p.head].edge[(self.p.headpos+2)%4]]

cdef class Component(_Disownable):
    cdef pd_component_t* p
    cdef readonly PlanarDiagram parent
    """The :py:class:`PlanarDiagram` to which this component belongs."""

    property nedges:
        def __get__(self):
            return self.p.nedges

    def __init__(self, PlanarDiagram parent):
        self.parent = parent

    def __len__(self):
        return self.p.nedges

    cdef disown(self):
        """Disconnect this data from its parent and copy the data."""
        self.p = <pd_component_t *>memcpy(PyMem_Malloc(sizeof(pd_component_t)),
                                     self.p, sizeof(pd_component_t))
        self.parent = None

    def __getitem__(self, i):
        if self.parent is None:
            raise NoParentException("This Component is not owned by any PlanarDiagram.")

        i = -i if i < 0 else i
        if i >= len(self):
            raise IndexError("Out of bounds error on component access")
        return self.parent.edges[self.p.edge[i]]

cdef class Face(_Disownable):
    cdef pd_face_t* p
    cdef readonly PlanarDiagram parent
    """The :py:class:`PlanarDiagram` to which this face belongs."""

    property nedges:
        def __get__(self):
            return self.p.nedges

    def __init__(self, PlanarDiagram parent):
        self.parent = parent

    cdef disown(self):
        """Disconnect this data from its parent and copy the data."""
        self.p = <pd_face_t *>memcpy(PyMem_Malloc(sizeof(pd_face_t)),
                                     self.p, sizeof(pd_face_t))
        self.parent = None


cdef class Crossing(_Disownable):
    cdef pd_crossing_t* p
    cdef readonly PlanarDiagram parent
    """The :py:class:`PlanarDiagram` to which this crossing belongs."""

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

    cdef disown(self):
        """Disconnect this data from its parent and copy the data."""
        self.p = <pd_crossing_t *>memcpy(PyMem_Malloc(sizeof(pd_crossing_t)),
                                     self.p, sizeof(pd_crossing_t))
        self.parent = None


    def __cmp__(self, Edge edge):
        return pd_cross_cmp(self.p, edge.p)

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
            for j in xrange(*i.indices(len(self))):
                (<_Disownable>self[j]).disown()
            super(_OwnedObjectList, self).__delitem__(i)
        else:
            (<_Disownable>self[i]).disown()
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
    """Class which represents a PD code (planar diagram). There are
    a few different ways to create ``PlanarDiagram``\ s:

    The constructor will create a new empty diagram with `max_verts`
    blank crossings.

    There are two class methods for reading diagrams from files, both
    from files written in ``PlanarDiagram``'s own format, and from files
    written in Mathematica's KnotTheory format.

    There are multiple class methods which create new ``PlanarDiagram``\ s
    of basic examples; their docstrings contain some specific information."""
    cdef pd_code_t *p
    cdef readonly _EdgeList edges
    """A sequence of :py:class:`Edge` which belong to this diagram"""
    cdef readonly _CrossingList crossings
    """A sequence of :py:class:`Crossing` which belong to this diagram"""
    cdef readonly _ComponentList components
    """A sequence of :py:class:`Component` which belong to this diagram"""
    cdef readonly _FaceList faces
    """A sequence of :py:class:`Face` which belong to this diagram"""

    cdef bool hashed

    property ncomps:
        """Number of components in this PlanarDiagram"""
        def __get__(self):
            return self.p.ncomps
    property nedges:
        """Number of edges in this PlanarDiagram"""
        def __get__(self):
            return self.p.nedges
    property ncross:
        """Number of crossings in this PlanarDiagram"""
        def __get__(self):
            return self.p.ncross
    property nfaces:
        """Number of faces in this PlanarDiagram"""
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
        """__init__([max_verts=15])

        Create a new PlanarDiagram. You may set the initial available number
        of vertices by passing a number to `max_verts`."""
        self.p = pd_code_new(max_verts)

    def copy(self):
        """copy() -> PlanarDiagram

        Returns a memory deepcopy of this PlanarDiagram
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
        """read(file_obj) -> PlanarDiagram

        Read the next pdcode from a file object.

        :return: A new :py:class:`PlanarDiagram`, or ``None`` on failure.
        """
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
        """read_knot_theory(file_obj) -> PlanarDiagram

        This function reads a pdcode from the Mathematica package KnotTheory,
        exported as text with something like:

        ``Export["7_2.txt",PD[Knot[7,2]]]``

        These PD codes don't have component or face information, so that is
        all regenerated once the crossings have been loaded from the file.
        This will only read one PD code per file.
        """
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

    # Methods which build standard PlanarDiagrams
    @classmethod
    def twist_knot(cls, n_twists):
        """twist_knot(n_twists) -> PlanarDiagram

        Create a new :py:class:`PlanarDiagram` which represents a twist knot with \\\\(n\\\\) twists.
        """
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        newobj.p = pd_build_twist_knot(n_twists)
        newobj.regenerate_py_os()
        return newobj
    from_twist_knot = twist_knot # Deprecated

    @classmethod
    def torus_knot(cls, p, q):
        """torus_knot(p, q) -> PlanarDiagram

        Create a new :py:class:`PlanarDiagram` which represents a \\\\((p,q)\\\\)-torus knot.

        *Caveat* Only implemented for \\\\(p=2\\\\)."""
        cdef PlanarDiagram newobj
        if p != 2:
            raise(Exception("torus_knot only implemented for p=2"))

        newobj = PlanarDiagram.__new__(cls)
        newobj.p = pd_build_torus_knot(p,q)
        newobj.regenerate_py_os()
        return newobj
    from_torus_knot = torus_knot # Deprecated

    @classmethod
    def simple_chain(cls, n_links):
        """simple_chain(n_links) -> PlanarDiagram

        Create a new :py:class:`PlanarDiagram` which represents an \\\\(n\\\\)-link chain."""
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        newobj.p = pd_build_simple_chain(n_links)
        newobj.regenerate_py_os()
        return newobj
    from_simple_chain = simple_chain # Deprecated

    @classmethod
    def unknot(cls, n_crossings):
        """unknot(n_crossings) -> PlanarDiagram

        Create a new :py:class:`PlanarDiagram` which represents an \\\\(n\\\\)-crossing unlink diagram"""
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        newobj.p = pd_build_unknot(n_crossings)
        newobj.regenerate_py_os()
        return newobj
    from_unknot = unknot # Deprecated

    @classmethod
    def unknot_wye(cls, a,b,c):
        """unknot_wye(a,b,c) -> PlanarDiagram

        Create a new :py:class:`PlanarDiagram` representing an unknot which is
        designed for hash collisions."""
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        newobj.p = pd_build_unknot_wye(a,b,c)
        newobj.regenerate_py_os()
        return newobj
    from_unknot_wye = unknot_wye # Deprecated

    def reorient_component(self, pd_idx_t component, pd_or_t sign):
        """reorient_component(component, sign)

        Reverse the orientation of component iff sign is :py:const:`PD_NEG_ORIENTATION`
        """
        pd_reorient_component(self.p, component, sign)

    # Regenerate pd code data
    def regenerate_crossings(self):
        """regenerate_crossings()

        Reorders the crossings cyclically to put the
        lowest index edge first and sorts crossings
        in dictionary order based on this reordering.
        """
        pd_regenerate_crossings(self.p)
    def regenerate_components(self):
        """regenerate_components()

        Generates randomly oriented and numbered
        edges from crossings, then strings them
        into components, sorts components by
        size, renumbers edges and orients them along
        components. Updates and regenerates crossings."""
        pd_regenerate_comps(self.p)
    def regenerate_faces(self):
        """regenerate_faces()

        Fills in faces from crossing, edge
        information, including orientation of
        each edge along the face."""
        pd_regenerate_faces(self.p)
    def regenerate_hash(self):
        """regenerate_hash()

        Fill in (printable) hash value from
        comps, faces, and crossings"""
        pd_regenerate_hash(self.p)
    def regenerate(self):
        """regenerate()

        Regenerates everything from cross data."""
        pd_regenerate(self.p)
        self.regenerate_py_os()

    # Validity checks
    def crossings_ok(self):
        """crossings_ok() -> bool

        Crossings are the fundamental data from which
        PD codes are regenerated. This checks that the crossing
        data is sane."""
        return pd_cross_ok(self.p)
    def edges_ok(self):
        """edges_ok() -> bool

        Check that the edge data agrees with crossing data."""
        return pd_edges_ok(self.p)
    def faces_ok(self):
        """faces_ok() -> bool

        Check that the face data agrees with crossing data."""
        return pd_faces_ok(self.p)
    def components_ok(self):
        """components_ok() -> bool

        Check that the component data agrees with crossing data."""
        return pd_comps_ok(self.p)
    def is_ok(self):
        """is_ok() -> bool

        Check that all data contained in this PD code is sane."""
        return pd_ok(self.p)

    def homfly(self):
        """homfly() -> str

        Compute the HOMFLY polynomial for this diagram (returned as string)."""
        return pd_homfly(self.p)

    cdef regenerate_py_os(self):
        cdef Edge e
        cdef Crossing c
        cdef Face f
        cdef Component m
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
