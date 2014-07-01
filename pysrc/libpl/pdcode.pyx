from pdcode cimport *
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.string cimport memcpy
import sys
from cython.operator cimport dereference as deref
cimport cython

cdef class Edge
cdef class Crossing
cdef class Face
cdef class Component
cdef class PlanarDiagram

ctypedef fused Edge_or_idx:
    short
    int
    long
    pd_idx_t
    Edge

cdef extern from "stdio.h":
    ctypedef struct FILE:
        pass

cdef extern from "Python.h":
    cdef FILE* PyFile_AsFile(file_obj)

def ori_char(pd_or_t ori):
    return chr(pd_print_or(ori))


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
    cdef readonly int index
    """The index by which this edge is known in the parent."""

    property head:
        """The index of the crossing towards which this edge points."""
        def __get__(self):
            return self.p.head
    property headpos:
        """The position in the crossing towards which this edge points."""
        def __get__(self):
            return self.p.headpos
    property tail:
        """The index of the crossing from which this edge leaves."""
        def __get__(self):
            return self.p.tail
    property tailpos:
        """The position in the crossing from which this edge leaves."""
        def __get__(self):
            return self.p.tailpos

    property head_tuple:
        """A tuple (vertex, pos) describing where this edge points.

        The position is an index [0..3] in crossing record of head vertex."""
        def __get__(self):
            return (self.p.head, self.p.headpos)
    property tail_tuple:
        """A tuple (vertex, pos) describing where this edge leaves.

        The position is an index [0..3] in crossing record of tail vertex."""
        def __get__(self):
            return (self.p.tail, self.p.tailpos)

    def __cinit__(self):
        self.p = NULL
        self.parent = None

    def __init__(self, PlanarDiagram parent):
        self.parent = parent

    def __dealloc__(self):
        if self.parent is None and self.p is not NULL:
            PyMem_Free(self.p)

    cdef disown(self):
        """Disconnect this data from its parent and copy the data."""
        self.p = <pd_edge_t *>memcpy(PyMem_Malloc(sizeof(pd_edge_t)),
                                     self.p, sizeof(pd_edge_t))
        self.parent = None

    cpdef Edge oriented(self, pd_or_t sign):
        """oriented(sign) -> Edge

        Returns a copy of original edge if ``sign ==
        PD_POS_EDGE_ORIENTATION``, or a new reversed edge if ``sign ==
        PD_NEG_EDGE_ORIENTATION``.

        The ``Edge`` returned is not connected to any
        :py:class:`PlanarDiagram` until it is attached.
        """
        cdef Edge ret = self.__new__(self.__class__)
        cdef pd_edge_t cedge = pd_oriented_edge(deref(self.p), sign)
        ret.p = <pd_edge_t*>memcpy(PyMem_Malloc(sizeof(pd_edge_t)),
                                   &cedge,
                                   sizeof(pd_edge_t))
        ret.parent = None
        return ret

    cpdef component_index_pos(self):
        """component_index_pos() -> (component index, edge pos)

        Returns the component on which this edge resides and this
        edge's position on that component.
        """
        cdef pd_idx_t comp, comp_pos
        if self.parent is None:
            raise NoParentException("This edge is not part of a diagram.")
        pd_component_and_pos(self.parent.p, self.index,
                             &comp, &comp_pos)
        return comp, comp_pos

    cpdef component_pos(self):
        """component_pos() -> (Component, edge pos)

        Returns the component on which this edge resides and this
        edge's position on that component. If you just want the index
        of the component, use :py:meth:`component_index_pos`
        """
        cdef pd_idx_t comp, comp_pos
        comp, comp_pos = self.component_index_pos()
        return self.parent.components[comp], comp_pos

    cpdef face_index_pos(self):
        """face_index_pos() -> ((pos_face, pos), (neg_face, pos))

        Returns the faces on which this edge resides and this
        edge's position on those faces.
        """
        cdef pd_idx_t plus, plus_pos, minus, minus_pos
        if self.parent is None:
            raise NoParentException("This edge is not part of a diagram.")
        pd_face_and_pos(self.parent.p, self.index,
                             &plus, &plus_pos,
                             &minus, &minus_pos)
        return (plus, plus_pos), (minus, minus_pos)

    cpdef face_pos(self):
        """face_pos() -> ((Face plus, pos), (Face minus, pos))

        Returns the Faces on which this edge resides and this
        edge's position on that face. If you just want the indices
        of the faces, use :py:meth:`face_index_pos`
        """
        cdef pd_idx_t plus, plus_pos, minus, minus_pos
        (plus,plus_pos),(minus,minus_pos) = self.face_index_pos()
        return ((self.parent.faces[plus], plus_pos),
                (self.parent.faces[minus], minus_pos))

    cpdef reorient(self, pd_or_t sign):
        """reorient(sign)

        Flips the edge if ``sign == PD_NEG_ORIENTATION``.
        """
        if self.parent is None:
            raise NoParentException("Edge does not belong to a pdcode")
        pd_reorient_edge(self.parent.p, self.index, sign)

    def __str__(self):
        return "%d_%d->%d_%d"%(
            self.p.tail, self.p.tailpos, self.p.head, self.p.headpos)
    def __repr__(self):
        return "Edge(%d_%d->%d_%d)"%(
            self.p.tail, self.p.tailpos, self.p.head, self.p.headpos)

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
            self.parent.p.cross[self.p.tail].edge[(self.p.tailpos+2)%4]]
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
    cdef readonly int index
    """The index by which this component is known in the parent."""

    property nedges:
        """The number of edges in this component. Equivalent to ``len(cmp)``."""
        def __get__(self):
            return self.p.nedges
    property edge:
        """A tuple of edge indices of which this component consists"""
        def __get__(self):
            cdef int i
            return tuple(self.p.edge[i] for i in range(self.p.nedges))

    def __init__(self, PlanarDiagram parent):
        self.parent = parent

    def __len__(self):
        return self.p.nedges

    def __cmp__(self, Component comp):
        return pd_component_cmp(self.p, comp.p)

    def __str__(self):
        return "Component(%s)"%("->".join(str(e) for e in self.edge))
    def __repr__(self):
        return "Component(%s)"%("->".join(str(e) for e in self.edge))

    cdef disown(self):
        """Disconnect this data from its parent and copy the data."""
        self.p = <pd_component_t *>memcpy(PyMem_Malloc(sizeof(pd_component_t)),
                                     self.p, sizeof(pd_component_t))
        self.parent = None

    def __getitem__(self, i):
        if self.parent is None:
            raise NoParentException("This Component is not owned by any PlanarDiagram.")

        if isinstance(i, slice):
            return tuple(self.p.edge[j] for j in range(*i.indices(self.p.nedges)))
        else:
            i = -i if i < 0 else i
            if i >= len(self):
                raise IndexError("Out of bounds error on component access")
            return self.parent.edges[self.p.edge[i]]

cdef class Face(_Disownable):
    cdef pd_face_t* p
    cdef readonly PlanarDiagram parent
    """The :py:class:`PlanarDiagram` to which this face belongs."""
    cdef readonly int index
    """The index by which this face is known in the parent."""

    property nedges:
        """The number of edges in this face"""
        def __get__(self):
            return self.p.nedges
    property edges:
        """A tuple of the *indices* of the edges around this face. The edges
        are in counterclockwise order around the face.

        If you want to access the Python objects, use this Face object
        like a sequence instead (e.g. ``face[4]`` or ``for edge in face:``).

        """
        def __get__(self):
            cdef int i
            return tuple(self.p.edge[i] for i in range(self.p.nedges))
    property signs:
        """A tuple of the signs of the edges around this face. The edges
        are in counterclockwise order around the face.
        """
        def __get__(self):
            cdef int i
            return tuple(self.p.ori[i] for i in range(self.p.nedges))

    def __init__(self, PlanarDiagram parent):
        self.parent = parent

    cdef disown(self):
        """Disconnect this data from its parent and copy the data."""
        self.p = <pd_face_t *>memcpy(PyMem_Malloc(sizeof(pd_face_t)),
                                     self.p, sizeof(pd_face_t))
        self.parent = None

    def __str__(self):
        return "Face(%s)"%("->".join(str(e) for e in self.edges))
    def __repr__(self):
        return "Face(%s)"%("->".join(str(e) for e in self.edges))

    def __cmp__(self, Face face):
        return pd_face_cmp(self.p, face.p)

    def __len__(self):
        return self.p.nedges

    def __getitem__(self, i):
        if self.parent is None:
            raise NoParentException("This Face is not owned by any PlanarDiagram.")

        if isinstance(i, slice):
            return tuple(self.p.edge[j] for j in range(*i.indices(self.p.nedges)))
        else:
            i = -i if i < 0 else i
            if i >= len(self):
                raise IndexError("Index is out of bounds")
            return self.parent.edges[self.p.edge[i]]


cdef class Crossing(_Disownable):
    cdef pd_crossing_t* p
    cdef readonly PlanarDiagram parent
    """The :py:class:`PlanarDiagram` to which this crossing belongs."""
    cdef readonly int index
    """The index by which this crossing is known in the parent."""

    property edges:
        """A tuple of the indices of the edges which connect to this crossing.

        If you want objects instead of indices, access this ``Crossing``
        like a sequence instead.
        """
        def __get__(self):
            return (self.p.edge[0],
                    self.p.edge[1],
                    self.p.edge[2],
                    self.p.edge[3])
    property sign:
        """The orientation of this crossing. The convention used to determine
        sign is this:

        **Positive crossing**: (upper tangent vector) \\\\(\\\\times\\\\)
        (lower tangent vector) points OUT of screen::

                ^
                |
           ----------->
                |
                |

        **Negative crossing**: (upper tangent vector) \\\\(\\\\times\\\\)
        (lower tangent vector) points INTO screen::

                ^
                |
           -----|----->
                |
                |

        You often simply want to know which of the strands (0-2) or
        (1-3) is on top. There are several cases, because it depends
        BOTH on the sign of the crossing and the orientation of the
        strands. It's recommended to use the methods
        :py:meth:`overstrand` and :py:meth:`understrand` to determine
        which is which. Also, the methods :py:meth:`overstrand_pos`
        and :py:meth:`understrand_pos` return the position in this
        crossing (that is, a number in 0..3) of the incoming and
        outgoing edges of the strand going over at this crossing.

        """
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

    def __len__(self):
        return 4
    def __cmp__(self, Crossing cross):
        return pd_cross_cmp(self.p, cross.p)

    cdef PlanarDiagram parent_x(self):
        if self.parent is None:
            raise NoParentException("This Crossing is not owned by any PlanarDiagram")
        return self.parent

    def overstrand(self):
        """overstrand(self) -> (Edge incoming, Edge outgoing)

        Return Edge objects which are the incoming and outgoing edges
        for the strand which passes *over* this crossing.

        If you want indices instead of objects, use
        :py:meth:`overstrand_indices`.
        """
        inc, out = self.overstrand_indices()
        return (self.parent_x().edges[inc], self.parent_x().edges[out])

    cpdef overstrand_indices(self):
        """overstrand_indices(self) -> (in_index, out_index)

        Return indices corresponding to the incoming and outgoing
        edges for the strand which passes *over* this crossing.

        If you want objects instead of indices, use
        :py:meth:`overstrand`.
        """
        cdef pd_idx_t inc,out
        pd_overstrand(self.parent_x().p, self.index, &inc, &out)
        return inc, out

    cpdef overstrand_pos(self):
        """overstrand_pos(self) -> (in_pos, out_pos)

        Return indices around this crossing corresponding to the
        incoming and outgoing edges for the strand which passes *over*
        this crossing.
        """
        cdef pd_idx_t incp,outp
        pd_overstrand_pos(self.parent_x().p, self.index, &incp, &outp)
        return incp, outp

    def understrand(self):
        """understrand(self) -> (Edge incoming, Edge outgoing)

        Return Edge objects which are the incoming and outgoing edges
        for the strand which passes *under* this crossing.

        If you want indices instead of objects, use
        :py:meth:`understrand_indices`.
        """
        inc, out = self.understrand_indices()
        return (self.parent_x().edges[inc], self.parent_x().edges[out])

    cpdef understrand_indices(self):
        """understrand_indices(self) -> (in_index, out_index)

        Return indices corresponding to the incoming and outgoing
        edges for the strand which passes *under* this crossing.

        If you want objects instead of indices, use
        :py:meth:`understrand`.
        """
        cdef pd_idx_t inc,out
        pd_understrand(self.parent_x().p, self.index, &inc, &out)
        return inc, out

    cpdef understrand_pos(self):
        """understrand_pos(self) -> (in_pos, out_pos)

        Return indices around this crossing corresponding to the
        incoming and outgoing edges for the strand which passes *under*
        this crossing.
        """
        cdef pd_idx_t incp,outp
        pd_understrand_pos(self.parent_x().p, self.index, &incp, &outp)
        return incp, outp

    def __str__(self):
        cdef int i
        return "%d.%d.%d.%d%s"%(
            tuple(self.p.edge[i] for i in range(4)) +
            (ori_char(self.p.sign),))
    def __repr__(self):
        cdef int i
        return "Crossing(%d.%d.%d.%d %s)"%(
            tuple(self.p.edge[i] for i in range(4)) +
            (ori_char(self.p.sign),))

# List subclasses which know how to properly delete/set elements within
cdef class _OwnedObjectList(list):
    def __delitem__(self, i):
        if len(self) == 0:
            return

        if isinstance(i, slice):
            for j in range(*i.indices(len(self))):
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

    def __richcmp__(PlanarDiagram self, PlanarDiagram other_pd, int op):
        if op == 2:
            return pd_isomorphic(self.p, other_pd.p)
        elif op == 3:
            return not pd_isomorphic(self.p, other_pd.p)
        else:
            raise NotImplementedError(
                "PlanarDiagrams do not support relative comparisons.")

    def isomorphic(self, PlanarDiagram other_pd):
        """isomorphic(PlanarDiagram other_pd) -> bool

        Returns whether or not this pdcode is isomorphic to the input
        pdcode ``other_pd``. This is equivalent to ``pd == other`` or
        ``not pd != other``.
        """
        return pd_isomorphic(self.p, other_pd.p)

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

    def write(self, f):
        """write(file f)

        Write this pdcode to the open file object ``f``.
        """
        if not ("w" in f.mode or "a" in f.mode or "+" in f.mode):
            raise IOError("File must be opened in a writable mode.")
        pd_write(PyFile_AsFile(f), self.p)

    @classmethod
    def read(cls, f):
        """read(file f) -> PlanarDiagram

        Read the next pdcode from a file object.

        :return: A new :py:class:`PlanarDiagram`, or ``None`` on failure.
        """
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        newobj.p = pd_read(PyFile_AsFile(f))
        if newobj.p is NULL:
            return None
        newobj.regenerate_py_os()
        return newobj

    @classmethod
    def read_all(cls, f):
        """read_all(file f) -> PlanarDiagram

        Returns a generator that iterates through the pdcodes stored
        in ``f``.
        """
        new_pd = cls.read(f)
        while new_pd is not None:
            yield new_pd
            new_pd = cls.read(f)

    @classmethod
    def read_knot_theory(cls, f):
        """read_knot_theory(file f) -> PlanarDiagram

        This function reads a pdcode which was exported from the
        Mathematica package KnotTheory, exported as text with
        something like:

        ``Export["7_2.txt",PD[Knot[7,2]]]``

        These PD codes don't have component or face information, so that is
        all regenerated once the crossings have been loaded from the file.
        This will only read one PD code from ``f``.
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

        Create a new :py:class:`PlanarDiagram` which represents a
        twist knot with \\\\(n\\\\) twists.
        """
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        newobj.p = pd_build_twist_knot(n_twists)
        newobj.regenerate_py_os()
        return newobj
    from_twist_knot = twist_knot # Deprecated

    @classmethod
    def torus_knot(cls, p, q):
        """torus_knot(p, q) -> PlanarDiagram

        Create a new :py:class:`PlanarDiagram` which represents a
        \\\\((p,q)\\\\)-torus knot. Only implemented for \\\\(p=2\\\\).
        """
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

        Create a new :py:class:`PlanarDiagram` which represents an
        \\\\(n\\\\)-link chain.
        """
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        newobj.p = pd_build_simple_chain(n_links)
        newobj.regenerate_py_os()
        return newobj
    from_simple_chain = simple_chain # Deprecated

    @classmethod
    def unknot(cls, n_crossings):
        """unknot(n_crossings) -> PlanarDiagram

        Create a new :py:class:`PlanarDiagram` which represents an
        \\\\(n\\\\)-crossing unlink diagram
        """
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
        cdef int i
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
                e.index = i
                self.edges.append(e)
            for i in range(self.p.ncross):
                c = Crossing(self)
                c.p = &(self.p.cross[i])
                c.index = i
                self.crossings.append(c)
            for i in range(self.p.ncomps):
                m = Component(self)
                m.p = &(self.p.comp[i])
                m.index = i
                self.components.append(m)
            for i in range(self.p.nfaces):
                f = Face(self)
                f.p = &(self.p.face[i])
                f.index = i
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
