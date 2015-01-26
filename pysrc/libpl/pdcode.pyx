from pdcode cimport *
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.string cimport memcpy
from libc.stdlib cimport malloc, free
import sys
from cython.operator cimport dereference as deref
import random
import re
from operator import itemgetter, mul
from itertools import islice, izip, cycle
import os
import libpl.data
from libpl.graphs import PlanarSignedFaceDigraph
from libc.stdlib cimport free
from collections import defaultdict

from pdisomorphism cimport pd_iso_t, pd_build_diagram_isotopies
from pdisomorphism cimport PlanarIsomorphism

#from cython.view cimport array
cimport cython

PD_VERBOSE = 10
#PD_LIVE_ON_ERROR = 1

cdef class Edge
cdef class Crossing
cdef class Face
cdef class Component

DEFAULT_PATH=os.path.join("data","pdstors")
SOURCE_DIR=libpl.data.dir

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

cdef bytes copy_and_free(char* charp):
    cdef bytes ret
    try:
        ret = <bytes>charp
    finally:
        free(charp)
    return ret



def ori_char(pd_or_t ori):
    return chr(pd_print_or(ori))

cdef sgn(x):
    if x == 0: return 0
    return x / abs(x)

cdef class _Disownable:
    cdef disown(self):
        raise NotImplementedError

class NoParentException(Exception):
    pass

cdef class Edge(_Disownable):
    """An oriented edge, joining two verts tail -> head

    Once you are finished setting any of the properties ``head``,
    ``headpos``, ``tail``, and ``tailpos``, you **must** call
    :py:meth:`PlanarDiagram.regenerate` on the parent so that sanity
    can be checked and faces and components regenerated.
    """
    cdef pd_edge_t *p
    cdef readonly PlanarDiagram parent
    """The :py:class:`PlanarDiagram` to which this edge belongs."""
    cdef readonly int index
    """The index by which this edge is known in the parent."""

    property head:
        """The index of the crossing towards which this edge points."""
        def __get__(self):
            return self.p.head
        def __set__(self, pd_idx_t i):
            self.p.head = i
    property headpos:
        """The position in the crossing towards which this edge points."""
        def __get__(self):
            return self.p.headpos
        def __set__(self, pd_idx_t i):
            self.p.headpos = i
    property tail:
        """The index of the crossing from which this edge leaves."""
        def __get__(self):
            return self.p.tail
        def __set__(self, pd_idx_t i):
            self.p.tail = i
    property tailpos:
        """The position in the crossing from which this edge leaves."""
        def __get__(self):
            return self.p.tailpos
        def __set__(self, pd_idx_t i):
            self.p.tailpos = i

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
        #return
        if self.parent is None and self.p is not NULL:
            pass
            #PyMem_Free(self.p)

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

    cpdef plug_head(self, Crossing x, pd_idx_t pos):
        if self.parent is not None and self.parent.p != x.parent.p:
            raise Exception("Cannot plug this edge into a different diagram")
        elif self.parent is None:
            raise Exception("Edge is not part of a PlanarDiagram")
        x.edges[pos] = self.index
        self.head = x.index
        self.headpos = pos

    cpdef plug_tail(self, Crossing x, pd_idx_t pos):
        if self.parent is not None and self.parent.p != x.parent.p:
            raise Exception("Cannot plug this edge into a different diagram")
        elif self.parent is None:
            raise Exception("Edge is not part of a PlanarDiagram")
        x.edges[pos] = self.index
        self.tail = x.index
        self.tailpos = pos

    cpdef swap_head(self, Edge other, pd_or_t head_or_tail):
        cdef Crossing x
        cdef pd_idx_t pos
        if head_or_tail == PD_NEG_ORIENTATION:
            # Plug into other's tail crossing
            pos = self.headpos
            x = self.next_crossing()
            self.plug_head(other.prev_crossing(), other.tailpos)
            other.plug_tail(x, pos)
        elif head_or_tail == PD_POS_ORIENTATION:
            # Plug into other's head crossing
            pos = self.headpos
            x = self.next_crossing()
            self.plug_head(other.next_crossing(), other.headpos)
            other.plug_head(x, pos)
        else:
            raise Exception("Invalid orientation describing head-or-tail")

    cpdef swap_tail(self, Edge other, pd_or_t head_or_tail):
        cdef Crossing x
        cdef pd_idx_t pos
        if head_or_tail == PD_NEG_ORIENTATION:
            # Plug into other's tail crossing
            pos = self.tailpos
            x = self.prev_crossing()
            self.plug_tail(other.prev_crossing(), other.tailpos)
            other.plug_tail(x, pos)
        elif head_or_tail == PD_POS_ORIENTATION:
            # Plug into other's head crossing
            pos = self.tailpos
            x = self.prev_crossing()
            self.plug_tail(other.next_crossing(), other.headpos)
            other.plug_head(x, pos)
        else:
            raise Exception("Invalid orientation describing head-or-tail")

    def __str__(self):
        return "%d_%d->%d_%d"%(
            self.p.tail, self.p.tailpos, self.p.head, self.p.headpos)
    def __repr__(self):
        return "Edge(\"%d_%d\",\"%d_%d\")"%(
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
        return sgn(pd_component_cmp(self.p, comp.p))

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
    property vertices:
        """A tuple of the crossing *indices* around this face in counterclockwise
        order.
        """
        def __get__(self):
            return tuple((edge.head if sign == PD_NEG_ORIENTATION else
                          edge.tail) for
                         edge, sign in zip(self, self.signs))

    def __init__(self, PlanarDiagram parent):
        self.parent = parent

    cdef disown(self):
        """Disconnect this data from its parent and copy the data."""
        self.p = <pd_face_t *>memcpy(PyMem_Malloc(sizeof(pd_face_t)),
                                     self.p, sizeof(pd_face_t))
        self.parent = None

    def get_vertices(self):
        return tuple(self.parent.crossings[vert] for vert in self.vertices)

    def __str__(self):
        return "Face(%s)"%("->".join(str(e) for e in self.edges))
    def __repr__(self):
        return "Face(%s)"%("->".join(str(e) for e in self.edges))

    def __cmp__(self, Face face):
        return sgn(pd_face_cmp(self.p, face.p))

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

    def __hash__(self):
        return self.index

    def __richcmp__(self, Face otherface, int op):
        if op == 2:
            if len(self.parent.crossings) == len(otherface.parent.crossings):
               truth_table = [self.parent.crossings[i].edges == \
               otherface.parent.crossings[i].edges \
               for i in range(0,len(self.parent.crossings))]
               if self.index == otherface.index and not False in truth_table:
                    return True

            return False
        elif op == 3:
            if len(self.parent.crossings) == len(otherface.parent.crossings):
                truth_table = [self.parent.crossings[i].edges == \
                otherface.parent.crossings[i].edges \
                for i in range(0,len(self.parent.crossings))]

                if self.index == otherface.index and not False in truth_table:
                    return False

            return True

        else:
            raise NotImplementedError(
                "PlanarDiagram Faces do not support relative comparisons.")


cdef class Crossing(_Disownable):
    """A crossing, which holds four edges in positions and an orientation.

    Crossing data is the **most important** piece of data in a
    PlanarDiagram. A full diagram object can be built from crossing
    data alone (although orientations of components may differ).

    Once you are finished setting any of the properties ``head``,
    ``headpos``, ``tail``, and ``tailpos``, you **must** call
    :py:meth:`PlanarDiagram.regenerate` on the parent so that sanity
    can be checked and faces and components regenerated.

    """

    cdef pd_crossing_t* p
    cdef readonly PlanarDiagram parent
    """The :py:class:`PlanarDiagram` to which this crossing belongs."""
    cdef readonly int index
    """The index by which this crossing is known in the parent."""

    cdef pd_idx_t[:] edgeview_get(self):
        return <pd_idx_t[:4]>self.p.edge
    property edges:
        """A tuple of the indices of the edges which connect to this crossing.

        If you want objects instead of indices, access this ``Crossing``
        like a sequence instead.
        """
        def __get__(self):
            return self.edgeview_get()
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

    def toggle_sign(self):
        if self.p.sign == PD_NEG_ORIENTATION:
            self.p.sign = PD_POS_ORIENTATION
        elif self.p.sign == PD_POS_ORIENTATION:
            self.p.sign = PD_NEG_ORIENTATION

    property faces:
        """A tuple of the indices of the faces around this vertex, in clockwise order.
        """
        def __get__(self):
            ret = [0]*4
            for i, edge in enumerate(self):
                (pos_face,_),(neg_face,_) = edge.face_index_pos()
                if edge.tail == self.index and edge.tailpos == i:
                    # Edge departed from this crossing
                    ret[i] = pos_face
                else:
                    ret[i] = neg_face
            return tuple(ret)

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
        return sgn(pd_cross_cmp(self.p, cross.p))
    def __getitem__(self, i):
        if self.parent is None:
            raise NoParentException("This Crossing is not owned by any PlanarDiagram.")

        if isinstance(i, slice):
            return tuple(self.p.edge[j] for j in range(*i.indices(4)))
        else:
            i = -i if i < 0 else i
            if i >= 4:
                raise IndexError("Index is out of bounds")
            return self.parent.edges[self.p.edge[i]]

    cdef PlanarDiagram parent_x(self):
        if self.parent is None:
            raise NoParentException("This Crossing is not owned by any PlanarDiagram")
        return self.parent

    def adjacent(self, pos):
        edge = self[pos]
        if (self.index, pos) == edge.head_tuple:
            return edge.prev_crossing(), edge.tailpos
        else:
            return edge.next_crossing(), edge.headpos

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

    def get_faces(self):
        """get_faces() -> tuple of Faces

        Return the faces on which this crossing is a vertex, in
        clockwise order.
        """
        return tuple(self.parent.faces[i] for i in self.faces)

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

cdef class _VirtualEdgeArray:
    """Discreetly wraps an edge array so that Python code can mutate C data."""
    pass

_latex_term_regex = re.compile("(-?[0-9]*)(a\\^\\{-?[0-9]+\\})?(z\\^\\{-?[0-9]+\\})?")
_latex_coef_regex = re.compile("[az]\\^\\{(-?[0-9]+)\\}")

cdef class HOMFLYTerm:
    cdef int C, alpha, zeta
    def __cinit__(self, int C, int alpha, int zeta):
        self.C = C
        self.alpha = alpha
        self.zeta = zeta

    def __invert__(self):
        return HOMFLYTerm(
            self.C,
            -self.alpha,
            self.zeta)
    def __mul__(HOMFLYTerm x, HOMFLYTerm y):
        return HOMFLYTerm(
            x.C * y.C,
            x.alpha + y.alpha,
            x.zeta + y.zeta)

    def __add__(HOMFLYTerm x, HOMFLYTerm y):
        if x.alpha != y.alpha or x.zeta != y.zeta:
            raise Exception("Terms do not match for summation")
        return HOMFLYTerm(
            x.C + y.C,
            x.alpha, x.zeta)

    def __nonzero__(self):
        if self.C == 0:
            return 0
        else:
            return 1
    def __cmp__(HOMFLYTerm x, HOMFLYTerm y):
        if x.alpha < y.alpha:
            return -1
        elif x.alpha == y.alpha:
            if x.zeta < y.zeta:
                return -1
            elif x.zeta == y.zeta:
                return 0
            else:
                return 1
        else:
            return 1
    cpdef bool equals(HOMFLYTerm x, HOMFLYTerm y):
        return x.alpha == y.alpha and x.zeta == y.zeta and x.C == y.C
    def __hash__(self):
        return hash((self.alpha, self.zeta))

    def __reduce__(self):
        """Required method to implement pickling of extension types.

        Returns a tuple: (Constructor-ish, Constructor args).
        """
        return (self.__class__, (self.C, self.alpha, self.zeta))


    def __str__(self):
        cdef str cterm = "-" if self.C == -1 else ("" if self.C == 1 else str(self.C))
        if self.alpha and self.zeta:
            return "%sa^{%s}z^{%s}"%(cterm, self.alpha, self.zeta)
        elif self.alpha:
            return "%sa^{%s}"%(cterm, self.alpha)
        elif self.zeta:
            return "%sz^{%s}"%(cterm, self.zeta)
        else:
            return "%s"%self.C
    def __repr__(self):
        cdef str cterm = "-" if self.C == -1 else ("" if self.C == 1 else str(self.C))
        if self.alpha and self.zeta:
            return "%sa^{%s}z^{%s}"%(cterm, self.alpha, self.zeta)
        elif self.alpha:
            return "%sa^{%s}"%(cterm, self.alpha)
        elif self.zeta:
            return "%sz^{%s}"%(cterm, self.zeta)
        else:
            return "%s"%self.C

cdef class HOMFLYPolynomial:
    """An [immutable] HOMFLY polynomial object."""
    cdef readonly tuple terms

    def __cinit__(self):
        self.terms = None
    def __init__(self, data, bool sort=True):
        """Create a new HOMFLY polynomial from a LaTeX representation"""
        # let's do this pythonically now and then improve it later
        cdef str latex

        if isinstance(data, basestring):
            latex = data
            self._load_from_latex(latex, sort=sort)
        elif getattr(data, '__iter__', False):
            self._load_from_terms(data, sort=sort)

    def _load_from_latex(self, str latex, bool sort=True):
        cdef object coeff, a_part, z_part, term
        cdef list result = []
        # change '- x' into '+ -x' for consistency
        latex = latex.replace("- ", "+ -")
        for term in latex.split(" + "):
            coeff, a_part, z_part = _latex_term_regex.match(term).groups()
            if coeff == "-":
                coeff = -1
            elif coeff == "":
                coeff = 1
            if a_part:
                a_part = _latex_coef_regex.match(a_part).group(1)
            else:
                a_part = 0
            if z_part:
                z_part = _latex_coef_regex.match(z_part).group(1)
            else:
                z_part = 0

            term = HOMFLYTerm(int(coeff),
                              int(a_part),
                              int(z_part))
            if term:
                result.append(term)
        # TODO: Raise error on parse failure or inappropriate latex rep

        if sort:
            self.terms = tuple(sorted(result))
        else:
            self.terms = tuple(result)

    def _load_from_terms(self, terms, bool sort=True):
        homfly_terms = []
        for term in terms:
            if not isinstance(term, HOMFLYTerm):
                term = HOMFLYTerm(*term)
            try:
                homfly_terms[homfly_terms.index(term)] += term
            except ValueError:
                homfly_terms.append(term)

        if sort:
            self.terms = tuple(sorted(homfly_terms))
        else:
            self.terms = tuple(homfly_terms)

    def __getitem__(self, x):
        try:
            return self.terms[x]
        except IndexError:
            raise
    def __len__(self):
        return len(self.terms)

    def __invert__(HOMFLYPolynomial self):
        cdef HOMFLYPolynomial newpoly = HOMFLYPolynomial.__new__(HOMFLYPolynomial)
        newpoly.terms = tuple(sorted((~A for A in self.terms)))
        return newpoly
    def __mul__(HOMFLYPolynomial x, HOMFLYPolynomial y):
        cdef list newterms = list()
        cdef HOMFLYPolynomial newpoly = HOMFLYPolynomial.__new__(HOMFLYPolynomial)
        cdef HOMFLYTerm A,B,C
        cdef dict found = dict()

        # There is room to replace this loop & sort with an O(n^2) solution
        for A in x.terms:
            for B in y.terms:
                C = A*B
                n = (C.alpha, C.zeta)
                if not C:
                    continue
                if n in found:
                    newterms[found[n]] = newterms[found[n]] + C
                else:
                    found[n] = len(newterms)
                    newterms.append(C)

        # remove any 0-coefficient terms... can possibly do this inline
        newpoly.terms = tuple(sorted(filter(None, newterms)))
        return newpoly
    def __pow__(HOMFLYPolynomial self, int k, z):
        if k == 0:
            return HOMFLYPolynomial('1')
        elif k < 0:
            raise Exception("Power must be positive")
        else:
            return reduce(mul, [self]*k)

    def __richcmp__(HOMFLYPolynomial x, HOMFLYPolynomial y, int op):
        cdef HOMFLYTerm A,B
        if op == 2:
            if not isinstance(x, HOMFLYPolynomial) or not isinstance(y, HOMFLYPolynomial):
                return False
            if len(x.terms) != len(y.terms):
                return False
            return not (False in (A.equals(B) for A,B in zip(x.terms, y.terms)))
        elif op == 3:
            if not isinstance(x, HOMFLYPolynomial) or not isinstance(y, HOMFLYPolynomial):
                return True
            if len(x.terms) != len(y.terms):
                return True
            return False in (A.equals(B) for A,B in zip(x.terms, y.terms))
        else:
            raise NotImplementedError(
                "Only equality checking implemented for HPs")

    def __reduce__(self):
        return (self.__class__, (self.terms,))

    def __hash__(self):
        return hash(self.terms)
    def __str__(self):
        return " + ".join(str(t) for t in self.terms)
    def __repr__(self):
        return " + ".join(repr(t) for t in self.terms)


cdef PlanarDiagram PlanarDiagram_wrap(pd_code_t *pd, bool thin=False):
    """Wraps a pd_code_t* in a PlanarDiagram"""
    cdef PlanarDiagram newobj = PlanarDiagram.__new__(PlanarDiagram)
    newobj.thin = thin
    if pd is NULL:
        return None
    newobj.p = pd
    newobj.regenerate_py_os()
    return newobj

cdef class PlanarDiagram:

    property ncomps:
        """Number of components in this PlanarDiagram"""
        def __get__(self):
            return self.p.ncomps
        def __set__(self, n):
            self.p.ncomps = n
            self.regenerate_py_os()
    property nedges:
        """Number of edges in this PlanarDiagram"""
        def __get__(self):
            return self.p.nedges
        def __set__(self, n):
            self.p.nedges = n
            self.regenerate_py_os()
    property ncross:
        """Number of crossings in this PlanarDiagram"""
        def __get__(self):
            return self.p.ncross
        def __set__(self, n):
            self.p.ncross = n
            self.regenerate_py_os()
    property nfaces:
        """Number of faces in this PlanarDiagram"""
        def __get__(self):
            return self.p.nfaces
        def __set__(self, n):
            self.p.nfaces = n
            self.regenerate_py_os()

    property hash:
        def __get__(self):
            if not self.hashed:
                self.hashed = True
                pd_regenerate_hash(self.p)
            return self.p.hash

    property _hash:
        def __get__(self):
            return self.p.hash

    def __cinit__(self):
        self.p = NULL
        self.thin = False
        self.hashed = False
        self._homfly = NULL
        self.edges = _EdgeList()
        self.crossings = _CrossingList()
        self.components = _ComponentList()
        self.faces = _FaceList()

    def __init__(self, max_verts=15):
        """__init__([max_verts=15])

        Create a new PlanarDiagram. You may set the initial available number
        of vertices by passing a number to `max_verts`."""
        self.p = pd_code_new(max_verts)

    @classmethod
    def _wrap(cls, object pdp, bool thin=False):
        """Should be cdef but can't because classmethod. If you must call from
        Python somehow, use this. Otherwise, use PlanarDiagram_wrap"""
        cdef pd_code_t *pd = <pd_code_t*>pdp
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        newobj.p = pd
        newobj.thin = thin
        newobj.regenerate_py_os()
        return newobj

    def __richcmp__(PlanarDiagram self, PlanarDiagram other_pd, int op):
        if op == 2:
            return pd_diagram_isotopic(self.p, other_pd.p)
        elif op == 3:
            return not pd_diagram_isotopic(self.p, other_pd.p)
        else:
            raise NotImplementedError(
                "PlanarDiagrams do not support relative comparisons.")

    def isotopic(self, PlanarDiagram other_pd):
        """isotopic(PlanarDiagram other_pd) -> bool

        Returns whether or not this pdcode is diagram-isotopic to the
        input pdcode ``other_pd``. This is equivalent to ``pd ==
        other`` or ``not pd != other``.
        """
        return pd_diagram_isotopic(self.p, other_pd.p)

    def build_isotopies(self, PlanarDiagram other_pd):
        cdef unsigned int nisos
        cdef pd_iso_t **isos

        isos = pd_build_diagram_isotopies(self.p, other_pd.p, &nisos)
        ret = tuple(PlanarIsomorphism() for _ in range(nisos))
        for i in range(nisos):
            (<PlanarIsomorphism>ret[i]).p = isos[i]

        return ret

    def isomorphic(self, PlanarDiagram other_pd):
        """isomorphic(PlanarDiagram other_pd) -> bool

        Returns whether or not this pdcode shadow is isomorphic to the
        input pdcode ``other_pd``.
        """
        return pd_isomorphic(self.p, other_pd.p)

    def copy(self, thin=False):
        """copy() -> PlanarDiagram

        Returns a memory deepcopy of this PlanarDiagram
        """
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(self.__class__)
        newobj.p = pd_copy(self.p)
        newobj.thin = thin
        newobj.regenerate_py_os()
        return newobj

    def assert_nonnull(self):
        if self.p is NULL:
            raise Exception("Initialization for this PlanarDiagram is incomplete.")

    def _bounds(self):
        print "MAXVERTS=%s" % self.p.MAXVERTS
        print "MAXEDGES=%s" % self.p.MAXEDGES
        print "MAXCOMPONENTS=%s" % self.p.MAXCOMPONENTS
        print "MAXFACES=%s" % self.p.MAXFACES

    def write(self, f):
        """write(file f)

        Write this pdcode to the open file object ``f``.
        """
        if not ("w" in f.mode or "a" in f.mode or "+" in f.mode):
            raise IOError("File must be opened in a writable mode.")
        pd_write(PyFile_AsFile(f), self.p)

    def write_knot_theory(self, f):
        """write_knot_theory(file f)

        Write this pdcode in KnotTheory format to the open file object ``f``.
        """
        if not ("w" in f.mode or "a" in f.mode or "+" in f.mode):
            raise IOError("File must be opened in a writable mode.")
        pd_write_KnotTheory(PyFile_AsFile(f), self.p)

    @classmethod
    def read(cls, f, read_header=False, thin=False):
        """read(file f) -> PlanarDiagram

        Read the next pdcode from a file object.

        :return: A new :py:class:`PlanarDiagram`, or ``None`` on failure.
        """
        cdef int err = 0
        if read_header == True:
            # Actually parse header and check validity
            f.readline()
            f.readline()
            f.readline()

        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        newobj.thin = thin
        newobj.p = pd_read_err(PyFile_AsFile(f), &err)
        if err:
            if err == PD_NOT_OK:
                raise Exception("Loaded PD code is not valid")
            elif err == PD_BAD_FORMAT:
                raise Exception("Error in pdcode file format data")
            elif err == PD_EOF:
                raise EOFError("Already at end of pdcode file")
        if newobj.p is NULL:
            return None
        newobj.regenerate_py_os()
        return newobj

    @classmethod
    def read_all(cls, f, read_header=False, thin=False):
        """read_all(file f) -> PlanarDiagram

        Returns a generator that iterates through the pdcodes stored
        in ``f``.
        """
        num_pds = float('inf')
        if read_header == True:
            # Actually parse header and check validity
            start = f.readline()
            if not start.strip() == "pdstor":
                raise Exception("Malformed pdstor header; expected `pdstor`, found %s"%start.strip())
            info = f.readline()
            args = info.strip().split(" ")
            num_pds = int(args[args.index("nelts")+1].split("/")[0])

            f.readline() # Blank

        new_pd = cls.read(f, read_header=False, thin=thin)
        i = 1
        while i < num_pds:
            yield new_pd
            try:
                new_pd = cls.read(f)
            except EOFError:
                return
            i += 1
        yield new_pd

    @classmethod
    def read_knot_theory(cls, f, thin=False):
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
        newobj.thin = thin
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
    def twist_knot(cls, n_twists, thin=False):
        """twist_knot(n_twists) -> PlanarDiagram

        Create a new :py:class:`PlanarDiagram` which represents a
        twist knot with \\\\(n\\\\) twists.
        """
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        newobj.p = pd_build_twist_knot(n_twists)
        newobj.thin = thin
        newobj.regenerate_py_os()
        return newobj
    from_twist_knot = twist_knot # Deprecated

    @classmethod
    def torus_knot(cls, p, q, thin=False):
        """torus_knot(p, q) -> PlanarDiagram

        Create a new :py:class:`PlanarDiagram` which represents a
        \\\\((p,q)\\\\)-torus knot. Only implemented for \\\\(p=2\\\\).
        """
        cdef PlanarDiagram newobj
        if p != 2:
            raise(Exception("torus_knot only implemented for p=2"))

        newobj = PlanarDiagram.__new__(cls)
        newobj.thin = thin
        newobj.p = pd_build_torus_knot(p,q)
        newobj.regenerate_py_os()
        return newobj
    from_torus_knot = torus_knot # Deprecated

    @classmethod
    def simple_chain(cls, n_links, thin=False):
        """simple_chain(n_links) -> PlanarDiagram

        Create a new :py:class:`PlanarDiagram` which represents an
        \\\\(n\\\\)-link chain.
        """
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        newobj.thin = thin
        newobj.p = pd_build_simple_chain(n_links)
        newobj.regenerate_py_os()
        return newobj
    from_simple_chain = simple_chain # Deprecated

    @classmethod
    def unknot(cls, n_crossings, thin=False):
        """unknot(n_crossings) -> PlanarDiagram

        Create a new :py:class:`PlanarDiagram` which represents an
        \\\\(n\\\\)-crossing unlink diagram
        """
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        newobj.p = pd_build_unknot(n_crossings)
        newobj.thin = thin
        newobj.regenerate_py_os()
        return newobj
    from_unknot = unknot # Deprecated

    @classmethod
    def unknot_wye(cls, a,b,c, thin=False):
        """unknot_wye(a,b,c) -> PlanarDiagram

        Create a new :py:class:`PlanarDiagram` representing an unknot which is
        designed for hash collisions."""
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        newobj.thin = thin
        newobj.p = pd_build_unknot_wye(a,b,c)
        newobj.regenerate_py_os()
        return newobj
    from_unknot_wye = unknot_wye # Deprecated

    def set_all_crossing_signs(self, signature):
        for i, sign in enumerate(signature):
            self.p.cross[i].sign = sign

    def reorient_component(self, pd_idx_t component, pd_or_t sign):
        """reorient_component(component, sign)

        Reverse the orientation of component iff sign is :py:const:`PD_NEG_ORIENTATION`
        """
        pd_reorient_component(self.p, component, sign)

    def toggle_crossing(self, pd_idx_t x_i):
        if self.p.cross[x_i].sign == PD_NEG_ORIENTATION:
            self.p.cross[x_i].sign = PD_POS_ORIENTATION
        elif self.p.cross[x_i].sign == PD_POS_ORIENTATION:
            self.p.cross[x_i].sign = PD_NEG_ORIENTATION


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

    # Knot-isomorphic modifications
    def R1_loop_deletion(self, pd_idx_t cr):
        """R1_loop_deletion(crossing_or_index) -> result pdcode

        Performs an (in-place) Reidemeister 1 loop deletion which
        removes the input crossing.
        """
        ret = PlanarDiagram_wrap(pd_R1_loopdeletion(self.p, cr))
        if ret is None:
            raise Exception("Error in R1 loopdeletion")
        return ret

    def R2_bigon_elimination(self, pd_idx_t cr1, pd_idx_t cr2):
        """R2_bigon_elimination(crossing_or_index, crossing_or_index)
        -> (Upper pdcode, [Lower pdcode])

        Performs a Reidemeister 2 bigon elimination removing the bigon
        with vertices the two inputs. This is **not an in-place operation**
        as there may be more than one PlanarDiagram as a result of the move.
        """
        cdef pd_idx_t *cr = [cr1, cr2]
        cdef pd_idx_t nout
        cdef pd_code_t **out_pds
        pd_R2_bigon_elimination(self.p, cr, &nout, &out_pds)
        if nout == 0:
            raise Exception("Error on R2 bigonelimination")
        return tuple(PlanarDiagram_wrap(out_pds[i]) for i in range(nout))

    def _R3_strand_swap(self, Face f, Edge e):
        # Stub for wrapping C implementation of tangle slide
        cdef pd_idx_t *strand_edges = <pd_idx_t*>malloc(3 * sizeof(pd_idx_t))
        cdef pd_idx_t *border_faces = <pd_idx_t*>malloc(3 * sizeof(pd_idx_t))
        cdef pd_idx_t npieces = 0
        cdef pd_code_t **pd_pieces
        cdef pd_tangle_t *tangle = pd_tangle_new(4)
        cdef Crossing x = [x for x in f.get_vertices() if
                           x.index != e.head and x.index != e.tail][0]
        for i, e_i in enumerate(x.edges):
            tangle.edge[i] = e_i
        pd_regenerate_tangle(self.p, tangle)

        strand_edges[0] = e.prev_edge().index
        strand_edges[1] = e.index
        strand_edges[2] = e.next_edge().index

        print "Strand.. ",
        print strand_edges[0], strand_edges[1], strand_edges[2]

        if e.face_index_pos()[0][0] == f.index:
            posneg = 0
        else:
            posneg = 1

        border_faces[0] = e.prev_edge().face_index_pos()[posneg][0]
        border_faces[1] = f.index
        border_faces[2] = e.next_edge().face_index_pos()[posneg][0]

        print "Face.. ",
        print border_faces[0], border_faces[1], border_faces[2]

        pd_tangle_slide(self.p, tangle,
                        3, strand_edges, border_faces,
                        &npieces, &pd_pieces)

    def R3_nice_swap(self, Face f, Edge e):
        fip = e.face_pos()
        pos = fip[0][1] if f == fip[0][0] else fip[1][1]
        a = f[(pos+1) %3]
        b = f[(pos+2) %3]
        x = a.next_crossing() if f.signs[(pos+1) %3] == 1 else a.prev_crossing()
        y = b.prev_crossing() if f.signs[(pos+2) %3] == 1 else b.next_crossing()
        a2 = a.next_edge() if f.signs[(pos+1) %3] == 1 else a.prev_edge()
        b2 = b.prev_edge() if f.signs[(pos+2) %3] == 1 else b.next_edge()
        return self.R3_strand_swap(x.index, a.index, b.index, a2.index, b2.index)


    def R3_strand_swap(self, pd_idx_t i_x,
                       pd_idx_t i_a_1, pd_idx_t i_b_1,
                       pd_idx_t i_a_2, pd_idx_t i_b_2):
        cdef PlanarDiagram result = self.copy()
        cdef Crossing x = result.crossings[i_x]
        cdef Edge a_1 = result.edges[i_a_1]
        cdef Edge a_2 = result.edges[i_a_2]
        cdef Edge b_1 = result.edges[i_b_1]
        cdef Edge b_2 = result.edges[i_b_2]
        cdef Edge a_0, b_0
        cdef bool a_1_in, b_1_in

        a_1_in = a_1.head == x.index
        b_1_in = b_1.head == x.index

        if a_1_in:
            a_0 = a_1.prev_edge()
        else:
            a_0 = a_1.next_edge()

        if b_1_in:
            b_0 = b_1.prev_edge()
        else:
            b_0 = b_1.next_edge()

        if a_1_in and b_1_in:
            #print "in in"
            a_0.swap_head(a_1, PD_POS_ORIENTATION)
            a_1.swap_tail(b_2, PD_NEG_ORIENTATION)
            b_0.swap_head(b_1, PD_POS_ORIENTATION)
            b_1.swap_tail(a_2, PD_NEG_ORIENTATION)

        elif a_1_in and not b_1_in:
            #print "in out"
            b_1.next_crossing().sign = (b_1.next_crossing().sign+1) % 2
            a_1.prev_crossing().sign = (a_1.prev_crossing().sign+1) % 2

            b_0.swap_tail(b_1, PD_NEG_ORIENTATION)
            b_1.swap_head(a_2, PD_NEG_ORIENTATION)
            b_1.swap_head(b_1, PD_NEG_ORIENTATION)

            a_0.swap_head(a_1, PD_POS_ORIENTATION)
            a_1.swap_tail(b_2, PD_POS_ORIENTATION)
            a_1.swap_tail(a_1, PD_POS_ORIENTATION)
        elif not a_1_in and not b_1_in:
            #print "out out"
            a_0.swap_tail(a_1, PD_NEG_ORIENTATION)
            a_1.swap_head(b_2, PD_POS_ORIENTATION)
            b_0.swap_tail(b_1, PD_NEG_ORIENTATION)
            b_1.swap_head(a_2, PD_POS_ORIENTATION)
        elif not a_1_in and b_1_in:
            #print "out in"
            b_1.prev_crossing().sign = (b_1.prev_crossing().sign+1) % 2
            a_1.next_crossing().sign = (a_1.next_crossing().sign+1) % 2

            a_0.swap_tail(a_1, PD_NEG_ORIENTATION)
            a_1.swap_head(b_2, PD_NEG_ORIENTATION)
            a_1.swap_head(a_1, PD_NEG_ORIENTATION)

            b_0.swap_head(b_1, PD_POS_ORIENTATION)
            b_1.swap_tail(a_2, PD_POS_ORIENTATION)
            b_1.swap_tail(b_1, PD_POS_ORIENTATION)

        result.regenerate()
        return result

    def tangle_slide(self, tangle, crs_edges, strand):
        pass

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

    def homfly(self, as_string=False):
        """homfly([as_string=False]) -> HOMFLYPolynomial

        Compute the HOMFLY polynomial for this diagram."""
        cdef bytes homflybytes = copy_and_free(pd_homfly(self.p))
        if as_string:
            return homflybytes
        else:
            return HOMFLYPolynomial(homflybytes)

    def unique_code(self):
        """unique_code() -> str

        Returns a 'pdcode-ish' string which uniquely identifies this
        planar diagram. If another planar diagram is not isomorphic to
        this, their ``unique_code()``\ s will differ, but two isomorphic
        planar diagrams may have different ``unique_code()``\ s.
        """
        pdstr = "PD[%s, CompSigns=[%s]]"%(
            ", ".join(("X{}[{},{},{},{}]".format(
                ori_char(crossing.sign),*crossing.edges) for
                       crossing in self.crossings)),
            ",".join(ori_char(component[0].prev_crossing().sign) for
                     component in self.components)
        )
        return pdstr

    @classmethod
    def _kt_db_search(cls, search_key, names_f, table_f):
        for line in names_f:
            if search_key in line:
                return cls.read_knot_theory(table_f)
            table_f.readline() # burn the line
        return None

    KNOT_NAMES = os.path.join(SOURCE_DIR,"rolfsennames.txt")
    KNOT_TABLE = os.path.join(SOURCE_DIR,"rolfsentable.txt")
    @classmethod
    def db_knot(cls, ncross, index, alternating=None):
        """db_knot(ncross, index, [alternating=None]) -> new PlanarDiagram

        Searches the Rolfsen table of knot types and returns a new
        PlanarDiagram which is isotopic. Raises KeyError if the specified
        knot is not found.

        The variable ``alternating`` defaults to None. If the crossing
        number is high enough that knots are distinguished by
        'Alternating' or 'NonAlternating,' be sure to set
        ``alternating`` to True or False respectively.

        """
        search_key = "Knot[%s, %s%s]"%(
            ncross,
            ("" if alternating is None else
             "Alternating, " if alternating else
             "NonAlternating, "),
            index)
        loaded_pd = None
        with open(cls.KNOT_NAMES) as names_f, open(cls.KNOT_TABLE) as table_f:
            loaded_pd = cls._kt_db_search(search_key, names_f, table_f)
        if loaded_pd is None:
            raise KeyError("%s not found in knots table"%search_key)
        return loaded_pd


    LINK_NAMES = os.path.join(SOURCE_DIR,"thistlethwaitenames.txt")
    LINK_TABLE = os.path.join(SOURCE_DIR,"thistlethwaitetable.txt")
    @classmethod
    def db_link(cls, ncross, index, alternating=True):
        """db_link(ncross, index, [alternating=True]) -> new PlanarDiagram

        Searches the Thistlethwaite table of link types and returns a new
        PlanarDiagram which is isotopic. Raises KeyError if the specified
        link is not found."""
        search_key = "Link[%s, %sAlternating, %s]"%(
            ncross,
            "" if alternating else "Non",
            index)
        loaded_pd = None
        with open(cls.LINK_NAMES) as names_f, open(cls.LINK_TABLE) as table_f:
            loaded_pd = cls._kt_db_search(search_key, names_f, table_f)
        if loaded_pd is None:
            raise KeyError("%s not found in links table"%search_key)
        return loaded_pd

    def get_R1_loops(self):
        """get_R1_loops() -> generator of Faces

        Returns a generator of loops which are viable candidates
        for Reidemeister 1 moves.
        """
        if self.ncross > 0:
            for face in self.faces:
                if len(face) == 1:
                    yield face

    def get_R2_bigons(self):
        """get_R2_bigons() -> generator of Faces

        Returns a generator of bigons which are viable candidates
        for Reidemeister 2 moves.
        """
        if self.ncross > 0:
            for face in self.faces:
                verts = face.get_vertices()
                if len(verts) == 2 and verts[0].sign != verts[1].sign:
                    yield face

    def get_R3_triangles(self):
        """get_R3_triangles() -> generator of (Face, Edge) pairs

        Returns a generator of trigons which are viable
        candidates for Reidemeister 3 moves, with the strand
        that would move."""
        if self.ncross > 2:
            for face in self.faces:
                verts = set(face.vertices)
                if len(face) == 3 and len(verts) == 3:
                    for edge in face:
                        if ((edge.index in edge.prev_crossing().overstrand_indices()) ==
                            (edge.index in edge.next_crossing().overstrand_indices())):
                            yield face, edge


    def neighbors(self):
        """neighbors() -> generator of tuples of PlanarDiagrams

        Returns a generator of knot-isomorphic PlanarDiagrams by
        identifying and performing knot moves which preserve
        isomorphism (i.e. Reidemeister moves). The generator yields
        *tuples* of PlanarDiagrams, in case the move untangled a link
        into two separate diagrams.
        """
        for loop in self.get_R1_loops():
            yield (self.R1_loop_deletion(loop.vertices[0]),)
        for bigon in self.get_R2_bigons():
            yield self.R2_bigon_elimination(*bigon.vertices)

    def simplify(self, seed=None):
        """simplify(seed=None) -> PlanarDiagram

        Returns an isomorphic PlanarDiagram whose crossing number
        software knows not how to decrease (i.e. this would
        **ideally** return a diagram of this knot with minimal number
        of crossings).

        If ``seed`` is not ``None``, pick a random move instead of
        using any order (using the integer seed provided).
        """
        if seed:
            return self._simplify_helper(random.Random(seed))
        else:
            return self._simplify_helper(None)

    def _simplify_helper(self, pyrng):
        if pyrng is None:
            try:
                if self.ncross == 0:
                    return (self.copy(),)
                if self.ncross == 1:
                    return (PlanarDiagram.unknot(0),)
                neighbor = self.neighbors().next()
                assert len(neighbor) > 0
                return sum((dia._simplify_helper(pyrng) for dia in neighbor),())

            except StopIteration:
                return (self.copy(),)

    def dual_graph(self):
        dual = PlanarSignedFaceDigraph()
        dual.verts = range(self.nfaces)
        for edge in self.edges:
            (pos_i,_),(neg_i,_) = edge.face_index_pos()
            dual.add_edge(neg_i, pos_i)
        for x in self.crossings:
            dual.add_face([dual.edges[e_i] for e_i in x.edges], x.sign)
        return dual

    def as_spherogram(self):
        """as_spherogram() -> spherogram.Link

        Returns a Spherogram Link object representing this PlanarDiagram

        Requires: ``spherogram``"""
        pdcode = self.pdcode()

        from spherogram import links
        sg_xings = [links.Crossing(x.index) for x in self.crossings]

        gluings = defaultdict(list)
        for x_i, (xing, pd_xing, sg_xing) in enumerate(
                zip(self.crossings, pdcode, sg_xings)):
            if xing.sign == PD_NEG_ORIENTATION:
                sg_xing.sign = -1
            elif xing.sign == PD_POS_ORIENTATION:
                sg_xing.sign = 1
            else:
                sg_xing.sign = 0
            for pos, e_i in enumerate(pd_xing):
                gluings[e_i].append((sg_xing, pos))
        for (x, xpos), (y, ypos) in gluings.itervalues():
            x[xpos] = y[ypos]

        return links.Link(sg_xings)

    def snappy_manifold(self):
        """snappy_manifold() -> snappy.Manifold

        Returns a SnapPy Manifold object which represents this knot/link's
        complement in $S^3$.

        Requires: `snappy`, `spherogram`"""
        from snappy import Manifold
        sgm_link = self.as_spherogram()
        return Manifold("DT:%s"%sgm_link.DT_code())

    @classmethod
    def from_dt_code(cls, dt_code):
        """from_dt_code() -> new PlanarDiagram

        Creates a new PlanarDiagram from a DT code (Dowker-Thistlethwaite).
        Requires Spherogram

        Warning: Currently depends on from_pdcode, which has bugs with
        components of length 2

        Requires: ``spherogram``"""
        import spherogram
        return cls.from_pdcode(spherogram.DTcodec(dt_code).PD_code())

    @classmethod
    def from_plink(cls, editor, thin=False):
        """from_plink(editor) -> new PlanarDiagram

        Creates a new PlanarDiagram from a plink LinkEditor object."""

        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        cdef pd_code_t* pd
        cdef pd_idx_t ncross = len(editor.Crossings)
        newobj.thin = thin

        components = editor.crossing_components()
        pd = pd_code_new(ncross+2)
        newobj.p = pd

        pd.ncross = ncross
        pd.nedges = ncross*2
        pd.nfaces = ncross+2
        pd.ncomps = len(components)

        i_e = 0
        for component in components:
            for etail, ehead in izip(component, islice(cycle(component), 1, None)):
                i_tail = editor.Crossings.index(etail.crossing)
                if etail.crossing.goes_over():
                    # TODO: ASSERT THAT THIS IS CORRECT INTERPRETATION
                    if etail.crossing.sign() == "RH":
                        pos = 1
                    elif etail.crossing.sign() == "LH":
                        pos = 3
                else:
                    pos = 2

                pd.cross[i_tail].edge[pos] = i_e
                pd.edge[i_e].tail = i_tail
                pd.edge[i_e].tailpos = pos

                i_head = editor.Crossings.index(ehead.crossing)
                pd.cross[i_head].edge[(pos+2)%4] = i_e
                pd.edge[i_e].head = i_head
                pd.edge[i_e].headpos = (pos+2)%4

                i_e += 1

        pd_regenerate(pd)
        return newobj

    def ccode(self):
        return copy_and_free(pdcode_to_ccode(self.p))

    def pdcode(self):
        cdef list pdcode = []
        cdef list x_temp
        cdef pd_idx_t o_in, o_out, u_in, u_out
        for x in self.crossings:
            x_temp = [-1,-1,-1,-1]
            o_in, o_out = x.overstrand_indices()
            o_ipos, o_opos = x.overstrand_pos()
            u_in, u_out = x.understrand_indices()
            u_ipos, u_opos = x.understrand_pos()

            x_temp[0] = u_in
            x_temp[2] = u_out

            if (u_ipos + 1) % 4 == o_ipos:
                x_temp[1] = o_in
                x_temp[3] = o_out
            else:
                x_temp[3] = o_in
                x_temp[1] = o_out

            pdcode.append(tuple(x_temp))
        return pdcode

    @classmethod
    def from_pdcode(cls, pdcode, thin=False):
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        cdef pd_code_t* pd
        cdef pd_idx_t ncross = len(pdcode)
        cdef pd_idx_t off = min([min(x) for x in pdcode])
        newobj.thin = thin
        pd = pd_code_new(ncross+2)
        newobj.p = pd
        pd.ncross = ncross
        pd.nedges = ncross*2
        pd.nfaces = ncross+2

        # tail dictionary for consistency checking
        tails = dict()

        # Set crossings and fill in some edges
        for i,x in enumerate(pdcode):
            for pos in (0,1,2,3):
                pd.cross[i].edge[pos] = x[pos] - off
        pd_regenerate_edges(pd)

        # Set signs on the crossings. It's important to do this
        # after regenerate_edges() as that removes sign info
        for i,x in enumerate(pdcode):
            tails[x[2]] = i
            if x[3] - x[1] == 1:
                if x[3] not in tails:
                    tails[x[3]] = i
                    pd.cross[i].sign = PD_NEG_ORIENTATION
                else:
                    tails[x[1]] = i
                    pd.cross[i].sign = PD_POS_ORIENTATION
            elif x[1] - x[3] == 1:
                if x[1] not in tails:
                    tails[x[1]] = i
                    pd.cross[i].sign = PD_POS_ORIENTATION
                else:
                    tails[x[3]] = i
                    pd.cross[i].sign = PD_NEG_ORIENTATION
            elif x[1] > x[3]:
                tails[x[3]] = i
                pd.cross[i].sign = PD_NEG_ORIENTATION
            elif x[3] > x[1]:
                tails[x[1]] = i
                pd.cross[i].sign = PD_POS_ORIENTATION
            else:
                print x
                raise Exception("Error in pdcode or pdcode reader")

        pd_regenerate(pd)
        newobj.regenerate_py_os()
        return newobj

    cdef regenerate_py_os(self):
        if self.thin:
            return
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

    # Special methods which allow pdcodes to be pickled.
    def __reduce__(self):
        """Required method to implement pickling of extension types.

        Returns a tuple: (Constructor-ish, Constructor args, __getstate__ output).
        """
        return (self.__class__,
                (self.p.ncross,),
                self.__getstate__())

    def __hash__(self):
        return hash((tuple(self.p.face[i].nedges for i in range(self.nfaces)), self.ncross,
                     tuple(self.p.comp[i].nedges for i in range(self.ncomps)), self.nedges))
        #return hash(self.__getstate__())

    def __getstate__(self):
        cdef pd_crossing_t* x
        cdef pd_edge_t* e
        cross = []
        edges = []
        for i in range(self.p.ncross):
            x = &(self.p.cross[i])
            cross.append((x.edge[0], x.edge[1], x.edge[2], x.edge[3],
                          x.sign))
        for i in range(self.p.nedges):
            e = &(self.p.edge[i])
            edges.append((e.head, e.headpos,
                          e.tail, e.tailpos))
        return (tuple(cross), tuple(edges))

    def __setstate__(self, state):
        cdef pd_crossing_t* x
        cdef pd_edge_t* e
        cross, edges = state # See __getstate__
        for i,data in enumerate(cross):
            x = &(self.p.cross[i])
            x.edge[0] = data[0]
            x.edge[1] = data[1]
            x.edge[2] = data[2]
            x.edge[3] = data[3]
            x.sign = data[4]
        self.p.ncross = len(cross)
        for i,data in enumerate(edges):
            e = &(self.p.edge[i])
            e.head = data[0]
            e.headpos = data[1]
            e.tail = data[2]
            e.tailpos = data[3]
        self.p.nedges = len(edges)
        self.regenerate()

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

cdef api pd_code_t **pd_simplify(pd_code_t *pd, int *ndias):
    cdef PlanarDiagram py_pd = PlanarDiagram_wrap(pd)
    cdef tuple simp_dias = py_pd.simplify()
    cdef PlanarDiagram simp
    cdef pd_code_t **ret
    py_pd.p = NULL # Cede the original diagram pointer to C
    ndias[0] = len(simp_dias)
    ret = <pd_code_t **>malloc(ndias[0] * sizeof(pd_code_t*))
    for i,simp in enumerate(simp_dias):
        ret[i] = simp.p
        simp.p = NULL # It's in C's hands now
    return ret
