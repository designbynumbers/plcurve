from .plctopology cimport *
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

from .components cimport *

#from cython.view cimport array
cimport cython

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

cdef Edge Edge_FromParent(PlanarDiagram parent, pd_idx_t index):
    cdef Edge new_edge = Edge.__new__(Edge)
    new_edge.parent = parent
    new_edge.p = &parent.p.edge[index]
    new_edge.index = index
    return new_edge

cdef class Component(_Disownable):
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
cdef Component Component_FromParent(PlanarDiagram parent, pd_idx_t index):
    cdef Component new_component = Component.__new__(Component)
    new_component.parent = parent
    new_component.p = &parent.p.comp[index]
    new_component.index = index
    return new_component

cdef class Face(_Disownable):
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
        """get_vertices() -> (Crossing, [Crossing, ...])

        A tuple of Crossings around this face in counterclockwise
        order. The first crossing is the one immediately before the first edge
        in this face.
        """
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

cdef Face Face_FromParent(PlanarDiagram parent, pd_idx_t index):
    cdef Face new_face = Face.__new__(Face)
    new_face.parent = parent
    new_face.p = &parent.p.face[index]
    new_face.index = index
    return new_face


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
        """toggle_sign(self)

        Switches the sign of the crossing to be positive if
        it was negative, or negative if it was positive. Has no
        effect on a crossing which is unsigned."""
        if self.p.sign == PD_NEG_ORIENTATION:
            self.p.sign = PD_POS_ORIENTATION
        elif self.p.sign == PD_POS_ORIENTATION:
            self.p.sign = PD_NEG_ORIENTATION

    property faces:
        """A tuple of the indices of the faces around this vertex,
        in counterclockwise order.
        """
        def __get__(self):
            ret = [0]*4
            for i, edge in enumerate(self):
                (pos_face,_),(neg_face,_) = edge.face_index_pos()
                if edge.tail == self.index and edge.tailpos == i:
                    # Edge departed from this crossing
                    ret[i] = neg_face
                else:
                    ret[i] = pos_face
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

cdef Crossing Crossing_FromParent(PlanarDiagram parent, pd_idx_t index):
    cdef Crossing new_crossing = Crossing.__new__(Crossing)
    new_crossing.parent = parent
    new_crossing.p = &parent.p.cross[index]
    new_crossing.index = index
    return new_crossing
