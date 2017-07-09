from .plctopology cimport (
    pd_idx_t, pd_or_t, pd_edge_t, pd_crossing_t, pd_face_t, pd_component_t,
    bool)
from .diagram cimport PlanarDiagram

cdef class _Disownable:
    cdef disown(self)

cdef class Edge(_Disownable):
    """
    An oriented edge, joining two vertices tail -> head.

    Once you are finished setting any of the properties `head`, `headpos`,
    `tail`, and `tailpos`, you **must** call
    `~libpl.pdcode.diagram.PlanarDiagram.regenerate()` on the parent so that
    sanity can be checked and faces and components regenerated.
    """
    cdef pd_edge_t *p
    cdef readonly PlanarDiagram parent
    """The `~.diagram.PlanarDiagram` to which this edge belongs."""
    cdef readonly int index
    """The index by which this edge is known in the parent."""

    cdef bool owned

    cpdef Edge oriented(self, pd_or_t sign)
    cpdef component_index_pos(self)
    cpdef component_pos(self)
    cpdef face_index_pos(self)
    cpdef face_pos(self)
    cpdef reorient(self, pd_or_t sign)
    cpdef plug_head(self, Crossing x, pd_idx_t pos)
    cpdef plug_tail(self, Crossing x, pd_idx_t pos)
    cpdef swap_head(self, Edge other, pd_or_t head_or_tail)
    cpdef swap_tail(self, Edge other, pd_or_t head_or_tail)

cdef Edge Edge_FromParent(PlanarDiagram parent, pd_idx_t index)

cdef class Component(_Disownable):
    """
    A component of a diagram.

    A component is an oriented cycle of edges in its parent diagram, formed by
    proceeding directly across each crossing to the following edge.  If a
    diagram has only one Component, it is a knot diagram.
    """
    cdef pd_component_t* p
    cdef readonly PlanarDiagram parent
    """The :py:class:`PlanarDiagram` to which this component belongs."""
    cdef readonly int index
    """The index by which this component is known in the parent."""
cdef Component Component_FromParent(PlanarDiagram parent, pd_idx_t index)

cdef class Face(_Disownable):
    """
    A face of a diagram.

    A face is an oriented cycle of edges in its parent diagram, formed by
    turning at each corner of a vertex.  Each face has the topology of a disk.
    """

    cdef pd_face_t* p
    cdef readonly PlanarDiagram parent
    """The :py:class:`PlanarDiagram` to which this face belongs."""
    cdef readonly int index
    """The index by which this face is known in the parent."""
cdef Face Face_FromParent(PlanarDiagram parent, pd_idx_t index)

cdef class Crossing(_Disownable):
    """
    A crossing, which holds four edges in positions and an orientation.

    Crossing data is the **most important** piece of data in a PlanarDiagram.
    A full diagram object can be built from crossing data alone (although
    orientations of components may differ).

    Once you are finished setting any of the edges, you **must** call
    `~.diagram.PlanarDiagram.regenerate()` on the parent so that sanity can be
    checked and faces and components regenerated.

    **Accessing crossing information**

    Crossings support two methods of viewing their child edge data;

        - Crossings ``X`` support being accessed as if they were lists;
          ``X[i]`` retrieves the `Edge` attached at position ``i``.  This
          method takes into account that vertices are cyclical, but does not
          support assignment.

        - The property `edges` is a direct-memory view of edge indices attached
          to this crossing.  It *is* possible to change these indices, but the
          parent must be regenerated and validity checked manually after doing
          this.
    """

    cdef pd_crossing_t* p
    cdef readonly PlanarDiagram parent
    """The :py:class:`PlanarDiagram` to which this crossing belongs."""
    cdef readonly int index
    """The index by which this crossing is known in the parent."""

    cdef pd_idx_t[:] edgeview_get(self)
    cdef PlanarDiagram parent_x(self)
    cpdef overstrand_indices(self)
    cpdef overstrand_pos(self)
    cpdef understrand_indices(self)
    cpdef understrand_pos(self)
cdef Crossing Crossing_FromParent(PlanarDiagram parent, pd_idx_t index)
