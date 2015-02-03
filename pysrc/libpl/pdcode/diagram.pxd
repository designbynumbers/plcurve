from .plctopology cimport pd_code_t, bool

cdef class _OwnedObjectList(list):
    cdef delitem(self, i)
    cdef clearout(self)

cdef class _EdgeList(_OwnedObjectList):
    pass
cdef class _CrossingList(_OwnedObjectList):
    pass
cdef class _ComponentList(_OwnedObjectList):
    pass
cdef class _FaceList(_OwnedObjectList):
    pass

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

    cdef readonly bool thin
    cdef bool hashed
    cdef char* _homfly

    cdef regenerate_py_os(self)

cdef PlanarDiagram PlanarDiagram_wrap(pd_code_t *pd, bool thin=*)
