from .plctopology cimport pd_code_t, bool

cdef class PlanarDiagram

cdef class _OwnedObjectList:
    cdef PlanarDiagram parent
    cdef object objs
    cdef object _new_child_object(self, long i)

cdef class _EdgeList(_OwnedObjectList):
    pass

cdef class _CrossingList(_OwnedObjectList):
    pass

cdef class _ComponentList(_OwnedObjectList):
    pass

cdef class _FaceList(_OwnedObjectList):
    pass

cdef class PlanarDiagram:
    """
    Class which represents a PD code (planar diagram).  There are a few
    different ways to create `PlanarDiagram`\ s:

    The constructor will create a new empty diagram with `max_verts` blank
    crossings.

    There are two class methods for reading diagrams from files, both from
    files written in `PlanarDiagram`'s own format, and from files written in
    Mathematica's KnotTheory format.

    There are multiple class methods which create new `PlanarDiagram`\ s of
    basic examples; their docstrings contain some specific information.

    `PlanarDiagram`\ s support the following operations; given ``K`` and ``L``
    two objects;

        - ``K == L`` checks whether `K` and `L` are `isotopic()` if
          `equiv_type` is `ISOTOPY`, or if they are `map_isomorphic()` if
          `equiv_type` is `UO_ISOMORPHISM`.

        - ``K != L`` checks whether `K` and `L` are not `isotopic()` if
          `equiv_type` is `ISOTOPY`, or if they are not `map_isomorphic()` if
          `equiv_type` is `UO_ISOMORPHISM`.

        - ``len(K)`` returns the number of components of ``K``, `ncomps`.
    """

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
    cdef readonly bool hashed
    cdef char* _homfly

    cdef public int equiv_type
    """
    What type of equivalence to use for the ``==`` operation. Valid values are,

        - `ISOTOPY`, the default; equality is dictated by `isotopic()`

        - `UO_ISOMORPHISM`; equality is dictated by `map_isomorphic()`
    """

    cdef regenerate_py_os(self)

cdef PlanarDiagram PlanarDiagram_wrap(pd_code_t *pd, bool thin=*)

cdef public api:
    pd_code_t **pd_simplify(pd_code_t *pd, int *ndias)
