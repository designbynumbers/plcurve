cimport cython
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.string cimport memcpy, strncpy
from libc.stdlib cimport calloc, malloc, free
import sys
from cython.operator cimport dereference as deref
import random
import re
from operator import itemgetter, mul
from itertools import islice, izip, cycle
import os
import libpl.data
from ..graphs import PlanarSignedFaceDigraph
from libc.stdlib cimport free
from collections import defaultdict
import weakref

from .pdisomorphism cimport (
    pd_iso_t, pd_build_diagram_isotopies,
    pd_build_isos, pd_build_map_isomorphisms)
from .isomorphism cimport PlanarIsomorphism
from .pd_invariants cimport *

from .plctopology cimport *

from .components cimport *
from .homfly cimport *
from .planarmap cimport *

#from cython.view cimport array
cimport cython
cdef extern char* pdcode_to_ccode(pd_code_t *pdC)
PD_VERBOSE = 10
#PD_LIVE_ON_ERROR = 1

DEFAULT_PATH=os.path.join("data","pdstors")
SOURCE_DIR=libpl.data.dir

POS_ORIENTATION = PD_POS_ORIENTATION
NEG_ORIENTATION = PD_NEG_ORIENTATION
UNSET_ORIENTATION = PD_UNSET_ORIENTATION

def pd_debug_off():
    """Disable debug output from libplCurve"""
    global PD_VERBOSE
    PD_VERBOSE = 9

def pd_debug_on():
    """Enable debug output from libplCurve"""
    global PD_VERBOSE
    PD_VERBOSE = 10

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

class NoParentException(Exception):
    pass

# List subclasses which know how to properly delete/set elements within
cdef class _OwnedObjectList:
    def __cinit__(self, PlanarDiagram parent):
        self.parent = parent
        self.objs = dict()

    def __len__(self):
        return 0

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def __repr__(self):
        return ("[" +
                ", ".join(repr(o) for o in self) +
                "]")

    def __getitem__(self, i):
        if not (isinstance(i, int) or isinstance(i, long)):
            raise TypeError("Index must be an integer.")

        if i < 0 or i > len(self):
            raise IndexError("Index %s is out of range"%i)

        if i not in self.objs:
            self.objs[i] = self._new_child_object(i)
        return self.objs[i]

    cdef object _new_child_object(self, long i):
        return None

    def _disown(self):
        cdef _Disownable o
        for o in self.objs.itervalues():
            o.disown()
        self.parent = None

cdef class _EdgeList:
    def __len__(self):
        return self.parent.nedges

    cdef object _new_child_object(self, long i):
        return Edge_FromParent(self.parent, i)

cdef class _CrossingList(_OwnedObjectList):
    def __len__(self):
        return self.parent.ncross

    cdef object _new_child_object(self, long i):
        return Crossing_FromParent(self.parent, i)

cdef class _ComponentList(_OwnedObjectList):
    def __len__(self):
        return self.parent.ncomps

    cdef object _new_child_object(self, long i):
        return Component_FromParent(self.parent, i)

cdef class _FaceList(_OwnedObjectList):
    def __len__(self):
        return self.parent.nfaces

    cdef object _new_child_object(self, long i):
        return Face_FromParent(self.parent, i)

cdef PlanarDiagram PlanarDiagram_wrap(pd_code_t *pd, bool thin=False):
    """Wraps a pd_code_t* in a PlanarDiagram"""
    cdef PlanarDiagram newobj = PlanarDiagram.__new__(PlanarDiagram)
    newobj.thin = thin
    if pd is NULL:
        return None
    newobj.p = pd
    newobj.regenerate_py_os()
    return newobj

ISOTOPY = 0
UO_ISOMORPHISM = 1

cdef class PlanarDiagram:

    property ncomps:
        """Number of components in this :class:`PlanarDiagram`"""
        def __get__(self):
            return self.p.ncomps
        def __set__(self, n):
            self.p.ncomps = n
            self.regenerate_py_os()
    property nedges:
        """Number of edges in this :class:`PlanarDiagram`"""
        def __get__(self):
            return self.p.nedges
        def __set__(self, n):
            self.p.nedges = n
            self.regenerate_py_os()
    property ncross:
        """Number of crossings in this :class:`PlanarDiagram`"""
        def __get__(self):
            return self.p.ncross
        def __set__(self, n):
            self.p.ncross = n
            self.regenerate_py_os()
    property nfaces:
        """Number of faces in this :class:`PlanarDiagram`"""
        def __get__(self):
            return self.p.nfaces
        def __set__(self, n):
            self.p.nfaces = n
            self.regenerate_py_os()

    property hash:
        """Hash of this of this :class:`PlanarDiagram`"""
        def __get__(self):
            if not self.hashed:
                self.hashed = True
                pd_regenerate_hash(self.p)
            return self.p.hash

    property uid:
        """UID of this :class:`PlanarDiagram`"""
        def __get__(self):
            return self.p.uid

    property _hash:
        def __get__(self):
            return self.p.hash

    def __cinit__(self):
        self.p = NULL
        self.thin = False
        self.hashed = False
        self._homfly = NULL
        self.edges = _EdgeList(self)
        self.crossings = _CrossingList(self)
        self.components = _ComponentList(self)
        self.faces = _FaceList(self)
        self.equiv_type = ISOTOPY

    def __init__(self, max_verts=15):
        """__init__([max_verts=15])

        Create a new PlanarDiagram. You may set the initial available number
        of vertices by passing a number to :obj:`max_verts`.

        :param int max_verts: The initial size (in vertices) of the object
        """
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
            if self.equiv_type == ISOTOPY:
                return pd_diagram_isotopic(self.p, other_pd.p)
            elif self.equiv_type == UO_ISOMORPHISM:
                return pd_map_isomorphic(self.p, other_pd.p)
        elif op == 3:
            if self.equiv_type == ISOTOPY:
                return not pd_diagram_isotopic(self.p, other_pd.p)
            elif self.equiv_type == UO_ISOMORPHISM:
                return not pd_map_isomorphic(self.p, other_pd.p)
        else:
            raise NotImplementedError(
                "PlanarDiagrams do not support relative comparisons.")

    def isotopic(self, PlanarDiagram other_pd):
        """
        isotopic(PlanarDiagram other_pd) -> bool

        Returns whether or not this :class:`PlanarDiagram` is
        diagram-isotopic to the input :class:`PlanarDiagram`
        :obj:`other_pd`. This is equivalent to ``pd == other`` or
        ``not pd != other``.

        :param PlanarDiagram other_pd: Diagram to check isotopy to
        :return: Whether or not the diagrams are isotopic
        :rtype: bool
        """
        return pd_diagram_isotopic(self.p, other_pd.p)

    def map_isomorphic(self, PlanarDiagram other_pd):
        """
        map_isomorphic(PlanarDiagram other_pd) -> bool

        Returns whether or not this pdcode is
        oriented-\\\\(S^2\\\\)-isotopic to the input
        :class:`PlanarDiagram` :obj:`other_pd`.

        :param PlanarDiagram other_pd: Diagram to check oriented-isomorphism
          to
        :return: Whether or not the diagrams are oriented-isomorphic
        :rtype: bool
        """
        return pd_map_isomorphic(self.p, other_pd.p)

    def build_isotopies(self, PlanarDiagram other_pd):
        """
        build_isotopies(PlanarDiagram other_pd) -> tuple of PlanarIsomorphism

        Compute all diagram isotopies to :obj:`other_pd`.

        :param PlanarDiagram other_pd: Diagram to compute diagram-isotopy
          to
        :return: Tuple of isotopies
        :rtype: tuple of PlanarIsomorphism
        """
        cdef unsigned int nisos
        cdef pd_iso_t **isos

        isos = pd_build_diagram_isotopies(self.p, other_pd.p, &nisos)
        ret = tuple(PlanarIsomorphism() for _ in range(nisos))
        for i in range(nisos):
            (<PlanarIsomorphism>ret[i]).p = isos[i]

        return ret

    def build_map_isomorphisms(self, PlanarDiagram other_pd):
        """
        build_map_isomorphisms(PlanarDiagram other_pd) -> tuple of PlanarIsomorphism

        Compute all oriented-sphere isomorphisms to :obj:`other_pd`.

        :param PlanarDiagram other_pd: Diagram to compute oriented-isomorphism
          to
        :return: Tuple of isotopies
        :rtype: tuple of PlanarIsomorphism
        """

        cdef unsigned int nisos
        cdef pd_iso_t **isos

        isos = pd_build_map_isomorphisms(self.p, other_pd.p, &nisos)
        ret = tuple(PlanarIsomorphism() for _ in range(nisos))
        for i in range(nisos):
            (<PlanarIsomorphism>ret[i]).p = isos[i]

        return ret

    def build_autoisotopies(self):
        """
        build_autoisotopies() -> tuple of PlanarIsomorphism

        Compute all diagram isotopies to itself

        :return: Tuple of isotopies
        :rtype: tuple of PlanarIsomorphism
        """

        return self.build_isotopies(self)

    def build_isomorphisms(self, PlanarDiagram other_pd):
        """
        build_isomorphisms(PlanarDiagram other_pd) -> tuple of PlanarIsomorphism

        Compute all unoriented-sphere isomorphisms to :obj:`other_pd`.

        :param PlanarDiagram other_pd: Diagram to compute unoriented-isomorphism
          to
        :return: Tuple of isomorphisms
        :rtype: tuple of PlanarIsomorphism
        """

        cdef unsigned int nisos
        cdef pd_iso_t **isos

        isos = pd_build_isos(self.p, other_pd.p, &nisos)
        ret = tuple(PlanarIsomorphism() for _ in range(nisos))
        for i in range(nisos):
            (<PlanarIsomorphism>ret[i]).p = isos[i]

        return ret

    def build_automorphisms(self):
        """
        build_autoisotopies() -> tuple of PlanarIsomorphism

        Compute all unoriented isomorphisms to itself

        :return: Tuple of isomorphisms
        :rtype: tuple of PlanarIsomorphism
        """
        return self.build_isomorphisms(self)

    def isomorphic(self, PlanarDiagram other_pd):
        """isomorphic(PlanarDiagram other_pd) -> bool

        Returns whether or not this pdcode is shadow isomorphic to the
        input pdcode ``other_pd``.

        :param PlanarDiagram other_pd: Diagram to check shadow isomorphism to
        :return: Whether or not the diagrams are shadow isomorphic
        :rtype: bool
        """
        return pd_isomorphic(self.p, other_pd.p)

    def identical(self, PlanarDiagram other_pd):
        """identical(PlanarDiagram other_pd) -> bool

        Returns whether or not this diagram is identical to :obj:`other_pd`,
        i.e., whether or not it has the exact same crossings, edges,
        faces, and components. Does not care about memory allocated,
        just the advertised structure.

        :param PlanarDiagram other_pd: Diagram to compare
        :return: Whether or not the diagrams are identical
        :rtype: bool
        """
        if (self.p.uid != other_pd.p.uid or
            <bytes>self.p.hash != <bytes>other_pd.p.hash):
            return False

        if (self.p.ncross != other_pd.p.ncross or
            self.p.nedges != other_pd.p.nedges or
            self.p.ncomps != other_pd.p.ncomps or
            self.p.nfaces != other_pd.p.nfaces):
            return False

        for i in range(self.p.ncross):
            if self.p.cross[i].sign != other_pd.p.cross[i].sign:
                return False
            for pos in range(4):
                if self.p.cross[i].edge[pos] != other_pd.p.cross[i].edge[pos]:
                    return False

        for i in range(self.p.nedges):
            if (self.p.edge[i].head != other_pd.p.edge[i].head or
                self.p.edge[i].headpos != other_pd.p.edge[i].headpos or
                self.p.edge[i].tail != other_pd.p.edge[i].tail or
                self.p.edge[i].tailpos != other_pd.p.edge[i].tailpos):
                return False

        for i in range(self.p.ncomps):
            if self.p.comp[i].tag != other_pd.p.comp[i].tag:
                return False
            if self.p.comp[i].nedges != other_pd.p.comp[i].nedges:
                return False
            for j in range(self.p.comp[i].nedges):
                if self.p.comp[i].edge[j] != other_pd.p.comp[i].edge[j]:
                    return False

        for i in range(self.p.nfaces):
            if self.p.face[i].nedges != other_pd.p.face[i].nedges:
                return False
            for j in range(self.p.face[i].nedges):
                if (self.p.face[i].edge[j] != other_pd.p.face[i].edge[j] or
                    self.p.face[i].ori[j] != other_pd.p.face[i].ori[j]):
                    return False

        return True

    def resize(self, pd_idx_t ncross):
        """resize(ncross)

        Resizes the :class:`PlanarDiagram` to fit more (or less, down to
        :obj:`ncross`) vertices
        """

        cdef pd_code_t *oldpd = self.p
        self.p = pd_copy_newsize(oldpd, ncross)
        pd_code_free(&oldpd)

    def add_crossing(self, N=1):
        cdef pd_idx_t old_n = self.p.ncross

        self.p.ncross += N

        if self.p.ncross > self.p.MAXVERTS:
            self.resize(self.p.ncross)

        for i in range(old_n, self.p.ncross):
            for pos in range(4):
                self.p.cross[i].edge[pos] = 0
            self.p.cross[i].sign = 0
        self.regenerate_py_os()

    def copy(self, thin=False):
        """copy() -> PlanarDiagram

        Returns a memory deepcopy of this PlanarDiagram

        :return: A deep memory copy of this PlanarDiagram
        :rtype: :class:`PlanarDiagram`
        """
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(self.__class__)
        newobj.p = pd_copy(self.p)
        newobj.thin = thin
        newobj.equiv_type = self.equiv_type
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
        """read(file f) -> new PlanarDiagram

        Read the next pdcode from a file object.

        :return: A new :class:`PlanarDiagram`, or :obj:`None` on failure.
        """

        if not isinstance(f, file):
            if isinstance(f, basestring):
                with open(f, "r") as real_file:
                    return cls.read(real_file, read_header=read_header, thin=thin)
            else:
                raise TypeError("read() first argument requires a file object or a path to a valid file")

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
        if newobj._hash:
            newobj.hashed = True
        newobj.regenerate_py_os()
        return newobj

    @classmethod
    def read_all(cls, f, read_header=False, thin=False):
        """read_all(file f) -> PlanarDiagram

        Returns a generator that iterates through the pdcodes stored
        in ``f``.
        """

        if not isinstance(f, file):
            if isinstance(f, basestring):
                with open(f, "r") as real_file:
                    return cls.read_all(real_file, read_header=read_header, thin=thin)
            else:
                raise TypeError("read_all() first argument requires a file object or a path to a valid file")

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
        """
        read_knot_theory(file f) -> PlanarDiagram

        This function reads a pdcode which was exported from the
        Mathematica package KnotTheory, exported as text with
        something like:

        ``Export["7_2.txt",PD[Knot[7,2]]]``

        These PD codes don't have component or face information, so
        that is all regenerated once the crossings have been loaded
        from the file.  This will only read one PD code from ``f``.
        """
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)

        if not isinstance(f, file):
            if isinstance(f, basestring):
                with open(f, "r") as real_file:
                    return cls.read_knot_theory(real_file, thin=thin)
            else:
                raise TypeError("read_knot_theory() first argument requires a file object or a path to a valid file")

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

        :return: A new diagram which representing the twist knot
        :rtype: PlanarDiagram
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
        """
        set_all_crossing_signs(signature)

        Set all signs of the crossings in this diagram to match
        signature; i.e. set ``self.crossing.sign[i] = signature[i]``.
        """
        for i, sign in enumerate(signature):
            self.p.cross[i].sign = sign

    def unset_crossing_signs(self):
        """
        unset_crossing_signs()

        Unset all signs on this PlanarDiagram to
        ``PD_UNSET_ORIENTATION``, i.e. make this diagram a shadow.
        """
        for x in self.crossings:
            x.sign = PD_UNSET_ORIENTATION

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

        Performs a Reidemeister 1 loop deletion which
        removes the input crossing.
        """
        ret = PlanarDiagram_wrap(pd_R1_loopdeletion(self.p, cr))
        if ret is None:
            raise Exception("Error in R1 loopdeletion")
        return ret

    def R1_loop_addition(self, pd_idx_t f, pd_idx_t e_on_f):
        return PlanarDiagram_wrap(pd_R1_loop_addition(self.p, f, e_on_f))

    def R2_bigon_elimination(self, Face f):
        """R2_bigon_elimination(crossing_or_index, crossing_or_index)
        -> (Upper pdcode, [Lower pdcode])

        Performs a Reidemeister 2 bigon elimination removing the bigon
        with vertices the two inputs. This is **not an in-place operation**
        as there may be more than one PlanarDiagram as a result of the move.
        """
        if len(f) != 2:
            raise Exception("Must R2 eliminate on a bigon, but f not a bigon")
        return self.R2_bigon_elimination_vertices(f[0].tail, f[0].head)

    def R2_bigon_addition(self, pd_idx_t f, pd_idx_t e1_on_f,
                          pd_idx_t e2_on_f, pd_or_t e1_over_e2_or):
        return PlanarDiagram_wrap(pd_R2_bigon_addition(self.p, f, e1_on_f, e2_on_f, e1_over_e2_or))

    def R2_bigon_elimination_vertices(self, pd_idx_t cr1, pd_idx_t cr2):
        cdef pd_idx_t *cr = [cr1, cr2]
        cdef pd_idx_t nout
        cdef pd_code_t **out_pds
        pd_R2_bigon_elimination(self.p, cr, &nout, &out_pds)
        if nout == 0:
            raise Exception("Error on R2 bigonelimination")
        return tuple(PlanarDiagram_wrap(out_pds[i]) for i in range(nout))

    def R3_triangle_flip(self, pd_idx_t f):
        return PlanarDiagram_wrap(pd_R3_triangle_flip(self.p, f))

    def connect_sum(self, Edge edge, PlanarDiagram other_pd, Edge other_edge):
        """
        connect_sum(Edge edge, PlanarDiagram other_pd, Edge other_edge) -> PlanarDiagram
        
        This function gives the connect sum of two PlanarDiagram objects at two different edges.
        """
        cdef pd_code_t* new_pd_p = pd_connect_sum(self.p, edge.index, other_pd.p, other_edge.index)
        return PlanarDiagram_wrap(new_pd_p)

    def _R3_strand_swap(self, Face f, Edge e):
        # Stub for wrapping C implementation of tangle slide
        cdef pd_idx_t *strand_edges = <pd_idx_t*>malloc(3 * sizeof(pd_idx_t))
        cdef pd_idx_t *border_faces = <pd_idx_t*>malloc(3 * sizeof(pd_idx_t))
        cdef pd_idx_t npieces = 0
        cdef pd_code_t **pd_pieces
        cdef pd_tangle_t *tangle = pd_tangle_new(4)
        cdef Crossing x = [x for x in f.get_vertices() if
                           x.index != e.head and x.index != e.tail][0]
        cdef int err = 0

        for i, e_i in enumerate(x.edges):
            tangle.edge[i] = e_i

        for i, k in enumerate(x.faces):
            tangle.face[i] = k
        pd_regenerate_tangle_err(self.p, tangle, &err)
        if err != 0:
            raise Exception("Error in creating tangle data structure")

        strand_edges[0] = e.prev_edge().index
        strand_edges[1] = e.index
        strand_edges[2] = e.next_edge().index

        if e.face_index_pos()[0][0] == f.index:
            posneg = 0
        else:
            posneg = 1

        border_faces[0] = e.prev_edge().face_index_pos()[posneg][0]
        border_faces[1] = f.index
        border_faces[2] = e.next_edge().face_index_pos()[posneg][0]

        #print "Face.. ",
        #print border_faces[0], border_faces[1], border_faces[2]

        err = 0
        pd_tangle_slide_err(self.p, tangle,
                            3, strand_edges, border_faces,
                            &npieces, &pd_pieces, &err)
        if err != 0:
            raise Exception("Error with tangle slide")

        ret = []
        for i in range(npieces):
            ret.append(PlanarDiagram_wrap(pd_pieces[i]))

        return ret

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

    def euler_characteristic(self):
        """euler_characteristic() -> int

        Return the Euler characteristic of the diagram."""
        return self.nfaces + self.ncross - self.nedges

    def homfly(self, as_string=False, timeout=None):
        """homfly([as_string=False, timeout=None]) -> HOMFLYPolynomial

        Compute the HOMFLY polynomial for this diagram. Returns None
        on error.
        """
        cdef char* homfly_from_pd
        cdef bytes homflybytes
        if timeout is None:
            homfly_from_pd = pd_homfly(self.p)
        else:
            homfly_from_pd = pd_homfly_timeout(self.p, timeout)
        if homfly_from_pd == NULL:
            return None
        else:
            homflybytes = copy_and_free(homfly_from_pd)

        if as_string:
            return homflybytes
        else:
            return HOMFLYPolynomial(homflybytes)

    def linking_number(self, unsigned int c1=0, unsigned int c2=0):
        """linking_number([pd_idx_t c1, [pd_idx_t c2]]) -> int

        Return the linking number of the two numbered components. If no
        arguments are supplied, checks the first component with
        itself. If one argument is supplied, the other defaults to the
        first component.

        """
        return pd_linking_number(self.p, c1, c2)

    def unsigned_linking_number(self, unsigned int c1=0, unsigned int c2=0):
        """linking_number([pd_idx_t c1, [pd_idx_t c2]]) -> int

        Return the linking number of the two numbered components. If no
        arguments are supplied, checks the first component with
        itself. If one argument is supplied, the other defaults to the
        first component.

        """
        return pd_unsigned_linking_number(self.p, c1, c2)

    def interlaced_crossings(self):
        """
        interlaced_crossings() -> tuple of int

        Calculate the interlaced crossings of this diagram.

        For a knot shadow, one can compute the average (over all possible
        crossing assignments) v2 invariant by v2=-1/4*C, where C is the number
        of interlaced crossings.

        :return: Tuple of :variable:`ncomps` ints
        """
        cdef int i
        cdef int *c_ret = pd_interlaced_crossings(self.p)
        ret = tuple(c_ret[i] for i in range(self.p.ncomps))
        free(c_ret)
        return ret

    def interlaced_crossings_unsigned(self):
        """
        interlaced_crossings() -> tuple of int

        Calculate the unsigned interlaced crossings of this diagram.

        :return: Tuple of :variable:`ncomps` ints
        """

        cdef int i
        cdef unsigned int *c_ret = pd_interlaced_crossings_unsigned(self.p)
        ret = tuple(c_ret[i] for i in range(self.p.ncomps))
        free(c_ret)
        return ret

    def unique_code(self):
        """unique_code() -> str

        Returns a 'pdcode-ish' string which uniquely identifies this
        planar diagram. If another planar diagram is not isomorphic to
        this, their :meth:`unique_code`\ s will differ, but two isomorphic
        planar diagrams may have different :meth:`unique_code`\ s.
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
        """db_knot(ncross, index, alternating=None) -> new PlanarDiagram

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
        """db_link(ncross, index, alternating=True) -> new PlanarDiagram

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
                if (len(verts) == 2
                    and ((verts[0].sign != verts[1].sign) or
                         (verts[0].sign == verts[1].sign == PD_UNSET_ORIENTATION))):
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
            yield self.R2_bigon_elimination_vertices(*bigon.vertices)

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

        Requires: :mod:`spherogram`"""
        pdcode = self.pdcode()

        from spherogram import links
        sg_xings = [links.Crossing(x.index) for x in self.crossings]

        gluings = defaultdict(list)
        cross_oriented = True
        for x_i, (xing, pd_xing, sg_xing) in enumerate(
                zip(self.crossings, pdcode, sg_xings)):
            if xing.sign == PD_NEG_ORIENTATION:
                sg_xing.sign = -1
            elif xing.sign == PD_POS_ORIENTATION:
                sg_xing.sign = 1
            else:
                sg_xing.sign = 0
                cross_oriented = False
            for pos, e_i in enumerate(pd_xing):
                gluings[e_i].append((sg_xing, pos))
        for (x, xpos), (y, ypos) in gluings.itervalues():
            x[xpos] = y[ypos]

        link = links.Link(sg_xings, check_planarity=False, build=False)

        # This was in Eric's original email and so I do not know if it needs
        # to be implemented here since I have passed build=False above
        if not cross_oriented:
            link._orient_crossings()

        component_starts = []
        for comp in self.components:

            # We need to find which CEP our first edge is; so let's take the crossing
            # at the head of the edge and then look at which entry point our first
            # edge is in our crossing's list of CEP's
            head_index = comp[0].head
            head = self.crossings[head_index]
            if comp[0].index == head.understrand_indices()[0]:
                component_starts.append( link.crossings[head_index].entry_points()[0] )

            elif comp[0].index == head.overstrand_indices()[0]:
                component_starts.append( link.crossings[head_index].entry_points()[1] )

            else:
                raise ValueError("Error in finding first edge in component {}".format(
                    self.components.index(comp)))

        # Now build components starting at the corresponding edges
        link._build_components(component_starts)

        # Lastly, check that our graph is indeed planar.
        if not link.is_planar():
            raise ValueError("Diagram not planar.")

        return link

    def snappy_manifold(self):
        """snappy_manifold() -> snappy.Manifold

        Returns a SnapPy Manifold object which represents this knot/link's
        complement in $S^3$.

        Requires: ``snappy``, ``spherogram``"""
        from snappy import Manifold
        sgm_link = self.as_spherogram()
        return Manifold("DT:%s"%sgm_link.DT_code())

    def edit_copy(self):
        """edit_copy() -> new PlanarDiagram

        Open a plink editor for this PlanarDiagram. On
        completion, returns a new PlanarDiagram object.

        Requires package ``plink`` and TKinter interface
        """
        from Tkinter import mainloop as _tk_mainloop
        editor = self.as_spherogram().view()
        _tk_mainloop()

        return self.from_plink(editor)

    def pdstor_index(self, pdstor_f, isotopy=True, read_header=False):
        """pdstor_index() -> int index

        Returns the index of this diagram's isotopy (or isomorphism, if isotopy=False)
        class in the opened, readable pdstor file pdstor_f.

        Treats pdstor_f's current position as beginning of file, and will re-seek
        file pointer back on completion. It may help to open the file as 'rb' to
        avoid inconsistencies with Windows.
        """
        bof = pdstor_f.tell()
        for i, pd in enumerate(
                self.__class__.read_all(pdstor_f, read_header=read_header)):
            if self.isomorphic(pd):
                pdstor_f.seek(bof)
                return i
        pdstor_f.seek(bof)
        return None

    @classmethod
    def random_diagram(cls, n_crossings, n_components=None, max_att=50, dia_type=None):
        u"""
        random_diagram(n_crossings, n_components=None, max_att=50,
        dia_type=None) -> new PlanarDiagram

        Create a new PlanarDiagram with n_crossings crossings by
        uniformly sampling a rooted 4-regular map using Giles
        Schaeffer's PlanarMap and uniformly sampling a sign at each
        crossing. If n_components is a positive integer rather
        than ``None``, then we will attempt up to max_att times to
        create a diagram with precisely n_components, and ``None``
        otherwise. If ``max_att=0`` then there is no limit on the
        number of attempts to satisfy n_components.

        Presently, dia_type is one of 'all', 'prime',
        '6conn', or 'biquart'. (``dia_type=None`` defaults to
        all)

        :param int n_crossings: The number of crossings of the result diagram
        :param n_components: The number of components of the result
           diagram, or None if no constraint
        :type n_components: int or None
        :param int max_att: The max number of attempts to create a
           diagram of :obj:`n_components` components.
        :param dia_type: The type of diagram to produce, or None for
           any. Valid types are 'all', 'prime', '6conn', and 'biquart',
           corresponding to 2-edge-connected (i.e. any) diagrams,
           4-connected (i.e. prime) diagrams, 6-connected diagrams,
           or bipartite diagrams. Notice that bipartite diagrams
           never have precisely 1 component.
        :type dia_type: str or None
        """

        import os
        cdef pmMap plmap
        cdef pmSize size
        cdef pmMethod meth
        cdef pmMemory mem
        cdef pm_edge *cur_e
        cdef pm_vertex *cur_v
        cdef PlanarDiagram newobj = PlanarDiagram.__new__(cls)
        cdef pd_code_t *pd = NULL

        ## Initialize "size", which determines what type of map to sample
        size.n_edges = 0
        size.n_verts = n_crossings
        size.n_faces = 0

        # We might actually want to use these now?
        size.min_loop_comps = 0
        size.max_loop_comps = 0

        # We don't currently do anything with these parameters
        size.r = 0 # #"red"
        size.g = 0 # #"green"
        size.d = 0 # something to do with degree on other types of maps
        size.t = 0 # allowed error on size

        size.dgArr = NULL

        if dia_type is None or dia_type == 'all':
            # 2-edge-connected quartic map
            size.map_type = PM_MAP_TYPE_QUART_2C
            size.basic_type = PM_BASIC_TYPE_QUART_2C

        elif dia_type == 'prime':
            # 4-edge-connected quartic map
            size.map_type = 5
            size.basic_type = 5

        elif dia_type == '6conn':
            # 6-edge-connected quartic map
            size.map_type = 6
            size.basic_type = 5

        elif dia_type == 'biquart':
            # bi-quartic map
            size.map_type = 9
            size.basic_type = 9

        else:
            raise Exception("dia_type is not an expected value")

        ## Initialize some parameters of the sampling method (namely, seed)
        meth.core = 0
        meth.pic = 0
        meth.seed = int(os.urandom(5).encode('hex'), 16) # TODO: Make randomish, 0 for testing
        meth.verbose = 0

        if not pmInitRND(&meth):
            raise Exception("Failure during init RND")
        if not pmSetParameters(&size, &meth):
            raise Exception("Failure during set size")

        if not pmMemoryInit(&size, &meth, &mem):
            raise Exception("Failure during memory init")
        if not pmExtendMemory(&size, &meth, &mem, 0):
            raise Exception("Failure during memory extend")

        att_N = 0
        while (pd == NULL or
               (pd != NULL and (n_components is not None and n_components != pd.ncomps))):

            if not pmPlanMap(&size, &meth, &mem, &plmap):
                raise Exception("Failure during map generation")

            map_nv = plmap.v

            if pd != NULL:
                pd_code_free(&pd)
            pd = pd_code_new(map_nv+2)

            pd.ncross = map_nv
            pd.nedges = map_nv*2
            pd.nfaces = map_nv+2

            cur_v = plmap.root.from_v
            v_idx = cur_v.label-1
            #print cur_v.label,
            cur_e = cur_v.root
            pos = 0
            while cur_e != cur_v.root.prev_e:
                e_idx = abs(cur_e.label)-1
                pd.cross[v_idx].edge[pos] = e_idx
                #print cur_e.label,
                cur_e = cur_e.next_e
                pos += 1
            e_idx = abs(cur_e.label)-1
            pd.cross[v_idx].edge[pos] = e_idx
            #print cur_e.label

            while cur_v.next_v != NULL:
                cur_v = cur_v.next_v
                v_idx = cur_v.label-1
                #print cur_v.label,
                cur_e = cur_v.root
                pos = 0
                while cur_e != cur_v.root.prev_e:
                    e_idx = abs(cur_e.label)-1
                    pd.cross[v_idx].edge[pos] = e_idx
                    #print cur_e.label,
                    cur_e = cur_e.next_e
                    pos += 1
                e_idx = abs(cur_e.label)-1
                pd.cross[v_idx].edge[pos] = e_idx
                #print cur_e.label
            pmFreeMap(&plmap)

            pd_regenerate_edges(pd)
            pd_regenerate(pd)
            att_N += 1
            meth.seed += 1
            if max_att and att_N > max_att:
                pd_code_free(&pd)
                return None

        uniform_mask = [random.randint(0,1) for _ in range(map_nv)]
        newobj.p = pd
        newobj.set_all_crossing_signs(uniform_mask)

        return newobj

    def randomly_assign_crossings(self):
        """
        randomly_assigned_crossings()

        Randomly assign a sign in {+,-} to each of the crossings of
        this PlanarDiagram.
        """
        uniform_mask = [random.randint(0,1) for _ in range(self.ncross)]
        self.set_all_crossing_signs(uniform_mask)

    @classmethod
    def from_dt_code(cls, dt_code):
        """from_dt_code() -> new PlanarDiagram

        Creates a new PlanarDiagram from a DT code (Dowker-Thistlethwaite).
        Requires :mod:`Spherogram` package

        Warning: Currently depends on :meth:`from_pdcode`, which has bugs with
        components of length 2

        Requires: :mod:`spherogram`"""
        import spherogram
        return cls.from_pdcode(spherogram.DTcodec(dt_code).PD_code())

    @classmethod
    def from_editor(cls, **kwargs):
        """
        from_editor() -> new PlanarDiagram

        Creates a new PlanarDiagram using a new :mod:`plink`
        :class:`LinkEditor`

        Requires package :mod:`plink` and Tkinter interface
        """
        from plink import LinkEditor
        from Tkinter import mainloop as _tk_mainloop

        editor = LinkEditor(**kwargs)
        _tk_mainloop()

        return cls.from_plink(editor)

    @classmethod
    def from_plink(cls, editor, thin=False):
        """
        from_plink(editor) -> new PlanarDiagram

        Creates a new :class:`PlanarDiagram` from a :mod:`plink`
        :class:`LinkManager` object.

        Requires package :mod:`plink`
        """

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

        crossing_indices = dict()
        cross_idx = 0
        edge_idx = 0
        for component in components:
            # Store this edge index for when we come back around
            first_edge = edge_idx

            # Iterate over "ECrossings"; an ECrossing is an instance of
            # when an arrow runs through a crossing, either over or under
            # as we go along the component
            for k, ec in enumerate(component):
                # Either add a new crossing index if we haven't run over
                # this crossing, or get the stored value
                if ec.crossing not in crossing_indices:
                    crossing_indices[ec.crossing] = cross_idx
                    cross_idx += 1
                x_i = crossing_indices[ec.crossing]

                # Previous edge can always be the 'edge_idx', but next edge
                # has to be careful; if we are at the end of the component
                # loop back to the first edge of the component
                prev_edge = edge_idx
                edge_idx += 1
                next_edge = edge_idx if k != len(component) - 1 else first_edge

                # Remember that PLink uses the following sign convention
                #
                # RH  ^       LH  ^
                #     |           |
                #  ------>     <------
                #     |           |
                # (+) |       (-) |
                #
                # For now, be PDCode-esque: Under strand goes 0-->2, over strand
                # depends on sign
                if ec.goes_over():
                    if ec.crossing.sign() == "RH":
                        in_pos, out_pos = 3, 1
                        pd.cross[x_i].sign = PD_POS_ORIENTATION
                    else:
                        in_pos, out_pos = 1, 3
                        pd.cross[x_i].sign = PD_NEG_ORIENTATION
                else:
                    in_pos, out_pos = 0, 2

                # Plug in edges and crossings; it is very important to actually
                # set the edge here so that pd_regenerate_edges doesn't get
                # called and discard all of our precious sign information
                pd.cross[x_i].edge[in_pos] = prev_edge
                pd.edge[prev_edge].head = x_i
                pd.edge[prev_edge].headpos = in_pos

                pd.cross[x_i].edge[out_pos] = next_edge
                pd.edge[next_edge].tail = x_i
                pd.edge[next_edge].tailpos = out_pos

        # Regenerate the object to load in faces and components
        # Notice that we could in theory just regenerate faces, and
        # use component information from the Plink editor for components
        newobj.regenerate()
        return newobj

    def ccode(self):
        """
        ccode() -> str

        Convert this PlanarDiagram to its Ewing-Millet ccode
        format. Wraps :c:func:`pdcode_to_ccode`.
        """

        return copy_and_free(pdcode_to_ccode(self.p))

    def serialize(self):
        """serialize() -> str

        Serialize this pdcode into a string matching the PDCode spec
        found in src/plcTopology.proto. Requires protobuf-python.
        """
        import plcTopology_pb2
        pb_code = plcTopology_pb2.PDCode()

        pb_code.uid = self.p.uid
        if <bytes>self.p.hash:
            pb_code.hash = <bytes>self.p.hash

        for i in range(self.p.ncross):
            pb_cross = pb_code.crossings.add()
            for pos in range(4):
                pb_cross.edges.append(self.p.cross[i].edge[pos])
            pb_cross.sign = self.p.cross[i].sign

        for i in range(self.p.nedges):
            pb_edge = pb_code.edges.add()
            pb_edge.head = self.p.edge[i].head
            pb_edge.headpos = self.p.edge[i].headpos
            pb_edge.tail = self.p.edge[i].tail
            pb_edge.tailpos = self.p.edge[i].tailpos

        for i in range(self.p.ncomps):
            pb_comp = pb_code.components.add()
            pb_comp.tag = self.p.comp[i].tag
            for j in range(self.p.comp[i].nedges):
                pb_comp.edges.append(self.p.comp[i].edge[j])

        for i in range(self.p.nfaces):
            pb_face = pb_code.faces.add()
            for j in range(self.p.face[i].nedges):
                pb_face_edge = pb_face.edges.add()
                pb_face_edge.edge = self.p.face[i].edge[j]
                pb_face_edge.sign = self.p.face[i].ori[j]

        return pb_code.SerializeToString()

    @classmethod
    def deserialize(cls, basestring buff):
        """deserialize(buffer) -> PlanarDiagram

        Deserialize this pdcode from a string matching the PDCode spec
        found in src/plcTopology.proto. Requires protobuf-python.

        Caveats: tries to reproduce a PlanarDiagram exactly from the
        serialization, hence will not regenerate() itself.
        """
        cdef PlanarDiagram pd = cls.__new__(cls)

        import plcTopology_pb2
        pb_code = plcTopology_pb2.PDCode.FromString(buff)

        pd.p = pd_code_new(len(pb_code.crossings)+2)
        pd.p.uid = pb_code.uid
        if pb_code.hash:
            strncpy(pd.p.hash, <bytes>pb_code.hash, PD_HASHSIZE)

        pd.p.ncross = len(pb_code.crossings)
        for i,pb_cross in enumerate(pb_code.crossings):
            for pos,edge in enumerate(pb_cross.edges):
                pd.p.cross[i].edge[pos] = edge
            pd.p.cross[i].sign = pb_cross.sign

        pd.p.nedges = len(pb_code.edges)
        for i,pb_edge in enumerate(pb_code.edges):
            pd.p.edge[i].head = pb_edge.head
            pd.p.edge[i].headpos = pb_edge.headpos
            pd.p.edge[i].tail = pb_edge.tail
            pd.p.edge[i].tailpos = pb_edge.tailpos

        pd.p.ncomps = len(pb_code.components)
        for i,pb_comp in enumerate(pb_code.components):
            pd.p.comp[i].tag = pb_comp.tag
            pd.p.comp[i].nedges = len(pb_comp.edges)
            pd.p.comp[i].edge = <pd_idx_t*>calloc(len(pb_comp.edges), sizeof(pd_idx_t))
            for j,pb_edge_idx in enumerate(pb_comp.edges):
                pd.p.comp[i].edge[j] = pb_edge_idx

        pd.p.nfaces = len(pb_code.faces)
        for i,pb_face in enumerate(pb_code.faces):
            pd.p.face[i].nedges = len(pb_face.edges)
            pd.p.face[i].edge = <pd_idx_t*>calloc(len(pb_face.edges), sizeof(pd_idx_t))
            pd.p.face[i].ori = <pd_or_t*>calloc(len(pb_face.edges), sizeof(pd_or_t))
            for j,pb_face_edge in enumerate(pb_face.edges):
                pd.p.face[i].edge[j] = pb_face_edge.edge
                pd.p.face[i].ori[j] = pb_face_edge.sign

        return pd

    def pdcode(self):
        """
        pdcode() -> list of ints

        Create a planar diagram code as a python list.
        """
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
        """
        from_pdcode(pdcode) -> new PlanarDiagram

        Create a PlanarDiagram from a planar diagram code list.
        """
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
        self.edges._disown()
        self.edges = _EdgeList(self)
        self.components._disown()
        self.components = _ComponentList(self)
        self.crossings._disown()
        self.crossings = _CrossingList(self)
        self.faces._disown()
        self.faces = _FaceList(self)

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

cdef pd_code_t **pd_simplify(pd_code_t *pd, int *ndias):
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
