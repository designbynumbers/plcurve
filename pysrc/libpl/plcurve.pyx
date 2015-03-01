from __future__ import division
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy
import numpy as np
cimport numpy as np
np.import_array()

from cython cimport view
from cython cimport numeric, floating, integral
from plcurve cimport *
from .pdcode.diagram cimport PlanarDiagram, PlanarDiagram_wrap
# if you MUST use Cython < 0.21, comment the above line and uncomment this.
#from libpl.pdcode.diagram cimport PlanarDiagram, PlanarDiagram_wrap

cdef plc_vector new_plc_vector(list l):
    return plc_build_vect(l[0], l[1], l[2])
cdef np.ndarray plc_vector_as_pyo(plc_vector *v):
    cdef np.npy_intp dim = 3
    return np.PyArray_SimpleNewFromData(1, &dim, np.NPY_DOUBLE, v.c)

cdef np.ndarray plc_vector_as_pyo_copy(plc_vector *v):
    cdef plc_vector *vc = <plc_vector*>memcpy(malloc(sizeof(plc_vector)), v, sizeof(plc_vector))
    cdef np.ndarray ret = plc_vector_as_pyo(vc)
    np.PyArray_UpdateFlags(ret, np.NPY_OWNDATA)
    return ret

cdef class RandomGenerator:
    def __cinit__(self):
        gsl_rng_env_setup()
        self.p = gsl_rng_alloc(gsl_rng_default)

    def set(self, seed):
        gsl_rng_set(self.p, seed)

cdef class Component:
    cdef update_array(self):
        self.vertex_array = view.array((self.p.nv,3), itemsize=sizeof(double), format="d", allocate_buffer=False)
        self.vertex_array.data = <char*>self.p.vt

    property vertices:
        def __get__(self):
            return self.vertex_array.memview

cdef Component Component_FromPointer(plc_strand *ptr):
    cdef Component ret = Component.__new__(Component)
    ret.p = ptr
    ret.update_array()
    return ret

cdef class PlCurve:
    """
    It's a PlCurve.
    """

    def __cinit__(self):
        self.own = False
    def __init__(self, component_spec=None):
        """Create a new PlCurve, optionally with empty components described by `component_spec`.

        If it is not `None`, `component_spec` is a sequence of the form
        (int num_vertices, bool open, int num_colors)."""
        self.own = True
        if component_spec is None:
            self.p = <plCurve*>malloc(sizeof(plCurve))
            self.p.nc = 0;
            self.p.cp = NULL;
            self.p.cst = NULL;
            self.p.quant = NULL;
            self.p.G = NULL;
        else:
            raise NotImplementedError("Oops")

    def __dealloc__(self):
        if self.p is not NULL and self.own:
            # really augmented free
            plc_free(self.p)

    def __len__(self):
        return self.p.nc

    property components:
        def __get__(self):
            if self._components is None:
                self.generate_py_components()
            return tuple(self._components)

    cdef void generate_py_components(self):
        self._components = []
        for i in range(self.p.nc):
            self._components.append(Component_FromPointer(&self.p.cp[i]))

    # Data Operations
    def add_component(self, component):
        """add_component(component)

        Add a component to this PlCurve."""
        cdef plc_vector *verts = <plc_vector*>malloc(len(component) * sizeof(plc_vector))
        for i, vert in enumerate(component):
            verts[i].c[0], verts[i].c[1], verts[i].c[2] = vert[0], vert[1], vert[2]
        plc_add_component(<plCurve*>self.p, 0, len(component), 0, 0, verts, NULL)

    def fix_wrap(self):
        """fix_wrap()

        """
        plc_fix_wrap(self.p)

    # Spline package

    # Geometric information
    property num_edges:
        """Number of edges in this PlCurve."""
        def __get__(self):
            return plc_num_edges(self.p)
    property num_vertices:
        """Number of vertices in this PlCurve."""
        def __get__(self):
            return plc_num_verts(self.p)

    def turning_angle(self, const int component, const int vertex):
        cdef bool no_error = True
        cdef double ret = plc_turning_angle(self.p, component, vertex, &no_error)
        if no_error:
            return ret
        else:
            raise Exception("Error in turning_angle!")
    def MR_curvature(self, const int component, const int vertex):
        return plc_MR_curvature(self.p, component, vertex)
    def total_curvature(self, get_component_tc=False):
        cdef double *component_tc
        if get_component_tc:
            raise NotImplementedError("Oops")
        else:
            return plc_totalcurvature(self.p, NULL)
    def total_torsion(self, get_component_tt=False):
        cdef double *component_tt
        if get_component_tt:
            raise NotImplementedError("Oops")
        else:
            return plc_totaltorsion(self.p, NULL)

    property pointset_diameter:
        """Diameter of the points"""
        def __get__(self):
            return plc_pointset_diameter(self.p)
    property center_of_mass:
        """Center of mass."""
        def __get__(self):
            cdef plc_vector v = plc_center_of_mass(self.p)
            return plc_vector_as_pyo_copy(&v)
    property gyradius:
        """Gyradius of this PlCurve."""
        def __get__(self):
            return plc_gyradius(self.p)

    # Geometric operations
    def perturb(self, floating radius):
        """perturb(radius)

        Randomly perturb this PlCurve"""
        plc_perturb(self.p, radius)

    # Random polygon library
    @classmethod
    def random_closed_polygon(cls, n_edges, RandomGenerator rng):
        """
        random_closed_polygon(int n_edges, RandomGenerator rng) -> PlCurve

        Create a random closed polygon.
        """
        cdef PlCurve plc = PlCurve.__new__(cls)
        plc.own = True
        plc.p = plc_random_closed_polygon(rng.p, n_edges)
        return plc
    @classmethod
    def random_open_polygon(cls, n_edges, RandomGenerator rng):
        """
        random_open_polygon(int n_edges, RandomGenerator rng) -> PlCurve

        Create a random open polygon.
        """
        cdef PlCurve plc = PlCurve.__new__(cls)
        plc.own = True
        plc.p = plc_random_open_polygon(rng.p, n_edges)
        return plc

    @classmethod
    def random_closed_plane_polygon(cls, n_edges, RandomGenerator rng):
        """
        random_closed_plane_polygon(int n_edges, RandomGenerator rng) -> PlCurve

        Create a random closed plane polygon.
        """
        cdef PlCurve plc = PlCurve.__new__(cls)
        plc.own = True
        plc.p = plc_random_closed_plane_polygon(rng.p, n_edges)
        return plc
    @classmethod
    def random_open_plane_polygon(cls, n_edges, RandomGenerator rng):
        """
        random_open_plane_polygon(int n_edges, RandomGenerator rng) -> PlCurve

        Create a random open plane polygon.
        """
        cdef PlCurve plc = PlCurve.__new__(cls)
        plc.own = True
        plc.p = plc_random_open_plane_polygon(rng.p, n_edges)
        return plc

    @classmethod
    def random_equilateral_closed_polygon(cls, n_edges, RandomGenerator rng):
        """
        random_equilateral_closed_polygon(int n_edges, RandomGenerator rng) -> PlCurve

        Generate a random equilateral polygon using the CSU algorithm.
        """
        cdef PlCurve plc = PlCurve.__new__(cls)
        plc.own = True
        plc.p = plc_random_equilateral_closed_polygon(rng.p, n_edges)
        return plc

    @classmethod
    def random_equilateral_open_polygon(cls, n_edges, RandomGenerator rng):
        """random_equilateral_open_polygon(int n_edges, RandomGenerator rng) -> PlCurve

        Create a random equilateral open polygon.
        """
        cdef PlCurve plc = PlCurve.__new__(cls)
        plc.own = True
        plc.p = plc_random_equilateral_open_polygon(rng.p, n_edges)
        return plc

    # Symmetry Functions

    # Topology Functions
    def projection_to_PlanarDiagram(self, RandomGenerator rng=None):
        """
        projection_to_PlanarDiagram(RandomGenerator rng) -> new PlanarDiagram

        Creates a new PlanarDiagram from a projection of this PlCurve.
        """
        if rng is None:
            rng = RandomGenerator()
        return PlanarDiagram_wrap(pd_code_from_plCurve(rng.p, self.p))
