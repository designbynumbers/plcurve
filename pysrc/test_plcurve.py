import unittest
from itertools import product, compress, izip
from math import *
from fractions import gcd
from libpl.plcurve import PlCurve, RandomGenerator

import numpy as np
from numpy.testing import assert_array_almost_equal

class ArcPresentation(object):
    def __init__(self):
        pass

    class Component(object):
        def __init__(self):
            pass

def plCurve_from_arcpres(AP):
    L = PlCurve()
    for component in AP.cp:
        verts = []
        for i, page in enumerate(component.page):
            verts.extend((
                (5*cos(page*2*pi/AP.pages),
                 5*sin(page*2*pi/AP.pages),
                 component.startlevel[i]),
                (5*cos(page*2*pi/AP.pages),
                 5*sin(page*2*pi/AP.pages),
                 component.startlevel[(i+1)%component.narcs]),
                (0,0,component.startlevel[(i+1)%component.narcs])
            ))
        L.add_component(verts)
    return L

def torus_knot(verts, p, q, major_r, minor_r):
    tstep = 2*pi/verts
    cp = gcd(p,q)
    knot = PlCurve()
    for j in range(cp):
        pofs = j*(2*pi/cp)
        vertices = []
        for vtx in range(verts):
            theta = vtx*tstep
            pangle = pofs + (p*1.0/cp)*theta
            qangle = (q*1.0/cp)*theta
            vertices.append((
                major_r*cos(qangle)*(1+(minor_r*1.0/major_r)*cos(pangle)),
                major_r*sin(qangle)*(1+(minor_r*1.0/major_r)*cos(pangle)),
                minor_r*sin(pangle)))
        knot.add_component(vertices)
    knot.fix_wrap()
    return knot

def equilateral_ngon(n):
    L = PlCurve()
    theta = pi*2/n
    r = (1.0/2)/sin(theta/2.0)
    vertices = []
    for i in range(n):
        vertices.append((r*cos(i*theta), r*sin(i*theta), 0))
    L.add_component(vertices)
    return L

def basic_hopf():
    AP = ArcPresentation()
    AP.pages = 4

    AP.cp = [ArcPresentation.Component(), ArcPresentation.Component()]
    AP.cp[0].narcs = 2
    AP.cp[0].page = [0, 2]
    AP.cp[0].startlevel = [0,2]

    AP.cp[1].narcs = 2
    AP.cp[1].page = [1, 3]
    AP.cp[1].startlevel = [1, 3]
    return plCurve_from_arcpres(AP)


class TestPlCurve(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.eq_ngon_nv = 10
        cls.eq_ngon = equilateral_ngon(cls.eq_ngon_nv)

        cls.hopf = basic_hopf()
        cls.torusknot = torus_knot(30, 2, 6, 2, 3)

    def setUp(self):
        pass

    def _rcp(self, seed, N):
        r = RandomGenerator()
        r.set(seed)
        return PlCurve.random_closed_polygon(r, N)

    def test_center_of_mass(self):
        # Here's a dumb test for center of mass
        L = self._rcp(3456, 50)
        assert_array_almost_equal(
            L.center_of_mass,
            np.array([-0.00134734, -0.02689932,  0.07941123]))

    @unittest.skip("test stub")
    def test_from_file(self):
        pass

    @unittest.skip("test stub")
    def test_write(self):
        pass

    def test_random_generators(self):
        r = RandomGenerator()
        gens = (
            PlCurve.random_closed_polygon,
            PlCurve.random_open_polygon,
            PlCurve.random_closed_plane_polygon,
            PlCurve.random_open_plane_polygon,
            PlCurve.random_equilateral_closed_polygon,
            PlCurve.random_equilateral_open_polygon
        )
        for pl_gen in gens:
            self.assertEqual(len(pl_gen(r, 50)),
                             1)

    def test_index_lengths(self):
        self.assertEqual(self.eq_ngon.num_vertices,
                         10)
        self.assertEqual(self.eq_ngon.num_vertices,
                         self.eq_ngon_nv)
        self.assertEqual(self.eq_ngon.num_edges,
                         self.eq_ngon_nv)
        self.assertEqual(len(self.eq_ngon),
                         1)
        self.assertEqual(len(self.eq_ngon.components),
                         1)

    def test_geometric_properties(self):
        self.assertAlmostEqual(self.eq_ngon.turning_angle(0, self.eq_ngon_nv//2),
                               0.6283185307179586)
        self.assertAlmostEqual(self.eq_ngon.MR_curvature(0, self.eq_ngon_nv//2),
                               0.6498393924658127)
        self.assertAlmostEqual(self.eq_ngon.total_curvature(),
                               6.283185307179586)
        self.assertAlmostEqual(self.torusknot.total_torsion(),
                               19.654663669424494)


if __name__=="__main__":
    unittest.main()
