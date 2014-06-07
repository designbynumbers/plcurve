import unittest
from itertools import product, compress, izip
from libplcurve.pd import PlanarDiagram
import libplcurve.plcurve as pl
from fractions import gcd
from math import pi,cos,sin

class ArcPresentation(object):
    def __init__(self):
        pass

    class Component(object):
        def __init__(self):
            pass

def plCurve_from_arcpres(AP):
    L = pl.PlCurve()
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
    knot = pl.PlCurve()
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
    L = pl.PlCurve()
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

class TestPlCurveTopology(unittest.TestCase):
    def setUp(self):
        self.rng = pl.RandomGenerator()

    def check_torus_knot(self, verts, q):
        # Check that a torus knot plCurve produces a diagram
        # isomorphic to an actual torus knot
        L = torus_knot(verts, q, 2, 5, 2)
        proj_pd = L.as_pd(self.rng)
        expc_pd = PlanarDiagram.from_torus_knot(2,q)
        self.assertTrue(proj_pd.isomorphic(expc_pd))

    def check_torus_knot_rotation(self, verts, q):
        # Check that a torus knot plCurve produces a diagram
        # isomorphic to an actual torus knot, independent of
        # vertex numbering shifts
        L = torus_knot(verts, q, 2, 5, 2)
        for ofs in range(verts):
            for component in L:
                for j in range(len(component)):
                    component[j] = component[j+1]
            L.fix_wrap()
            proj_pd = L.as_pd(self.rng)
            expc_pd = PlanarDiagram.from_torus_knot(2,q)
            self.assertTrue(proj_pd.isomorphic(expc_pd))

    def test_torus_knot(self):
        self.check_torus_knot(150, 4)
        self.check_torus_knot(550, 8)

    def test_torus_knot_rotation(self):
        self.check_torus_knot_rotation(150, 3)
        self.check_torus_knot_rotation(150, 4)

    def test_unknot_component(self):
        L = equilateral_ngon(50)
        proj_pd = L.as_pd(self.rng)
        self.assertTrue(proj_pd.is_ok())
        self.assertEqual(proj_pd.ncross, 1)

    def test_hopf_link(self):
        L = basic_hopf()
        proj_pd = L.as_pd(self.rng)
        self.assertTrue(proj_pd.is_ok())

if __name__=="__main__":
    unittest.main()
