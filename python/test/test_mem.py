import unittest
from itertools import product, compress, izip
from libplcurve.pd import PlanarDiagram
import libplcurve.plcurve as pl
from fractions import gcd
from math import pi,cos,sin

class TestPlCurveMemoryManagement(unittest.TestCase):
    def setUp(self):
        self.rng = pl.RandomGenerator()

    def test_add_drop_cmp_ownership(self):
        L = pl.PlCurve.random_closed_polygon(self.rng, 500)
        C = L.components[0]
        print
        C[-1][1] = 4
        print C.vertices
        self.assertFalse(C.own())
        L.drop_component(0)
        print C.vertices
        print L
        self.assertTrue(C.own())

    def test_del_plcurve_cmp_ownership(self):
        L = pl.PlCurve.random_closed_polygon(self.rng, 500)
        C = L.components[0]
        print
        C[-1][1] = 4
        print C.vertices
        self.assertFalse(C.own())
        del L
        print C.vertices
        self.assertTrue(C.own())

    def test_append_drop_components(self):
        L = pl.PlCurve.random_closed_polygon(self.rng, 50)
        L.append(L.components[0][10:40])
        L.append(L.components[0][10:30])
        #L.append(L.components[0][10:20])
        L.drop_component(2)
        #L.drop_component(1)
        #L.drop_component(0)
        #L.drop_component(-1)

    def test_append_components_2(self):
        L = pl.PlCurve.random_closed_polygon(self.rng, 60)
        L.append(L[0][1:5])
        L.append(L[0][2:4])
        L.append(L[0][3:4])
        L.drop_component(0)
        L.drop_component(2)
        L.drop_component(1)

    def test_append_components(self):
        L = pl.PlCurve.random_closed_polygon(self.rng, 50)
        L.append(L.components[0][10:40])
        L.append(L.components[0][10:30])
        L.append(L.components[0][10:20])
        L.append(L.components[0][10:12])
        L.append(L.components[0][10:15])

if __name__=="__main__":
    unittest.main()
