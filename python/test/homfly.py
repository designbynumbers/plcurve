import unittest
from itertools import product, compress
from libplcurve.pd import PlanarDiagram

class TestCCode(unittest.TestCase):
    def setUp(self):
        pass

    def check_all_signs(self, pd):
        # Check that ccode does not croak on any crossing signatures for pd
        for mask in product([0,1], repeat=len(pd.crossings)):
            temp_pd = pd.copy()
            for crossing in compress(temp_pd.crossings, mask):
                crossing.sign = (crossing.sign+1)%2
            self.assertIsNotNone(temp_pd.ccode())

    def test_trefoil_ccode(self):
        # Check that a trefoil diagram produces the correct ccode
        trefoil_pd = PlanarDiagram.from_torus_knot(2,3)
        self.assertEqual(trefoil_pd.ccode(),
                         "1+2b3a3d2c\n"
                         "2+3b1a1d3c\n"
                         "3+1b2a2d1c\n"
                         "\n\n")

    def test_unknot_generation(self):
        # Generate codes for unknots with 2-10 + crossings
        for pd in (PlanarDiagram.from_unknot(k) for k in range(2,11)):
            self.assertIsNotNone(pd.ccode())

    def test_unknot_signs(self):
        # Test all crossing signs for 5-crossing diagram
        pd = PlanarDiagram.from_unknot(5)
        self.check_all_signs(pd)

        # Test all crossing signs for 4-crossing diagram
        pd = PlanarDiagram.from_unknot(4)
        self.check_all_signs(pd)

class TestHomfly(unittest.TestCase):
    def setUp(self):
        pass

    def check_all_signs_unknot(self, pd):
        # Check that homfly polynomial is 1 regardless of signs for unknot pd
        for mask in product([0,1], repeat=len(pd.crossings)):
            temp_pd = pd.copy()
            for crossing in compress(temp_pd.crossings, mask):
                crossing.sign = (crossing.sign+1)%2
            self.assertEqual(temp_pd.homfly(), "1")

    def test_unknot_homfly(self):
        # Test unknot with simple twists
        for pd in (PlanarDiagram.from_unknot(k) for k in range(2,7)):
            self.check_all_signs_unknot(pd)

        # Test unknot (wye/plectonemes)
        for pd in (PlanarDiagram.from_unknot_wye(*p) for
                   p in product(range(2,4), repeat=3)):
            self.check_all_signs_unknot(pd)

if __name__=="__main__":
    unittest.main()
