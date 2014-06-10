import unittest
from itertools import product, compress, izip
from libplcurve.plcurve import PlCurve, RandomGenerator

class TestPlCurve(unittest.TestCase):
    def setUp(self):
        pass

    def _rcp(self, seed, N):
        r = RandomGenerator()
        r.set(seed)
        return PlCurve.random_closed_polygon(r, N)

    def test_center_of_mass(self):
        L = self._rcp(3456, 50)
        com = L.center_of_mass
        print com.__repr__()
        del L
        print com

if __name__=="__main__":
    unittest.main()
