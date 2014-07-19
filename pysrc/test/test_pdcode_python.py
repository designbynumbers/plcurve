import unittest
from itertools import product, compress, izip
from libpl.pdcode import PlanarDiagram

class TestPDCode(unittest.TestCase):
    def setUp(self):
        pass

    def test_pickle(self):
        import cPickle
        pdc = PlanarDiagram.from_torus_knot(2, 11)
        pdc2 = PlanarDiagram.from_torus_knot(2, 7)
        pdc_new = cPickle.loads(cPickle.dumps(pdc))
        self.assertTrue(pdc.isomorphic(pdc_new))
        self.assertFalse(pdc2.isomorphic(pdc_new))
        self.assertEqual(pdc.homfly(), pdc_new.homfly())

if __name__=="__main__":
    unittest.main()
