import unittest
from itertools import product, compress, izip
from libpl.pdcode import PlanarDiagram
from test_pdcode import PlanarDiagramAssertMixin

class TestPDCode(PlanarDiagramAssertMixin, unittest.TestCase):
    def setUp(self):
        pass

    @unittest.skip("Requires protobuf from google")
    def test_serialize(self):
        test_cases = [
            PlanarDiagram.torus_knot(2,7),
            PlanarDiagram.unknot_wye(2,3,8),
            PlanarDiagram.unknot(0)
        ]
        for K in test_cases:
            K.serialize()
            L = PlanarDiagram.deserialize(K.serialize())
            L.regenerate_hash()
            print L.hash == K.hash
            print repr(L), repr(K)
            print K.isotopic(L)

            #L = PlanarDiagram.unknot(0)
            #self.assertEqual(K, L)

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
