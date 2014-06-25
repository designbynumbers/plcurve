import unittest
from itertools import product, compress, izip
from libpl.pdcode import PlanarDiagram

class TestPDCode(unittest.TestCase):
    def setUp(self):
        pass

    def check_num_components(self, pd, check_len):
        # Ensure the number of components in a pd is check_len
        self.assertEqual(pd.ncomps, check_len)

        # Ensure that all different component lengths agree
        self.assertEqual(pd.ncomps, len(pd.components))

    def check_component_lengths(self, pd, lengths):
        lengths = reversed(sorted(lengths))
        for component, length in izip(pd.components, lengths):
            self.assertEqual(len(component), length)
            self.assertEqual(len(component), component.nedges)

    def check_face_lengths(self, pd, lengths):
        lengths = reversed(sorted(lengths))
        for face, length in izip(pd.faces, lengths):
            self.assertEqual(face.nedges, length)

    def check_twist_knot(self, n):
        # Generate a twist knot and make some sanity checks
        pd = PlanarDiagram.from_twist_knot(n)
        self.assertIsNotNone(pd)
        self.assertTrue(pd.is_ok())
        self.check_num_components(pd, 1)

        # Test the number of edges in each component
        self.check_component_lengths(pd, [pd.nedges])
        # Test the number of edges on each face
        face_lengths = [2, 3]
        face_lengths.extend([2]*(n-1))
        face_lengths.extend([3,n+1,n+1])
        self.check_face_lengths(pd, face_lengths)

    def check_torus_2q(self, q):
        pd = PlanarDiagram.from_torus_knot(2, q)
        # Sanity check
        self.assertIsNotNone(pd)
        self.assertTrue(pd.is_ok())
        # Number of components check
        even = q%2 == 0
        self.check_num_components(pd, 2 if even else 1)
        # Edges in components check
        self.check_component_lengths(pd, [pd.nedges/2]*2 if even else [pd.nedges])
        # Edges in faces check
        face_lengths = [2]*q
        face_lengths.extend([q, q])
        self.check_face_lengths(pd, face_lengths)

    def check_simple_chain(self, n):
        # Simple chain tests
        pd = PlanarDiagram.from_simple_chain(n)
        # Sanity check
        self.assertIsNotNone(pd)
        self.assertTrue(pd.is_ok())
        # Number of components check
        self.check_num_components(pd, n)
        # Edges in components check
        self.check_component_lengths(pd, [2, 2]+([4]*(n-2)))
        # Edges in faces check
        self.check_face_lengths(pd, [2,2]+[2]*(n-1)+[4]*(n-2)+[2*(n-2)+2])

    def check_unknot(self, n):
        # Generate an unknot and make some sanity checks
        pd = PlanarDiagram.from_unknot(n)
        self.assertIsNotNone(pd)
        self.assertTrue(pd.is_ok())
        self.check_num_components(pd, 1)

        # Test the number of edges in each component
        self.check_component_lengths(pd, [2*n])
        # Test the number of edges on each face
        self.check_face_lengths(pd, [1,1]+[2]*(n-1)+[2*(n-1)+2])

    def check_unknot_wye_abcsum(self, n):
        # Test unknot wyes with a+b+c = n
        pd = PlanarDiagram.from_unknot_wye(n, 0,0)
        gold_hash = pd.hash

        for A in range(n+1):
            for B in range(A,n+1):
                a,b,c = (A, B-A, n-B)
                pd = PlanarDiagram.from_unknot_wye(a,b,c)
                # Sanity check
                self.assertIsNotNone(pd)
                self.assertTrue(pd.is_ok())
                # Number of components check
                self.check_num_components(pd, 1)
                # Edges in components check
                self.check_component_lengths(pd, [pd.nedges])
                # Edges in faces check
                face_lengths = [1,1,1]
                face_lengths.extend([2]*(n))
                face_lengths.extend([3, 2*(n+3)])
                self.check_face_lengths(pd, face_lengths)
                # Check hash
                self.assertEqual(pd.hash, gold_hash)

    def test_twist(self):
        for n in range(1, 7):
            self.check_twist_knot(n)

    def test_torus(self):
        for n in range(2, 7):
            self.check_torus_2q(n)

    def test_simple_chain(self):
        for n in range(3, 7):
            self.check_simple_chain(n)

    def test_unknot(self):
        for n in range(2, 7):
            self.check_unknot(n)

    def test_unknotwye(self):
        for n in range(1, 7):
            self.check_unknot_wye_abcsum(n)

if __name__=="__main__":
    unittest.main()
