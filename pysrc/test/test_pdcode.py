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

    def test_r1_move(self):
        pd_fnames = ("data/r1_testB_before.pdstor", "data/r1_testB_after.pdstor")

        files = tuple(open(fname) for fname in pd_fnames)
        for f in files:
            f.readline()
            f.readline()
            f.readline()
        f_before, f_after = files
        pd_before = PlanarDiagram.read(f_before)
        self.assertIsNotNone(pd_before)
        pd_before.regenerate()
        f_before.close()
        pd_after = PlanarDiagram.read(f_after)
        self.assertIsNotNone(pd_after)
        pd_after.regenerate()
        f_after.close()

        print pd_before.faces
        print pd_before.crossings
        print pd_before.edges
        print pd_before.components
        print pd_after.faces
        print pd_after.crossings
        print pd_after.edges
        print pd_after.components
        pd_result = pd_before.R1_loop_deletion(2)
        self.assertEqual(pd_result, pd_after)

    def test_r2_move(self):
        case_faces = {
            "A": 6,
            "B": 6,
            "C": 2,
            "D": 10,
            "E": 1,
            "F": 4,
            "G": 9,
            "H": 8,
            "I": 8,
            "J": 12,
            "K": 5,
            "L": 5,
            "M": 8,
            "N": 5,
            "O": 10,
            "P": 10,
            "R": 3,
            "S": 7,
        }
        failed_cases = "JMP"
        working_cases = "ABCDEFGHIKLNO"
        for case in working_cases + failed_cases:
            print
            print "Case %s:"%case
            pd_fnames = ("data/r2_test%s_before.pdstor"%case,
                         "data/r2_test%s_after.pdstor"%case)

            try:
                files = tuple(open(fname) for fname in pd_fnames)
                f_before, f_after = files

                pd_before = PlanarDiagram.read(f_before, read_header=True)
                self.assertIsNotNone(pd_before)
                #pd_before.regenerate()
                f_before.close()

                pds_after = tuple(PlanarDiagram.read_all(f_after, read_header=True))
                self.assertIsNotNone(pds_after)
                for pdcode in pds_after:
                    #pdcode.regenerate()
                    pass
                f_after.close()
            except Exception as e:
                print "Error: %s"%e
                continue

            face = pd_before.faces[case_faces[case]]
            edge = face[0]
            print face
            print edge.head, edge.tail
            pd_results = pd_before.R2_bigon_elimination(edge.head,edge.tail)
            self.assertEqual(len(pd_results), len(pds_after))
            for result, check in zip(pd_results, pds_after):
                print result
                #result.regenerate()
                self.assertEqual(result, check)


if __name__=="__main__":
    unittest.main()
