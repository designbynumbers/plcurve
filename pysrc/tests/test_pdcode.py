import unittest
from itertools import product, compress, izip
from libpl.pdcode import PlanarDiagram, pd_debug_off
import os, os.path
from libpl.pdcode import Crossing, Edge
from suite import ROOT_DIR

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
        pd.regenerate_hash()
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
                #pd.regenerate_hash()
                #print repr(pd)
                #self.assertEqual(pd.hash, gold_hash)

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

    def test_getcrossing(self):
        pd = PlanarDiagram.torus_knot(2,7)
        self.assertEqual(len(pd.crossings), 7)
        self.assertIsInstance(pd.crossings[3],Crossing)
        pd1 = PlanarDiagram.db_knot(8,3)
        self.assertEqual(len(pd1.crossings), 8)
        self.assertIsInstance(pd1.crossings[3],Crossing)

class TestCrossing(unittest.TestCase):
    def test_length(self):
        pd = PlanarDiagram.simple_chain(3)
        x = pd.crossings[1]
        self.assertEqual(len(x),4)

    def test_togglesign(self):
        pd = PlanarDiagram.simple_chain(3)
        x = pd.crossings[1]
        self.assertEqual(x.sign,1)
        x.toggle_sign()
        self.assertEqual(x.sign,0)
        self.assertEqual(x.sign,pd.crossings[1].sign)
        simp_result = pd.simplify()
        self.assertEqual(len(simp_result),2)
        x.toggle_sign()
        self.assertEqual(x.sign,1)
        x.sign = 2
        self.assertEqual(x.sign,2)
        x.toggle_sign()
        self.assertEqual(x.sign,2)

    def test_iterate(self):
        pd = PlanarDiagram.simple_chain(3)
        n = 0
        for e in pd.crossings[1]:
            self.assertIsInstance(e,Edge)
            n = n + 1
        self.assertEqual(n,4)

class PlanarDiagramAssertMixin(object):
    def __init__(self, *args, **kwargs):
        super(PlanarDiagramAssertMixin, self).__init__(*args, **kwargs)
        self.addTypeEqualityFunc(PlanarDiagram, 'assertDiagramEqual')

    def assertDiagramEqual(self, pd1, pd2, msg=None):
        for cross1, cross2 in zip(pd1.crossings, pd2.crossings):
            if not cross1.index == cross2.index:
                self.fail("Crossing mismatch: Crossing %s is not %s"%(cross2.index, cross1.index))

        for edge1, edge2 in zip(pd1.edges, pd2.edges):
            if not edge1.index == edge2.index:
                self.fail("Edge mismatch: Edge %s is not %s"%(edge2.index, edge1.index))

        for component1, component2 in zip(pd1.components, pd2.components):
            if not component1.index == component2.index:
                self.fail("Component mismatch: Component %s is not %s"%(component2.index, component1.index))

        for face1, face2 in zip(pd1.faces, pd2.faces):
            if not face1.index == face2.index:
                self.fail("Face mismatch: Face %s is not %s"%(face2.index, face1.index))

        if not pd1 == pd2:
            self.fail(self._formatMessage(msg, "%s is not equal to %s." % (repr(pd1), repr(pd2))))

class BeforeAfterFileMixin(object):
    def _get_before_fname(self, tag):
        return os.path.join(self.DATA_DIR,
                            self.DATA_BEFORE_TEMPLATE % tag)

    def _get_after_fname(self, tag):
        return os.path.join(self.DATA_DIR,
                            self.DATA_AFTER_TEMPLATE % tag)

class TestPDCodeMethods(PlanarDiagramAssertMixin, unittest.TestCase):
    def setUp(self):
        self.link = PlanarDiagram.torus_knot(2,12)
        self.knot = PlanarDiagram.torus_knot(2,11)

    def test_ccode(self):
        self.knot.ccode()

class TestR1LoopDeletion(BeforeAfterFileMixin, PlanarDiagramAssertMixin, unittest.TestCase):
    DATA_DIR = os.path.join(ROOT_DIR, "data")
    DATA_BEFORE_TEMPLATE = "r1_test%s_before.pdstor"
    DATA_AFTER_TEMPLATE = "r1_test%s_after.pdstor"

    def checkLoopDeletion(self, tag, face_n):
        with open(self._get_before_fname(tag)) as before_f:
            before_pd = PlanarDiagram.read(before_f, read_header=True)
        self.assertIsNotNone(before_pd)

        with open(self._get_after_fname(tag)) as after_f:
            after_pd = PlanarDiagram.read(after_f, read_header=True)
            after_pd.regenerate_hash()
        self.assertIsNotNone(after_pd)

        before_pd_copy = before_pd.copy()
        result_pd = before_pd.R1_loop_deletion(face_n)

        print result_pd._hash, after_pd._hash
        self.assertEqual(before_pd, before_pd_copy)
        self.assertDiagramEqual(result_pd, after_pd)

    def test_B(self):
        self.checkLoopDeletion("B", 2)

class TestSimplify(BeforeAfterFileMixin, PlanarDiagramAssertMixin, unittest.TestCase):
    DATA_DIR = os.path.join(ROOT_DIR, "data", "simplify")
    DATA_BEFORE_TEMPLATE = "%s_before.pdstor"
    DATA_AFTER_TEMPLATE = "%s_after.pdstor"

    def checkSimplify(self, tag):
        with open(self._get_before_fname(tag)) as before_f:
            before_pd = PlanarDiagram.read(before_f, read_header=True)
        self.assertIsNotNone(before_pd)

        with open(self._get_after_fname(tag)) as after_f:
            after_pds = tuple(PlanarDiagram.read_all(after_f, read_header=True))
        self.assertIsNotNone(after_pds)
        for after_pd in after_pds:
            after_pd.regenerate_hash()
            self.assertIsNotNone(after_pd)

        before_pd_copy = before_pd.copy()
        result_pds = before_pd.simplify()

        self.assertDiagramEqual(before_pd, before_pd_copy)
        self.assertSequenceEqual(result_pds, after_pds)

    @unittest.skip("Getting segfault 8/24.. used to work")
    def test_A_trivial(self):
        self.checkSimplify("A")

    @unittest.skip("Getting segfault 8/24.. used to work")
    def test_3_trivial(self):
        self.checkSimplify("3")

    @unittest.skip("Getting segfault 8/24.. used to work")
    def test_4_trivial(self):
        self.checkSimplify("4")

    @unittest.skip("Crashes (assertions.. meh)")
    def test_5_trivial(self):
        self.checkSimplify("5")

    @unittest.skip("Doesn't simplify; needs R3")
    def test_6_split_link(self):
        self.checkSimplify("6")

    @unittest.skip("Bad _after; Doesn't simplify, needs R3")
    def test_7_split_link(self):
        self.checkSimplify("7")

    @unittest.skip("No _after, needs R3 to simplify")
    def test_deftref1(self):
        self.checkSimplify("deftref1")

    def test_deftref2(self):
        self.checkSimplify("deftref2")

    @unittest.skip("No _after, PDcode not valid.")
    def test_deftref3(self):
        self.checkSimplify("deftref3")

class TestR2BigonElimination(BeforeAfterFileMixin, PlanarDiagramAssertMixin, unittest.TestCase):
    DATA_DIR = os.path.join(ROOT_DIR, "data")
    DATA_BEFORE_TEMPLATE = "r2_test%s_before.pdstor"
    DATA_AFTER_TEMPLATE = "r2_test%s_after.pdstor"

    def checkBigonElimination(self, tag, face_n):
        with open(self._get_before_fname(tag)) as before_f:
            before_pd = PlanarDiagram.read(before_f, read_header=True)
        self.assertIsNotNone(before_pd)

        with open(self._get_after_fname(tag)) as after_f:
            after_pds = tuple(PlanarDiagram.read_all(after_f, read_header=True))
        self.assertIsNotNone(after_pds)
        for after_pd in after_pds:
            self.assertIsNotNone(after_pd)

        before_pd_copy = before_pd.copy()
        face = before_pd.faces[face_n]
        edge = face[0]
        result_pds = before_pd.R2_bigon_elimination_vertices(edge.head, edge.tail)

        self.assertEqual(before_pd, before_pd_copy)
        self.assertEqual(len(result_pds), len(after_pds))
        for result_pd, after_pd in zip(result_pds, after_pds):
            self.assertEqual(result_pd, after_pd)

    @unittest.skip("Bad data file (crossing sign error)")
    def test_A(self):
        self.checkBigonElimination("A", 6)

    def test_B(self):
        self.checkBigonElimination("B", 6)

    def test_C(self):
        self.checkBigonElimination("C", 2)

    def test_D(self):
        self.checkBigonElimination("D", 10)

    @unittest.skip("Getting segfault 8/24.. used to work")
    def test_E(self):
        self.checkBigonElimination("E", 1)

    def test_F(self):
        self.checkBigonElimination("F", 4)

    @unittest.skip("Error with data")
    def test_G(self):
        self.checkBigonElimination("G", 9)

    def test_H(self):
        self.checkBigonElimination("H", 8)

    @unittest.skip("Error with data")
    def test_I(self):
        self.checkBigonElimination("I", 8)

    @unittest.skip("Error with data")
    def test_J(self):
        self.checkBigonElimination("J", 12)

    def test_K(self):
        self.checkBigonElimination("K", 5)

    def test_L(self):
        self.checkBigonElimination("L", 5)

    @unittest.skip("Bad data file (bad results)")
    def test_M(self):
        self.checkBigonElimination("M", 8)

    def test_N(self):
        self.checkBigonElimination("N", 5)

    def test_O(self):
        self.checkBigonElimination("O", 10)

    @unittest.skip("No such face")
    def test_P(self):
        self.checkBigonElimination("P", 10)

    @unittest.skip("No file")
    def test_R(self):
        self.checkBigonElimination("R", 3)

    @unittest.skip("No file")
    def test_S(self):
        self.checkBigonElimination("S", 7)

class TestR3StrandSwap(BeforeAfterFileMixin, PlanarDiagramAssertMixin, unittest.TestCase):
    DATA_DIR = os.path.join(ROOT_DIR, "data", "r3")
    DATA_BEFORE_TEMPLATE = "%s_before.pdstor"
    DATA_AFTER_TEMPLATE = "%s_after.pdstor"

    def setUp(self):
        pass#pd_debug_on()

    def checkStrandSwap(self, tag, face_n, edge_n):
        with open(self._get_before_fname(tag)) as before_f:
            before_pd = PlanarDiagram.read(before_f, read_header=False)
        self.assertIsNotNone(before_pd)

        # with open(self._get_after_fname(tag)) as after_f:
        #     after_pds = tuple(PlanarDiagram.read_all(after_f, read_header=True))
        # self.assertIsNotNone(after_pds)
        # for after_pd in after_pds:
        #     self.assertIsNotNone(after_pd)

        before_pd_copy = before_pd.copy()
        result_pd = before_pd._R3_strand_swap(before_pd.faces[face_n],
                                              before_pd.edges[edge_n])
        #print result_pd[0]

        # self.assertEqual(before_pd, before_pd_copy)
        # self.assertEqual(len(result_pds), len(after_pds))
        # for result_pd, after_pd in zip(result_pds, after_pds):
        #     self.assertEqual(result_pd, after_pd)

    @unittest.skip("We apparently can't handle tangles with the same edge entering twice")
    def test_A_simple(self):
        self.checkStrandSwap("A", 1, 4)

    @unittest.skip("Lacks an appropriately exterior anchor point")
    def test_B_noloop(self):
        self.checkStrandSwap("B", 3, 9)

    def test_c_noloop_anchored(self):
        self.checkStrandSwap("C", 3, 11)


if __name__=="__main__":
    unittest.main()
