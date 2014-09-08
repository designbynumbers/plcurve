import unittest
from itertools import product, compress, izip
from libpl.pdcode import PlanarDiagram, HOMFLYPolynomial, HOMFLYTerm

class TestHOMFLYTerm(unittest.TestCase):
    def setUp(self):
        pass

    def test_init(self):
        A = HOMFLYTerm(3, 14, 159)

    def test_invert(self):
        A = HOMFLYTerm(3, 14, 159)
        B = HOMFLYTerm(3, -14, 159)
        self.assertTrue((~A).equals(B))
        self.assertTrue(A.equals(~B))

    def test_mul(self):
        A = HOMFLYTerm(5, 7, 11)
        B = HOMFLYTerm(13, 19, 29)
        C = HOMFLYTerm(-3, 19, 29)
        D = HOMFLYTerm(3, -19, 29)
        E = HOMFLYTerm(3, 19, -29)
        const = HOMFLYTerm(9, 0, 0)

        F = HOMFLYTerm(65, 26, 40)
        G = HOMFLYTerm(-15, 26, 40)
        H = HOMFLYTerm(15, -12, 40)
        I = HOMFLYTerm(15, 26, -18)

        self.assertTrue((A*B).equals(F))
        self.assertTrue((A*C).equals(G))
        self.assertTrue((A*D).equals(H))
        self.assertTrue((A*E).equals(I))
        self.assertTrue((D*E).equals(const))

    def test_add(self):
        A = HOMFLYTerm(13, 19, 29)
        B = HOMFLYTerm(-3, 19, 29)
        C = HOMFLYTerm(3, -19, 29)
        D = HOMFLYTerm(3, 9, 29)

        E = HOMFLYTerm(10, 19, 29)

        self.assertTrue((A+B).equals(E))

        with self.assertRaises(Exception):
            A+C
        with self.assertRaises(Exception):
            A+D

    def test_cmp_equals(self):
        A = HOMFLYTerm(13, 19, 29)
        B = HOMFLYTerm(-3, 19, 29)
        C = HOMFLYTerm(3, -19, 29)
        D = HOMFLYTerm(3, 9, 29)
        E = HOMFLYTerm(13, 19, 29)
        F = HOMFLYTerm(13, 19, 10)
        G = HOMFLYTerm(13, 18, 30)

        self.assertEquals(A, B)
        self.assertTrue(A == B)
        self.assertFalse(A.equals(B))

        self.assertEquals(A, E)
        self.assertTrue(A == E)
        self.assertTrue(A.equals(E))

        self.assertNotEquals(A, C)
        self.assertNotEquals(A, D)

        self.assertTrue(A > C)
        self.assertTrue(A > D)
        self.assertTrue(A > F)
        self.assertTrue(A > G)

        self.assertTrue(C < A)
        self.assertTrue(C < D)

    def test_hash(self):
        A = HOMFLYTerm(13, 19, 10)
        B = HOMFLYTerm(13, -19, 10)
        C = HOMFLYTerm(13, -19, -10)

        table = dict()
        table[A] = 10
        table[B] = 9

        self.assertEquals(table[A], 10)
        self.assertEquals(table[B], 9)
        with self.assertRaises(KeyError):
            table[C]

    def test_str_repr(self):
        A = HOMFLYTerm(2, 3, 4)
        B = HOMFLYTerm(1, 3, 4)
        C = HOMFLYTerm(-1, 3, 4)
        D = HOMFLYTerm(-2, 3, 4)

        self.assertEquals(str(A), "2a^{3}z^{4}")
        self.assertEquals(str(B), "a^{3}z^{4}")
        self.assertEquals(str(C), "-a^{3}z^{4}")
        self.assertEquals(str(D), "-2a^{3}z^{4}")

        self.assertEquals(repr(A), "2a^{3}z^{4}")
        self.assertEquals(repr(B), "a^{3}z^{4}")
        self.assertEquals(repr(C), "-a^{3}z^{4}")
        self.assertEquals(repr(D), "-2a^{3}z^{4}")

        A = HOMFLYTerm(2, 0, 4)
        B = HOMFLYTerm(1, 0, 4)
        C = HOMFLYTerm(2, 3, 0)
        D = HOMFLYTerm(1, 3, 0)
        E = HOMFLYTerm(2, 0, 0)
        F = HOMFLYTerm(1, 0, 0)
        G = HOMFLYTerm(-1, 0, 0)

        self.assertEquals(str(A), "2z^{4}")
        self.assertEquals(str(B), "z^{4}")
        self.assertEquals(str(C), "2a^{3}")
        self.assertEquals(str(D), "a^{3}")
        self.assertEquals(str(E), "2")
        self.assertEquals(str(F), "1")
        self.assertEquals(str(G), "-1")

        self.assertEquals(repr(A), "2z^{4}")
        self.assertEquals(repr(B), "z^{4}")
        self.assertEquals(repr(C), "2a^{3}")
        self.assertEquals(repr(D), "a^{3}")
        self.assertEquals(repr(E), "2")
        self.assertEquals(repr(F), "1")
        self.assertEquals(repr(G), "-1")

class TestHOMFLYPolynomial(unittest.TestCase):
    def setUp(self):
        pass

    def test_init(self):
        P_trivial = HOMFLYPolynomial("-a^{1}z^{-1} - a^{-1}z^{-1}")

        P_noa = HOMFLYPolynomial("-a^{1}z^{-1} - z^{-1}")
        P_nila = HOMFLYPolynomial("-a^{1}z^{-1} - a^{0}z^{-1}")
        self.assertTrue(P_noa == P_nila)

        P_noz = HOMFLYPolynomial("-a^{1}z^{-1} - a^{-1}")
        P_nilz = HOMFLYPolynomial("-a^{1}z^{-1} - a^{-1}z^{0}")
        self.assertTrue(P_noz == P_nilz)

        P_minus = HOMFLYPolynomial("-a^{1}z^{-1} - a^{-1}")
        P_noone = HOMFLYPolynomial("-a^{1}z^{-1} + -a^{-1}")
        P_one = HOMFLYPolynomial("-a^{1}z^{-1} + -1a^{-1}")
        self.assertTrue(P_minus == P_noone)
        self.assertTrue(P_one == P_noone)
        self.assertTrue(P_one == P_minus)

        # TODO: Write more tests

    def test_getitem(self):
        P_trivial = HOMFLYPolynomial("-a^{1}z^{-1} - a^{-1}z^{-1}")
        A = P_trivial[0]
        B = P_trivial[1]
        C = P_trivial[-1]
        self.assertTrue(B.equals(C))
        with self.assertRaises(IndexError):
            P_trivial[2]

    def test_invert(self):
        P       = HOMFLYPolynomial("-a^{-4} - 2a^{-2}z^{1} + a^{-2}z^{2}")
        P_star  = HOMFLYPolynomial("-a^{4} - 2a^{2}z^{1} + a^{2}z^{2}")
        self.assertTrue(P != P_star)
        self.assertTrue(P == ~P_star)
        self.assertTrue(~P == P_star)
        self.assertTrue(~P != ~P_star)

    def test_mul(self):
        P        = HOMFLYPolynomial("-a^{1}z^{-1} - a^{-1}z^{-1}")
        Psquared = HOMFLYPolynomial("a^{2}z^{-2} + 2z^{-2} + a^{-2}z^{-2}")

        self.assertTrue(P*P == Psquared)

        # TODO: check more; check ordering...

    def test_richcmp(self):
        A = HOMFLYPolynomial("-a^{-4} - 2a^{-2}z^{1} + a^{-2}z^{2}")
        B = HOMFLYPolynomial("-a^{-4} - 2a^{-2}z^{1} + a^{-2}z^{2}")
        C = HOMFLYPolynomial("-a^{-4} - 2a^{-2}z^{1}")

        # == check
        self.assertTrue(A == A)
        self.assertTrue(A == B)
        self.assertTrue(C == C)
        self.assertFalse(A == C)
        self.assertFalse(C == A)

        # != check
        self.assertFalse(A != A)
        self.assertFalse(A != B)
        self.assertFalse(C != C)
        self.assertTrue(A != C)
        self.assertTrue(C != A)

        with self.assertRaises(NotImplementedError):
            A < B
        with self.assertRaises(NotImplementedError):
            A <= B
        with self.assertRaises(NotImplementedError):
            A > B
        with self.assertRaises(NotImplementedError):
            A <= B

    def test_hash(self):
        A = HOMFLYPolynomial("-a^{-4} - 2a^{-2}z^{1} + a^{-2}z^{2}")
        B = HOMFLYPolynomial("-a^{1}z^{-1} - a^{-1}z^{-1}")
        C = HOMFLYPolynomial("-a^{-4} - 2a^{-2}z^{1}")

        table = dict()
        table[A] = 5
        table[C] = 10

        self.assertEqual(table[A], 5)
        self.assertEqual(table[C], 10)

        with self.assertRaises(KeyError):
            table[B]

    def test_str_repr(self):
        A = HOMFLYPolynomial("-a^{-4} - 2a^{-2}z^{1} + a^{-2}z^{2}")

        self.assertEqual(str(A),  "-a^{-4} + -2a^{-2}z^{1} + a^{-2}z^{2}")
        self.assertEqual(repr(A), "-a^{-4} + -2a^{-2}z^{1} + a^{-2}z^{2}")

if __name__=="__main__":
    unittest.main()
