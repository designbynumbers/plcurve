import unittest
from itertools import product, compress, izip
from libpl.pdcode import PlanarDiagram, HOMFLYPolynomial, HOMFLYTerm

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

    def test_richcmp(self):
        A = HOMFLYPolynomial("-a^{-4} - 2a^{-2}z^{1} + a^{-2}z^{2}")
        B = HOMFLYPolynomial("-a^{-4} - 2a^{-2}z^{1} + a^{-2}z^{2}")
        C = HOMFLYPolynomial("-a^{-4} - 2a^{-2}z^{1}")

        self.assertTrue(A == A)
        self.assertTrue(A == B)
        self.assertTrue(C == C)
        self.assertFalse(A == C)
        self.assertFalse(C == A)

if __name__=="__main__":
    unittest.main()
