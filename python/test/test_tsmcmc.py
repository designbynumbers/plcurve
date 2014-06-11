import unittest
from itertools import product, compress, izip
from libplcurve.plcurve import PlCurve, RandomGenerator
import libplcurve.tsmcmc as ts

class TestTSMCMCPythonCallbackTypes(unittest.TestCase):
    def setUp(self):
        self.rp = ts.RunParams.default_unconfined()
        self.r = RandomGenerator()
        self.r.set(500)
        self.T = ts.Triangulation.new_fan(50)

    def test_noargs(self):
        def cb(L):
            return 42
        self.assertEqual(
            ts.equilateral_expectation(self.r, cb, 200, 10,
                                       self.T, self.rp)[0],
            42)

    def test_args(self):
        def cb(L, retval):
            return retval
        self.assertEqual(
            ts.equilateral_expectation(self.r, (cb, (42,)), 200, 10,
                                       self.T, self.rp)[0],
            42)

    def test_args_kwargs(self):
        def cb(L, delta, base=2):
            return base + delta
        self.assertEqual(
            ts.equilateral_expectation(self.r, (cb, (2,), {"base": 40}), 200, 10,
                                       self.T, self.rp)[0],
            42)

    def test_bound_method(self):
        class Callbacks(object):
            def cb(self, L):
                return self.retval
        C = Callbacks()
        C.retval = 42
        self.assertEqual(
            ts.equilateral_expectation(self.r, C.cb, 200, 10,
                                       self.T, self.rp)[0],
            42)

if __name__=="__main__":
    unittest.main()
