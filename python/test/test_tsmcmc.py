import unittest
from itertools import product, compress, izip
from libplcurve.plcurve import PlCurve, RandomGenerator
import libplcurve.tsmcmc as ts
import copy

class TestTSMCMCPythonCallbackTypes(unittest.TestCase):
    '''
    TestCase which tests the various types of callback objects which
    can be passed to C, and that the callbacks receive actual plCurve
    data [i.e., without segfaults]
    '''
    def setUp(self):
        # Set up parameters for the checks
        self.rp = ts.RunParams.default_unconfined()
        self.r = RandomGenerator()
        self.r.set(500)
        self.nv = 50
        self.T = ts.Triangulation.new_fan(self.nv)

        self.max_steps = 200
        self.max_secs = 10
        self.confinement_radius = 5
        self.close_failure = 1

        # Here's all of the expectation callables and arguments.
        #  the integrand `df` is always in the second position of args.
        self.exp_calls = (
            (ts.equilateral_expectation,
             [self.r, None, self.max_steps, self.max_secs, self.T, self.rp]),
            (ts.confined_equilateral_expectation,
             [self.r, None, self.confinement_radius, self.nv, self.max_steps,
              self.max_secs, self.rp]),
            (ts.fixed_ftc_expectation,
             [self.r, None, self.close_failure,
              self.max_steps, self.max_secs, self.T, self.rp]),
        )

    def run_integrand_checks(self, exp_calls, to_check):
        for integral, args in exp_calls:
            for df, exp_val in to_check:
                real_args = copy.copy(args)
                real_args[1] = df
                self.assertAlmostEqual(integral(*real_args)[0], exp_val)

    def test_noargs(self):
        # If you'd like to add new callbacks and expected results, put it here
        def cb(L):
            return 42
        def cb_2(L):
            return L[0][1][0]
        to_check = (
            (cb, 42), (cb_2, 1)
        )

        # Run the checks.
        self.run_integrand_checks(self.exp_calls, to_check)

    def test_args(self):
        # If you'd like to add new callbacks and expected results, put it here
        def cb(L, retval):
            return retval
        def cb_2(L, retval):
            return L[0][1][0] + retval
        to_check = (
            ((cb, (42,)), 42),
            ((cb_2, (42,)), 1+42)
        )

        # Run the checks.
        self.run_integrand_checks(self.exp_calls, to_check)

    def test_args_kwargs(self):
        # If you'd like to add new callbacks and expected results, put it here
        def cb(L, delta, base=2):
            return base + delta
        def cb_2(L, delta, base=2):
            return L[0][1][0] * base + delta
        to_check = (
            ((cb, (2,), {"base": 40}), 42),
            ((cb, (2,)), 4),
            ((cb, (2,), {"base": 40}), 40+2)
        )

        # Run the checks.
        self.run_integrand_checks(self.exp_calls, to_check)

    def test_bound_method(self):
        # If you'd like to add new callbacks and expected results, put it here
        class Callbacks(object):
            def cb(self, L):
                return self.retval
            def cb_2(self, L):
                return self.retval - L[0][1][0]
        C = Callbacks()
        C.retval = 42
        to_check = (
            (C.cb, C.retval),
            (C.cb_2, C.retval-1)
        )

        # Run the checks.
        self.run_integrand_checks(self.exp_calls, to_check)

if __name__=="__main__":
    unittest.main()
