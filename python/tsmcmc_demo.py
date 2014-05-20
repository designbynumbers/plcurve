import libplcurve.tsmcmc as ts
import libplcurve.plcurve as pl
from IPython import start_ipython

def f(plc, args, newarg=50):
    return plc.total_curvature()

def f_noargs(plc):
    return plc.total_curvature()

class MethHolder(object):
    def m(self, plc):
        return 5.0

r = pl.RandomGenerator()
r.set(345)
T = ts.Triangulation.new_fan(50)
rp = ts.RunParams.default_unconfined()

import sys

dummy = ("no one's ever typed this before", "nope nope")
fakedict = {"newarg": 40}

for i in range(1):
    r.set(345)
    print sys.getrefcount(f), sys.getrefcount(dummy), sys.getrefcount(fakedict)
    res,s,e = ts.equilateral_expectation(r, (f, (dummy,), fakedict), 200, 10, T, rp)
    print ts.confined_equilateral_expectation(r, f_noargs, 4.0, 40, 200, 10, rp)

    print s.max_lagged_covariance_used
    print s.lagged_covariances_available

    print s.dihedral_steps
    print s.mp_steps
    print s.permute_steps

    print s.total_seconds
    print s.geyer_ips_seconds
