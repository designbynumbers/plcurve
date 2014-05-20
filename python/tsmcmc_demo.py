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

for i in range(100):
    r.set(345)
    print sys.getrefcount(f), sys.getrefcount(dummy), sys.getrefcount(fakedict)
    print ts.equilateral_expectation(r, (f, (dummy,), fakedict), 200, 10, T, rp, None, None)
    print ts.confined_equilateral_expectation(r, f_noargs, 4.0, 40, 200, 10, rp, None, None)
print sys.getrefcount(dummy)
