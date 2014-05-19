import libplcurve.tsmcmc as ts
import libplcurve.plcurve as pl

def f(plc, args):
    return plc.total_curvature()

r = pl.RandomGenerator()
r.set(345)
T = ts.Triangulation.new_fan(50)
rp = ts.RunParams.default_unconfined()

print ts.equilateral_expectation(r, f, None, 200, 10, T, rp, None, None)
print ts.equilateral_expectation(r, f, None, 200, 10, T, rp, None, None)
