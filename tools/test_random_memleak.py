from libpl.plcurve import *

r = RandomGenerator()
for i in range(100000):
    K = PlCurve.random_equilateral_closed_polygon(50, r)
