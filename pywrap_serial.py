import libplcurve.plcurve as pl
r = pl.RandomGenerator()
p = pl.PlCurve.random_closed_polygon(r, 44)
serial = p.serialize()
print(":".join("{0:02x}".format(ord(c)) for c in serial))