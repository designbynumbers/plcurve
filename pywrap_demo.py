from libplcurve import plCurve

v1 = plCurve.plc_random_vect()
print v1
for i in range(3):
    print "{}: {}".format(i, plCurve.darray_get(v1.c, i))

r = plCurve.make_gsl_rng()
random_poly = plCurve.plc_random_closed_polygon(r, 10)

print plCurve.plc_pointset_diameter(random_poly)
print plCurve.plc_arclength(random_poly, None)
