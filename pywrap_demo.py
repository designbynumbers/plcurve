#!/usr/bin/python

from libplcurve import plCurve
import random

# Set up some functions we use a lot
random_closed_polygon = plCurve.plc_random_closed_polygon
classify = plCurve.plc_classify
iget = plCurve.iarray_get
free_knottype = plCurve.free_knottype_struct
free_plcurve = plCurve.plc_free

def generate_and_classify(rng, n, _=None):
    poly = random_closed_polygon(rng, n)

    knottype, nposs = classify(poly)

    cross, idx = None, None
    if knottype is not None and nposs != 0:
        cross = iget(knottype.cr, 0)
        idx = iget(knottype.ind, 0)
        if cross == 3:
            print(cross, idx)

    free_knottype(knottype)
    free_plcurve(poly)

    return cross,idx

def run_sequential(num_tests=1000, seed=None):
    num_trefoil = 0

    if seed is None:
        seed = random.getrandbits(32)
    print ("Running {} trials, with seed: {}".format(num_tests, seed))

    rng = plCurve.make_gsl_rng()
    plCurve.gsl_rng_set(rng, seed)

    for _ in xrange(num_tests):
        generate_and_classify(rng, 6)

    plCurve.gsl_rng_free(rng)
    print "{} trefoils in {} trials".format(num_trefoil, num_tests)
    print "Fraction is {}".format(num_trefoil/num_tests)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run experiments on random n-gons")
    parser.add_argument(
        "num_trials", metavar="N", type=int, help="Number of trials to run")
    parser.add_argument("--seed", type=int, help="Seed to use, if predetermined")

    args = parser.parse_args()
    run_sequential(args.num_trials, seed=args.seed)


# v1 = plCurve.plc_random_vect()
# print v1
# for i in range(3):
#     print "{}: {}".format(i, plCurve.darray_get(v1.c, i))

# for test in range(num_tests):
#     random_poly = plCurve.plc_random_closed_polygon(r, 6)

#     #print plCurve.plc_pointset_diameter(random_poly)
#     #print plCurve.plc_arclength(random_poly, None)

#     knottype, nposs = plCurve.plc_classify(random_poly)
#     #print knottype, nposs
#     if knottype is not None and nposs != 0:
#         cross = plCurve.iarray_get(knottype.cr, 0)
#         idx = plCurve.iarray_get(knottype.ind, 0)
#         if cross == 3:
#             num_trefoil += 1
#             print cross, idx

#     plCurve.plc_free(random_poly)
#     #if (knottype is not None and knottype.nf > 1) or nposs > 1:
#     #    print knottype.nf, nposs

# print num_trefoil
# print num_trefoil*1.0 / num_tests
