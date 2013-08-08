from libplcurve import plCurve
import random
import multiprocessing

def random_polygons(n_trials, rng):
    for _ in range(n_trials):
        yield plCurve.plc_random_closed_polygon(rng, 6)

def classify_polygons():
    poly = None
    while True:
        poly = (yield poly)

        knottype, nposs = plCurve.plc_classify(poly)
        #print knottype, nposs
        if knottype is not None and nposs != 0:
            cross = plCurve.iarray_get(knottype.cr, 0)
            idx = plCurve.iarray_get(knottype.ind, 0)
            if cross == 3:
                print(cross, idx)
        plCurve.free_knottype_struct(knottype)

def free_polygons():
    while True:
        poly = (yield)
        plCurve.plc_free(poly)

def run_pipeline():
    r = plCurve.make_gsl_rng()
    plCurve.gsl_rng_set(r, random.getrandbits(32))

    num_tests = 1000
    num_trefoil = 0

    stages = [random_polygons(num_tests, r),
              classify_polygons(),
              free_polygons()]

    pl = pipeline.Pipeline(stages)
    pl.run_parallel(128)

    plCurve.gsl_rng_free(r)

from functools import partial

def knot_process(_, n=6):
    rng = plCurve.make_gsl_rng()
    plCurve.gsl_rng_set(rng, random.getrandbits(32))

    poly = plCurve.plc_random_closed_polygon(rng, n)

    knottype, nposs = plCurve.plc_classify(poly)

    cross, idx = None, None
    if knottype is not None and nposs != 0:
        cross = plCurve.iarray_get(knottype.cr, 0)
        idx = plCurve.iarray_get(knottype.ind, 0)
        if cross == 3:
            print(cross, idx)

    plCurve.free_knottype_struct(knottype)
    plCurve.plc_free(poly)
    plCurve.gsl_rng_free(rng)

    return cross, idx

def run_parallel():
    num_tests = 1000
    num_trefoil = 0

    pool = multiprocessing.Pool(processes=128)
    for cross, idx in pool.imap_unordered(knot_process, range(num_tests)):
        if cross:
            print("{}_{}".format(cross, idx))

def run_sequential(num_tests=1000, seed=None):
    num_trefoil = 0

    if seed is None:
        seed = random.getrandbits(32)

    rng = plCurve.make_gsl_rng()
    plCurve.gsl_rng_set(rng, seed)

    for _ in range(num_tests):
        poly = plCurve.plc_random_closed_polygon(rng, 6)

        knottype, nposs = plCurve.plc_classify(poly)

        cross, idx = None, None
        if knottype is not None and nposs != 0:
            cross = plCurve.iarray_get(knottype.cr, 0)
            idx = plCurve.iarray_get(knottype.ind, 0)
            if cross == 3:
                print(cross, idx)
                num_trefoil += 1

        plCurve.free_knottype_struct(knottype)
        plCurve.plc_free(poly)

    plCurve.gsl_rng_free(rng)
    print "{} trefoils in {} trials".format(num_trefoil, num_tests)
    print "Fraction is {}".format(num_trefoil/num_tests)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run experiments on random n-gons")
    parser.add_argument(
        "num_trials", metavar="N", type=int, help="Number of trials to run")

    args = parser.parse_args()
    run_sequential(args.num_trials)


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
