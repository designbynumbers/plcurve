#!/usr/bin/python

from libplcurve.plcurve import PlCurve, RandomGenerator
import random

def generate_and_classify(rng_or_state, n, _=None):
    if isinstance(rng_or_state, RandomGenerator):
        rng = rng_or_state
    elif isinstance(rng_or_state, str):
        try:
            rng = RandomGenerator.from_state(rng_or_state)
        except Exception:
            raise Exception("Not a valid RNG state string")
    else:
        raise Exception("Must pass a RNG instance or state string")

    # Create a random closed polygon
    poly = PlCurve.random_closed_polygon(rng, n)

    knottype, nposs = poly.classify()

    factorization = []
    if knottype is not None and nposs > 0:
        #for factor in knottype.factors:
        print knottype.factors
        kt_str = "#".join("{}_{}".format(f.cr, f.ind) for f in knottype.factors)
        print "Knot is probably {}".format(kt_str)

    return factorization
    # Note that poly, knottype (and if necessary) rng are free'd

def run_sequential(num_tests=1000, seed=None):
    num_trefoil = 0

    if seed is None:
        seed = random.getrandbits(32)
    print ("Running {} trials, with seed: {}".format(num_tests, seed))

    rng = RandomGenerator()
    rng.set(seed)

    results = []
    for _ in xrange(num_tests):
        results.append(generate_and_classify(rng, 6))

    knots = [knotsum for knotsum in results if knotsum[0][0] != 0]
    print knots

    # We no longer need the RNG, so...
    del rng

    num_trefoil = len(knots)
    print " === RESULTS FOR SEED {} ===".format(seed)
    print "{} knots in {} trials".format(num_trefoil, num_tests)
    print "Fraction is {}".format(1.0*num_trefoil/num_tests)

    return num_trefoil

## The following command execution should guarantee that the very
## 377th knot is a trefoil:
#    $  python pywrap_demo.py --seed 3630318575 377

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run experiments on random n-gons")
    parser.add_argument(
        "num_trials", metavar="N", type=int, help="Number of trials to run")
    parser.add_argument("--seed", type=int, help="Seed to use, if predetermined")

    args = parser.parse_args()
    print run_sequential(args.num_trials, seed=args.seed)
