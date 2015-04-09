from database import *
from libpl.pdstor import *
from collections import defaultdict
P_unk = HOMFLYPolynomial('1')
from collections import Counter
from operator import mul
from itertools import combinations_with_replacement, groupby

def knot_partition(n, k):
    """knot_partition(n, k) -> [a_1, ...]

    Partition the number n into a sum, with each summand >= k."""
    if n == 0:
        yield []
        return
    if n == k:
        yield [n]
        return

    for i in range(k, n+1):
        for subpart in knot_partition(n-i, i):
            yield [i] + subpart

def knot_rle_partition(n):
    for p in knot_partition(n, 3):
        counts = Counter(p)
        yield sorted(counts.items(), key=lambda x:x[0])

def all_knotgraded_sums_for(rp, L):
    if not rp:
        yield tuple()
        return

    x, count = rp[0]
    for comb in combinations_with_replacement(L[x], count):
        for tail_sum in all_knotgraded_sums_for(rp[1:], L):
            yield comb + tail_sum

def all_knotgraded_sums(n, L):
    """Where L is a set graded from 3...n, take all "knot" partitions of n
    and pick summands from L_i for each term in the partition
    """
    for rp in knot_rle_partition(n):
        for kgsum in all_knotgraded_sums_for(rp, L):
            yield kgsum

if __name__ == "__main__":
    session = Session()
    n_cross = 10
    cdb = ClassifyDatabase()
    cdb.load()

    primes_in_db = defaultdict(list)
    for prime in session.query(LinkFactor.name, LinkFactor.comp_mask, LinkFactor.mirrored).all():
        primes_in_db[prime.name].append( (prime.comp_mask, prime.mirrored,) )
    print primes_in_db

    for homfly, factors in cdb.cls_dict.iteritems():
        for factor in factors:
            # We only care about knots for now
            if factor.ncomps != 1:
                continue

            # Add in the factor to the DB
            if not (factor.ktname in primes_in_db and
                     (bin_list_to_int(factor.compmask), factor.mirror) in primes_in_db[factor.ktname]):

                db_factor = LinkFactor(
                    homfly = homfly,
                    name = factor.ktname,
                    n_cross = factor.ncross,
                    n_comps = factor.ncomps,
                    pd = factor.pdcode,
                    mirrored = factor.mirror,
                    comp_mask = factor.compmask
                )
                session.add(db_factor)

            # Add in the mirrored factor to the DB, if it differs
            if (not (factor.ktname in primes_in_db and
                     (bin_list_to_int(factor.compmask), not factor.mirror) in
                     primes_in_db[factor.ktname]) and
                homfly != ~homfly):

                db_factor = LinkFactor(
                    homfly = ~homfly,
                    name = factor.ktname,
                    n_cross = factor.ncross,
                    n_comps = factor.ncomps,
                    pd = factor.pdcode,
                    mirrored = not factor.mirror,
                    comp_mask = factor.compmask
                )
                session.add(db_factor)
    session.commit()

    MAX_N = n_cross
    MIN_N = 3 # knots
    prime_knots = defaultdict(list)
    for lf in session.query(LinkFactor).filter(LinkFactor.n_comps==1):
        prime_knots[lf.n_cross].append(lf)
    for N in range(MIN_N, MAX_N+1):
        for fzn in all_knotgraded_sums(N, prime_knots):
            factors = []
            for k,g in groupby(fzn):
                fassoc = FactorizationFactor(multiplicity=len(list(g)))
                fassoc.factor = k
                factors.append(fassoc)

            db_fzn = LinkFactorization(
                homfly = reduce(mul, (k.homfly for k in fzn)),
                n_splits = 0,
                n_cross = N,
                n_comps = 1,
                factors = factors,
            )
            session.add(db_fzn)

    session.commit()
    session.close()
