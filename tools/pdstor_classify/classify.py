#!/usr/bin/python
from functools import reduce
from operator import mul
from collections import defaultdict

from libpl.pdcode import *
from libpl.pdstor import *

class LinkFactor(object):
    def __init__(self,
                 name, mirrored, pd):

        self.name = name
        self.mirrored = mirrored

        self._pd = pd
        self._homfly = self.pd.homfly()
        self._n_cross = pd.ncross
        self._n_comps = pd.ncomps

    @property
    def pd(self):
        return self._pd
    @property
    def n_cross(self):
        return self._n_cross
    @property
    def n_comps(self):
        return self._n_comps
    @property
    def homfly(self):
        return self.homfly

    def __str__(self):
        return "%s%s"%(self.name, "m" if self.mirrored else "")
    def __repr__(self):
        return self.__str__()

class LinkFactorization(object):
    def __init__(self, factors, n_splits=0):
        self._factors = factors
        self._homfly = reduce(mul, (K.homfly for K in self.factors), 1)
        self.n_splits = n_splits

    @property
    def factors(self):
        return self._factors
    @property
    def homfly(self):
        return self._homfly

    def __str__(self):
        return "#".join(self.factors)

class PDStorClassifier(object):
    def __init__(self):
        self.factors_by_name = dict()
        self.factors_by_homfly = defaultdict(list)

        self.load_prime_factors()

    def load_prime_factors(self):
        cdb = ClassifyDatabase()
        cdb.load(load_links=False)
        for homfly, results in cdb.cls_dict.iteritems():
            for kt in results:
                factor = LinkFactor(kt.ktname, kt.mirror, kt.pdcode)
                self.factors_by_homfly[homfly].append(factor)
                self.factors_by_name[kt.concise()] = factor

        print self.factors_by_name

    def load_factorizations(self):
        pass

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Classify knots in a pdstor file.')
    parser.add_argument(
        'pdstor', type=argparse.FileType('r'),
        help='.pdstor file to classify')

    args = parser.parse_args()

    classifier = PDStorClassifier()
    classifier.classify(args.pdstor)
