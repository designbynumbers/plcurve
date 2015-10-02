#!/usr/bin/python
from functools import reduce
from operator import mul
from collections import defaultdict
from os import path
import csv
import itertools

from libpl.pdcode import *
from libpl.pdstor import *
from libpl.pdcode.pdstor import *
from libpl import data

MAX_CROSS = 10

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
        return self._homfly

    def __str__(self):
        return "%s%s"%(self.name, "m" if self.mirrored else "")
    def __repr__(self):
        return self.__str__()

class LinkFactorization(object):
    def __init__(self, factors, n_splits=0):
        self._factors = factors
        self._homfly = reduce(mul, (K.homfly for K in self.factors), HOMFLYPolynomial('1'))
        self.n_splits = n_splits

    @property
    def factors(self):
        return self._factors
    @property
    def n_cross(self):
        return sum(f.n_cross for f in self.factors)
    @property
    def homfly(self):
        return self._homfly

    def __str__(self):
        return "#".join(str(fc) for fc in self.factors)
    def __repr__(self):
        return self.__str__()

class PDStorClassifier(object):
    def __init__(self):
        self.factors_by_name = dict()
        self.factors_by_homfly = defaultdict(list)

        self.kt_by_homfly = defaultdict(list)
        self.ambiguous_homflys = dict()

        self.load_prime_factors()
        self.load_factorizations()

    def load_prime_factors(self):
        # Load in prime factors (don't worry about mirrors yet)
        cdb = ClassifyDatabase()
        cdb.load(load_links=False)
        for homfly, results in cdb.cls_dict.iteritems():
            for kt in results:
                if kt.ncross > MAX_CROSS:
                    continue
                factor = LinkFactor(kt.ktname, kt.mirror, kt.pdcode)
                self.factors_by_homfly[homfly].append(factor)
                self.factors_by_name[kt.concise(mirror_char="m")] = factor

        # Add in mirrors based on data
        with open(path.join(data.dir, "knotchart-prime.txt"), 'rb') as prime_tsv_f:
            prime_reader = csv.reader(prime_tsv_f, delimiter="\t")
            for shortkt, symtype in prime_reader:
                if ("m" not in shortkt and
                    "r" not in shortkt and
                    symtype in "reversible, none" and
                    shortkt in self.factors_by_name):
                    factor = self.factors_by_name[shortkt]
                    name, pd = factor.name, factor.pd
                    mirrored_pd = pd.copy()
                    for crossing in mirrored_pd.crossings:
                        crossing.toggle_sign()

                    mirror_homfly = mirrored_pd.homfly()
                    assert(mirror_homfly == ~pd.homfly())
                    mirrored_factor = LinkFactor(name, True, mirrored_pd)

                    self.factors_by_homfly[mirror_homfly].append(mirrored_factor)
                    self.factors_by_name["%sm"%shortkt] = mirrored_factor

    def load_factorizations(self):
        # Push prime kts
        for factor in self.factors_by_name.itervalues():
            self.kt_by_homfly[factor.homfly].append(LinkFactorization([factor]))

        # Add in link factorizations
        with open(path.join(data.dir, "knotchart-composite-only.txt"), 'rb') as cmp_tsv_f:
            cmp_reader = csv.reader(cmp_tsv_f, delimiter="\t")
            for ktype_str, symtype in cmp_reader:
                fc_names = ktype_str.split(" ")
                if not all([name in self.factors_by_name for name in fc_names]):
                    continue
                new_fzn = LinkFactorization(
                    [self.factors_by_name[name] for name in fc_names])
                if new_fzn.n_cross > MAX_CROSS:
                    continue

                self.kt_by_homfly[new_fzn.homfly].append(new_fzn)

        print
        print "Ambiguous knot types (by HOMFLY):"
        for homfly, kts in self.kt_by_homfly.iteritems():
            n = len(kts)
            if n > 1:
                N = max([kt.n_cross for kt in kts])
                print "%s: %d, %d %s"%(kts, N, n,homfly)
                self.ambiguous_homflys[homfly] = N

        print self.ambiguous_homflys
        print "~~~"

    def classify(self, pdstor_f, read_header=True):
        ambiguous = defaultdict(list)
        kt_counts = defaultdict(lambda: 0)
        P_unk = HOMFLYPolynomial('1')

        print
        n_dia_total = 0
        i = 0

        fname_prefix = os.path.basename(pdstor_f.name)
        with open("%s.ambig.pdstor"%fname_prefix, "w") as ambig_out:
            PDStorage.start_incremental(ambig_out)
            n_amb_hash, n_amb_elt = 0, 0

            for shadow in PlanarDiagram.read_all(pdstor_f, read_header=read_header):
                i += 1
                if (i-1)%100 == 0:
                    print "Now classifying shadow %d"%(i-1)
                if shadow.ncomps > 1:
                    continue

                iso_classes = defaultdict(set)
                for comp_pd, cmask in PDStoreExpander.component_combinations(
                        shadow, amortize=True, thin=False):
                    for pd, xmask in PDStoreExpander.crossing_combinations(comp_pd):
                        sign_count = sum(x*2-1 for x in xmask)
                        if not (sign_count in iso_classes and pd in iso_classes[sign_count]):
                            # Do classify magic
                            homfly = pd.homfly()
                            if (homfly in self.ambiguous_homflys and
                                pd.ncross >= self.ambiguous_homflys[homfly]):
                                ambiguous[homfly].append(pd)
                                ambig_stor = PDStorage()
                                ambig_stor.add(pd)
                                n_amb_hash, n_amb_elt = ambig_stor.add_to_incremental(
                                    ambig_out, n_amb_hash, n_amb_elt)

                            elif homfly == P_unk:
                                kt_counts["Unknot[]"] += 1
                            else:
                                kts = self.kt_by_homfly[homfly]
                                #print kts
                                if len(kts) == 1:
                                    kt, = kts
                                else:
                                    kts_nx = [k.n_cross for k in kts]
                                    kt = kts[kts_nx.index(min(kts_nx))]

                                #print kt
                                kt_counts[str(kt)] += 1

                            iso_classes[sign_count].add(pd)

                new_n = sum(len(sign_iso) for sign_iso in iso_classes.itervalues())
                n_dia_total += new_n

            with open("%s.counts.tsv"%fname_prefix, "wb") as counts_out:
                countwriter = csv.writer(counts_out, delimiter="\t")
                for row in kt_counts.iteritems():
                    countwriter.writerow(row)

            PDStorage.finish_incremental(ambig_out, n_amb_hash, n_amb_elt)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Classify knots in a pdstor file.')
    parser.add_argument(
        'pdstor', type=argparse.FileType('r'),
        help='.pdstor file to classify')

    args = parser.parse_args()

    classifier = PDStorClassifier()
    classifier.classify(args.pdstor)
