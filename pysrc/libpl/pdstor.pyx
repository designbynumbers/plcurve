from libpl.pdcode import PlanarDiagram, HOMFLYPolynomial
from itertools import combinations, chain, compress, product
from cpython cimport bool
import pprint
pp = pprint.PrettyPrinter(indent=2)
import os
import cPickle
import libpl.data
import re
from itertools import combinations_with_replacement as combs
from collections import namedtuple, defaultdict
from operator import itemgetter, mul

class PrimeFactor(object):
    def __init__(self, ncross, ncomps, ktname, filepos, pdcode, compmask,
                 mirror=False):
        self.ncross = ncross
        self.ncomps = ncomps
        self.ktname = ktname
        self.filepos = filepos
        self.pdcode = pdcode
        self.compmask = compmask
        self.mirror = mirror
    def __invert__(self):
        return PrimeFactor(self.ncross, self.ncomps, self.ktname,
                           self.filepos, self.pdcode, self.compmask,
                           not self.mirror)

    def concise(self):
        if self.ktname == "Unknot[]":
            return "0_1"

        i,j = self.ktname.index("["), self.ktname.index("]")
        uid = self.ktname[i+1:j].split(", ")
        if len(uid) == 3 and uid[1] == "Alternating":
                index_name = "a%s"%uid[2]
        else:
            index_name = "%s"%uid[-1]

        return "%s%s%s%s"%(
            "L" if self.ncomps > 1 else "",
            "%s_%s"%(self.ncross, index_name),
            "*" if self.mirror else "",
            "_%s"%"".join(str(i) for i in self.compmask) if self.ncomps > 1 else ""
        )

    def __str__(self):
        return "%s%s%s"%(
            self.ktname,
            "*" if self.mirror else "",
            "_%s"%"".join(str(i) for i in self.compmask) if self.ncomps > 1 else ""
        )
    def __repr__(self):
        return ("PrimeFactor(ncross=%s, ncomps=%s, ktname=%s, filepos=%s, "
                "pdcode=%s, compmask=%s, mirror=%s)"%(
                    repr(self.ncross), repr(self.ncomps), repr(self.ktname),
                    repr(self.filepos), repr(self.pdcode), repr(self.compmask),
                    self.mirror))

UnknotFactor = PrimeFactor(0, 1, "Unknot[]", (0,-1), PlanarDiagram.unknot(0), ())

class ClassifyResult(object):
    def __init__(self, factor_list, n_splits=0):
        if n_splits >= len(factor_list):
            import ipdb; ipdb.set_trace()
            raise Exception("A link can not be split more times than it has factors")
        self.factor_list = factor_list
        self.n_splits = n_splits

    @classmethod
    def all_sums(cls, factor_results, n_splits):
        for knot_sum in product(*factor_results):
            factor_list = [prime.factor_list[0] for prime in knot_sum]
            #print "################" + str(factor_list)
            yield cls(factor_list, n_splits=n_splits)

    def __getitem__(self, i):
        return self.factor_list[i]
    def __len__(self):
        return len(self.factor_list)

    def __add__(self, other):
        """Knot connect sum operation"""
        return ClassifyResult(self.factor_list + other.factor_list,
                              self.n_splits + other.n_splits)

    def __and__(self, other):
        """Split link operation"""
        return ClassifyResult(self.factor_list + other.factor_list,
                              self.n_splits + other.n_splits + 1)

    def concise(self):
        if self.n_splits == 0:
            return "#".join(factor.concise() for
                              factor in self.factor_list)
        elif self.n_splits == len(self.factor_list) - 1:
            return "U".join(factor.concise() for
                              factor in self.factor_list)
        else:
            return (",".join(factor.concise() for
                              factor in self.factor_list) +
                    "(%sU)"%self.n_splits)


    def __str__(self):
        if self.n_splits == 0:
            return " # ".join(str(factor) for
                              factor in self.factor_list)
        elif self.n_splits == len(self.factor_list) - 1:
            return " U ".join(str(factor) for
                              factor in self.factor_list)
        else:
            return (", ".join(str(factor) for
                             factor in self.factor_list) +
                    "split %s times"%self.n_splits)

    def __repr__(self):
        return "ClassifyResult(factor_list=%s, n_splits=%s)"%(
            self.factor_list,
            self.n_splits)

DEFAULT_PATH=os.path.join("data","pdstors")
SOURCE_DIR=libpl.data.dir

#def bin_list_to_int(blist):
#    return sum((i*2)**n for i,n in enumerate(reversed(blist)))

# Some constant, common HOMFLY polynomials
P_unknot = HOMFLYPolynomial("1")
Ptlink = HOMFLYPolynomial("-a^{1}z^{-1} + -a^{-1}z^{-1}")

class ClassifyDatabase(object):
    KNOT=0
    LINK=1
    def __init__(self):
        self.cls_dict = dict()
        self.by_type = dict()
        self.by_ncross = dict()
        self.prod_max_x = 0
        self.prod_dict = dict()
        self.prod_nx = dict()
        self.trim_dict = dict()
        self.trim_nx = dict()

    @classmethod
    def _iter_name_and_pd(self, names, table):
        name = names.readline().strip()
        pd = PlanarDiagram.read_knot_theory(table)
        while pd is not None and name is not None:
            yield name, pd
            del pd
            name = names.readline().strip()
            pd = PlanarDiagram.read_knot_theory(table)

    @staticmethod
    def component_sign_mask(pdc, signature, fix_first=True):
        new_pd = pdc.copy()
        start = 1 if fix_first else 0
        for comp in compress(range(start, len(new_pd.components)), signature):
            new_pd.reorient_component(comp, 0)
        return new_pd

    comb_table = dict()
    @classmethod
    def _bin_strings(cls, n, amortize=True):
        if amortize: # performance > memory
            if n not in cls.comb_table:
                cls.comb_table[n] = tuple(
                    product([0,1], repeat=n))
            return cls.comb_table[n]
        else:
            return product([0,1], repeat=n)

    @classmethod
    def component_combinations(cls, pdc, amortize=True, fix_first=True):
        """component_combinations(PlanarDiagram) -> generator(new
        PlanarDiagrams, masks)

        Generate all possible component direction permutations.
        """
        mask_len = len(pdc.components) - 1 if fix_first else len(pdc.components)
        for comp_set in cls._bin_strings(mask_len, amortize):
            yield cls.component_sign_mask(pdc, comp_set, fix_first=fix_first), comp_set


    def _load_names_and_table(self, names_fname, table_fname, filetype):
        with open(names_fname) as names_f, open(table_fname) as table_f:
            names_head = names_f.readline()
            names_type, names_count = names_head.split(" ")
            table_head = table_f.readline()
            table_type, table_count = table_head.split(" ")
            count = int(names_count)
            if count != int(table_count):
                raise Exception("Two files differ in number of rows")
            for i, (name, pd_undir) in enumerate(self._iter_name_and_pd(names_f, table_f)):
                #print "\nLine: %s"%i,
                for pd,mask in self.component_combinations(pd_undir):
                    homfly = pd.homfly()
                    result = PrimeFactor(pd.ncross, pd.ncomps, name, (filetype, i), pd, mask)
                    if homfly in self.cls_dict:
                        self.cls_dict[homfly].append(result)#(name, (filetype, i), pd))
                    else:
                        self.cls_dict[homfly] = [result]#(name, (filetype, i), pd),]
                    if i >= count-1:
                        break

    def load(self, load_knots=True, load_links=True):
        if load_knots:
            self.load_rolfsen(
                os.path.join(SOURCE_DIR,"rolfsennames.txt"),
                os.path.join(SOURCE_DIR,"rolfsentable.txt"))
        if load_links:
            self.load_thistlethwaite(
                os.path.join(SOURCE_DIR,"thistlethwaitenames.txt"),
                os.path.join(SOURCE_DIR,"thistlethwaitetable.txt"))

    def load_rolfsen(self, names_fname, table_fname):
        self._load_names_and_table(names_fname, table_fname, self.KNOT)
    def load_thistlethwaite(self, names_fname, table_fname):
        self._load_names_and_table(names_fname, table_fname, self.LINK)

    def calculate_composites(self, max_n_cross, num_components=float('inf')):
        """Use facts about HOMFLY polynomials to populate a new search table
        which lists HOMFLYs for composite knots and split links

        max_n_cross: Max number of crossings of the composite knots in the table
        (warning: if this is too large, memory WILL be devoured without compassion)
        """
        self.prod_max_x = max_n_cross

        # prepare trim_by_nx as list of buckets of crossing-limited entries
        self.trim_by_nx = [[] for _ in range(1 + max_n_cross)]
        self.trim_by_nx[0] = [P_unknot]

        # We only want to worry about polynomials which have a knot with a low
        # enough crossing number in their result bucket
        for P, bucket in self.cls_dict.iteritems():
            minx = max_n_cross + 1 # Basically, infinity
            for res in bucket:
                nx = res.ncross

                if nx > max_n_cross:
                    continue

                minx = min(minx, nx)
            if minx <= max_n_cross:
                self.trim_by_nx[minx].append(P)

        self.prod_dict = defaultdict(list)

        # The most times it makes sense to split link (each adds 2 crossings)
        max_split_links = min(
            (max_n_cross//2) + 1,
            num_components-1)

        # Smallest nontrivial crossing number:
        #  2 for links (e.g. hopf)
        #  3 for knots (trefoil)
        if num_components > 1:
            MIN_NONTRIVIAL_NX = 2
        else:
            MIN_NONTRIVIAL_NX = 3

        # The most times it makes sense to split link OR connect sum
        max_joins = (max_n_cross//MIN_NONTRIVIAL_NX) + 1

        unknot_polys = (P_unknot,) * (max_split_links + 1)
        tlink_pows = [Ptlink**(n_unions) for n_unions in range(max_split_links+1)]
        # Trivial split links
        for k in range(1, max_split_links):
            self.prod_dict[tlink_pows[k]].append((k, 2*k, unknot_polys[:k+1]))

        # Explode trim_by_nx into pairs (nx, poly)
        nx_poly_pairs = []
        for nx, bucket in enumerate(self.trim_by_nx[2:], start=2):
            for poly in bucket:
                nx_poly_pairs.append((nx, poly))

        # Nontrivial connect sums and split links
        for k in range(1, max_joins):
            for k_fold_comb in combs(nx_poly_pairs, k):
                # Smallest number of crossings a connect sum yields
                n_cs = sum((n for n,_ in k_fold_comb))
                if n_cs > max_n_cross:
                    continue # Too many crossings, move on

                # List of polynomials
                print n_cs, k_fold_comb
                P     = [p for _,p in k_fold_comb]
                print P
                Pstar = [~p for p in P]

                P_cmp = [] # List of composite knot homflys & factors
                for mask in product((False, True), repeat=k):
                    chain = tuple(B if mirrored else A for
                                  A,B,mirrored in zip(P,Pstar,mask))
                    P_cmp.append( (reduce(mul, chain), chain) )

                # Each split link adds 2 crossings, how many can we have?
                # We can also only be split fewer times than we have factors
                assert(n_cs <= max_n_cross)
                max_n_unions = min(
                    (num_components-1,
                     (max_n_cross - n_cs)//2,
                     k-1))

                # All possible number of split links for these factors
                min_n_trivlks = 1 if k == 1 else 0
                # If we're a 1-fold-combination, need at least 1 trivial link

                for n_unions in range(max_n_unions + 1):
                    nx_with_split = n_unions*2+n_cs
                    assert(nx_with_split <= max_n_cross)

                    # For each mirror-masked homfly...
                    for Q, factors in P_cmp:
                        max_n_trivlks = min(
                            (max_n_cross - nx_with_split) // 2,
                            num_components - n_unions - 1)
                        # For each possible number of split-linked unknots
                        for i in range(min_n_trivlks, max_n_trivlks+1):
                            if n_unions >= len(factors):
                                print "Error:"
                                print n_unions
                                print factors
                                assert False

                            print max_n_unions
                            print max_n_trivlks
                            print n_unions, i
                            total_times_split = n_unions + i
                            total_nx = nx_with_split + (2*i)

                            assert total_times_split < len(factors)+len(unknot_polys[:i])
                            # Factor from split links
                            S = tlink_pows[total_times_split]

                            self.prod_dict[S * Q].append(
                                (total_times_split,
                                 total_nx,
                                 factors + unknot_polys[:i]))

    def classify_prime_homfly(self, hf, ncross=float("inf")):
        mirror_hf = ~hf
        results = []

        if hf == P_unknot:
            results.append(ClassifyResult([UnknotFactor]))

        if hf in self.cls_dict:
            results.extend(ClassifyResult([prime_factor]) for
                           prime_factor in self.cls_dict[hf] if
                           prime_factor.ncross <= ncross)

        if mirror_hf != hf and mirror_hf in self.cls_dict:
            results.extend(ClassifyResult([~prime_factor]) for
                           prime_factor in self.cls_dict[mirror_hf] if
                           prime_factor.ncross <= ncross)

        return results

    def classify(self, pd):
        return self.classify_homfly(pd.homfly(), ncross=pd.ncross)

    def classify_homfly(self, hf, ncross=float("inf")):
        results = []

        # Search for prime matches
        results.extend(self.classify_prime_homfly(hf, ncross))

        # Search for composite and split matches
        if hf in self.prod_dict:
            n_splits, min_nx, factors = self.prod_dict[hf][0]

            if min_nx <= ncross:
                results.extend(
                    ClassifyResult.all_sums(
                        [self.classify_prime_homfly(hf_factor) for
                         hf_factor in factors],
                        n_splits=n_splits))

        return results

class PDStoreExpander(object):
    def __init__(self, dirloc=DEFAULT_PATH,
                 amortize=True, debug=False):
        self.amortize = amortize
        self.dirloc = dirloc

    def open(self, crossings_list, debug=False, orient_all=True,
             thin=False, homflys=True, max_components=None):
        pd_files = (self._db_file(x, self.dirloc) for x in crossings_list)
        iters = []
        for f in pd_files:
            _,npds,_ = self.parse_header(f, debug=debug)
            iters.append(self.read_pdstor(f, debug=debug, num_pds=int(npds),
                                          thin=thin, max_components=max_components,
                                          homflys=homflys, orient_all=orient_all))
        return chain(*iters)

    def open_shadows(self, crossings_list, debug=False, orient_all=True,
             thin=False, homflys=True, max_components=None):
        pd_files = (self._db_file(x, self.dirloc) for x in crossings_list)
        iters = []
        for f in pd_files:
            _,npds,_ = self.parse_header(f, debug=debug)
            iters.append(self.read_pdstor_shadows(f, debug=debug, num_pds=int(npds),
                                                  thin=thin, max_components=max_components,
                                                  homflys=homflys, orient_all=orient_all))
        return chain(*iters)

    @staticmethod
    def _db_file(cross_count, dirloc):
        return file(os.path.join(dirloc,"%d.pdstor"%cross_count))

    @staticmethod
    def uid(pd, pos_in_stor, crs_sgn, cmp_sgn):
        """
        Present uid spec is:
        (number of crossings,
         position in the .pdstor,
         cross sign mapping,
         component sign mapping).

        The crossing number is extraneous but should eventually become
        a filename descriptor instead; crossing and component count info
        is sufficiently encoded in the sign mapping (length of the tuple)
        """
        return (pd.ncross,
                pos_in_stor,
                crs_sgn,
                cmp_sgn)

    @staticmethod
    def _pds_in_file(f):
        pd = PlanarDiagram.read(f)
        while pd is not None:
            yield pd
            del pd
            pd = PlanarDiagram.read(f)

    @staticmethod
    def crossing_sign_mask(pdc, signature, thin=False):
        new_pd = pdc.copy(thin=thin)
        new_pd.set_all_crossing_signs(signature)
        return new_pd

    @staticmethod
    def component_sign_mask(pdc, signature, thin=False):
        new_pd = pdc.copy(thin=thin)
        unsignature = (1-x for x in signature) # Make 0's 1s and  1's 0s
        # We want to flip only components with 0 (i.e. negative) sign
        for comp in compress(range(new_pd.ncomps), unsignature):
            new_pd.reorient_component(comp, 0)
        return new_pd

    comb_table = dict()
    @classmethod
    def _bin_strings(cls, n, amortize=True):
        if amortize: # performance > memory
            if n not in cls.comb_table:
                cls.comb_table[n] = tuple(
                    product([0,1], repeat=n))
            return cls.comb_table[n]
        else:
            return product([0,1], repeat=n)
    @classmethod
    def crossing_combinations(cls, pdc, bool amortize=True, thin=False):
        """crossing_combinations(PlanarDiagram) -> generator(new
        PlanarDiagrams, masks)

        Generate all possible crossing sign permutations.  First crossing
        is fixed as (+) so as to avoid mirror images
        """
        for cross_set in cls._bin_strings(pdc.ncross, amortize):
            yield cls.crossing_sign_mask(pdc, cross_set, thin=thin), cross_set

    @classmethod
    def component_combinations(cls, pdc, amortize=True, thin=False):
        """component_combinations(PlanarDiagram) -> generator(new
        PlanarDiagrams, masks)

        Generate all possible component permutations.
        """
        for comp_set in cls._bin_strings(pdc.ncomps, amortize):
            yield cls.component_sign_mask(pdc, comp_set, thin=thin), comp_set

    def add_data(self, cross_count, dirloc=DEFAULT_PATH):
        f = self._db_file(cross_count, dirloc)
        self.parse_header(f)
        self.read_pdstor(f)

    def add_db(self, other_db):
        for homfly, pds in other_db.pd_dict.iteritems():
            if homfly in self.pd_dict:
                self.pd_dict[homfly] |= pds
            else:
                self.pd_dict[homfly] = pds

    @classmethod
    def pd_by_filepos(cls, cross_count, pos, dirloc=DEFAULT_PATH):
        f = cls._db_file(cross_count, dirloc)
        cls.parse_header(f)
        # there is a better way to do this using file spec but.. eh
        for i, pd in enumerate(cls._pds_in_file(f)):
            if i == pos:
                return pd

    @classmethod
    def pd_by_id(cls, uid, dirloc=DEFAULT_PATH):
        cnum, pos, cross_sgn, comp_sgn = uid
        return (cls.crossing_sign_mask(
            cls.component_sign_mask(
                cls.pd_by_filepos(cnum, pos, dirloc),
                comp_sgn), cross_sgn))

    @classmethod
    def parse_header(cls, f, debug=False):
        # TODO: Actually parse the header
        header = f.readline()
        if "pdstor" not in header:
            return None # something bad actually happened
        # try to parse the data line
        header = f.readline().split()
        try:
            claimed, actual = header[1].split("/")
            nhashes = header[4]
        except:
            return None # actually throw an error here too
        return (claimed,actual,nhashes)

    def shadow_iter(self, pd_id, pd, debug=False, orient_all=True,
                    thin=False, homflys=True):
        if orient_all:
            cmp_gen = self.component_combinations(pd, self.amortize, thin=thin)
        else:
            cmp_gen = ((pd, (1,)),)

        for cmp_pd, cmp_sgn in cmp_gen:
            for crs_pd, crs_sgn in self.crossing_combinations(cmp_pd, self.amortize, thin=thin):
                if homflys:
                    homfly = crs_pd.homfly()
                else:
                    homfly = None
                pduid = self.uid(crs_pd, pd_id, crs_sgn, cmp_sgn)
                yield (crs_pd, homfly, pduid)

    def read_pdstor(self, f, debug=False, num_pds=None, orient_all=True,
                    thin=False, homflys=True, max_components=None):
        # Read in the pd codes one-by-one.
        if debug: print "Reading file %s"%f.name
        for pd_id, pd in enumerate(PlanarDiagram.read_all(f, read_header=False, thin=thin)):
            if debug: print "> Reading %s of %s (%0.1f%%)"%(
                    pd_id+1,
                    "Unknown" if num_pds is None else num_pds,
                    "??" if num_pds is None else 100.0*pd_id/num_pds),
            if max_components is not None and pd.ncomps > max_components:
                if debug: print "... skipped (too many components)"
                continue
            else:
                if debug: print

            for oriented_pd in self.shadow_iter(pd_id, pd, debug=debug,
                                                orient_all=orient_all, thin=thin, homflys=homflys):
                yield oriented_pd
        f.close()

    def read_pdstor_shadows(self, f, debug=False, num_pds=None, orient_all=True,
                            thin=False, homflys=True, max_components=None):
        # Read in the pd codes one-by-one.
        if debug: print "Reading file %s"%f.name
        for pd_id, pd in enumerate(PlanarDiagram.read_all(f, read_header=False, thin=thin)):
            if debug: print "> Reading %s of %s (%0.1f%%)"%(
                    pd_id+1,
                    "Unknown" if num_pds is None else num_pds,
                    "??" if num_pds is None else 100.0*pd_id/num_pds),
            if max_components is not None and pd.ncomps > max_components:
                if debug: print "... skipped (too many components)"
                continue
            else:
                if debug: print
            yield self.shadow_iter(pd_id, pd, debug=debug,
                                    orient_all=orient_all, thin=thin, homflys=homflys)
        f.close()

class PDStorage:
    def __init__(self, isotopy=False):
        """PDStorage(isotopy=False)

        Create a new PDStorage object. If isotopy is True,
        then equivalence among diagrams will be diagram isotopy,
        rather than diagram isomorphism."""
        
        self.hashes = defaultdict(list)
        self.isotopy = isotopy

    def copy(self):
        """Return a shallow copy of self."""
        pdstor = self.__class__(isotopy=self.isotopy)
        for hash_, bucket in self.hashes.iteritems():
            pdstor.hashes[hash_] = bucket.copy()
        return pdstor
            
    @classmethod
    def read(cls, f, read_header=False, isotopy=False):
        """read(isotopy=False) -> new PDStorage

        Read a PDStorage object from a file. If isotopy is True,
        then equivalence among diagrams will be diagram isotopy,
        rather than diagram isomorphism."""
        
        ret = cls(isotopy=isotopy)
        for pdcode in PlanarDiagram.read_all(f, read_header=read_header):
            ret.hashes[pdcode.hash].append(pdcode)
        return ret

    def __le__(self, other):
        if not isinstance(other, PDStorage):
            return False

        if len(self.hashes) > len(other.hashes):
            return False
        
        for hash_, bucket in self.hashes.iteritems():
            if hash_ in other.hashes:
                if len(bucket) > len(other.hashes[hash_]):
                    return False
                
                for pdcode in bucket:
                    if not True in (pdcode.isomorphic(opd) for
                                    opd in other.hashes[hash_]):
                        
                        return False
            else:
                return False
        return True

    def __ge__(self, other):
        return NotImplemented

    ## Really lazy quadruple of set inclusion/equality operators
    def __lt__(self, other):
        return self <= other and not self >= other
    def __gt__(self, other):
        return NotImplemented
    def __eq__(self, other):
        return self <= other and self >= other
    def __ne__(self, other):
        return not self == other

    def __sub__(self, other):
        if not isinstance(other, PDStorage):
            raise TypeError

        pdstor = self.__class__(isotopic=self.isotopic)
        
        for hash_, bucket in self.hashes.iteritems():
            newbucket = bucket.copy()
            if hash_ in other.hashes:
                for pdcode in other.hashes[hash_]:
                    newbucket = [subpd for subpd in newbucket if
                                 (self.isotopic and subpd.isotopic(pdcode)) or
                                 (not self.isotopic and subpd.isomorphic(pdcode))]
            pdstor.hashes[hash_] = newbucket
        return pdstor

class PDDatabase(PDStoreExpander):
    def save(self, fname):
        f = open(fname, "wb")
        cPickle.dump(self, f, cPickle.HIGHEST_PROTOCOL)
        f.close()

    def __init__(self, crossings_list=[3,4,5], dirloc=DEFAULT_PATH,
                 amortize=True, debug=False, callbacks=None, keep_pd=False):
        super(PDDatabase,self).__init__(dirloc=dirloc,
                                        amortize=amortize,
                                        debug=debug)
        self.callbacks = [] if callbacks is None else callbacks
        self.pd_dict = {}
        self.all_pds = {}
        for pd, homfly, uid in self.open(crossings_list):
            if keep_pd:
                self.add_pd(homfly, pd)
            else:
                self.add_pd(homfly, uid)
            for cb in self.callbacks:
                cb(pd, homfly, uid)

    @staticmethod
    def load(fname):
        f = open(fname, "rb")
        ret = cPickle.load(f)
        f.close()
        return ret

    def add_pd(self, homfly, pduid):
        if homfly in self.pd_dict:
            self.pd_dict[homfly].add(pduid)
        else:
            self.pd_dict[homfly] = set([pduid])


    def subdb(self, select):
        """Returns a new database which only contains diagrams held in
        files contained in the list select"""
        newdb = PDDatabase([])
        for homfly, pds in self.pd_dict.iteritems():
            for pduid in pds:
                if pduid[0] in select:
                    newdb.add_pd(homfly, pduid)
        return newdb

    def classify(self, clsdb):
        """Separates the database by key"""
        known = dict()
        unknown = dict()
        for homfly, pdids in self.pd_dict.iteritems():
            if homfly in clsdb.cls_dict:
                known[homfly] = (clsdb.cls_dict[homfly], pdids)
            elif homfly == "1": # Unknot
                known[homfly] = ("Knot[0, 0]", pdids)
            else:
                unknown[homfly] = (None, pdids)

        return known, unknown

if __name__ == "__main__":
    pddb = PDDatabase([3], dirloc="../data/pdstors")
    clsdb = ClassifyDatabase()
    clsdb.load_rolfsen(
        "../test/rolfsennames.txt",
        "../test/rolfsentable.txt")
    clsdb.load_thistlethwaite(
        "../test/thistlethwaitenames.txt",
        "../test/thistlethwaitetable.txt")

    pp.pprint(clsdb.cls_dict)
    pp.pprint(pddb.classify(clsdb))
    #pp.pprint(pddb.pd_dict)
    #pddb.save("database_main.pystor")
