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
from collections import namedtuple
from operator import itemgetter

Result = namedtuple('Result',
                    'ncross ncomps ktname filepos pdcode')

DEFAULT_PATH=os.path.join("data","pdstors")
SOURCE_DIR=libpl.data.dir

def bin_list_to_int(blist):
    return sum((i*2)**n for i,n in enumerate(reversed(blist)))

# Some constant, common HOMFLY polynomials
Punk = HOMFLYPolynomial("1")
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
    def component_sign_mask(pdc, signature):
        new_pd = pdc.copy()
        for comp in compress(range(len(new_pd.components)), signature):
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
    def component_combinations(cls, pdc, amortize=True):
        """component_combinations(PlanarDiagram) -> generator(new
        PlanarDiagrams, masks)

        Generate all possible component direction permutations.
        """
        for comp_set in cls._bin_strings(len(pdc.components), amortize):
            yield cls.component_sign_mask(pdc, comp_set), comp_set


    def _load_names_and_table(self, names_fname, table_fname, filetype):
        names_f = None
        table_f = None
        try:
            names_f = open(names_fname)
            table_f = open(table_fname)

            names_head = names_f.readline()
            names_type, names_count = names_head.split(" ")
            table_head = table_f.readline()
            table_type, table_count = table_head.split(" ")
            count = int(names_count)
            if count != int(table_count):
                raise Exception("Two files differ in number of rows")
            for i, (name, pd_undir) in enumerate(self._iter_name_and_pd(names_f, table_f)):
                for pd,_ in self.component_combinations(pd_undir):
                    homfly = pd.homfly()
                    if homfly in self.cls_dict:
                        self.cls_dict[homfly].append((name, (filetype, i), pd))
                    else:
                        self.cls_dict[homfly] = [(name, (filetype, i), pd),]
                    if i >= count-1:
                        break
        finally:
            if names_f:
                names_f.close()
            if table_f:
                table_f.close()

    def load(self):
        self.load_rolfsen(
            os.path.join(SOURCE_DIR,"rolfsennames.txt"),
            os.path.join(SOURCE_DIR,"rolfsentable.txt"))
        self.load_thistlethwaite(
            os.path.join(SOURCE_DIR,"thistlethwaitenames.txt"),
            os.path.join(SOURCE_DIR,"thistlethwaitetable.txt"))

    def load_rolfsen(self, names_fname, table_fname):
        self._load_names_and_table(names_fname, table_fname, self.KNOT)
    def load_thistlethwaite(self, names_fname, table_fname):
        self._load_names_and_table(names_fname, table_fname, self.LINK)

    def calculate_composites(self, max_n_cross):
        """Use facts about HOMFLY polynomials to populate a new search table
        which lists HOMFLYs for composite knots and split links

        max_n_cross: Max number of crossings of the composite knots in the table
        (warning: if this is too large, memory WILL be devoured without compassion)
        """
        self.prod_max_x = max_n_cross

        # prepare trim_dict of crossing-limited entries
        self.trim_dict = dict()
        self.trim_nx = dict()
        for P, bucket in self.cls_dict.iteritems():
            nbuk = dict()
            minx = max_n_cross
            for res in bucket:
                nx = res[2].ncross
                if nx > max_n_cross:
                    continue

                minx = min(minx, nx)
                nres = Result(nx, res[2].ncomps, res[0], res[1], res[2])
                if nx in nbuk:
                    nbuk[nx].append(nres)
                else:
                    nbuk[nx] = [nres]
            if nbuk:
                self.trim_dict[P] = nbuk
                self.trim_nx[P] = minx

        self.prod_dict = dict()
        for Pa, Pb in combs(self.trim_dict.iterkeys(), 2):
            # Smallest number of crossings which this prod can represent:
            n_cs = self.trim_nx[Pa] + self.trim_nx[Pb] # As connected sum
            n_sl = n_cs + 2 # As split link

            # HOMFLY of mirrored diagrams
            Pastar = ~Pa
            Pbstar = ~Pb

            # Knot connected sum
            Pab     = (Pa * Pb)
            Pabstar = (Pa * Pbstar)
            self.prod_dict[Pab]     = (1, Pa, Pb)
            self.prod_dict[Pabstar] = (1, Pa, Pbstar)
            self.prod_nx[Pab]     = n_cs
            self.prod_nx[Pabstar] = n_cs

            # Split links
            Pabsplit     = (Ptlink * Pab)
            Pabstarsplit = (Ptlink * Pabstar)
            self.prod_dict[Pabsplit]     = (2, Pa, Pb)
            self.prod_dict[Pabstarsplit] = (2, Pa, Pbstar)
            self.prod_nx[Pabsplit]     = n_sl
            self.prod_nx[Pabstarsplit] = n_sl

            self.prod_dict[Ptlink] = (2, Punk, Punk) # Trivial split link
            self.prod_nx[Ptlink] = 2
            # WARNING: lmpoly is not returning the correct homfly for this
        for P in self.trim_dict.iterkeys():
            # Split link with unknot
            Psplit = (Ptlink * P)
            self.prod_dict[Psplit] = (2, P, Punk)
            self.prod_nx[Psplit]   = self.trim_nx[Pa] + 2

    def classify_prime_homfly(self, homfly, ncross=1000):
        mhomfly = ~homfly
        if homfly == Punk:
            return (False, {0: [Result(0, 1, "Unknot[]", (0,-1), PlanarDiagram.unknot(0))]})
        if homfly in self.trim_dict and self.trim_nx[homfly] <= ncross:
            return (False, {nx: h for nx, h in self.trim_dict[homfly].iteritems()
                            if nx <= ncross})
        elif mhomfly in self.trim_dict and self.trim_nx[mhomfly] <= ncross:
            return (True, {nx: h for nx, h in self.trim_dict[mhomfly].iteritems()
                           if nx <= ncross})

        return (None, dict())

    def classify(self, pd):
        return self.classify_homfly(pd.homfly(), ncross=pd.ncross)

    def classify_homfly(self, homfly, ncross=1000):
        mhomfly = ~homfly

        p_res = (None, dict())
        c_res, c_mres = [], []

        # Search for prime matches
        p_res = self.classify_prime_homfly(homfly, ncross)
        #p_mres = self.classify_prime_homfly(mhomfly, ncross)

        # Search for composite and split matches
        if homfly in self.prod_dict and self.prod_nx[homfly] <= ncross:
            sumkind, Pa, Pb = self.prod_dict[homfly]
            c_res = (sumkind,
                     self.classify_prime_homfly(Pa, ncross),
                     self.classify_prime_homfly(Pb, ncross))
        if mhomfly in self.prod_dict and self.prod_nx[mhomfly] <= ncross:
            sumkind, Pa, Pb = self.prod_dict[mhomfly]
            c_mres = (sumkind,
                     self.classify_prime_homfly(Pa, ncross),
                     self.classify_prime_homfly(Pb, ncross))

        return p_res, c_res, c_mres


class PDStoreExpander(object):
    def __init__(self, dirloc=DEFAULT_PATH,
                 amortize=True, debug=False):
        self.amortize = amortize
        self.dirloc = dirloc

    def open(self, crossings_list, debug=False):
        pd_files = (self._db_file(x, self.dirloc) for x in crossings_list)
        iters = []
        for f in pd_files:
            _,npds,_ = self.parse_header(f, debug=debug)
            iters.append(self.read_pdstor(f, debug=debug, num_pds=int(npds)))
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
        return (len(pd.crossings),
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
    def crossing_sign_mask(pdc, signature):
        new_pd = pdc.copy()
        for crs, sign in zip(new_pd.crossings, signature):
            crs.sign = sign
        return new_pd

    @staticmethod
    def component_sign_mask(pdc, signature):
        new_pd = pdc.copy()
        for comp in compress(range(len(new_pd.components)), signature):
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
    def crossing_combinations(cls, pdc, bool amortize=True):
        """crossing_combinations(PlanarDiagram) -> generator(new
        PlanarDiagrams, masks)

        Generate all possible crossing sign permutations.  First crossing
        is fixed as (+) so as to avoid mirror images
        """
        for cross_set in cls._bin_strings(len(pdc.crossings), amortize):
            yield cls.crossing_sign_mask(pdc, cross_set), cross_set

    @classmethod
    def component_combinations(cls, pdc, amortize=True):
        """component_combinations(PlanarDiagram) -> generator(new
        PlanarDiagrams, masks)

        Generate all possible component permutations.
        """
        for comp_set in cls._bin_strings(len(pdc.components), amortize):
            yield cls.component_sign_mask(pdc, comp_set), comp_set

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
        return (cls.component_sign_mask(
            cls.crossing_sign_mask(
                cls.pd_by_filepos(cnum, pos, dirloc),
                cross_sgn), comp_sgn))

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

    def read_pdstor(self, f, debug=False, num_pds=None):
        # Read in the pd codes one-by-one.
        if debug: print "Reading file %s"%f.names
        for pd_id, pd in enumerate(self._pds_in_file(f)):
            if debug: print "> Reading %s of %s (%0.1f%%)"%(
                    pd_id+1,
                    "Unknown" if num_pds is None else num_pds,
                    "??" if num_pds is None else 100.0*pd_id/num_pds)
            for crs_pd, crs_sgn in self.crossing_combinations(pd, self.amortize):
                for cmp_pd, cmp_sgn in self.component_combinations(crs_pd, self.amortize):
                    homfly = cmp_pd.homfly()
                    pduid = self.uid(cmp_pd, pd_id, crs_sgn, cmp_sgn)
                    yield (cmp_pd, homfly, pduid)
        f.close()

class PDDatabase(PDStoreExpander):
    def save(self, fname):
        f = open(fname, "wb")
        cPickle.dump(self, f, cPickle.HIGHEST_PROTOCOL)
        f.close()

    def __init__(self, crossings_list=[3,4,5], dirloc=DEFAULT_PATH,
                 amortize=True, debug=False, callbacks=None):
        super(PDDatabase,self).__init__(dirloc=dirloc,
                                        amortize=amortize,
                                        debug=debug)
        self.callbacks = [] if callbacks is None else callbacks
        self.pd_dict = {}
        self.all_pds = {}
        for pd, homfly, uid in self.open(crossings_list):
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
