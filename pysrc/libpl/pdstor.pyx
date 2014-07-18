from libpl.pdcode import PlanarDiagram
from itertools import combinations, chain, compress, product
from cpython cimport bool
import pprint
pp = pprint.PrettyPrinter(indent=2)
import os
import cPickle

DEFAULT_PATH=os.path.join("data","pdstors")

def bin_list_to_int(blist):
    return sum((i*2)**n for i,n in enumerate(reversed(blist)))

class ClassifyDatabase(object):
    KNOT=0
    LINK=1
    def __init__(self):
        self.cls_dict = dict()

    @classmethod
    def _iter_name_and_pd(self, names, table):
        name = names.readline().strip()
        pd = PlanarDiagram.read_knot_theory(table)
        while pd is not None and name is not None:
            yield name, pd
            del pd
            name = names.readline().strip()
            pd = PlanarDiagram.read_knot_theory(table)

    def _load_names_and_table(self, names_fname, table_fname, filetype):
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
            for i, (name, pd) in enumerate(self._iter_name_and_pd(names_f, table_f)):
                homfly = pd.homfly()
                if homfly in self.cls_dict:
                    self.cls_dict[homfly].append((name, (filetype, i)))
                else:
                    self.cls_dict[homfly] = [(name, (filetype, i)),]
                if i >= count-1:
                    break
        finally:
            if names_f:
                names_f.close()
            if table_f:
                table_f.close()

    def load_rolfsen(self, names_fname, table_fname):
        self._load_names_and_table(names_fname, table_fname, self.KNOT)
    def load_thistlethwaite(self, names_fname, table_fname):
        self._load_names_and_table(names_fname, table_fname, self.LINK)

class PDDatabase(object):
    def save(self, fname):
        f = open(fname, "wb")
        cPickle.dump(self, f, cPickle.HIGHEST_PROTOCOL)
        f.close()

    def __init__(self, crossings_list=[3,4,5,6], dirloc=DEFAULT_PATH,
                 amortize=True, debug=False):
        self.amortize = amortize
        pd_files = (self._db_file(ncrs, dirloc) for ncrs in crossings_list)
        self.pd_dict = {}
        for f in pd_files:
            _,npds,_ = self.parse_header(f, debug=debug)
            self.read_pdstor(f, debug=debug, num_pds=int(npds))
            f.close()

    @staticmethod
    def _db_file(cross_count, dirloc):
        return file(os.path.join(dirloc,"%d.pdstor"%cross_count))

    @staticmethod
    def load(fname):
        f = open(fname, "rb")
        ret = cPickle.load(f)
        f.close()
        return ret

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
        for crs in compress(new_pd.crossings, signature):
            crs.sign = (crs.sign+1)%2
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
        f.close()

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

    def add_pd(self, homfly, pduid):
        if homfly in self.pd_dict:
            self.pd_dict[homfly].add(pduid)
        else:
            self.pd_dict[homfly] = set([pduid])


    def read_pdstor(self, f, debug=False, num_pds=None):
        # Read in the pd codes one-by-one.
        if debug: print "Reading file %s"%f.name
        for pd_id, pd in enumerate(self._pds_in_file(f)):
            if debug: print "> Reading %s of %s (%0.1f%%)"%(
                    pd_id+1,
                    "Unknown" if num_pds is None else num_pds,
                    "??" if num_pds is None else 100.0*pd_id/num_pds)
            for crs_pd, crs_sgn in self.crossing_combinations(pd, self.amortize):
                for cmp_pd, cmp_sgn in self.component_combinations(crs_pd, self.amortize):
                    self.add_pd(cmp_pd.homfly(), self.uid(cmp_pd, pd_id, crs_sgn, cmp_sgn))
                    del cmp_pd
                del crs_pd
            del pd

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
