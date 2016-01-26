import unittest
from libpl.pdstor import *
import os.path
import gzip
from suite import ROOT_DIR

class TestPDStorFiles(unittest.TestCase):
    DATADIR = os.path.join(ROOT_DIR, "..", "..", "data", "pdstors")
    
    def check_pdstor_isomorphism(self, n_cross):
        with open(os.path.join(self.DATADIR, "%s.pdstor"%n_cross), "rb") as pdstor_f:
            #pdstor_f.readline()
            #pdstor_f.readline()
            #pdstor_f.readline()
            for pdcode in PlanarDiagram.read_all(pdstor_f, read_header=True):
                try:
                    self.assertIsNone(pdcode.pdstor_index(pdstor_f))
                except EOFError:
                    pass

    @unittest.skip("Depends on pdstor data files, which are not included in dist")
    def test_pdstor_isomorphism(self):
        for i in range(3,2):
            self.check_pdstor_isomorphism(i)

    @unittest.skip("Depends on pdstor data files, which are not included in dist")
    def assert_no_isomorph_doubles(self, shadows):
        for i, shadow in enumerate(shadows):
            for other_shadow in shadows[i+1:]:
                self.assertFalse(shadow is other_shadow)
                self.assertNotEqual(shadow.uid, other_shadow.uid)
                self.assertFalse(shadow.isomorphic(other_shadow))

    @unittest.skip("Depends on pdstor data files, which are not included in dist")
    def check_pdstor_hashclasses(self, n_cross):
        with open(os.path.join(self.DATADIR, "%s.pdstor"%n_cross), "rb") as pdstor_f:
            pdstor = PDStorage.read(pdstor_f, read_header=True)
            print " PDstor read in."
            for hash_, bucket in pdstor.hashes.iteritems():
                self.assert_no_isomorph_doubles(bucket)

    @unittest.skip("Depends on pdstor data files, which are not included in dist")
    def test_pdstor_hashclasses(self):
        # Test passed for me but takes a long time on 8+
        for i in range(3,5):
            print "Checking hashclasses in %s.pdstor"%i
            self.check_pdstor_hashclasses(i)

    @unittest.skip("Depends on pdstor data files, which are not included in dist")
    def check_pdstor_hashes(self, n_cross):
        with open(os.path.join(self.DATADIR, "%s.pdstor"%n_cross), "rb") as pdstor_f:
            for pdcode in PlanarDiagram.read_all(pdstor_f, read_header=True):
                self.assertTrue(pdcode.hashed)

                oldhash = pdcode._hash
                pdcode.regenerate_hash()
                newhash = pdcode._hash
                self.assertEqual(oldhash, newhash)

    @unittest.skip("Depends on pdstor data files, which are not included in dist")
    def test_pdstor_hashes(self):
        for i in range(3,5):
            self.check_pdstor_hashes(i)

    def check_prime_subpdstor(self, n_cross):
        with open(os.path.join(self.DATADIR, "%s.pdstor"%n_cross), "rb") as pdstor_f, \
             open(os.path.join(self.DATADIR, "%sprime.pdstor"%n_cross), "rb") as primes_f:
            pdstor = PDStorage.read(pdstor_f, read_header=True)
            primes = PDStorage.read(primes_f, read_header=True)
            print " %s PDstors read in."%n_cross
            self.assertTrue(primes <= pdstor)

    @unittest.skip("Depends on pdstor data files, which are not included in dist")
    def test_prime_subpdstor(self):
        for i in range(3,8):
            self.check_prime_subpdstor(i)
            
class TestPDDatabase(unittest.TestCase):
    DATADIR = os.path.join(ROOT_DIR, "..", "..", "data", "pdstors")
    def setUp(self):
        pass
    
    @unittest.skip("Test takes a long time and requires pdstor files which are not in dist")
    def test_pd_dict_homflys(self):
        def assert_homfly(pd, homfly, uid):
            uid_pd = PDDatabase.pd_by_id(uid, dirloc=self.DATADIR)
            self.assertTrue(pd.isotopic(uid_pd))
            self.assertEqual(pd.homfly(), homfly)
            self.assertEqual(uid_pd.homfly(), homfly)
        pdb = PDDatabase(crossings_list=[3], dirloc=self.DATADIR, callbacks=[assert_homfly])
        

if __name__=="__main__":
    unittest.main()
