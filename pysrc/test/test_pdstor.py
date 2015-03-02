import unittest
from libpl.pdstor import *
import os.path

class TestPDStorFiles(unittest.TestCase):
    DATADIR = os.path.join("..", "..", "data", "pdstors")
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
                
    def test_pdstor_isomorphism(self):
        for i in range(3,2):
            self.check_pdstor_isomorphism(i)

    def assert_no_isomorph_doubles(self, shadows):
        for i, shadow in enumerate(shadows):
            for other_shadow in shadows[i+1:]:
                self.assertFalse(shadow is other_shadow)
                self.assertNotEqual(shadow.uid, other_shadow.uid)
                self.assertFalse(shadow.isomorphic(other_shadow))

    def check_pdstor_hashclasses(self, n_cross):
        with open(os.path.join(self.DATADIR, "%s.pdstor"%n_cross), "rb") as pdstor_f:
            pdstor = PDStorage.read(pdstor_f, read_header=True)
            print " PDstor read in."
            for hash_, bucket in pdstor.hashes.iteritems():
                self.assert_no_isomorph_doubles(bucket)
                
    def test_pdstor_hashclasses(self):
        # Test passed for me but takes a long time on 8+
        for i in range(3,8):
            print "Checking hashclasses in %s.pdstor"%i
            self.check_pdstor_hashclasses(i)

    def check_pdstor_hashes(self, n_cross):
        with open(os.path.join(self.DATADIR, "%s.pdstor"%n_cross), "rb") as pdstor_f:
            for pdcode in PlanarDiagram.read_all(pdstor_f, read_header=True):
                self.assertTrue(pdcode.hashed)

                oldhash = pdcode._hash
                pdcode.regenerate_hash()
                newhash = pdcode._hash
                self.assertEqual(oldhash, newhash)
                    
    def test_pdstor_hashes(self):
        for i in range(3,6):
            self.check_pdstor_hashes(i)
            
class TestPDDatabase(unittest.TestCase):
    DATADIR = os.path.join("..", "..", "data", "pdstors")
    def setUp(self):
        pass
    
    #@unittest.skip("Test takes a long time")
    def test_pd_dict_homflys(self):
        def assert_homfly(pd, homfly, uid):
            uid_pd = PDDatabase.pd_by_id(uid, dirloc=self.DATADIR)
            self.assertTrue(pd.isotopic(uid_pd))
            self.assertEqual(pd.homfly(), homfly)
            self.assertEqual(uid_pd.homfly(), homfly)
        pdb = PDDatabase(crossings_list=[3], dirloc=self.DATADIR, callbacks=[assert_homfly])
        

if __name__=="__main__":
    unittest.main()
