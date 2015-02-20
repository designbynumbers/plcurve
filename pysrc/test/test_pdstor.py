import unittest
from libpl.pdstor import *
import os.path

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
