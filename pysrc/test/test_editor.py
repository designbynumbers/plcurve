from libpl.pdcode import PlanarDiagram, HOMFLYPolynomial
from test_pdcode import PlanarDiagramAssertMixin
import unittest
from plink import *
        
class TestPDCodeFromPlink(PlanarDiagramAssertMixin, unittest.TestCase):
    DATA_DIR = os.path.join("data", "plink")
    def _get_data_file_path(self, file_name):
        return os.path.join(self.DATA_DIR, file_name)

    def setUp(self):
        pass

    def test_figure_8(self):
        editor = LinkEditor(file_name=self._get_data_file_path("unknot_1x.lnk"))
        K = PlanarDiagram.from_plink(editor)
        self.assertEqual(K.homfly(), HOMFLYPolynomial("1"))

    def test_trefoil(self):
        editor = LinkEditor(file_name=self._get_data_file_path("3_1_simple.lnk"))
        K = PlanarDiagram.from_plink(editor)
        self.assertEqual(K.homfly(), HOMFLYPolynomial("-2a^{2} + a^{2}z^{2} + -a^{4}"))

if __name__=="__main__":
    unittest.main()
