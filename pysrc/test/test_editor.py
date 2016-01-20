from libpl.pdcode import PlanarDiagram, HOMFLYPolynomial
from test_pdcode import PlanarDiagramAssertMixin
import unittest
from plink import *

class TestPDCodeFromPlink(PlanarDiagramAssertMixin, unittest.TestCase):
    DATA_DIR = os.path.join("data", "plink")
    def _get_data_file_path(self, file_name):
        return os.path.join(self.DATA_DIR, file_name)

    def plink_file_to_Diagram(self, file_name):
        manager = LinkManager()
        with open(self._get_data_file_path(file_name)) as f:
            manager._from_string(f.read())
        return PlanarDiagram.from_plink(manager)

    def setUp(self):
        pass

    @unittest.skip("Requires the ability to make a TK window")
    def test_plink_LinkManager(self):
        """ test_plink_LinkManager()

        Asserts that a Link produced (really loaded) by a LinkEditor
        graphical editor is the same as the (silent) LinkManager loading
        method used to run the primary tests from this class.
        """
        FILENAME = "3_1_simple.lnk"
        editor = LinkEditor(file_name=self._get_data_file_path(FILENAME))
        K = PlanarDiagram.from_plink(editor)
        L = self.plink_file_to_Diagram(FILENAME)

        self.assertEqual(K, L)

    def test_figure_8(self):
        K = self.plink_file_to_Diagram("unknot_1x.lnk")
        self.assertEqual(K.homfly(), HOMFLYPolynomial("1"))

    def test_trefoil(self):
        K = self.plink_file_to_Diagram("3_1_simple.lnk")
        self.assertEqual(K.homfly(), HOMFLYPolynomial("-2a^{2} + a^{2}z^{2} + -a^{4}"))

    def test_unknot_trefoil_shadow(self):
        K = self.plink_file_to_Diagram("3_1_shadow_unknot.lnk")
        self.assertEqual(K.homfly(), HOMFLYPolynomial("1"))

    def test_hopf_link(self):
        links = {
            (1,1): self.plink_file_to_Diagram("hopf_link_pospos.lnk"),
            (1,0): self.plink_file_to_Diagram("hopf_link_posneg.lnk"),
            (0,1): self.plink_file_to_Diagram("hopf_link_negpos.lnk"),
            (0,0): self.plink_file_to_Diagram("hopf_link_negneg.lnk"),
        }
        self.assertEqual(
            links[(1,1)].homfly(),
            HOMFLYPolynomial("a^{-3}z^{-1} + a^{-1}z^{-1} + -a^{-1}z^{1}"))
        self.assertEqual(
            links[(0,0)].homfly(),
            HOMFLYPolynomial("a^{-3}z^{-1} + a^{-1}z^{-1} + -a^{-1}z^{1}"))
        self.assertEqual(
            links[(0,1)].homfly(),
            HOMFLYPolynomial("a^{1}z^{-1} + -a^{1}z^{1} + a^{3}z^{-1}"))
        self.assertEqual(
            links[(1,0)].homfly(),
            HOMFLYPolynomial("a^{1}z^{-1} + -a^{1}z^{1} + a^{3}z^{-1}"))

    @unittest.skip("Not implemented")
    def test_split_link(self):
        pass

    @unittest.skip("Not implemented")
    def test_hopf_trefoil_unknot(self):
        pass

    @unittest.skip("Not implemented")
    def test_split_trefoil_unknot(self):
        pass

if __name__=="__main__":
    unittest.main()
