import unittest
from libpl.pdcode import *
from test_pdcode import BeforeAfterFileMixin, PlanarDiagramAssertMixin
from suite import ROOT_DIR
import os

class TestR2BigonElimination(BeforeAfterFileMixin, PlanarDiagramAssertMixin, unittest.TestCase):
    DATA_DIR = os.path.join(ROOT_DIR, "data")
    DATA_BEFORE_TEMPLATE = "sr2_test%s_before.pdstor"
    DATA_AFTER_TEMPLATE = "sr2_test%s_after.pdstor"

    def checkBigonElimination(self, tag, face_n):
        with open(self._get_before_fname(tag)) as before_f:
            before_pd = PlanarDiagram.read(before_f, read_header=True)
        self.assertIsNotNone(before_pd)

        with open(self._get_after_fname(tag)) as after_f:
            after_pds = tuple(PlanarDiagram.read_all(after_f, read_header=True))
        self.assertIsNotNone(after_pds)
        for after_pd in after_pds:
            self.assertIsNotNone(after_pd)

        before_pd_copy = before_pd.copy()
        face = before_pd.faces[face_n]
        result_pds = before_pd.R2_bigon_elimination(face)

        self.assertEqual(before_pd, before_pd_copy)
        self.assertEqual(len(result_pds), len(after_pds))
        for result_pd, after_pd in zip(result_pds, after_pds):
            self.assertTrue(result_pd.isomorphic(after_pd))

    def test_A(self):
        self.checkBigonElimination("A", 18)
