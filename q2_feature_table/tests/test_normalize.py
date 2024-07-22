# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
from unittest import TestCase, main
from unittest.mock import MagicMock, patch

import numpy as np
import numpy.testing as npt
import pandas as pd
from biom.table import Table
from pandas._testing import assert_series_equal
from q2_types.feature_data import SequenceCharacteristicsDirectoryFormat

from q2_feature_table import rarefy
from q2_feature_table._normalize import (_validate_parameters,
                                         _convert_lengths, normalize)


class RarefyTests(TestCase):

    def test_rarefy(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        a = rarefy(t, 2)
        self.assertEqual(a.shape, (2, 2))
        self.assertEqual(set(a.ids(axis='sample')), set(['S2', 'S3']))
        self.assertEqual(set(a.ids(axis='observation')), set(['O1', 'O2']))
        npt.assert_array_equal(a.sum(axis='sample'), np.array([2., 2.]))

    def test_rarefy_replacement(self):
        t = Table(np.array([[0, 10, 30], [10, 10, 20]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        rt = rarefy(t, 3, with_replacement=True)
        self.assertEqual(rt.shape, (2, 3))

        # IMPORTANT: samples below subsample depth should be removed
        for n_draws in range(11, 21):
            rt = rarefy(t, n_draws, with_replacement=True)
            npt.assert_array_equal(rt.sum('sample'),
                                   np.array([n_draws] * 2))
        for n_draws in range(21, 50):
            rt = rarefy(t, n_draws, with_replacement=True)
            npt.assert_array_equal(rt.sum('sample'),
                                   np.array([n_draws] * 1))

    def test_rarefy_depth_error(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])

        with self.assertRaisesRegex(ValueError, 'shallow enough'):
            rarefy(t, 50)


class NormalizeTests(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.lengths = pd.Series(
            {
                "ARO1": 1356.0,
                "ARO2": 1173.0,
            },
            name="values",
        )
        cls.lengths.index.name = "index"
        cls.table = pd.DataFrame({
            'ID': ['sample1', 'sample2'],
            'ARO1': [2.0, 2.0],
            'ARO2': [0.0, 0.0]
        }).set_index('ID')

    def test_validate_parameters_uq_with_m_a_trim(self):
        # Test Error raised if gene-length is given with UQ method
        with self.assertRaisesRegex(
                ValueError,
                "Parameters m-trim and a-trim can only "
                "be used with methods TMM and CTF.",
        ):
            _validate_parameters("uq", 0.2, 0.05, None)

    def test_validate_parameters_tpm_missing_gene_length(self):
        # Test Error raised if gene-length is missing with TPM method
        with self.assertRaisesRegex(
                ValueError, "gene-length input is missing."):
            _validate_parameters("tpm", None, None, None)

    def test_validate_parameters_tmm_gene_length(self):
        # Test Error raised if gene-length is given with TMM method
        with self.assertRaisesRegex(
                ValueError,
                "gene-length input can only be used with FPKM and "
                "TPM methods."
        ):
            _validate_parameters(
                "tmm", None, None, gene_length=MagicMock())

    def test_validate_parameters_default_m_a_trim(self):
        # Test if m_trim and a_trim get set to default values if None
        m_trim, a_trim = _validate_parameters("tmm", None, None, None)
        self.assertEqual(m_trim, 0.3)
        self.assertEqual(a_trim, 0.05)

    def test_validate_parameters_m_a_trim(self):
        # Test if m_trim and a_trim are not modified if not None
        m_trim, a_trim = _validate_parameters("tmm", 0.1, 0.06, None)
        self.assertEqual(m_trim, 0.1)
        self.assertEqual(a_trim, 0.06)

    def test_convert_lengths_gene_length(self):
        # Test _convert_lengths
        gene_length = SequenceCharacteristicsDirectoryFormat()
        with open(os.path.join(
                str(gene_length), "sequence_characteristics.tsv"),
                  'w') as file:
            file.write("id\tlength\nARO1\t1356.0\nARO2\t1173.0")

        obs = _convert_lengths(self.table, gene_length=gene_length)
        assert_series_equal(obs, self.lengths)

    def test_convert_lengths_short_gene_length(self):
        # Test Error raised if gene-length is missing genes
        gene_length = SequenceCharacteristicsDirectoryFormat()
        with open(os.path.join(
                str(gene_length),
                "sequence_characteristics.tsv"), 'w') as file:
            file.write("id\tlength\nARO1\t1356.0")
        with self.assertRaisesRegex(
                ValueError,
                "There are genes present in the FeatureTable that are "
                "not present in the gene-length input. Missing lengths "
                "for genes: {'ARO2'}",
        ):
            _convert_lengths(self.table, gene_length=gene_length)

    @patch("q2_feature_table._normalize.TPM")
    def test_tpm_fpkm_with_valid_inputs(self, mock_tpm):
        # Test valid inputs for TPM method
        gene_length = SequenceCharacteristicsDirectoryFormat()
        with open(os.path.join(
                str(gene_length), "sequence_characteristics.tsv"),
                  'w') as file:
            file.write("id\tlength\nARO1\t1356.0\nARO2\t1173.0")
        normalize(table=self.table, gene_length=gene_length, method="tpm")

    @patch("q2_feature_table._normalize.TMM")
    def test_tmm_uq_cuf_ctf_with_valid_inputs(self, mock_tmm):
        # Test valid inputs for TMM method
        gene_length = SequenceCharacteristicsDirectoryFormat()
        with open(os.path.join(
                str(gene_length), "sequence_characteristics.tsv"),
                  'w') as file:
            file.write("id\tlength\nARO1\t1356.0\nARO2\t1173.0")
        normalize(table=self.table, method="tmm", a_trim=0.06, m_trim=0.4)


if __name__ == "__main__":
    main()
