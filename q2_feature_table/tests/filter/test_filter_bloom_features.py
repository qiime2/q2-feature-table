# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import qiime2
import numpy as np
import pandas as pd
from biom.table import Table

from q2_feature_table import filter_bloom_features


class FilterBloomFeaturesTests(unittest.TestCase):
    def setUp(self):
        self.seqs = pd.Series(['ACGT', 'GCTA'],
                              index=['O1', 'O2'])
        self.drop_seqs = pd.Series(['GCTAT', 'TGACGAC'])
        self.table = Table(np.array([[0, 1, 1], [1, 1, 2]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])

    def test_filter_blooms(self):
        actual = filter_bloom_features(self.table, self.seqs,
                                       bloom_sequences=self.drop_seqs)
        expected = Table(np.array([[0, 1, 1]]),
                         ['O1'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

    def test_no_bloom_sequences_specified(self):
        seq1 = "TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTT"\
               "AAGTTGAATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCCAAAACTGGCAAGCTAGA"\
               "GTATGGTAGAGGGTAGTGGAATTTCCTG"
        seq2 = "TGACAC"
        seq3 = "TACG"
        seqs = pd.Series([seq1, seq2, seq3], index=['O1', 'O2', 'O3'])
        actual = filter_bloom_features(self.table, seqs)
        expected = Table(np.array([[1, 1, 2]]),
                         ['O2'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

    def test_table_seqs_have_different_lengths(self):
        seqs = pd.Series(['ACGT', 'GCTAT'],
                         index=['O1', 'O2'])

        actual = filter_bloom_features(self.table, seqs,
                                       bloom_sequences=self.drop_seqs)
        expected = Table(np.array([[0, 1, 1]]),
                         ['O1'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

    def test_table_has_sequences_longer_than_bloom_seqs(self):
        table = Table(np.array([[0, 1, 1], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        seqs = pd.Series(['ACGTCGAGCT', 'AGCT'], index=['O1', 'O2'])
        with self.assertRaisesRegex(ValueError, 'must be at least as long as'):
            filter_bloom_features(table, seqs, bloom_sequences=self.drop_seqs)

    def test_table_has_sequences_not_in_feature_data(self):
        table = Table(np.array([[0, 1, 1], [1, 1, 2]]),
                      ['O1', 'O3'],
                      ['S1', 'S2', 'S3'])
        with self.assertRaisesRegex(ValueError,
                                    'Table features must all be present in '
                                    'the Sequence Data.'):
            filter_bloom_features(table, self.seqs,
                                  bloom_sequences=self.drop_seqs)

    def test_all_features_filtered(self):
        table = Table(np.array([[0, 1, 1], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        drop_seqs = pd.Series(['GCTA', 'ACGT'])

        with self.assertRaisesRegex(ValueError, 'All features were filtered '
                                                'out of the data'):
            filter_bloom_features(table, self.seqs, bloom_sequences=drop_seqs)


if __name__ == "__main__":
    unittest.main()
