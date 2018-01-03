# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import skbio
import numpy as np
from biom.table import Table
import biom.exception
import pandas as pd
import pandas.util.testing as pdt

from q2_feature_table import merge, merge_seqs, merge_taxa
from q2_feature_table._merge import _merge_feature_data


class MergeTableTests(unittest.TestCase):
    def test_single_table(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        obs = merge([t])

        self.assertEqual(t, obs)

    def test_valid_overlapping_feature_ids(self):
        t1 = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                   ['O1', 'O2'],
                   ['S1', 'S2', 'S3'])
        t2 = Table(np.array([[0, 2, 6], [2, 2, 4]]),
                   ['O1', 'O3'],
                   ['S4', 'S5', 'S6'])
        obs = merge([t1, t2])
        exp = Table(np.array([[0, 1, 3, 0, 2, 6], [1, 1, 2, 0, 0, 0],
                              [0, 0, 0, 2, 2, 4]]),
                    ['O1', 'O2', 'O3'],
                    ['S1', 'S2', 'S3', 'S4', 'S5', 'S6'])
        self.assertEqual(obs, exp)

    def test_valid_non_overlapping_feature_ids(self):
        t1 = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                   ['O1', 'O2'],
                   ['S1', 'S2', 'S3'])
        t2 = Table(np.array([[0, 2, 6], [2, 2, 4]]),
                   ['O3', 'O4'],
                   ['S4', 'S5', 'S6'])
        obs = merge([t1, t2])
        exp = Table(np.array([[0, 1, 3, 0, 0, 0], [1, 1, 2, 0, 0, 0],
                              [0, 0, 0, 0, 2, 6], [0, 0, 0, 2, 2, 4]]),
                    ['O1', 'O2', 'O3', 'O4'],
                    ['S1', 'S2', 'S3', 'S4', 'S5', 'S6'])
        self.assertEqual(obs, exp)

    def test_invalid_overlapping_feature_ids(self):
        t1 = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                   ['O1', 'O2'],
                   ['S1', 'S2', 'S3'])
        t2 = Table(np.array([[0, 2, 6], [2, 2, 4]]),
                   ['O1', 'O3'],
                   ['S4', 'S5', 'S6'])
        with self.assertRaisesRegex(biom.exception.DisjointIDError,
                                    'IDs are not disjoint'):
            merge([t1, t2], 'error_on_overlapping_feature')

    def test_valid_overlapping_sample_ids(self):
        t1 = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                   ['O1', 'O2'],
                   ['S1', 'S2', 'S3'])
        t2 = Table(np.array([[0, 2, 6], [2, 2, 4]]),
                   ['O3', 'O4'],
                   ['S1', 'S5', 'S6'])
        obs = merge([t1, t2], 'error_on_overlapping_feature')
        exp = Table(np.array([[0, 1, 3, 0, 0], [1, 1, 2, 0, 0],
                              [0, 0, 0, 2, 6], [2, 0, 0, 2, 4]]),
                    ['O1', 'O2', 'O3', 'O4'],
                    ['S1', 'S2', 'S3', 'S5', 'S6'])
        self.assertEqual(obs, exp)

    def test_invalid_overlapping_sample_ids(self):
        t1 = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                   ['O1', 'O2'],
                   ['S1', 'S2', 'S3'])
        t2 = Table(np.array([[0, 2, 6], [2, 2, 4]]),
                   ['O1', 'O3'],
                   ['S1', 'S5', 'S6'])
        with self.assertRaisesRegex(biom.exception.DisjointIDError,
                                    'IDs are not disjoint'):
            merge([t1, t2])

    def test_invalid_overlap_method(self):
        t1 = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                   ['O1', 'O2'],
                   ['S1', 'S2', 'S3'])
        t2 = Table(np.array([[0, 2, 6], [2, 2, 4]]),
                   ['O1', 'O3'],
                   ['S1', 'S5', 'S6'])
        with self.assertRaisesRegex(ValueError, 'overlap method'):
            merge([t1, t2], 'peanut')

    def test_sum_full_overlap(self):
        t1 = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                   ['O1', 'O2'],
                   ['S1', 'S2', 'S3'])
        t2 = Table(np.array([[0, 2, 6], [2, 2, 4]]),
                   ['O1', 'O2'],
                   ['S1', 'S2', 'S3'])
        obs = merge([t1, t2], 'sum')
        exp = Table(np.array([[0, 3, 9], [3, 3, 6]]),
                    ['O1', 'O2'],
                    ['S1', 'S2', 'S3'])
        self.assertEqual(obs, exp)

    def test_sum_triple_overlap(self):
        t1 = Table(np.array([[1, 1, 1], [1, 1, 1]]),
                   ['O1', 'O2'],
                   ['S1', 'S2', 'S3'])
        obs = merge([t1] * 3, 'sum')
        exp = Table(np.array([[3, 3, 3], [3, 3, 3]]),
                    ['O1', 'O2'],
                    ['S1', 'S2', 'S3'])
        self.assertEqual(obs, exp)

    def test_sum_some_overlap(self):
        # Did I stutter?
        t1 = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                   ['O1', 'O2'],
                   ['S1', 'S2', 'S3'])
        t2 = Table(np.array([[0, 2, 6], [2, 2, 4]]),
                   ['O1', 'O3'],
                   ['S4', 'S2', 'S5'])
        obs = merge([t1, t2], 'sum')
        exp = Table(np.array([[0, 3, 3, 0, 6], [1, 1, 2, 0, 0],
                              [0, 2, 0, 2, 4]]),
                    ['O1', 'O2', 'O3'],
                    ['S1', 'S2', 'S3', 'S4', 'S5'])
        self.assertEqual(obs, exp)

    def test_sum_overlapping_sample_ids(self):
        # This should produce the same result as `error_on_overlapping_feature`
        t1 = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                   ['O1', 'O2'],
                   ['S1', 'S2', 'S3'])
        t2 = Table(np.array([[0, 2, 6], [2, 2, 4]]),
                   ['O3', 'O4'],
                   ['S1', 'S5', 'S6'])
        obs = merge([t1, t2], 'sum')
        exp = Table(np.array([[0, 1, 3, 0, 0], [1, 1, 2, 0, 0],
                              [0, 0, 0, 2, 6], [2, 0, 0, 2, 4]]),
                    ['O1', 'O2', 'O3', 'O4'],
                    ['S1', 'S2', 'S3', 'S5', 'S6'])
        self.assertEqual(obs, exp)

    def test_sum_overlapping_feature_ids(self):
        # This should produce the same result as `error_on_overlapping_sample`
        t1 = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                   ['O1', 'O2'],
                   ['S1', 'S2', 'S3'])
        t2 = Table(np.array([[0, 2, 6], [2, 2, 4]]),
                   ['O1', 'O3'],
                   ['S4', 'S5', 'S6'])
        obs = merge([t1, t2], 'sum')
        exp = Table(np.array([[0, 1, 3, 0, 2, 6], [1, 1, 2, 0, 0, 0],
                              [0, 0, 0, 2, 2, 4]]),
                    ['O1', 'O2', 'O3'],
                    ['S1', 'S2', 'S3', 'S4', 'S5', 'S6'])
        self.assertEqual(obs, exp)


class MergeFeatureDataTests(unittest.TestCase):
    def test_merge_single(self):
        d = pd.Series(['ACGT', 'ACCT'], index=['f1', 'f2'])
        obs = _merge_feature_data([d])
        pdt.assert_series_equal(obs, d)

    def test_valid_overlapping_feature_ids(self):
        d1 = pd.Series(['ACGT', 'ACCT'], index=['f1', 'f2'])
        d2 = pd.Series(['ACGT', 'ACCA'], index=['f1', 'f3'])
        obs = _merge_feature_data([d1, d2])
        exp = pd.Series(['ACGT', 'ACCT', 'ACCA'], index=['f1', 'f2', 'f3'])
        pdt.assert_series_equal(obs, exp)

    def test_first_feature_data_retained(self):
        d1 = pd.Series(['ACGT', 'ACCT'], index=['f1', 'f2'])
        d2 = pd.Series(['ACGAAA', 'ACCA'], index=['f1', 'f3'])

        obs = _merge_feature_data([d1, d2])
        exp = pd.Series(['ACGT', 'ACCT', 'ACCA'], index=['f1', 'f2', 'f3'])
        pdt.assert_series_equal(obs, exp)

        # swapping input order changes f1 data
        obs = _merge_feature_data([d2, d1])
        exp = pd.Series(['ACGAAA', 'ACCT', 'ACCA'], index=['f1', 'f2', 'f3'])
        pdt.assert_series_equal(obs, exp)

    def test_multiple_overlapping_feature_ids_order_maintained(self):
        d1 = pd.Series(['ACGT', 'ACCT'], index=['f1', 'f2'])
        d2 = pd.Series(['ACGAAA', 'ACCA'], index=['f1', 'f3'])
        d3 = pd.Series(['AGGA', 'ATAT'], index=['f3', 'f4'])

        obs = _merge_feature_data([d1, d2, d3])
        exp = pd.Series(['ACGT', 'ACCT', 'ACCA', 'ATAT'],
                        index=['f1', 'f2', 'f3', 'f4'])
        pdt.assert_series_equal(obs, exp)

        # swapping input order changes f1 and f3
        obs = _merge_feature_data([d3, d2, d1])
        exp = pd.Series(['ACGAAA', 'ACCT', 'AGGA', 'ATAT'],
                        index=['f1', 'f2', 'f3', 'f4'])
        pdt.assert_series_equal(obs, exp)

    def test_valid_non_overlapping_feature_ids(self):
        d1 = pd.Series(['ACGT', 'ACCT'], index=['f1', 'f2'])
        d2 = pd.Series(['ACGT', 'ACCA'], index=['f3', 'f4'])
        obs = _merge_feature_data([d1, d2])
        exp = pd.Series(['ACGT', 'ACCT', 'ACGT', 'ACCA'],
                        index=['f1', 'f2', 'f3', 'f4'])
        pdt.assert_series_equal(obs, exp)


class MergeFeatureSequenceTests(unittest.TestCase):
    # More extensive testing is performed in MergeFeatureDataTests, which
    # tests the shared private API.

    def test_merge_seqs(self):
        d1 = pd.Series([skbio.DNA('ACGT', metadata={'id': 'abc'}),
                        skbio.DNA('ACCT', metadata={'id': 'xyz'})],
                       index=['f1', 'f2'])
        d2 = pd.Series([skbio.DNA('ACGT', metadata={'id': 'abc'}),
                        skbio.DNA('ACCA', metadata={'id': 'wxy'})],
                       index=['f1', 'f3'])
        obs = merge_seqs([d1, d2])
        exp = pd.Series([skbio.DNA('ACGT', metadata={'id': 'abc'}),
                         skbio.DNA('ACCT', metadata={'id': 'xyz'}),
                         skbio.DNA('ACCA', metadata={'id': 'wxy'})],
                        index=['f1', 'f2', 'f3'])
        pdt.assert_series_equal(obs, exp)


class MergeFeatureTaxonomyTests(unittest.TestCase):
    # More extensive testing is performed in MergeFeatureDataTests, which
    # tests the shared private API.

    def test_merge_taxa(self):
        # this test calls the public API directly
        d1 = pd.Series(['a;b;c;d', 'a;b;c;e'], index=['f1', 'f2'])
        d2 = pd.Series(['a;b;c;d', 'a;b;c;e'], index=['f1', 'f3'])
        obs = merge_taxa([d1, d2])
        exp = pd.Series(['a;b;c;d', 'a;b;c;e', 'a;b;c;e'],
                        index=['f1', 'f2', 'f3'])
        pdt.assert_series_equal(obs, exp)


if __name__ == "__main__":
    unittest.main()
