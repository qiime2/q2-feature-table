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
import pandas as pd
import pandas.util.testing as pdt

from q2_feature_table import (merge, merge_seq_data,
                              merge_taxa_data)
from q2_feature_table._merge import _merge_feature_data


class MergeTableTests(unittest.TestCase):

    def test_valid_overlapping_feature_ids(self):
        t1 = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                   ['O1', 'O2'],
                   ['S1', 'S2', 'S3'])
        t2 = Table(np.array([[0, 2, 6], [2, 2, 4]]),
                   ['O1', 'O3'],
                   ['S4', 'S5', 'S6'])
        obs = merge(t1, t2)
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
        obs = merge(t1, t2)
        exp = Table(np.array([[0, 1, 3, 0, 0, 0], [1, 1, 2, 0, 0, 0],
                              [0, 0, 0, 0, 2, 6], [0, 0, 0, 2, 2, 4]]),
                    ['O1', 'O2', 'O3', 'O4'],
                    ['S1', 'S2', 'S3', 'S4', 'S5', 'S6'])
        self.assertEqual(obs, exp)

    def test_invalid_overlapping_sample_ids(self):
        t1 = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                   ['O1', 'O2'],
                   ['S1', 'S2', 'S3'])
        t2 = Table(np.array([[0, 2, 6], [2, 2, 4]]),
                   ['O1', 'O3'],
                   ['S1', 'S5', 'S6'])
        with self.assertRaises(ValueError):
            merge(t1, t2)


class MergeFeatureDataTests(unittest.TestCase):

    def test_valid_overlapping_feature_ids(self):
        d1 = pd.Series(['ACGT', 'ACCT'], index=['f1', 'f2'])
        d2 = pd.Series(['ACGT', 'ACCA'], index=['f1', 'f3'])
        obs = _merge_feature_data(d1, d2)
        exp = pd.Series(['ACGT', 'ACCT', 'ACCA'], index=['f1', 'f2', 'f3'])
        pdt.assert_series_equal(obs, exp)

    def test_first_feature_data_retained(self):
        d1 = pd.Series(['ACGT', 'ACCT'], index=['f1', 'f2'])
        d2 = pd.Series(['ACGAAA', 'ACCA'], index=['f1', 'f3'])

        obs = _merge_feature_data(d1, d2)
        exp = pd.Series(['ACGT', 'ACCT', 'ACCA'], index=['f1', 'f2', 'f3'])
        pdt.assert_series_equal(obs, exp)

        # swapping input order changes f1 data
        obs = _merge_feature_data(d2, d1)
        exp = pd.Series(['ACGAAA', 'ACCT', 'ACCA'], index=['f1', 'f2', 'f3'])
        pdt.assert_series_equal(obs, exp)

    def test_valid_non_overlapping_feature_ids(self):
        d1 = pd.Series(['ACGT', 'ACCT'], index=['f1', 'f2'])
        d2 = pd.Series(['ACGT', 'ACCA'], index=['f3', 'f4'])
        obs = _merge_feature_data(d1, d2)
        exp = pd.Series(['ACGT', 'ACCT', 'ACGT', 'ACCA'],
                        index=['f1', 'f2', 'f3', 'f4'])
        pdt.assert_series_equal(obs, exp)


class MergeFeatureSequenceTests(unittest.TestCase):
    # More extensive testing is performed in MergeFeatureDataTests, which
    # tests the shared private API.

    def test_merge_seq_data(self):
        d1 = pd.Series([skbio.DNA('ACGT', metadata={'id': 'abc'}),
                        skbio.DNA('ACCT', metadata={'id': 'xyz'})],
                       index=['f1', 'f2'])
        d2 = pd.Series([skbio.DNA('ACGT', metadata={'id': 'abc'}),
                        skbio.DNA('ACCA', metadata={'id': 'wxy'})],
                       index=['f1', 'f3'])
        obs = merge_seq_data(d1, d2)
        exp = pd.Series([skbio.DNA('ACGT', metadata={'id': 'abc'}),
                         skbio.DNA('ACCT', metadata={'id': 'xyz'}),
                         skbio.DNA('ACCA', metadata={'id': 'wxy'})],
                        index=['f1', 'f2', 'f3'])
        pdt.assert_series_equal(obs, exp)


class MergeFeatureTaxonomyTests(unittest.TestCase):
    # More extensive testing is performed in MergeFeatureDataTests, which
    # tests the shared private API.

    def test_merge_taxa_data(self):
        # this test calls the public API directly
        d1 = pd.Series(['a;b;c;d', 'a;b;c;e'], index=['f1', 'f2'])
        d2 = pd.Series(['a;b;c;d', 'a;b;c;e'], index=['f1', 'f3'])
        obs = merge_taxa_data(d1, d2)
        exp = pd.Series(['a;b;c;d', 'a;b;c;e', 'a;b;c;e'],
                        index=['f1', 'f2', 'f3'])
        pdt.assert_series_equal(obs, exp)


if __name__ == "__main__":
    unittest.main()
