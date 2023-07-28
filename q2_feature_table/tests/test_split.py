# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
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

from q2_feature_table import split


class SplitTests(unittest.TestCase):

    def test_one_split(self):
        md_column = qiime2.CategoricalMetadataColumn(
            pd.Series(['a', 'a', 'a'], name='foo',
                      index=pd.Index(['S1', 'S2', 'S3'], name='id')))
        table = Table(np.array([[5, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = split(table, metadata=md_column)
        expected1 = Table(np.array([[5, 1, 3], [1, 1, 2]]),
                          ['O1', 'O2'],
                          ['S1', 'S2', 'S3'])
        self.assertEqual(actual, {'a': expected1})

    def test_two_splits(self):
        md_column = qiime2.CategoricalMetadataColumn(
            pd.Series(['a', 'b', 'a'], name='foo',
                      index=pd.Index(['S1', 'S2', 'S3'], name='id')))
        table = Table(np.array([[5, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = split(table, metadata=md_column)
        expected1 = Table(np.array([[5, 3], [1, 2]]),
                          ['O1', 'O2'],
                          ['S1', 'S3'])
        expected2 = Table(np.array([[1], [1]]),
                          ['O1', 'O2'],
                          ['S2',])
        self.assertEqual(actual, {'a': expected1, 'b': expected2})

    def test_three_splits(self):
        md_column = qiime2.CategoricalMetadataColumn(
            pd.Series(['a', 'b', 'c'], name='foo',
                      index=pd.Index(['S1', 'S2', 'S3'], name='id')))
        table = Table(np.array([[5, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = split(table, metadata=md_column)
        expected1 = Table(np.array([[5], [1]]),
                          ['O1', 'O2'],
                          ['S1'])
        expected2 = Table(np.array([[1], [1]]),
                          ['O1', 'O2'],
                          ['S2',])
        expected3 = Table(np.array([[3], [2]]),
                          ['O1', 'O2'],
                          ['S3'])
        self.assertEqual(actual,
                         {'a': expected1, 'b': expected2, 'c': expected3})

    def test_invalid_values(self):
        table = Table(np.array([[5, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])

        md_column = qiime2.CategoricalMetadataColumn(
            pd.Series(['a', 'b c', 'a'], name='foo',
                      index=pd.Index(['S1', 'S2', 'S3'], name='id')))
        with self.assertRaisesRegex(KeyError, 'invalid.*: b c'):
            split(table, metadata=md_column)

        md_column = qiime2.CategoricalMetadataColumn(
            pd.Series(['@a', '-a', '!a'], name='foo',
                      index=pd.Index(['S1', 'S2', 'S3'], name='id')))
        with self.assertRaisesRegex(KeyError, 'invalid.*: !a, @a'):
            split(table, metadata=md_column)

    def test_missing_data_samples_dropped(self):
        md_column = qiime2.CategoricalMetadataColumn(
            pd.Series(['a', 'b', np.nan], name='foo',
                      index=pd.Index(['S1', 'S2', 'S3'], name='id')))
        table = Table(np.array([[5, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = split(table, metadata=md_column)
        expected1 = Table(np.array([[5], [1]]),
                          ['O1', 'O2'],
                          ['S1'])
        expected2 = Table(np.array([[1], [1]]),
                          ['O1', 'O2'],
                          ['S2',])
        self.assertEqual(actual, {'a': expected1, 'b': expected2})

    def test_drop_features_w_zero_count(self):
        md_column = qiime2.CategoricalMetadataColumn(
            pd.Series(['a', 'a', 'b'], name='foo',
                      index=pd.Index(['S1', 'S2', 'S3'], name='id')))
        table = Table(np.array([[0, 0, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = split(table, metadata=md_column, filter_empty_features=True)
        expected1 = Table(np.array([[1, 1]]),
                          ['O2'],
                          ['S1', 'S2'])
        expected2 = Table(np.array([[3], [2]]),
                          ['O1', 'O2'],
                          ['S3',])
        self.assertEqual(actual, {'a': expected1, 'b': expected2})

        md_column = qiime2.CategoricalMetadataColumn(
            pd.Series(['a', 'a', 'b'], name='foo',
                      index=pd.Index(['S1', 'S2', 'S3'], name='id')))
        table = Table(np.array([[0, 0, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = split(table, metadata=md_column, filter_empty_features=False)
        expected1 = Table(np.array([[0, 0], [1, 1]]),
                          ['O1', 'O2'],
                          ['S1', 'S2'])
        expected2 = Table(np.array([[3], [2]]),
                          ['O1', 'O2'],
                          ['S3',])
        self.assertEqual(actual, {'a': expected1, 'b': expected2})
