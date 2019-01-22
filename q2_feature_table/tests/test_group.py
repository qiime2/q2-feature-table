# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import biom
import qiime2
import pandas as pd
import numpy as np

from q2_feature_table import group


class TestGroup(unittest.TestCase):
    def test_identity_groups(self):
        # These map to the same values as before
        sample_mc = qiime2.CategoricalMetadataColumn(
            pd.Series(['a', 'b', 'c'], name='foo',
                      index=pd.Index(['a', 'b', 'c'], name='sampleid')))
        feature_mc = qiime2.CategoricalMetadataColumn(
            pd.Series(['x', 'y'], name='foo',
                      index=pd.Index(['x', 'y'], name='featureid')))
        table = biom.Table(np.array([[1, 2, 3], [30, 20, 10]]),
                           sample_ids=sample_mc.to_series().index,
                           observation_ids=feature_mc.to_series().index)

        # Sample x Sum
        result = group(table, axis='sample', metadata=sample_mc, mode='sum')
        self.assertEqual(table, result)

        # Sample x Mean
        result = group(table, axis='sample', metadata=sample_mc,
                       mode='mean-ceiling')
        self.assertEqual(table, result)

        # Sample x Median
        result = group(table, axis='sample', metadata=sample_mc,
                       mode='median-ceiling')
        self.assertEqual(table, result)

        # Feature x Sum
        result = group(table, axis='feature', metadata=feature_mc, mode='sum')
        self.assertEqual(table, result)

        # Feature x Mean
        result = group(table, axis='feature', metadata=feature_mc,
                       mode='mean-ceiling')
        self.assertEqual(table, result)

        # Feature x Median
        result = group(table, axis='feature', metadata=feature_mc,
                       mode='median-ceiling')
        self.assertEqual(table, result)

    def test_one_to_one_rename(self):
        sample_mc = qiime2.CategoricalMetadataColumn(
            pd.Series(['a_new', 'b_new', 'c_new'], name='foo',
                      index=pd.Index(['a', 'b', 'c'], name='sampleid')))
        original_sample_ids = sample_mc.to_series().index
        new_sample_ids = list(sample_mc.to_series())

        feature_mc = qiime2.CategoricalMetadataColumn(
            pd.Series(['x_new', 'y_new'], name='foo',
                      index=pd.Index(['x', 'y'], name='featureid')))
        original_feature_ids = feature_mc.to_series().index
        new_feature_ids = list(feature_mc.to_series())

        data = np.array([[1, 2, 3], [30, 20, 10]])
        table = biom.Table(data, sample_ids=original_sample_ids,
                           observation_ids=original_feature_ids)

        # Sample renames
        expected = biom.Table(data, sample_ids=new_sample_ids,
                              observation_ids=original_feature_ids)

        # Sample x Sum
        result = group(table, axis='sample', metadata=sample_mc, mode='sum')
        self.assertEqual(expected, result)

        # Sample X Mean
        result = group(table, axis='sample', metadata=sample_mc,
                       mode='mean-ceiling')
        self.assertEqual(expected, result)

        # Sample X Mean
        result = group(table, axis='sample', metadata=sample_mc,
                       mode='median-ceiling')
        self.assertEqual(expected, result)

        # Feature renames
        expected = biom.Table(data, sample_ids=original_sample_ids,
                              observation_ids=new_feature_ids)

        # Feature X Sum
        result = group(table, axis='feature', metadata=feature_mc, mode='sum')
        self.assertEqual(expected, result)

        # Feature X Mean
        result = group(table, axis='feature', metadata=feature_mc,
                       mode='mean-ceiling')
        self.assertEqual(expected, result)

        # Feature X Median
        result = group(table, axis='feature', metadata=feature_mc,
                       mode='median-ceiling')
        self.assertEqual(expected, result)

    def test_superset_feature_group(self):
        feature_mc = qiime2.CategoricalMetadataColumn(
            pd.Series(
                ['g0', 'g0', 'g1', 'g2', 'g1', 'g2', 'extra'], name='foo',
                index=pd.Index(['a', 'b', 'c', 'd', 'e', 'f', 'g'],
                               name='featureid')))
        data = np.array([
            [1, 0, 0],
            [1, 10, 10],
            [0, 0, 100],
            [5, 5, 5],
            [0, 1, 100],
            [7, 8, 9]])
        # g is missing on purpose
        table = biom.Table(data, sample_ids=['s1', 's2', 's3'],
                           observation_ids=['a', 'b', 'c', 'd', 'e', 'f'])

        expected = biom.Table(
                np.array([[2, 10, 10], [0, 1, 200], [12, 13, 14]]),
                sample_ids=['s1', 's2', 's3'],
                observation_ids=['g0', 'g1', 'g2'])
        result = group(table, axis='feature', metadata=feature_mc, mode='sum')
        self.assertEqual(expected, result)

    def test_superset_sample_group(self):
        sample_mc = qiime2.CategoricalMetadataColumn(
            pd.Series(['g0', 'g1', 'g2', 'g0', 'g1', 'g2'], name='foo',
                      index=pd.Index(['s1', 's2', 's3', 's4', 's5', 's6'],
                                     name='sampleid')))
        data = np.array([
            [0, 1, 2, 3],
            [10, 11, 12, 13],
            [100, 110, 120, 130]])
        table = biom.Table(data, sample_ids=['s1', 's2', 's4', 's5'],
                           observation_ids=['x', 'y', 'z'])

        expected = biom.Table(
            np.array([[2, 4], [22, 24], [220, 240]]),
            sample_ids=['g0', 'g1'],
            observation_ids=['x', 'y', 'z'])

        result = group(table, axis='sample', metadata=sample_mc, mode='sum')
        self.assertEqual(expected, result)

    def test_missing_feature_ids(self):
        feature_mc = qiime2.CategoricalMetadataColumn(
            pd.Series(['g0', 'g1', 'g2', 'g1', 'g2', 'extra'], name='foo',
                      index=pd.Index(['a', 'c', 'd', 'e', 'f', 'g'],
                                     name='featureid')))
        data = np.array([
            [1, 0, 0],
            [1, 10, 10],
            [0, 0, 100],
            [5, 5, 5],
            [0, 1, 100],
            [7, 8, 9]])
        # g is missing on purpose
        table = biom.Table(data, sample_ids=['s1', 's2', 's3'],
                           observation_ids=['a', 'b', 'c', 'd', 'e', 'f'])

        with self.assertRaisesRegex(ValueError, "not present.*'b'"):
            group(table, axis='feature', metadata=feature_mc, mode='sum')

    def test_missing_sample_ids(self):
        sample_mc = qiime2.CategoricalMetadataColumn(
            pd.Series(['g0', 'g2', 'g0', 'g2'], name='foo',
                      index=pd.Index(['s1', 's3', 's4', 's6'],
                                     name='sampleid')))
        data = np.array([
            [0, 1, 2, 3],
            [10, 11, 12, 13],
            [100, 110, 120, 130]])
        table = biom.Table(data, sample_ids=['s1', 's2', 's4', 's5'],
                           observation_ids=['x', 'y', 'z'])

        with self.assertRaisesRegex(ValueError, 'not present.*s2.*s5') as e:
            group(table, axis='sample', metadata=sample_mc, mode='sum')

        self.assertIn('s2', str(e.exception))
        self.assertIn('s5', str(e.exception))

    def test_reorder(self):
        sample_mc = qiime2.CategoricalMetadataColumn(
            pd.Series(['c', 'b', 'a'], name='foo',
                      index=pd.Index(['c', 'b', 'a'], name='sampleid')))

        data = np.array([[1, 2, 3], [30, 20, 10]])
        table = biom.Table(data, sample_ids=['a', 'b', 'c'],
                           observation_ids=['x', 'y'])

        expected = biom.Table(np.array([[3, 2, 1], [10, 20, 30]]),
                              sample_ids=['c', 'b', 'a'],
                              observation_ids=['x', 'y'])
        result = group(table, axis='sample', metadata=sample_mc, mode='sum')
        self.assertEqual(expected, result)

    def test_group_to_single_id(self):
        sample_mc = qiime2.CategoricalMetadataColumn(
            pd.Series(['all_samples', 'all_samples', 'all_samples'],
                      name='foo', index=pd.Index(['a', 'b', 'c'],
                                                 name='sampleid')))

        data = np.array([[1, 2, 3], [30, 20, 10]])
        table = biom.Table(data, sample_ids=['a', 'b', 'c'],
                           observation_ids=['x', 'y'])

        expected = biom.Table(np.array([[6], [60]]),
                              sample_ids=['all_samples'],
                              observation_ids=['x', 'y'])
        result = group(table, axis='sample', metadata=sample_mc, mode='sum')
        self.assertEqual(expected, result)

    def test_empty_metadata_values(self):
        # Trusting that the code is sane enough to not invent a distinction
        # between feature and sample metadata where there is none
        sample_mc = qiime2.CategoricalMetadataColumn(
            pd.Series(['a_new', 'a_new', None], name='foo',
                      index=pd.Index(['a', 'b', 'c'], name='sampleid')))
        sample_ids = sample_mc.to_series().index

        data = np.array([[1, 2, 3], [30, 20, 10]])
        table = biom.Table(data, sample_ids=sample_ids,
                           observation_ids=['x', 'y'])

        with self.assertRaisesRegex(ValueError, "missing.*value.*'c'"):
            group(table, axis='sample', metadata=sample_mc, mode='sum')

        nan_mc = qiime2.CategoricalMetadataColumn(
            pd.Series(['a_new', float('nan'), 'a_new'], name='foo',
                      index=pd.Index(['a', 'b', 'c'], name='id')))

        with self.assertRaisesRegex(ValueError, "missing.*value.*'b'"):
            group(table, axis='sample', metadata=nan_mc, mode='sum')

    def test_empty_only_in_superset(self):
        # Trusting that the code is sane enough to not invent a distinction
        # between feature and sample metadata where there is none
        sample_mc = qiime2.CategoricalMetadataColumn(
            pd.Series(['a_new', 'a_new', 'b_new', None], name='foo',
                      index=pd.Index(['a', 'b', 'c', 'd'], name='sampleid')))

        data = np.array([[1, 2, 3], [30, 20, 10]])
        table = biom.Table(data, sample_ids=['a', 'b', 'c'],
                           observation_ids=['x', 'y'])
        expected = biom.Table(np.array([[2, 3], [25, 10]]),
                              sample_ids=['a_new', 'b_new'],
                              observation_ids=['x', 'y'])
        result = group(table, axis='sample', metadata=sample_mc,
                       mode='mean-ceiling')
        self.assertEqual(expected, result)

    def test_numeric_strings(self):
        data = np.array([[1, 2, 3], [30, 20, 10]])
        table = biom.Table(data, sample_ids=['a', 'b', 'c'],
                           observation_ids=['x', 'y'])

        sample_mc = qiime2.CategoricalMetadataColumn(
            pd.Series(['-4.2', '-4.2', '-4.2'], name='foo',
                      index=pd.Index(['a', 'b', 'c'], name='sampleid')))

        expected = biom.Table(np.array([[6], [60]]),
                              sample_ids=['-4.2'],
                              observation_ids=['x', 'y'])
        result = group(table, axis='sample', metadata=sample_mc, mode='sum')
        self.assertEqual(expected, result)

    def test_empty_table(self):
        mc = qiime2.CategoricalMetadataColumn(
            pd.Series(['a_new', 'b_new'], name='foo',
                      index=pd.Index(['a', 'b'], name='id')))

        table = biom.Table(np.array([[]]), sample_ids=[], observation_ids=[])

        with self.assertRaisesRegex(ValueError, 'empty table'):
            group(table, axis='sample', metadata=mc, mode='sum')

        with self.assertRaisesRegex(ValueError, 'empty table'):
            group(table, axis='feature', metadata=mc, mode='sum')

    def _shared_setup(self):
        sample_mc = qiime2.CategoricalMetadataColumn(
            pd.Series(['treatment', 'treatment', 'control', 'other',
                       'control', 'other', 'other'], name='foo',
                      index=pd.Index(['a', 'b', 'c', 'd', 'e', 'f', 'g'],
                                     name='sampleid')))

        feature_mc = qiime2.CategoricalMetadataColumn(
            pd.Series(['g0', 'g1', 'g1', 'g1', 'g0'], name='foo',
                      index=pd.Index(['v', 'w', 'x', 'y', 'z'],
                                     name='featureid')))

        data = np.array([
            # t  t   c   o    c     o    o
            # a  b   c   d    e     f    g
            [0,  0,  0,  0,   1,    0,   2],    # v  g0
            [10, 10, 10, 10,  10,   100, 1],    # w  g1
            [12, 3,  14, 0,   0,    3,   34],   # x  g1
            [1,  1,  1,  1,   1,    1,   1],    # y  g1
            [0,  1,  11, 111, 1111, 20,  20]])  # z  g0

        table = biom.Table(data, sample_ids=sample_mc.to_series().index,
                           observation_ids=feature_mc.to_series().index)

        return sample_mc, feature_mc, table

    def test_feature_sum(self):
        sample_mc, feature_mc, table = self._shared_setup()
        expected = biom.Table(
            np.array([[0, 1, 11, 111, 1112, 20, 22],
                      [23, 14, 25, 11, 11, 104, 36]]),
            sample_ids=sample_mc.to_series().index,
            observation_ids=['g0', 'g1'])

        result = group(table, axis='feature', metadata=feature_mc, mode='sum')
        self.assertEqual(expected, result)

    def test_feature_mean_ceiling(self):
        sample_mc, feature_mc, table = self._shared_setup()
        expected = biom.Table(
            np.array([[0, 1, 6, 56, 556, 10, 11],
                      [8, 5, 9, 4, 4, 35, 12]]),
            sample_ids=sample_mc.to_series().index,
            observation_ids=['g0', 'g1'])

        result = group(table, axis='feature', metadata=feature_mc,
                       mode='mean-ceiling')
        self.assertEqual(expected, result)

    def test_feature_median_ceiling(self):
        sample_mc, feature_mc, table = self._shared_setup()
        expected = biom.Table(
            np.array([[0, 1, 6, 56, 556, 10, 11],
                      [10, 3, 10, 1, 1, 3, 1]]),
            sample_ids=sample_mc.to_series().index,
            observation_ids=['g0', 'g1'])

        result = group(table, axis='feature', metadata=feature_mc,
                       mode='median-ceiling')
        self.assertEqual(expected, result)

    def test_sample_sum(self):
        sample_mc, feature_mc, table = self._shared_setup()
        expected = biom.Table(
            np.array([[0, 1, 2],
                      [20, 20, 111],
                      [15, 14, 37],
                      [2, 2, 3],
                      [1, 1122, 151]]),
            sample_ids=['treatment', 'control', 'other'],
            observation_ids=feature_mc.to_series().index)

        result = group(table, axis='sample', metadata=sample_mc, mode='sum')
        self.assertEqual(expected, result)

    def test_sample_mean_ceiling(self):
        sample_mc, feature_mc, table = self._shared_setup()
        expected = biom.Table(
            np.array([[0, 1, 1],
                      [10, 10, 37],
                      [8, 7, 13],
                      [1, 1, 1],
                      [1, 561, 51]]),
            sample_ids=['treatment', 'control', 'other'],
            observation_ids=feature_mc.to_series().index)

        result = group(table, axis='sample', metadata=sample_mc,
                       mode='mean-ceiling')
        self.assertEqual(expected, result)

    def test_sample_median_ceiling(self):
        sample_mc, feature_mc, table = self._shared_setup()
        expected = biom.Table(
            np.array([[0, 1, 0],
                      [10, 10, 10],
                      [8, 7, 3],
                      [1, 1, 1],
                      [1, 561, 20]]),
            sample_ids=['treatment', 'control', 'other'],
            observation_ids=feature_mc.to_series().index)

        result = group(table, axis='sample', metadata=sample_mc,
                       mode='median-ceiling')
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
