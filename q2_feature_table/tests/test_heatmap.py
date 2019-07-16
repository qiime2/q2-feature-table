# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import tempfile
import os.path

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
import qiime2

from q2_feature_table import heatmap
from q2_feature_table._heatmap._visualizer import (
    _munge_metadata, _munge_feature_metadata)


class TestHeatmap(unittest.TestCase):
    def setUp(self):
        self.table = pd.DataFrame(data=[[0, 10], [10, 12], [10, 11]],
                                  columns=['O1', 'O2'],
                                  index=['S1', 'S2', 'S3'])
        self.output_dir_obj = tempfile.TemporaryDirectory(
                prefix='q2-feature-table-test-temp-')
        self.output_dir = self.output_dir_obj.name

    def tearDown(self):
        self.output_dir_obj.cleanup()

    def assertBasicVizValidity(self, viz_dir, normalize=True):
        index_fp = os.path.join(viz_dir, 'index.html')
        self.assertTrue(os.path.exists(index_fp))

        with open(index_fp) as fh:
            index_html = fh.read()

        normalize_str = '(normalized)' if normalize else '(not normalized)'
        self.assertTrue(normalize_str in index_html)

        for ext in ['png', 'svg']:
            fp = os.path.join(viz_dir, 'feature-table-heatmap.%s' % ext)
            self.assertTrue(os.path.exists(fp))

    def test_defaults(self):
        heatmap(self.output_dir, self.table)

        self.assertBasicVizValidity(self.output_dir)

    def test_with_title(self):
        heatmap(self.output_dir, self.table, title='foo')

        self.assertBasicVizValidity(self.output_dir)

    def test_with_metadata(self):
        md = qiime2.CategoricalMetadataColumn(
            pd.Series(['milo', 'summer', 'russ'], name='pet',
                      index=pd.Index(['S1', 'S2', 'S3'], name='id')))
        heatmap(self.output_dir, self.table, sample_metadata=md)

        self.assertBasicVizValidity(self.output_dir)

    def test_with_feature_metadata(self):
        feature_md = qiime2.CategoricalMetadataColumn(
            pd.Series(['peanut', 'dog'], name='species',
                      index=pd.Index(['O1', 'O2'], name='id')))
        heatmap(self.output_dir, self.table, feature_metadata=feature_md)

        self.assertBasicVizValidity(self.output_dir)

    def test_with_sample_and_feature_metadata(self):
        md = qiime2.CategoricalMetadataColumn(
            pd.Series(['milo', 'summer', 'russ'], name='pet',
                      index=pd.Index(['S1', 'S2', 'S3'], name='id')))
        feature_md = qiime2.CategoricalMetadataColumn(
            pd.Series(['peanut', 'dog'], name='species',
                      index=pd.Index(['O1', 'O2'], name='id')))
        heatmap(self.output_dir, self.table, sample_metadata=md,
                feature_metadata=feature_md)

        self.assertBasicVizValidity(self.output_dir)

    def test_empty_table(self):
        empty_table = pd.DataFrame([], [], [])

        with self.assertRaisesRegex(ValueError, 'empty'):
            heatmap(self.output_dir, empty_table)

    def test_table_ids_are_subset_of_metadata_ids(self):
        md = qiime2.CategoricalMetadataColumn(
            pd.Series(['milo', 'russ'], name='pet',
                      index=pd.Index(['S1', 'S3'], name='id')))

        with self.assertRaisesRegex(ValueError, 'not present.*S2'):
            heatmap(self.output_dir, self.table, sample_metadata=md)

    def test_extra_metadata_ids(self):
        md = qiime2.CategoricalMetadataColumn(
            pd.Series(['milo', 'summer', 'russ', 'peanut'], name='pet',
                      index=pd.Index(['S1', 'S2', 'S3', 'S4'], name='id')))

        heatmap(self.output_dir, self.table, sample_metadata=md)

        self.assertBasicVizValidity(self.output_dir)

    def test_no_normalization(self):
        heatmap(self.output_dir, self.table, normalize=False)

        self.assertBasicVizValidity(self.output_dir, normalize=False)

    def test_no_sample_cluster(self):
        md = qiime2.CategoricalMetadataColumn(
            pd.Series(['milo', 'summer', 'russ'], name='pet',
                      index=pd.Index(['S1', 'S2', 'S3'], name='id')))

        heatmap(self.output_dir, self.table, sample_metadata=md,
                cluster='features')

        self.assertBasicVizValidity(self.output_dir)

    def test_no_cluster(self):
        md = qiime2.CategoricalMetadataColumn(
            pd.Series(['milo', 'summer', 'russ'], name='pet',
                      index=pd.Index(['S1', 'S2', 'S3'], name='id')))

        heatmap(self.output_dir, self.table, sample_metadata=md,
                cluster='none')

        self.assertBasicVizValidity(self.output_dir)


class TestPrivateHelpers(unittest.TestCase):
    def setUp(self):
        self.table = pd.DataFrame(data=[[0, 10], [10, 12], [10, 11]],
                                  columns=['O1', 'O2'],
                                  index=['S1', 'S2', 'S3'])

    def test_munge_metadata_simple(self):
        md = qiime2.CategoricalMetadataColumn(
            pd.Series(['milo', 'russ', 'russ'], name='pet',
                      index=pd.Index(['S1', 'S2', 'S3'], name='id')))
        obs = _munge_metadata(md, self.table, 'both')

        exp_idx = pd.Index(['milo | S1', 'russ | S2', 'russ | S3'],
                           name='pet | id')
        exp = pd.DataFrame([[0, 10], [10, 12], [10, 11]], columns=['O1', 'O2'],
                           index=exp_idx)
        assert_frame_equal(exp, obs)

    def test_munge_metadata_ids_different_order(self):
        md = qiime2.CategoricalMetadataColumn(
            pd.Series(['russ', 'milo', 'russ'], name='pet',
                      index=pd.Index(['S2', 'S1', 'S3'], name='id')))
        obs = _munge_metadata(md, self.table, 'both')

        exp_idx = pd.Index(['milo | S1', 'russ | S2', 'russ | S3'],
                           name='pet | id')
        exp = pd.DataFrame([[0, 10], [10, 12], [10, 11]], columns=['O1', 'O2'],
                           index=exp_idx)
        assert_frame_equal(exp, obs)

    def test_munge_metadata_missing_samples(self):
        md = qiime2.CategoricalMetadataColumn(
            pd.Series(['milo', 'russ'], name='pet',
                      index=pd.Index(['S1', 'S3'], name='id')))
        with self.assertRaisesRegex(ValueError, 'not present.*S2'):
            _munge_metadata(md, self.table, 'both')

    def test_munge_metadata_empty_values(self):
        md = qiime2.CategoricalMetadataColumn(
            pd.Series([None, 'russ', np.nan], name='pet',
                      index=pd.Index(['S1', 'S2', 'S3'], name='id')))
        obs = _munge_metadata(md, self.table, 'both')

        exp_idx = pd.Index(['[No Value] | S1', 'russ | S2', '[No Value] | S3'],
                           name='pet | id')
        exp = pd.DataFrame([[0, 10], [10, 12], [10, 11]], columns=['O1', 'O2'],
                           index=exp_idx)
        assert_frame_equal(exp, obs)

    def test_munge_metadata_sort_samples(self):
        md = qiime2.CategoricalMetadataColumn(
            pd.Series(['peanut', 'milo', 'russ'], name='pet',
                      index=pd.Index(['S1', 'S2', 'S3'], name='id')))
        obs = _munge_metadata(md, self.table, 'features')

        exp_idx = pd.Index(['milo | S2', 'peanut | S1', 'russ | S3'],
                           name='pet | id')
        exp = pd.DataFrame([[10, 12], [0, 10], [10, 11]], columns=['O1', 'O2'],
                           index=exp_idx)
        assert_frame_equal(exp, obs)

    def test_munge_feature_metadata_simple(self):
        feature_md = qiime2.CategoricalMetadataColumn(
            pd.Series(['peanut', 'dog'], name='species',
                      index=pd.Index(['O1', 'O2'], name='id')))
        obs = _munge_feature_metadata(feature_md, self.table, 'both')

        exp = pd.DataFrame(
            [[0, 10], [10, 12], [10, 11]], columns=['peanut', 'dog'],
            index=pd.Index(['S1', 'S2', 'S3']))
        assert_frame_equal(exp, obs)

    def test_munge_feature_metadata_ids_different_order(self):
        feature_md = qiime2.CategoricalMetadataColumn(
            pd.Series(['dog', 'peanut'], name='species',
                      index=pd.Index(['O2', 'O1'], name='id')))
        obs = _munge_feature_metadata(feature_md, self.table, 'both')

        exp = pd.DataFrame(
            [[0, 10], [10, 12], [10, 11]], columns=['peanut', 'dog'],
            index=pd.Index(['S1', 'S2', 'S3']))
        assert_frame_equal(exp, obs)

    def test_munge_feature_metadata_missing_features(self):
        feature_md = qiime2.CategoricalMetadataColumn(
            pd.Series(
                ['dog'], name='species', index=pd.Index(['O2'], name='id')))
        with self.assertRaisesRegex(ValueError, 'not present.*O1'):
            _munge_feature_metadata(feature_md, self.table, 'both')

    def test_munge_feature_metadata_is_superset(self):
        feature_md = qiime2.CategoricalMetadataColumn(
            pd.Series(['peanut', 'dog', 'cujo'], name='species',
                      index=pd.Index(['O1', 'O2', 'O3'], name='id')))
        obs = _munge_feature_metadata(feature_md, self.table, 'both')

        exp = pd.DataFrame(
            [[0, 10], [10, 12], [10, 11]], columns=['peanut', 'dog'],
            index=pd.Index(['S1', 'S2', 'S3']))
        assert_frame_equal(exp, obs)

    def test_munge_feature_metadata_sort_samples(self):
        feature_md = qiime2.CategoricalMetadataColumn(
            pd.Series(['peanut', 'dog'], name='species',
                      index=pd.Index(['O1', 'O2'], name='id')))
        obs = _munge_feature_metadata(feature_md, self.table, 'samples')

        exp = pd.DataFrame(
            [[10, 0], [12, 10], [11, 10]], columns=['dog', 'peanut'],
            index=pd.Index(['S1', 'S2', 'S3']))
        assert_frame_equal(exp, obs)


if __name__ == "__main__":
    unittest.main()
