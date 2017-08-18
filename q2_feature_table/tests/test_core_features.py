# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import unittest
import tempfile
import os.path

import biom
import numpy as np
import pandas as pd
import pandas.testing as pdt

from q2_feature_table import core_features
from q2_feature_table._core_features._visualizer import (
    _get_core_features, _seven_number_summary, _round_fractions)


class TestCoreFeatures(unittest.TestCase):
    def setUp(self):
        self.table = biom.Table(np.array([[0, 11, 11], [13, 11, 11]]),
                                ['O1', 'O2'], ['S1', 'S2', 'S3'])
        self.output_dir_obj = tempfile.TemporaryDirectory(
                prefix='q2-feature-table-test-temp-')
        self.output_dir = self.output_dir_obj.name

    def tearDown(self):
        self.output_dir_obj.cleanup()

    def assertBasicVizValidity(self, viz_dir, exp_sample_count=3,
                               exp_feature_count=2, exp_no_core=False):
        index_fp = os.path.join(viz_dir, 'index.html')
        self.assertTrue(os.path.exists(index_fp))
        with open(index_fp, 'r') as fh:
            index_contents = fh.read()

        self.assertIn('Sample count: %d' % exp_sample_count, index_contents)
        self.assertIn('Feature count: %d' % exp_feature_count, index_contents)
        if exp_no_core:
            self.assertIn('No core features', index_contents)

        svg_fp = os.path.join(viz_dir, 'core-feature-counts.svg')
        self.assertTrue(os.path.exists(svg_fp))

    def assertCoreFeaturesPresent(self, filepath, positive_ids, negative_ids):
        self.assertTrue(os.path.exists(filepath))

        with open(filepath, 'r') as fh:
            contents = fh.read()

        for positive_id in positive_ids:
            self.assertIn(positive_id, contents)

        for negative_id in negative_ids:
            self.assertNotIn(negative_id, contents)

    def test_defaults(self):
        core_features(self.output_dir, self.table)

        self.assertBasicVizValidity(self.output_dir)

        core_50_fp = os.path.join(self.output_dir, 'core-features-0.500.tsv')
        self.assertCoreFeaturesPresent(core_50_fp, ['O1', 'O2'], [])

        core_100_fp = os.path.join(self.output_dir, 'core-features-1.000.tsv')
        self.assertCoreFeaturesPresent(core_100_fp, ['O2'], ['O1'])

    def test_fraction_range(self):
        core_features(self.output_dir, self.table, min_fraction=0.55,
                      max_fraction=0.95, steps=9)

        self.assertBasicVizValidity(self.output_dir)

        core_50_fp = os.path.join(self.output_dir, 'core-features-0.500.tsv')
        self.assertFalse(os.path.exists(core_50_fp))

        core_55_fp = os.path.join(self.output_dir, 'core-features-0.550.tsv')
        self.assertCoreFeaturesPresent(core_55_fp, ['O1', 'O2'], [])

        core_95_fp = os.path.join(self.output_dir, 'core-features-0.950.tsv')
        self.assertCoreFeaturesPresent(core_95_fp, ['O2'], ['O1'])

        core_100_fp = os.path.join(self.output_dir, 'core-features-1.000.tsv')
        self.assertFalse(os.path.exists(core_100_fp))

    def test_minimum_number_of_steps(self):
        core_features(self.output_dir, self.table, min_fraction=0.55,
                      max_fraction=0.95, steps=2)

        self.assertBasicVizValidity(self.output_dir)

        core_55_fp = os.path.join(self.output_dir, 'core-features-0.550.tsv')
        core_95_fp = os.path.join(self.output_dir, 'core-features-0.950.tsv')
        tsv_files = sorted(glob.glob(os.path.join(self.output_dir, '*.tsv')))
        self.assertEqual(tsv_files, [core_55_fp, core_95_fp])

        self.assertCoreFeaturesPresent(core_55_fp, ['O1', 'O2'], [])
        self.assertCoreFeaturesPresent(core_95_fp, ['O2'], ['O1'])

    def test_equal_min_max_fractions(self):
        core_features(self.output_dir, self.table, min_fraction=0.55,
                      max_fraction=0.55)

        self.assertBasicVizValidity(self.output_dir)

        core_55_fp = os.path.join(self.output_dir, 'core-features-0.550.tsv')
        tsv_files = glob.glob(os.path.join(self.output_dir, '*.tsv'))
        self.assertEqual(tsv_files, [core_55_fp])

        self.assertCoreFeaturesPresent(core_55_fp, ['O1', 'O2'], [])

    def test_tiny_steps_precision(self):
        # Bug in a previous iteration of the code used a fixed number of digits
        # in the "feature list" TSV filenames, causing files to be overwritten
        # with the wrong data if a small enough step size was used.
        core_features(self.output_dir, self.table, min_fraction=0.1,
                      max_fraction=0.11, steps=101)

        self.assertBasicVizValidity(self.output_dir)

        fp = os.path.join(self.output_dir, 'core-features-0.1000.tsv')
        self.assertTrue(os.path.exists(fp))

        fp = os.path.join(self.output_dir, 'core-features-0.1001.tsv')
        self.assertTrue(os.path.exists(fp))

        fp = os.path.join(self.output_dir, 'core-features-0.1099.tsv')
        self.assertTrue(os.path.exists(fp))

        fp = os.path.join(self.output_dir, 'core-features-0.1100.tsv')
        self.assertTrue(os.path.exists(fp))

    def test_no_core_features(self):
        table = biom.Table(np.array([[0, 11, 11], [11, 11, 0]]), ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])

        core_features(self.output_dir, table)

        self.assertBasicVizValidity(self.output_dir, exp_no_core=True)

        core_50_fp = os.path.join(self.output_dir, 'core-features-0.500.tsv')
        self.assertCoreFeaturesPresent(core_50_fp, ['O1', 'O2'], [])

        # No core features exist at fraction=1.0
        core_100_fp = os.path.join(self.output_dir, 'core-features-1.000.tsv')
        self.assertFalse(os.path.exists(core_100_fp))

    def test_empty_table(self):
        table = biom.Table(np.array([[]]), [], [])

        core_features(self.output_dir, table)

        self.assertBasicVizValidity(self.output_dir, exp_sample_count=0,
                                    exp_feature_count=0, exp_no_core=True)

        self.assertFalse(glob.glob(os.path.join(self.output_dir, '*.tsv')))

    def test_invalid_parameters(self):
        with self.assertRaisesRegex(ValueError, 'fraction'):
            core_features(self.output_dir, self.table, min_fraction=0.75,
                          max_fraction=0.5)


class TestCoreFeaturesPrivateFunctions(unittest.TestCase):
    # Most of the work happens in private functions for this visualizer,
    # so it's important to have some tests of the private functionality

    def test_seven_number_summary(self):
        # validated against calling pd.Series.describe directly
        exp = pd.Series([11., 11., 11., 11., 12., 12.64, 12.92],
                        index=['2%', '9%', '25%', '50%', '75%', '91%', '98%'])
        pdt.assert_series_equal(_seven_number_summary([13, 11, 11]), exp)

        exp = pd.Series([1., 1., 1., 1., 1., 1., 1.],
                        index=['2%', '9%', '25%', '50%', '75%', '91%', '98%'])
        pdt.assert_series_equal(_seven_number_summary([1]), exp)

        exp = pd.Series([np.nan, np.nan, np.nan, np.nan, np.nan,
                         np.nan, np.nan],
                        index=['2%', '9%', '25%', '50%', '75%', '91%', '98%'])
        pdt.assert_series_equal(_seven_number_summary([]), exp)

    def test_get_core_features(self):
        table = biom.Table(np.array([[0, 11, 11], [13, 11, 11]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        o2_seven_num = pd.Series(
            [11., 11., 11., 11., 12., 12.64, 12.92],
            index=['2%', '9%', '25%', '50%', '75%', '91%', '98%'])
        exp = pd.DataFrame([o2_seven_num], index=['O2'])

        obs = _get_core_features(table, fraction=0.67)
        pdt.assert_frame_equal(obs, exp)

        obs = _get_core_features(table, fraction=1.0)
        pdt.assert_frame_equal(obs, exp)

    def test_get_core_features_all(self):
        table = biom.Table(np.array([[0, 11, 11], [13, 11, 11]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        o1_seven_num = pd.Series(
            [0.44, 1.98, 5.5, 11., 11., 11., 11.],
            index=['2%', '9%', '25%', '50%', '75%', '91%', '98%'])
        o2_seven_num = pd.Series(
            [11., 11., 11., 11., 12., 12.64, 12.92],
            index=['2%', '9%', '25%', '50%', '75%', '91%', '98%'])
        exp = pd.DataFrame([o1_seven_num, o2_seven_num], index=['O1', 'O2'])

        # fraction = 2/3 - 0.006
        obs = _get_core_features(table, fraction=0.660)
        pdt.assert_frame_equal(obs, exp)

        obs = _get_core_features(table, fraction=0.0)
        pdt.assert_frame_equal(obs, exp)

    def test_round_fractions_single(self):
        obs = _round_fractions([0.0])
        self.assertEqual(obs, ['0.000'])

        obs = _round_fractions([0.1236])
        self.assertEqual(obs, ['0.124'])

    def test_round_fractions_3_decimals(self):
        obs = _round_fractions([0.1236, 0.45, 0.123])
        self.assertEqual(obs, ['0.124', '0.450', '0.123'])

    def test_round_fractions_4_decimals(self):
        obs = _round_fractions([0.1234, 0.45, 0.0, 0.12346])
        self.assertEqual(obs, ['0.1234', '0.4500', '0.0000', '0.1235'])

    def test_round_fractions_full_precision(self):
        obs = _round_fractions([0.01234567890, 0.01234567891])
        self.assertEqual(obs, ['0.0123456789', '0.01234567891'])


if __name__ == "__main__":
    unittest.main()
