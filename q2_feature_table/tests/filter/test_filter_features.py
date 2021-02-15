# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
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

from q2_feature_table import filter_features


class FilterFeaturesTests(unittest.TestCase):
    """ These tests are minimal relative to FilterSamplesTests, since the
        two functions being tested using the same private function under the
        hood. These tests cover the two places where the axis parameter is
        passed, to ensure that the tests work on the 'observation' axis as
        well as the 'sample' axis.
    """

    def test_min_frequency(self):
        # no filtering
        table = Table(np.array([[0, 1, 1], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_features(table, min_frequency=2)
        expected = Table(np.array([[0, 1, 1], [1, 1, 2]]),
                         ['O1', 'O2'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter one
        table = Table(np.array([[0, 1, 1], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_features(table, min_frequency=3)
        expected = Table(np.array([[1, 1, 2]]),
                         ['O2'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter all
        table = Table(np.array([[0, 1, 1], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_features(table, min_frequency=5)
        expected = Table(np.array([]), [], [])
        self.assertEqual(actual, expected)

    def test_filter_empty_samples(self):
        # no filtering
        table = Table(np.array([[0, 1, 1], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_features(table, min_frequency=2,
                                 filter_empty_samples=False)
        expected = Table(np.array([[0, 1, 1], [1, 1, 2]]),
                         ['O1', 'O2'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter all
        table = Table(np.array([[0, 1, 1], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_features(table, min_frequency=5,
                                 filter_empty_samples=False)
        expected = Table(np.empty((0, 3)), [], ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

    def test_feature_metadata(self):
        # no filtering
        df = pd.DataFrame({'SequencedGenome': ['yes', 'yes']},
                          index=pd.Index(['O1', 'O2'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_features(table, metadata=metadata)
        expected = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                         ['O1', 'O2'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter one
        df = pd.DataFrame({'SequencedGenome': ['yes']},
                          index=pd.Index(['O1'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_features(table, metadata=metadata)
        expected = Table(np.array([[1, 3]]),
                         ['O1'],
                         ['S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter all
        df = pd.DataFrame({}, index=pd.Index(['foo'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_features(table, metadata=metadata)
        expected = Table(np.array([]), [], [])
        self.assertEqual(actual, expected)

        # exclude one
        df = pd.DataFrame({'SequencedGenome': ['yes']},
                          index=pd.Index(['O1'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_features(table, metadata=metadata,
                                 exclude_ids=True)
        expected = Table(np.array([[1, 1, 2]]),
                         ['O2'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

        # exclude all
        df = pd.DataFrame({'SequencedGenome': ['yes', 'yes']},
                          index=pd.Index(['O1', 'O2'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_features(table, metadata=metadata,
                                 exclude_ids=True)
        expected = Table(np.array([]), [], [])
        self.assertEqual(actual, expected)

    def test_where(self):
        # no filtering
        df = pd.DataFrame({'SequencedGenome': ['yes', 'no']},
                          index=pd.Index(['O1', 'O2'], name='feature-id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        where = "SequencedGenome='yes' OR SequencedGenome='no'"
        actual = filter_features(table, metadata=metadata, where=where)
        expected = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                         ['O1', 'O2'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter one
        df = pd.DataFrame({'SequencedGenome': ['yes', 'no']},
                          index=pd.Index(['O1', 'O2'], name='feature-id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        where = "SequencedGenome='yes'"
        actual = filter_features(table, metadata=metadata, where=where)
        expected = Table(np.array([[1, 3]]),
                         ['O1'],
                         ['S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter all
        df = pd.DataFrame({'SequencedGenome': ['yes', 'no']},
                          index=pd.Index(['O1', 'O2'], name='feature-id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        where = "SequencedGenome='yes' AND SequencedGenome='no'"
        actual = filter_features(table, metadata=metadata, where=where)
        expected = Table(np.array([]), [], [])
        self.assertEqual(actual, expected)

        # filter one -> exclude one
        df = pd.DataFrame({'SequencedGenome': ['yes', 'no']},
                          index=pd.Index(['O1', 'O2'], name='feature-id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        where = "SequencedGenome='yes'"
        actual = filter_features(table,
                                 exclude_ids=True,
                                 metadata=metadata,
                                 where=where)
        expected = Table(np.array([[1, 1, 2]]),
                         ['O2'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
