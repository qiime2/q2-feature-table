# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
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
from q2_feature_table._filter import check_relative_frequency


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

        # filter all raising ValueError
        table = Table(np.array([[0, 1, 1], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        with self.assertRaisesRegex(ValueError, 'table is empty'):
            filter_features(table, min_frequency=5)

        # filter all and allow empty table
        actual = filter_features(table, min_frequency=5,
                                 allow_empty_table=True)
        self.assertTrue(actual.is_empty())

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

        # filter all raising ValueError
        table = Table(np.array([[0, 1, 1], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        with self.assertRaisesRegex(ValueError, 'table is empty'):
            filter_features(table, min_frequency=5, filter_empty_samples=False)

        # filter all and allow empty table
        actual = filter_features(table, min_frequency=5,
                                 filter_empty_samples=False,
                                 allow_empty_table=True)
        self.assertTrue(actual.is_empty())

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

        # filter all and raise ValueError
        df = pd.DataFrame({}, index=pd.Index(['foo'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        with self.assertRaisesRegex(ValueError, 'table is empty'):
            filter_features(table, metadata=metadata)

        # filter all and allow empty table
        actual = filter_features(table, metadata=metadata,
                                 allow_empty_table=True)
        self.assertTrue(actual.is_empty())

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

        # exclude all and raise ValueError
        df = pd.DataFrame({'SequencedGenome': ['yes', 'yes']},
                          index=pd.Index(['O1', 'O2'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        with self.assertRaisesRegex(ValueError, 'table is empty'):
            filter_features(table, metadata=metadata, exclude_ids=True,
                            allow_empty_table=False)

        # exclude all and allow empty table
        actual = filter_features(table, metadata=metadata, exclude_ids=True,
                                 allow_empty_table=True)
        self.assertTrue(actual.is_empty())

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

        # filter all with ValueError
        df = pd.DataFrame({'SequencedGenome': ['yes', 'no']},
                          index=pd.Index(['O1', 'O2'], name='feature-id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        where = "SequencedGenome='yes' AND SequencedGenome='no'"

        with self.assertRaisesRegex(ValueError, 'table is empty'):
            filter_features(table, metadata=metadata, where=where,
                            allow_empty_table=False)

        # Filter all and allow empty table
        actual = filter_features(table, metadata=metadata, where=where,
                                 allow_empty_table=True)
        self.assertTrue(actual.is_empty())

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

    def test_check_relative_frequency_true(self):
        table = Table(np.array([[0.2, 0.5, 0.3],
                                [0.8, 0.5, 0.7],
                                [0.0, 0.0, 0.0]]),
                      ['O1', 'O2', 'O3'],
                      ['S1', 'S2', 'S3'])

        self.assertTrue(check_relative_frequency(table))

    def test_check_relative_frequency_false(self):
        table = Table(np.array([[0.2, 0.5, 0.3],
                                [0.3, 0.1, 0.6]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])

        self.assertFalse(check_relative_frequency(table))

    def test_relative_frequency(self):
        table = Table(np.array([[0.2, 0.5, 0.3, 0.0],
                                [0.3, 0.1, 0.6, 0.0],
                                [0.5, 0.4, 0.1, 0.0]]),
                      ['O1', 'O2', 'O3'],
                      ['S1', 'S2', 'S3', 'S4'])
        metadata = qiime2.Metadata(
            pd.DataFrame({'keep': ['yes', 'no', 'yes']},
                         index=pd.Index(['O1', 'O2', 'O3'], name='id')))
        actual = filter_features(table, metadata=metadata,
                                 where="keep='yes'",
                                 filter_empty_samples=False)

        self.assertEqual(actual.shape, (2, 4))
        self.assertEqual(set(actual.ids(axis='sample')),
                         set(['S1', 'S2', 'S3', 'S4']))
        self.assertEqual(set(actual.ids(axis='observation')),
                         set(['O1', 'O3']))
        np.testing.assert_allclose(
            actual.sum(axis='sample'), np.array([1., 1., 1., 0.])
        )
        np.testing.assert_allclose(
            actual.matrix_data.toarray(), np.array([[2/7, 5/9, 3/4, 0.],
                                                    [5/7, 4/9, 1/4, 0.]])
        )

    def test_relative_frequency_ids(self):
        table = Table(np.array([[0.2, 0.5, 0.3, 0.0],
                                [0.3, 0.1, 0.6, 0.0],
                                [0.5, 0.4, 0.1, 0.0]]),
                      ['O1', 'O2', 'O3'], ['S1', 'S2', 'S3', 'S4'])
        actual = filter_features(table, ids=['O1', 'O3'],
                                 filter_empty_samples=False)

        self.assertEqual(set(actual.ids(axis='sample')),
                         set(['S1', 'S2', 'S3', 'S4']))
        self.assertEqual(set(actual.ids(axis='observation')),
                         set(['O1', 'O3']))
        np.testing.assert_allclose(
            actual.sum(axis='sample'), np.array([1., 1., 1., 0.])
        )
        np.testing.assert_allclose(
            actual.matrix_data.toarray(), np.array([[2/7, 5/9, 3/4, 0.],
                                                    [5/7, 4/9, 1/4, 0.]])
        )

    def test_feature_ids(self):
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'], ['S1', 'S2', 'S3'])
        actual = filter_features(table, ids=['O1'])
        expected = Table(np.array([[1, 3]]), ['O1'], ['S2', 'S3'])
        self.assertEqual(actual, expected)

        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'], ['S1', 'S2', 'S3'])
        actual = filter_features(table, ids=['O1'], exclude_ids=True)
        expected = Table(np.array([[1, 1, 2]]), ['O2'], ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

    def test_feature_ids_and_metadata_are_combined(self):
        table = Table(np.array([[0, 1, 3], [1, 1, 2], [4, 0, 0]]),
                      ['O1', 'O2', 'O3'], ['S1', 'S2', 'S3'])
        metadata = qiime2.Metadata(pd.DataFrame(
            {'keep': ['yes', 'no']},
            index=pd.Index(['O2', 'O3'], name='id')))
        actual = filter_features(table, ids=['O1'], metadata=metadata,
                                 where="keep='yes'")
        expected = Table(np.array([[0, 1, 3], [1, 1, 2]]), ['O1', 'O2'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

        table = Table(np.array([[0, 1, 3], [1, 1, 2], [4, 0, 0]]),
                      ['O1', 'O2', 'O3'], ['S1', 'S2', 'S3'])
        actual = filter_features(table, ids=['O1'], metadata=metadata,
                                 where="keep='yes'", exclude_ids=True)
        expected = Table(np.array([[4]]), ['O3'], ['S1'])
        self.assertEqual(actual, expected)

    def test_feature_ids_with_sample_filter(self):
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'], ['S1', 'S2', 'S3'])
        actual = filter_features(table, ids=['O1', 'O2'], min_samples=3)
        expected = Table(np.array([[1, 1, 2]]), ['O2'], ['S1', 'S2', 'S3'])

        self.assertEqual(actual, expected)

    def test_feature_ids_without_empty_sample_filtering(self):
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'], ['S1', 'S2', 'S3'])
        actual = filter_features(table, ids=['O1'],
                                 filter_empty_samples=False)
        expected = Table(np.array([[0, 1, 3]]), ['O1'], ['S1', 'S2', 'S3'])

        self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
