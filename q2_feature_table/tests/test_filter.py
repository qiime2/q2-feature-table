# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
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

from q2_feature_table import filter_samples, filter_features


class FilterSamplesTests(unittest.TestCase):

    def test_invalid_args(self):
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        with self.assertRaises(ValueError):
            filter_samples(table)

        with self.assertRaises(ValueError):
            filter_samples(table, where="Subject='subject-1'")

    def test_min_frequency(self):
        # no filtering
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, min_frequency=1)
        expected = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                         ['O1', 'O2'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter one
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, min_frequency=2)
        expected = Table(np.array([[1, 3], [1, 2]]),
                         ['O1', 'O2'],
                         ['S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter two
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, min_frequency=3)
        expected = Table(np.array([[3], [2]]),
                         ['O1', 'O2'],
                         ['S3'])
        self.assertEqual(actual, expected)

        # filter all
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, min_frequency=42)
        expected = Table(np.array([]), [], [])
        self.assertEqual(actual, expected)

    def test_max_frequency(self):
        # no filtering
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, max_frequency=42)
        expected = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                         ['O1', 'O2'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter one
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, max_frequency=4)
        expected = Table(np.array([[0, 1], [1, 1]]),
                         ['O1', 'O2'],
                         ['S1', 'S2'])
        self.assertEqual(actual, expected)

        # filter two
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, max_frequency=1)
        expected = Table(np.array([[1]]),
                         ['O2'],
                         ['S1'])
        self.assertEqual(actual, expected)

        # filter all
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, max_frequency=0)
        expected = Table(np.array([]), [], [])
        self.assertEqual(actual, expected)

    def test_min_features(self):
        # no filtering
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, min_features=1)
        expected = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                         ['O1', 'O2'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter one
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, min_features=2)
        expected = Table(np.array([[1, 3], [1, 2]]),
                         ['O1', 'O2'],
                         ['S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter all
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, min_features=3)
        expected = Table(np.array([]), [], [])
        self.assertEqual(actual, expected)

    def test_max_features(self):
        # no filtering
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, max_features=2)
        expected = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                         ['O1', 'O2'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter two
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, max_features=1)
        expected = Table(np.array([[1]]),
                         ['O2'],
                         ['S1'])
        self.assertEqual(actual, expected)

        # filter all
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, max_features=0)
        expected = Table(np.array([]), [], [])
        self.assertEqual(actual, expected)

    def test_sample_metadata(self):
        # no filtering
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=['S1', 'S2', 'S3'])
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, metadata=metadata)
        expected = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                         ['O1', 'O2'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter one
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-2'],
                           'SampleType': ['tongue', 'gut']},
                          index=['S2', 'S3'])
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, metadata=metadata)
        expected = Table(np.array([[1, 3], [1, 2]]),
                         ['O1', 'O2'],
                         ['S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter all
        df = pd.DataFrame({})
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, metadata=metadata)
        expected = Table(np.array([]), [], [])
        self.assertEqual(actual, expected)

    def test_sample_metadata_extra_ids(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=['S-not-in-table', 'S2', 'S3'])
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, metadata=metadata)
        expected = Table(np.array([[1, 3], [1, 2]]),
                         ['O1', 'O2'],
                         ['S2', 'S3'])
        self.assertEqual(actual, expected)

    def test_where(self):
        # no filtering
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='#SampleID'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        where = "Subject='subject-1' OR Subject='subject-2'"
        actual = filter_samples(table, metadata=metadata, where=where)
        expected = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                         ['O1', 'O2'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter one
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='#SampleID'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        where = "Subject='subject-1'"
        actual = filter_samples(table, metadata=metadata, where=where)
        expected = Table(np.array([[0, 1], [1, 1]]),
                         ['O1', 'O2'],
                         ['S1', 'S2'])
        self.assertEqual(actual, expected)

        # filter two
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='#SampleID'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        where = "Subject='subject-1' AND SampleType='gut'"
        actual = filter_samples(table, metadata=metadata, where=where)
        expected = Table(np.array([[1]]),
                         ['O2'],
                         ['S1'])
        self.assertEqual(actual, expected)

        # filter all
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='#SampleID'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        where = "Subject='subject-1' AND Subject='subject-2'"
        actual = filter_samples(table, metadata=metadata, where=where)
        expected = Table(np.array([]), [], [])
        self.assertEqual(actual, expected)

    def test_combine_id_and_frequency_filters(self):
        # no filtering
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='#SampleID'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        where = "Subject='subject-1' OR Subject='subject-2'"
        actual = filter_samples(table, metadata=metadata, where=where,
                                min_frequency=1)
        expected = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                         ['O1', 'O2'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

        # id and frequency filters active
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='#SampleID'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        where = "Subject='subject-1'"
        actual = filter_samples(table, metadata=metadata, where=where,
                                min_frequency=2)
        expected = Table(np.array([[1], [1]]),
                         ['O1', 'O2'],
                         ['S2'])
        self.assertEqual(actual, expected)

    def test_combine_frequency_filters(self):
        # no filtering
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, min_frequency=1, max_frequency=5)
        expected = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                         ['O1', 'O2'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter two
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, min_frequency=2, max_frequency=2)
        expected = Table(np.array([[1], [1]]),
                         ['O1', 'O2'],
                         ['S2'])
        self.assertEqual(actual, expected)

        # filter two
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, max_frequency=2, min_features=2)
        expected = Table(np.array([[1], [1]]),
                         ['O1', 'O2'],
                         ['S2'])
        self.assertEqual(actual, expected)


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

    def test_feature_metadata(self):
        # no filtering
        df = pd.DataFrame({'SequencedGenome': ['yes', 'yes']},
                          index=['O1', 'O2'])
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
                          index=['O1'])
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
        df = pd.DataFrame({})
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_features(table, metadata=metadata)
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


if __name__ == "__main__":
    unittest.main()
