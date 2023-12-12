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

from q2_feature_table import filter_samples


class FilterSamplesTests(unittest.TestCase):

    def test_invalid_args(self):
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        with self.assertRaisesRegex(ValueError, "No filtering"):
            filter_samples(table)

        with self.assertRaisesRegex(ValueError,
                                    "'where' is specified."):
            filter_samples(table, where="Subject='subject-1'")

        with self.assertRaisesRegex(ValueError,
                                    "'exclude_ids' is True."):
            filter_samples(table, exclude_ids=True)

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

        # filter all raising ValueError
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        with self.assertRaisesRegex(ValueError, 'table is empty'):
            filter_samples(table, min_frequency=42)

        # filter all and allow empty table
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]), ['O1-alt', 'O2-alt'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, min_frequency=42,
                                allow_empty_table=True)
        self.assertTrue(actual.is_empty())

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

        # filter all raising ValueError
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        with self.assertRaisesRegex(ValueError, 'table is empty'):
            filter_samples(table, max_frequency=0)

        # filter all and allow empty table
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]), ['O1-alt', 'O2-alt'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, max_frequency=0, allow_empty_table=True)
        self.assertTrue(actual.is_empty())

    def test_filter_empty_features(self):
        # no filtering
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, max_frequency=42,
                                filter_empty_features=False)
        expected = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                         ['O1', 'O2'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter one
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, max_frequency=4,
                                filter_empty_features=False)
        expected = Table(np.array([[0, 1], [1, 1]]),
                         ['O1', 'O2'],
                         ['S1', 'S2'])
        self.assertEqual(actual, expected)

        # filter two
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, max_frequency=1,
                                filter_empty_features=False)
        expected = Table(np.array([[0], [1]]),
                         ['O1', 'O2'],
                         ['S1'])
        self.assertEqual(actual, expected)

        # filter all raising ValueError
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        with self.assertRaisesRegex(ValueError, 'table is empty'):
            filter_samples(table, max_frequency=0,
                           filter_empty_features=False)

        # filter all and allow empty table
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]), ['O1-alt', 'O2-alt'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, max_frequency=0,
                                filter_empty_features=False,
                                allow_empty_table=True)
        self.assertTrue(actual.is_empty())

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

        # filter all raising ValueError
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        with self.assertRaisesRegex(ValueError, 'table is empty'):
            filter_samples(table, min_features=3)

        # filter all and allow empty table
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]), ['O1-alt', 'O2-alt'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, min_features=3, allow_empty_table=True)
        self.assertTrue(actual.is_empty())

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

        # filter all raising ValueError
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        with self.assertRaisesRegex(ValueError, 'table is empty'):
            filter_samples(table, max_features=0)

        # filter all and allow empty table
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1-alt', 'O2-alt'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, max_features=0, allow_empty_table=True)
        self.assertTrue(actual.is_empty())

    def test_sample_metadata(self):
        # no filtering
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
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
                          index=pd.Index(['S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, metadata=metadata)
        expected = Table(np.array([[1, 3], [1, 2]]),
                         ['O1', 'O2'],
                         ['S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter all raising ValueError
        df = pd.DataFrame({}, index=pd.Index(['foo'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        with self.assertRaisesRegex(ValueError, 'table is empty'):
            filter_samples(table, metadata=metadata)

        # filter all and allow empty table
        df = pd.DataFrame({}, index=pd.Index(['foo'], name='id'))
        metadata = qiime2.Metadata(df)
        actual = filter_samples(table, metadata=metadata,
                                allow_empty_table=True)
        self.assertTrue(actual.is_empty())

        # exclude none
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['gut']},
                          index=pd.Index(['S90'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, metadata=metadata, exclude_ids=True)
        self.assertEqual(actual, table)

        # exclude one
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['gut']},
                          index=pd.Index(['S1'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, metadata=metadata, exclude_ids=True)
        expected = Table(np.array([[1, 3], [1, 2]]),
                         ['O1', 'O2'],
                         ['S2', 'S3'])
        self.assertEqual(actual, expected)

        # exclude two
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1'],
                           'SampleType': ['gut', 'tongue']},
                          index=pd.Index(['S1', 'S2'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table, metadata=metadata, exclude_ids=True)
        expected = Table(np.array([[3], [2]]),
                         ['O1', 'O2'],
                         ['S3'])
        self.assertEqual(actual, expected)

        # exclude all raising ValueError
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)
        with self.assertRaisesRegex(ValueError, 'table is empty'):
            filter_samples(table, metadata=metadata,
                           exclude_ids=True)

        # exclude all and allow empty table
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)
        actual = filter_samples(table, metadata=metadata,
                                exclude_ids=True,
                                allow_empty_table=True)
        self.assertTrue(actual.is_empty())

    def test_sample_metadata_extra_ids(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S-not-in-table', 'S2', 'S3'],
                                         name='id'))
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

        # filter all raising ValueError
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='#SampleID'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        where = "Subject='subject-1' AND Subject='subject-2'"
        with self.assertRaisesRegex(ValueError, 'table is empty'):
            filter_samples(table, metadata=metadata, where=where)

        # filter all allowing empty table
        actual = filter_samples(table, metadata=metadata, where=where,
                                allow_empty_table=True)
        self.assertTrue(actual.is_empty())

        # filter none -> exclude none
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='#SampleID'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        where = "Subject='subject-1' AND SampleType='elbow'"
        actual = filter_samples(table,
                                metadata=metadata,
                                where=where,
                                exclude_ids=True)
        self.assertEqual(actual, table)

        # filter one -> exclude one
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='#SampleID'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        where = "Subject='subject-1' AND SampleType='gut'"
        actual = filter_samples(table,
                                metadata=metadata,
                                where=where,
                                exclude_ids=True)
        expected = Table(np.array([[1, 3], [1, 2]]),
                         ['O1', 'O2'],
                         ['S2', 'S3'])
        self.assertEqual(actual, expected)

        # filter two -> exclude two
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='#SampleID'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        where = "Subject='subject-1'"
        actual = filter_samples(table,
                                metadata=metadata,
                                where=where,
                                exclude_ids=True)
        expected = Table(np.array([[3], [2]]),
                         ['O1', 'O2'],
                         ['S3'])
        self.assertEqual(actual, expected)

        # filter all -> exclude all raising ValueError
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='#SampleID'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        where = "Subject='subject-1' OR Subject='subject-2'"
        with self.assertRaisesRegex(ValueError, 'table is empty'):
            filter_samples(table,
                           metadata=metadata,
                           where=where,
                           exclude_ids=True)

        # exclude all and allow empty table
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='#SampleID'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        where = "Subject='subject-1' OR Subject='subject-2'"
        actual = filter_samples(table,
                                metadata=metadata,
                                where=where,
                                exclude_ids=True,
                                allow_empty_table=True)
        self.assertTrue(actual.is_empty())

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

    def test_combine_exclude_ids_and_frequency_filters(self):
        # exclude one, min_frequency filter none
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['gut']},
                          index=pd.Index(['S1'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table,
                                metadata=metadata,
                                exclude_ids=True,
                                min_frequency=2)
        expected = Table(np.array([[1, 3], [1, 2]]),
                         ['O1', 'O2'],
                         ['S2', 'S3'])
        self.assertEqual(actual, expected)

        # exclude one, min_frequency filter one
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['gut']},
                          index=pd.Index(['S1'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table,
                                metadata=metadata,
                                exclude_ids=True,
                                min_frequency=3)
        expected = Table(np.array([[3], [2]]),
                         ['O1', 'O2'],
                         ['S3'])
        self.assertEqual(actual, expected)

        # exclude one, min_frequency filter for same one
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['gut']},
                          index=pd.Index(['S1'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [0, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table,
                                metadata=metadata,
                                exclude_ids=True,
                                min_frequency=1)
        expected = Table(np.array([[1, 3], [1, 2]]),
                         ['O1', 'O2'],
                         ['S2', 'S3'])
        self.assertEqual(actual, expected)

        # exclude one, max_frequency filter none
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['gut']},
                          index=pd.Index(['S1'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table,
                                metadata=metadata,
                                exclude_ids=True,
                                max_frequency=42)
        expected = Table(np.array([[1, 3], [1, 2]]),
                         ['O1', 'O2'],
                         ['S2', 'S3'])
        self.assertEqual(actual, expected)

        # exclude one, max_frequency filter one
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['gut']},
                          index=pd.Index(['S1'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table,
                                metadata=metadata,
                                exclude_ids=True,
                                max_frequency=4)
        expected = Table(np.array([[1], [1]]),
                         ['O1', 'O2'],
                         ['S2'])
        self.assertEqual(actual, expected)

        # exclude one, max_frequency filter for same one
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['gut']},
                          index=pd.Index(['S1'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[5, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table,
                                metadata=metadata,
                                exclude_ids=True,
                                max_frequency=5)
        expected = Table(np.array([[1, 3], [1, 2]]),
                         ['O1', 'O2'],
                         ['S2', 'S3'])
        self.assertEqual(actual, expected)

        # exclude one, max & min_frequency filter none
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['gut']},
                          index=pd.Index(['S1'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table,
                                metadata=metadata,
                                exclude_ids=True,
                                max_frequency=42,
                                min_frequency=0)
        expected = Table(np.array([[1, 3], [1, 2]]),
                         ['O1', 'O2'],
                         ['S2', 'S3'])
        self.assertEqual(actual, expected)

        # exclude one, min_frequency filter one,
        # max_frequency filter one raising ValueError
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['gut']},
                          index=pd.Index(['S1'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        with self.assertRaisesRegex(ValueError, 'table is empty'):
            filter_samples(table,
                           metadata=metadata,
                           exclude_ids=True,
                           max_frequency=4,
                           min_frequency=3)

        # allow empty table
        actual = filter_samples(table,
                                metadata=metadata,
                                exclude_ids=True,
                                max_frequency=4,
                                min_frequency=3,
                                allow_empty_table=True)
        self.assertTrue(actual.is_empty())

        # where filter one -> exclude one,
        # min_frequency filter one,
        # max_frequency filter one raising ValueError
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue']},
                          index=pd.Index(['S1', 'S2'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        where = "Subject='subject-1'"
        with self.assertRaisesRegex(ValueError, 'table is empty'):
            filter_samples(table,
                           metadata=metadata,
                           exclude_ids=True,
                           where=where,
                           max_frequency=4,
                           min_frequency=3)

        # allow empty table
        actual = filter_samples(table,
                                metadata=metadata,
                                exclude_ids=True,
                                where=where,
                                max_frequency=4,
                                min_frequency=3,
                                allow_empty_table=True)
        self.assertTrue(actual.is_empty())

    def test_combine_exclude_ids_and_features_filters(self):
        # exclude one, min_features filter none
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['gut']},
                          index=pd.Index(['S1'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table,
                                metadata=metadata,
                                exclude_ids=True,
                                min_features=2)
        expected = Table(np.array([[1, 3], [1, 2]]),
                         ['O1', 'O2'],
                         ['S2', 'S3'])
        self.assertEqual(actual, expected)

        # exclude one, min_features filter one
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['gut']},
                          index=pd.Index(['S2'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table,
                                metadata=metadata,
                                exclude_ids=True,
                                min_features=2)
        expected = Table(np.array([[3], [2]]),
                         ['O1', 'O2'],
                         ['S3'])
        self.assertEqual(actual, expected)

        # exclude one, min_features filter for same one
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['gut']},
                          index=pd.Index(['S2'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[4, 1, 3], [6, 0, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table,
                                metadata=metadata,
                                exclude_ids=True,
                                min_features=1)
        expected = Table(np.array([[4, 3], [6, 2]]),
                         ['O1', 'O2'],
                         ['S1', 'S3'])
        self.assertEqual(actual, expected)

        # exclude one, max_features filter none
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['gut']},
                          index=pd.Index(['S1'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table,
                                metadata=metadata,
                                exclude_ids=True,
                                max_features=3)
        expected = Table(np.array([[1, 3], [1, 2]]),
                         ['O1', 'O2'],
                         ['S2', 'S3'])
        self.assertEqual(actual, expected)

        # exclude one, max_features filter one
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['gut']},
                          index=pd.Index(['S2'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table,
                                metadata=metadata,
                                exclude_ids=True,
                                max_features=1)
        expected = Table(np.array([[1]]), ['O2'], ['S1'])
        self.assertEqual(actual, expected)

        # exclude one, max_features filter for same one
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['gut']},
                          index=pd.Index(['S2'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [1, 10, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table,
                                metadata=metadata,
                                exclude_ids=True,
                                max_features=9)
        expected = Table(np.array([[0, 3], [1, 2]]),
                         ['O1', 'O2'],
                         ['S1', 'S3'])
        self.assertEqual(actual, expected)

        # exclude one, max_features filter none,
        # min_features filter none
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['gut']},
                          index=pd.Index(['S2'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [0, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_samples(table,
                                metadata=metadata,
                                exclude_ids=True,
                                min_features=0,
                                max_features=5)
        expected = Table(np.array([[0, 3], [0, 2]]),
                         ['O1', 'O2'],
                         ['S1', 'S3'])
        self.assertEqual(actual, expected)

        # exclude one, max_features filter one,
        # min_features filter one raising ValueError
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['gut']},
                          index=pd.Index(['S2'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [0, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        with self.assertRaisesRegex(ValueError, 'table is empty'):
            filter_samples(table,
                           metadata=metadata,
                           exclude_ids=True,
                           min_features=1,
                           max_features=1)
        # allow empty table
        actual = filter_samples(table,
                                metadata=metadata,
                                exclude_ids=True,
                                min_features=1,
                                max_features=1,
                                allow_empty_table=True)

        # where filter one -> exclude one,
        # max_features filter one,
        # min_features filter one raising ValueError
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue']},
                          index=pd.Index(['S1', 'S2'], name='id'))
        metadata = qiime2.Metadata(df)
        table = Table(np.array([[0, 1, 3], [0, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        where = "SampleType='tongue'"
        with self.assertRaisesRegex(ValueError, 'table is empty'):
            filter_samples(table,
                           metadata=metadata,
                           where=where,
                           exclude_ids=True,
                           min_features=1,
                           max_features=1)

        # allow empty table
        actual = filter_samples(table,
                                metadata=metadata,
                                where=where,
                                exclude_ids=True,
                                min_features=1,
                                max_features=1,
                                allow_empty_table=True)
        self.assertTrue(actual.is_empty())


if __name__ == "__main__":
    unittest.main()
